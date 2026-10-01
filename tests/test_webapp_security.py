from __future__ import annotations

import http.client
import socket
import tempfile
import threading
import unittest
from pathlib import Path
from unittest.mock import patch

from hcrprobeforge import webapp


class WebAppSecurityTests(unittest.TestCase):
    token = "local-test-token-1234567890"

    def _server(self, *, require_auth: bool = False) -> webapp.HCRProbeForgeHTTPServer:
        return webapp.HCRProbeForgeHTTPServer(
            ("127.0.0.1", 0),
            webapp.HCRProbeForgeHandler,
            access_token=self.token,
            require_auth=require_auth,
            allowed_hosts=webapp._allowed_hostnames("127.0.0.1"),
        )

    def _start(self, server: webapp.HCRProbeForgeHTTPServer) -> threading.Thread:
        thread = threading.Thread(target=server.serve_forever, daemon=True)
        thread.start()
        return thread

    def _request(
        self,
        server: webapp.HCRProbeForgeHTTPServer,
        method: str,
        path: str,
        *,
        headers: dict[str, str] | None = None,
        body: bytes | None = None,
    ) -> tuple[int, dict[str, str], bytes]:
        connection = http.client.HTTPConnection("127.0.0.1", server.server_port, timeout=3)
        try:
            connection.request(method, path, body=body, headers={"Connection": "close", **(headers or {})})
            response = connection.getresponse()
            return response.status, dict(response.getheaders()), response.read()
        finally:
            connection.close()

    def test_root_sets_capability_cookie_and_renders_form_token(self):
        server = self._server()
        self._start(server)
        try:
            status, headers, body = self._request(server, "GET", "/")
        finally:
            server.shutdown()
            server.server_close()

        self.assertEqual(status, 200)
        self.assertIn("HCRProbeForge-Token=local-test-token-1234567890", headers["Set-Cookie"])
        self.assertIn("no-store", headers["Cache-Control"])
        self.assertEqual(headers["Pragma"], "no-cache")
        self.assertIn(b"const hcrRequestToken=\"local-test-token-1234567890\"", body)
        self.assertNotIn(b"beforeunload", body)
        self.assertNotIn(b"pagehide", body)

    def test_setup_refreshes_status_when_restored_from_browser_history(self):
        html = webapp.WEBAPP_ENHANCEMENTS
        self.assertIn("window.addEventListener('pageshow',()=>refreshIndexStatus())", html)

    def test_custom_preset_badge_uses_registered_reference_metadata(self):
        preset = webapp.references.SpeciesPreset(
            key="c_elegans",
            display_name="Caenorhabditis elegans",
            scientific_name="Caenorhabditis elegans",
            assembly_name="WBcel235",
            assembly_accession="GCF_000002985.6",
        )
        with tempfile.TemporaryDirectory() as temporary:
            database = Path(temporary) / "annotation.sqlite"
            database.touch()
            metadata = {
                "index_species_alias": "WBcel235",
                "species": "Caenorhabditis elegans",
                "display_name": "Caenorhabditis elegans",
                "scientific_name": "Caenorhabditis elegans",
                "status": "ready",
                "reference_data_directory": temporary,
                "annotation_database_path": str(database),
                "annotation_database_status": "ready",
            }
            with (
                patch.object(
                    webapp.references,
                    "supported_species",
                    return_value=tuple(webapp.references.SPECIES_PRESETS.values()) + (preset,),
                ),
                patch.object(webapp.references, "custom_species_presets", return_value=(preset,)),
                patch.object(webapp.references, "list_installed_references", return_value=[metadata]),
                patch.object(webapp.references, "registered_index_is_ready", return_value=True),
                patch("hcrprobeforge.premrna.annotation_database_is_structurally_ready") as structural_check,
            ):
                html = webapp._render_form(self.token)
        self.assertIn('data-annotation-ready="true" value="c_elegans"', html)
        structural_check.assert_not_called()

    def test_installed_alternate_badge_uses_fast_metadata_check(self):
        with tempfile.TemporaryDirectory() as temporary:
            database = Path(temporary) / "annotation.sqlite"
            database.touch()
            metadata = {
                "index_species_alias": "WBcel235",
                "species": "Caenorhabditis elegans",
                "scientific_name": "Caenorhabditis elegans",
                "assembly": "WBcel235",
                "status": "ready",
                "annotation_database_path": str(database),
                "annotation_database_status": "ready",
            }
            with (
                patch.object(webapp.references, "supported_species", return_value=()),
                patch.object(webapp.references, "custom_species_presets", return_value=()),
                patch.object(webapp.references, "registered_index_is_ready", return_value=True),
                patch("hcrprobeforge.premrna.annotation_database_is_structurally_ready") as structural_check,
            ):
                html = webapp._species_options_html(
                    installed_references=[metadata],
                    annotation_structural_check=False,
                )
        self.assertIn('data-annotation-ready="true" value="WBcel235"', html)
        structural_check.assert_not_called()

    def test_legacy_xtr10_alias_is_not_a_species_option(self):
        self.assertIsNone(webapp.references.get_species_preset("xtr10"))
        self.assertIsNone(webapp.references.find_installed_reference("xtr10"))
        stale_builtin = {
            "index_species_alias": "xtr10",
            "species": "Xenopus tropicalis",
            "scientific_name": "Xenopus tropicalis",
            "assembly": "UCB_Xtro_10.0",
            "assembly_accession": "GCF_000004195.4",
            "status": "ready",
        }
        with (
            patch.object(webapp.references, "registered_index_is_ready", return_value=True),
            patch.object(webapp.references, "custom_species_presets", return_value=()),
        ):
            html = webapp._species_options_html(
                installed_references=[stale_builtin, {"index_species_alias": "xtr10", "status": "ready"}]
            )
        self.assertNotIn('value="xtr10"', html)
        self.assertNotIn("UCB_Xtro_10.0", html)

    def test_index_status_uses_fast_metadata_check(self):
        with tempfile.TemporaryDirectory() as temporary:
            database = Path(temporary) / "annotation.sqlite"
            database.touch()
            metadata = {
                "index_species_alias": "WBcel235",
                "species": "Caenorhabditis elegans",
                "status": "ready",
                "annotation_database_path": str(database),
                "annotation_database_status": "ready",
            }
            server = self._server()
            self._start(server)
            try:
                with (
                    patch.object(webapp.references, "get_species_preset", return_value=None),
                    patch.object(webapp.references, "find_installed_reference", return_value=metadata),
                    patch.object(webapp.references, "registered_index_is_ready", return_value=True),
                    patch("hcrprobeforge.premrna.annotation_database_is_structurally_ready") as structural_check,
                ):
                    status, _, body = self._request(server, "GET", "/index-status?species=WBcel235")
            finally:
                server.shutdown()
                server.server_close()
        self.assertEqual(status, 200)
        payload = __import__("json").loads(body)
        self.assertTrue(payload["ready"])
        self.assertTrue(payload["annotation_ready"])
        structural_check.assert_not_called()

    def test_mutating_request_requires_token_and_rejects_foreign_origin(self):
        server = self._server()
        self._start(server)
        try:
            status, headers, _ = self._request(server, "GET", "/")
            self.assertEqual(status, 200)
            cookie = headers["Set-Cookie"].split(";", 1)[0]

            with patch.object(webapp, "_choose_local_folder", return_value="/tmp/selected") as chooser:
                missing_status, _, _ = self._request(server, "POST", "/choose-folder")
                foreign_status, _, _ = self._request(
                    server,
                    "POST",
                    "/choose-folder",
                    headers={"Cookie": cookie, "Origin": "https://evil.example"},
                )
                valid_status, _, _ = self._request(
                    server,
                    "POST",
                    "/choose-folder",
                    headers={
                        "Cookie": cookie,
                        "Origin": f"http://127.0.0.1:{server.server_port}",
                    },
                )
        finally:
            server.shutdown()
            server.server_close()

        self.assertEqual(missing_status, 403)
        self.assertEqual(foreign_status, 403)
        self.assertEqual(valid_status, 200)
        chooser.assert_called_once_with()

    def test_invalid_host_is_rejected_even_with_token(self):
        server = self._server()
        self._start(server)
        try:
            status, _, _ = self._request(
                server,
                "POST",
                "/choose-folder",
                headers={
                    "Host": f"attacker.example:{server.server_port}",
                    webapp.WEB_TOKEN_HEADER: self.token,
                },
            )
        finally:
            server.shutdown()
            server.server_close()

        self.assertEqual(status, 400)

    def test_remote_mode_requires_token_and_supports_access_url_handshake(self):
        server = self._server(require_auth=True)
        self._start(server)
        try:
            denied_status, _, _ = self._request(server, "GET", "/")
            redirect_status, redirect_headers, _ = self._request(
                server,
                "GET",
                f"/?access_token={self.token}",
            )
            cookie = redirect_headers["Set-Cookie"].split(";", 1)[0]
            allowed_status, _, _ = self._request(
                server,
                "GET",
                "/",
                headers={"Cookie": cookie},
            )
        finally:
            server.shutdown()
            server.server_close()

        self.assertEqual(denied_status, 403)
        self.assertEqual(redirect_status, 303)
        self.assertEqual(allowed_status, 200)

    def test_bind_retries_when_requested_port_is_occupied(self):
        occupied = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        occupied.bind(("127.0.0.1", 0))
        occupied.listen(1)
        requested_port = occupied.getsockname()[1]
        server = webapp._bind_web_server(
            "127.0.0.1",
            requested_port,
            access_token=self.token,
            require_auth=False,
            allowed_hosts=webapp._allowed_hostnames("127.0.0.1"),
        )
        try:
            self.assertNotEqual(server.server_port, requested_port)
            self.assertGreater(server.server_port, requested_port)
        finally:
            server.server_close()
            occupied.close()

    def test_shutdown_and_species_changes_are_rejected_while_work_is_active(self):
        server = self._server()
        self._start(server)
        token = "queued-job"
        with webapp.JOBS_LOCK:
            webapp.JOBS[token] = {"status": "queued"}
        try:
            headers = {webapp.WEB_TOKEN_HEADER: self.token}
            shutdown_status, _, shutdown_body = self._request(
                server, "POST", "/shutdown", headers=headers
            )
            register_status, _, register_body = self._request(
                server, "POST", "/register-species", headers=headers
            )
            delete_status, _, delete_body = self._request(
                server, "POST", "/delete-species", headers=headers
            )
            preview_status, _, preview_body = self._request(
                server, "GET", "/delete-species-preview?species=custom_species"
            )
            alive_status, _, _ = self._request(server, "GET", "/")
        finally:
            with webapp.JOBS_LOCK:
                webapp.JOBS.pop(token, None)
            server.shutdown()
            server.server_close()

        for status in (shutdown_status, register_status, delete_status, preview_status):
            self.assertEqual(status, 409)
        for body in (shutdown_body, register_body, delete_body, preview_body):
            self.assertIn(b"queued or running", body)
        self.assertEqual(alive_status, 200)

    def test_request_staging_directories_are_removed_on_failure(self):
        captured: list[Path] = []

        def fail(*args, **kwargs):
            captured.append(Path(kwargs.get("job_dir", "")))
            raise RuntimeError("synthetic failure")

        with patch.object(webapp, "_run_submission_in_directory", side_effect=fail):
            with self.assertRaises(RuntimeError):
                webapp._run_submission_inner({}, [])

        self.assertEqual(len(captured), 1)
        self.assertFalse(captured[0].exists())

        before = set(Path(tempfile.gettempdir()).glob("hcrprobeforge-inspect-*"))
        with self.assertRaises(ValueError):
            webapp._inspect_submission({"mode": "not-a-workflow"}, [])
        after = set(Path(tempfile.gettempdir()).glob("hcrprobeforge-inspect-*"))
        self.assertEqual(after, before)

    def test_premrna_validation_is_not_species_alias_specific(self):
        fields = {
            "mode": "design",
            "species": "xla",
            "organism": "Xenopus laevis",
            "gene": "sox9",
            "target_type": "pre-mrna",
            "auto_curate": "",
        }
        with patch.object(webapp.references, "registered_index_is_ready", return_value=True):
            webapp._validate_form_values(fields)

    def test_index_annotation_database_is_an_explicit_opt_in(self):
        self.assertFalse(webapp._field_enabled({"index_annotation_database": "0"}, "index_annotation_database"))
        self.assertTrue(webapp._field_enabled({"index_annotation_database": "1"}, "index_annotation_database"))
        html = webapp._render_form(self.token)
        self.assertIn("Prepare this index for intron design?", html)
        self.assertIn("first Pre-mRNA run", html)
        self.assertIn('placeholder="Optional, versioned GCF_ or GCA_ accession"', html)
        self.assertNotIn('GCF_000002035.6', html)

        common = {
            "mode": "index",
            "species": "xtr",
            "index_threads": "1",
            "index_assembly_accession": "GCF_000004195.4",
            "index_assembly_name": "UCB_Xtro_10.0",
        }
        deferred = webapp._inspect_submission({**common, "index_annotation_database": "0"}, [])
        requested = webapp._inspect_submission({**common, "index_annotation_database": "1"}, [])
        self.assertIn("deferred until the first Pre-mRNA run", deferred["message"])
        self.assertIn("also be built now", requested["message"])

    def test_annotation_badge_recognizes_registered_alternate_alias(self):
        metadata = {
            "index_species_alias": "c_elegans",
            "species": "c_elegans",
            "status": "ready",
            "reference_data_directory": "/references/c_elegans/WBcel235",
            "annotation_database_path": "/references/c_elegans/WBcel235/annotation.sqlite",
        }
        with (
            patch.object(webapp.references, "get_species_preset", return_value=None),
            patch.object(webapp.references, "find_installed_reference", return_value=metadata),
            patch.object(webapp.references, "registered_index_is_ready", return_value=True),
            patch("hcrprobeforge.premrna.annotation_database_is_structurally_ready", return_value=True),
        ):
            self.assertTrue(webapp._annotation_database_ready("c_elegans"))


if __name__ == "__main__":
    unittest.main()
