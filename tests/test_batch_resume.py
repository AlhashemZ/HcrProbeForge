import json
import tempfile
import unittest
from pathlib import Path

from hcrprobeforge import batch


class ManifestResumeMarkerTests(unittest.TestCase):
    def test_marker_requires_matching_parameter_signature(self):
        payload = batch._manifest_resume_payload(
            extra_args=("--target-probes", "20", "--qc-stringency", "balanced"),
            channels=("B1", "B3"),
            smart_default=True,
            species="xtr",
        )
        signature = batch._manifest_resume_signature(payload)
        changed_payload = dict(payload)
        changed_payload["extra_args"] = ["--target-probes", "24", "--qc-stringency", "balanced"]
        changed_signature = batch._manifest_resume_signature(changed_payload)

        with tempfile.TemporaryDirectory() as temporary:
            marker = Path(temporary) / "target.done"
            marker.write_text(
                json.dumps({"signature": signature, "parameters": payload}) + "\n",
                encoding="utf-8",
            )
            self.assertTrue(batch._manifest_done_marker_matches(marker, signature))
            self.assertFalse(batch._manifest_done_marker_matches(marker, changed_signature))

    def test_legacy_empty_marker_is_not_reused(self):
        with tempfile.TemporaryDirectory() as temporary:
            marker = Path(temporary) / "target.done"
            marker.touch()
            self.assertFalse(batch._manifest_done_marker_matches(marker, "signature"))

    def test_only_true_species_mismatch_is_skipped(self):
        self.assertEqual(
            batch._classify_input_failure("record does not belong to the selected organism"),
            "skipped_species_mismatch",
        )
        self.assertEqual(
            batch._classify_input_failure("No NCBI gene match was found"),
            "failed",
        )
        self.assertIsNone(batch._classify_input_failure("temporary NCBI failure"))


if __name__ == "__main__":
    unittest.main()
