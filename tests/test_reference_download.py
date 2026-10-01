from __future__ import annotations

import gzip
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

import hcrprobeforge.references as references


class FakeResponse:
    def __init__(self, *, status_code=200, payload=None, content=b"", headers=None):
        self.status_code = status_code
        self._payload = payload
        self._content = content
        self.headers = headers or {}
        self.text = content.decode("utf-8", errors="replace")

    def raise_for_status(self):
        if self.status_code >= 400:
            raise RuntimeError(f"HTTP {self.status_code}")

    def json(self):
        return self._payload

    def iter_content(self, chunk_size=1024 * 1024):
        del chunk_size
        yield self._content


class ReferenceDownloadTests(unittest.TestCase):
    def test_assembly_filename_uses_ftp_directory_convention(self):
        record = {
            "AssemblyAccession": "GCF_000001215.4",
            "AssemblyName": "Release 6 plus ISO1 MT",
            "FtpPath_RefSeq": (
                "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/215/"
                "GCF_000001215.4_Release_6_plus_ISO1_MT"
            ),
        }

        urls = references._assembly_download_urls(
            record, "GCF_000001215.4", "Release 6 plus ISO1 MT"
        )

        self.assertTrue(
            urls[0].endswith(
                "GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna.gz"
            )
        )
        self.assertNotIn("%20", urls[0])

    def test_assembly_record_rejects_first_uid_when_summary_does_not_match(self):
        class Session:
            def get(self, endpoint, params, timeout):
                del timeout
                if endpoint.endswith("esearch.fcgi"):
                    return FakeResponse(payload={"esearchresult": {"idlist": ["1", "2"]}})
                if params["id"] == "1":
                    return FakeResponse(
                        payload={"result": {"uids": ["1"], "1": {"AssemblyAccession": "GCF_9.1"}}}
                    )
                return FakeResponse(
                    payload={
                        "result": {
                            "uids": ["2"],
                            "2": {
                                "AssemblyAccession": "GCF_000001215.4",
                                "AssemblyName": "Release 6 plus ISO1 MT",
                                "FtpPath_RefSeq": "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/215/GCF_000001215.4_Release_6_plus_ISO1_MT",
                            },
                        }
                    }
                )

        record = references._assembly_record(
            "gcf_000001215.4", session=Session(), email=None, api_key=None
        )

        self.assertEqual(record["AssemblyAccession"], "GCF_000001215.4")
        self.assertEqual(record["uid"], "2")

    def test_download_retries_transient_error_and_replaces_atomically(self):
        payload = gzip.compress(b">contig\nACGT\n")

        class Session:
            def __init__(self):
                self.calls = 0

            def get(self, url, stream, timeout):
                del url, stream, timeout
                self.calls += 1
                if self.calls == 1:
                    return FakeResponse(status_code=503)
                return FakeResponse(
                    content=payload,
                    headers={"Content-Length": str(len(payload))},
                )

        session = Session()
        with TemporaryDirectory() as directory, patch.object(references.time, "sleep", lambda seconds: None):
            destination = Path(directory) / "genome.fna.gz"
            references._download(
                "https://example.invalid/genome.fna.gz", destination, session=session
            )
            self.assertEqual(session.calls, 2)
            self.assertEqual(destination.read_bytes(), payload)
            self.assertEqual(list(Path(directory).glob("*.part.*")), [])

    def test_directory_listing_is_a_last_resort_genome_fallback(self):
        record = {
            "FtpPath_RefSeq": (
                "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/215/"
                "GCF_000001215.4_unknown"
            )
        }

        class Session:
            def get(self, url, timeout):
                self.url = url
                del timeout
                return FakeResponse(
                    content=b'<a href="actual_genomic.fna.gz">actual_genomic.fna.gz</a>'
                )

        session = Session()
        urls = references._discover_assembly_download_urls(record, session=session)

        self.assertEqual(
            urls,
            [
                "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/215/"
                "GCF_000001215.4_unknown/actual_genomic.fna.gz"
            ],
        )

    def test_fetch_passes_validated_genome_to_existing_index_builder(self):
        record = {
            "AssemblyAccession": "GCF_000001215.4",
            "AssemblyName": "Release 6 plus ISO1 MT",
            "FtpPath_RefSeq": (
                "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/215/"
                "GCF_000001215.4_Release_6_plus_ISO1_MT"
            ),
        }
        payload = gzip.compress(b">contig\nACGT\n")

        class RequestsModule:
            class Session:
                def __init__(self):
                    self.headers = {}

        def write_download(url, destination, **kwargs):
            del url, kwargs
            destination.write_bytes(payload)

        with TemporaryDirectory() as directory:
            root = Path(directory)
            with (
                patch.object(references.core, "requests", RequestsModule),
                patch.object(references, "species_data_root", return_value=root / "refs"),
                patch.object(references, "hcrprobedesign_data_root", return_value=root / "design"),
                patch.object(references, "_assembly_record", return_value=record),
                patch.object(references, "_download", side_effect=write_download),
                patch.object(references, "_run_build_genome_index", return_value=(0, "", "")) as build,
            ):
                result = references.fetch_and_build_index(
                    "xtr",
                    assembly_accession="gcf_000001215.4",
                    assembly_name="Release 6 plus ISO1 MT",
                )

            self.assertTrue(
                result["genome_url"].endswith(
                    "GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna.gz"
                )
            )
            self.assertEqual(build.call_args.kwargs["species"], "xtr_Release_6_plus_ISO1_MT")

    def test_local_build_persists_uncompressed_genome_annotation_and_lookup(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            fasta = root / "input.fna"
            fasta.write_text(">chr1\nACGTACGTACGT\n", encoding="utf-8")
            annotation = root / "input.gff3"
            annotation.write_text(
                "chr1\tncbi\tmRNA\t1\t12\t.\t+\t.\tID=rna1;transcript_id=NM_1.1;gene=sox9\n"
                "chr1\tncbi\texon\t1\t3\t.\t+\t.\tParent=rna1\n"
                "chr1\tncbi\texon\t10\t12\t.\t+\t.\tParent=rna1\n",
                encoding="utf-8",
            )
            with (
                patch.object(references, "species_data_root", return_value=root / "refs"),
                patch.object(references, "hcrprobedesign_data_root", return_value=root / "design"),
                patch.object(references, "_run_build_genome_index", return_value=(0, "", "")),
            ):
                result = references.build_local_index(
                    "xtr",
                    fasta,
                    annotation=annotation,
                    assembly="TestAssembly",
                )
            reference_dir = root / "refs" / "xtr" / "TestAssembly"
            self.assertEqual(result["status"], "ready")
            self.assertTrue((reference_dir / "genome.fna").is_file())
            self.assertTrue((reference_dir / "annotation.gff").is_file())
            self.assertFalse((reference_dir / "annotation.sqlite").exists())
            self.assertEqual(list(reference_dir.glob("*.tmp.*")), [])
            self.assertEqual(list(reference_dir.glob("*.gz")), [])

            with (
                patch.object(references, "species_data_root", return_value=root / "refs"),
                patch.object(references, "hcrprobedesign_data_root", return_value=root / "design"),
                patch.object(references, "_run_build_genome_index", return_value=(0, "", "")),
            ):
                references.build_local_index(
                    "xtr",
                    fasta,
                    annotation=annotation,
                    assembly="TestAssembly",
                    force=True,
                    build_annotation_database=True,
                )
            self.assertTrue((reference_dir / "annotation.sqlite").is_file())

    def test_local_build_restores_existing_reference_after_commit_failure(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            fasta = root / "new.fna"
            fasta.write_text(">chr1\nACGTACGT\n", encoding="utf-8")
            reference_dir = root / "refs" / "xtr" / "TestAssembly"
            reference_dir.mkdir(parents=True)
            existing = reference_dir / "genome.fna"
            existing.write_text(">chr1\nTTTT\n", encoding="utf-8")

            original_write_metadata = references._write_metadata
            calls = {"count": 0}

            def write_then_fail(path, payload):
                calls["count"] += 1
                if calls["count"] == 2:
                    raise OSError("simulated metadata commit failure")
                return original_write_metadata(path, payload)

            with (
                patch.object(references, "species_data_root", return_value=root / "refs"),
                patch.object(references, "hcrprobedesign_data_root", return_value=root / "design"),
                patch.object(references, "_run_build_genome_index", return_value=(0, "", "")),
                patch.object(references, "_write_metadata", side_effect=write_then_fail),
                self.assertRaises(OSError),
            ):
                references.build_local_index("xtr", fasta, assembly="TestAssembly", force=True)

            self.assertEqual(existing.read_text(encoding="utf-8"), ">chr1\nTTTT\n")
