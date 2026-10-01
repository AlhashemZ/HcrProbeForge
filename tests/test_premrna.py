from __future__ import annotations

import unittest
import sqlite3
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
from unittest.mock import patch

from hcrprobeforge import core
from hcrprobeforge import references
from hcrprobeforge import premrna
from hcrprobeforge.premrna import (
    _accession_compatible,
    _accession_matches,
    _decompress_once,
    annotation_database_is_valid,
    build_annotation_index,
    parse_gff3_models,
    parse_intron_selection,
    prepare_target,
    resolve_reference_metadata,
)


class PremrnaTests(unittest.TestCase):
    def test_gff3_links_exons_to_refseq_transcript(self):
        with TemporaryDirectory() as directory:
            path = Path(directory) / "annotation.gff"
            path.write_text(
                "##gff-version 3\n"
                "chr1\tncbi\tmRNA\t10\t100\t.\t+\t.\tID=rna1;transcript_id=XM_123.1;gene=sox10\n"
                "chr1\tncbi\texon\t10\t20\t.\t+\t.\tParent=rna1\n"
                "chr1\tncbi\texon\t40\t50\t.\t+\t.\tParent=rna1\n"
                "chr1\tncbi\texon\t80\t100\t.\t+\t.\tParent=rna1\n",
                encoding="utf-8",
            )
            models = parse_gff3_models(path)
        self.assertEqual(len(models), 1)
        self.assertEqual(models[0]["aliases"][0], "rna1")
        self.assertEqual(models[0]["exons"], [{"start": 10, "end": 20}, {"start": 40, "end": 50}, {"start": 80, "end": 100}])

    def test_intron_selection_is_positive_unique_and_optional(self):
        self.assertIsNone(parse_intron_selection(None))
        self.assertEqual(parse_intron_selection("1, 3, 1"), (1, 3))
        with self.assertRaises(ValueError):
            parse_intron_selection("0,2")

    def test_premrna_target_folder_names_identify_intron_scope(self):
        args = SimpleNamespace(gene="sox9")
        base = {
            "source": "ncbi_premrna",
            "accession": "NM_001016853.2",
        }
        all_introns = dict(base, premrna_target={
            "selection_mode": "intronic",
            "distribution": "even_by_intron",
            "selected_introns": [1, 2],
        })
        selected = dict(base, premrna_target={
            "selection_mode": "intronic",
            "distribution": "selected_introns",
            "selected_introns": [2, 3],
        })
        whole = dict(base, premrna_target={
            "selection_mode": "whole",
            "distribution": "whole_target",
            "selected_introns": [],
        })
        self.assertTrue(core.target_label_for_record(args, all_introns).endswith("_premrna_introns"))
        self.assertTrue(core.target_label_for_record(args, selected).endswith("_premrna_introns2,3"))
        self.assertTrue(core.target_label_for_record(args, whole).endswith("_premrna_whole"))

    def test_multi_intron_scope_is_sanitized_only_for_designprobes(self):
        args = SimpleNamespace(
            species="xtr",
            channel="B1",
            tile_size=52,
            min_gc=45,
            max_gc=55,
            min_gibbs=-70,
            max_gibbs=-50,
            target_gibbs=-60,
            max_run_mismatches=2,
            max_probes=30,
            num_hits_allowed=1,
            index=None,
            no_genomemask=False,
            dtm_filter=False,
            dtm_max=5,
        )
        target = "sox9_NM_001016853.2_premrna_introns1,2"
        command = core.build_designprobes_command(
            "designProbes",
            Path("target.fa"),
            target,
            Path("probes.tsv"),
            Path("idt.tsv"),
            args,
        )
        target_name = command[command.index("--targetName") + 1]
        self.assertEqual(target_name, "sox9_NM_001016853.2_premrna_introns1_2")
        self.assertEqual(core.target_label_for_record(SimpleNamespace(gene="sox9"), {
            "source": "ncbi_premrna",
            "accession": "NM_001016853.2",
            "premrna_target": {
                "selection_mode": "intronic",
                "distribution": "selected_introns",
                "selected_introns": [1, 2],
            },
        }), target)

    def test_transcript_model_cache_is_human_readable_json(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            gff = root / "annotation.gff3"
            gff.write_text(
                "chr1\tncbi\tmRNA\t1\t12\t.\t+\t.\tID=rna1;transcript_id=NM_123.1;gene=sox9\n"
                "chr1\tncbi\texon\t1\t3\t.\t+\t.\tParent=rna1\n"
                "chr1\tncbi\texon\t10\t12\t.\t+\t.\tParent=rna1\n",
                encoding="utf-8",
            )
            with patch.dict("os.environ", {"HCRPROBEFORGE_CACHE_DIR": str(root / "cache")}):
                from hcrprobeforge.premrna import _load_or_parse_target_models

                _load_or_parse_target_models(gff, {"accession": "NM_123.1"})
            cached = next((root / "cache").rglob("transcript_NM_123.1_model.json"))
            text = cached.read_text(encoding="utf-8")
        self.assertTrue(text.startswith("{\n"))
        self.assertIn("\n  \"models\": [", text)

    def test_intron_annotation_and_balanced_selector(self):
        annotation = {
            "regions": [
                {"kind": "intron", "intron_number": 1, "start": 1, "end": 100, "genomic_start": 21, "genomic_end": 40},
                {"kind": "separator", "start": 101, "end": 152},
                {"kind": "intron", "intron_number": 2, "start": 153, "end": 252, "genomic_start": 51, "genomic_end": 80},
            ]
        }
        rows = [
            {"candidate_id": "a", "start": 5, "end": 56, "length": 52, "auto_qc_reject_reason": ""},
            {"candidate_id": "b", "start": 160, "end": 211, "length": 52, "auto_qc_reject_reason": ""},
            {"candidate_id": "cross", "start": 105, "end": 156, "length": 52, "auto_qc_reject_reason": ""},
        ]
        core.annotate_probe_rows(rows, annotation)
        self.assertEqual([row["transcript_region"] for row in rows], ["intron_1", "intron_2", "intron_boundary"])
        selected = core.select_premrna_candidates(
            rows,
            2,
            set(),
            sequence_length=252,
            record={"transcript_annotation": annotation},
        )
        self.assertEqual([row["premrna_intron_number"] for row in selected], [1, 2])

    def test_intronic_selector_blocks_intervals_selected_by_quota_pass(self):
        annotation = {
            "regions": [
                {"kind": "intron", "intron_number": 1, "start": 1, "end": 100},
                {"kind": "separator", "start": 101, "end": 100},
                {"kind": "intron", "intron_number": 2, "start": 101, "end": 200},
            ]
        }
        rows = [
            {"candidate_id": "i1a", "start": 1, "end": 20, "length": 20, "auto_qc_reject_reason": ""},
            {"candidate_id": "i1b", "start": 30, "end": 49, "length": 20, "auto_qc_reject_reason": ""},
            {"candidate_id": "i1c", "start": 60, "end": 79, "length": 20, "auto_qc_reject_reason": ""},
            {"candidate_id": "i2a", "start": 101, "end": 120, "length": 20, "auto_qc_reject_reason": ""},
            {"candidate_id": "i2b", "start": 105, "end": 124, "length": 20, "auto_qc_reject_reason": ""},
        ]
        selected = core.select_premrna_candidates(
            rows,
            4,
            set(),
            sequence_length=200,
            record={"transcript_annotation": annotation},
        )
        intervals = sorted((int(row["start"]), int(row["end"])) for row in selected)
        self.assertEqual(len(selected), 4)
        self.assertEqual(intervals, [(1, 20), (30, 49), (60, 79), (101, 120)])
        self.assertTrue(
            all(left_end < right_start for (_, left_end), (right_start, _) in zip(intervals, intervals[1:]))
        )

    def test_probe_name_exclusive_coordinate_suffix_is_normalized(self):
        self.assertEqual(
            core.normalize_probe_name_coordinates("candidate_101-121", 101, 120, 20),
            "candidate_101-120",
        )
        self.assertEqual(
            core.normalize_probe_name_coordinates("user_label", 101, 120, 20),
            "user_label",
        )

    def test_one_pass_premrna_order_table_is_nonoverlapping(self):
        rows = [
            {"start": 1, "end": 20, "length": 20, "GC": 50, "dTm": 1, "GibbsFE": -60},
            {"start": 15, "end": 34, "length": 20, "GC": 50, "dTm": 1, "GibbsFE": -60},
            {"start": 40, "end": 59, "length": 20, "GC": 50, "dTm": 1, "GibbsFE": -60},
        ]
        selected = core.enforce_premrna_nonoverlap(rows, 100)
        self.assertEqual([(row["start"], row["end"]) for row in selected], [(1, 20), (40, 59)])

    def test_premrna_intron_capacity_counts_nonoverlapping_tiles(self):
        record = {
            "transcript_annotation": {
                "regions": [
                    {"kind": "exon", "start": 1, "end": 9},
                    {"kind": "intron", "intron_number": 1, "start": 10, "end": 69},
                    {"kind": "exon", "start": 70, "end": 79},
                    {"kind": "intron", "intron_number": 2, "start": 80, "end": 139},
                ]
            }
        }
        sequence = "A" * 139
        self.assertEqual(core.premrna_intron_capacity(record, sequence, 52), 2)
        self.assertEqual(core.premrna_intron_capacity(record, sequence, 36), 2)

    def test_premrna_intron_capacity_respects_ambiguous_segments(self):
        record = {
            "transcript_annotation": {
                "regions": [
                    {"kind": "intron", "intron_number": 1, "start": 1, "end": 100},
                ]
            }
        }
        sequence = "A" * 52 + "N" + "A" * 47
        self.assertEqual(core.premrna_intron_capacity(record, sequence, 52), 1)

    def test_ready_hcrprobedesign_index_supplies_missing_forge_metadata(self):
        with patch.object(references, "find_installed_reference", return_value=None), \
             patch.object(references, "registered_index_is_ready", return_value=True), \
             patch.object(references, "registered_index_prefix", return_value=Path("/tmp/xtr/xtr")):
            metadata = resolve_reference_metadata(references, "xtr")

        self.assertIsNotNone(metadata)
        self.assertTrue(metadata["metadata_fallback"])
        self.assertEqual(metadata["assembly"], "UCB_Xtro_10.0")
        self.assertEqual(metadata["assembly_accession"], "GCF_000004195.4")
        self.assertEqual(metadata["index_prefix"], "/tmp/xtr/xtr")

    def test_annotation_index_handles_exons_before_transcript_and_versioned_alias(self):
        with TemporaryDirectory() as directory:
            path = Path(directory) / "annotation.gff3"
            path.write_text(
                "##gff-version 3\n"
                "chr1\tncbi\texon\t1\t3\t.\t+\t.\tParent=rna1\n"
                "chr1\tncbi\texon\t10\t12\t.\t+\t.\tParent=rna1\n"
                "chr1\tncbi\tmRNA\t1\t12\t.\t+\t.\tID=rna1;transcript_id=NM_123.1;gene=sox9\n",
                encoding="utf-8",
            )
            database = build_annotation_index(path)
            from hcrprobeforge.premrna import parse_gff3_model_for_record

            models = parse_gff3_model_for_record(
                path,
                {"accession": "NM_123.1", "gene_symbol": "sox9"},
                database_path=database,
            )
        self.assertEqual(len(models), 1)
        self.assertEqual(models[0]["exons"], [{"start": 1, "end": 3}, {"start": 10, "end": 12}])

    def test_annotation_index_reports_byte_progress_and_cleans_sqlite_sidecars(self):
        with TemporaryDirectory() as directory:
            path = Path(directory) / "annotation.gff3"
            path.write_text(
                "chr1\tncbi\tmRNA\t1\t12\t.\t+\t.\tID=rna1;transcript_id=NM_123.1;gene=sox9\n"
                "chr1\tncbi\texon\t1\t3\t.\t+\t.\tParent=rna1\n"
                "chr1\tncbi\texon\t10\t12\t.\t+\t.\tParent=rna1\n",
                encoding="utf-8",
            )
            updates = []
            database = build_annotation_index(path, progress_callback=updates.append)
            self.assertTrue(database.is_file())
            self.assertTrue(annotation_database_is_valid(database))
            self.assertTrue(any(update.get("total") == path.stat().st_size for update in updates))
            self.assertEqual(list(Path(directory).glob("annotation.sqlite*")), [database])

    def test_version_matching_rejects_different_versions_but_accepts_unversioned_annotation(self):
        self.assertTrue(_accession_matches("NM_123.1", "NM_123.1"))
        self.assertFalse(_accession_matches("NM_123.2", "NM_123.1"))
        self.assertFalse(_accession_compatible("NM_123.1", "NM_123.2"))
        self.assertTrue(_accession_compatible("NM_123", "NM_123.2"))

    def test_gene_id_fallback_requires_an_exact_match(self):
        model = {"aliases": [], "gene_id": "12345", "gene_symbol": "gene1"}
        record = {"accession": "", "gene_id": "123", "gene_symbol": "other"}
        self.assertEqual(premrna._target_model_match(model, record), (False, False))

    def test_output_rollback_restores_overwritten_existing_files(self):
        with TemporaryDirectory() as directory:
            root = Path(directory) / "runs"
            existing = root / "target" / "summary.json"
            existing.parent.mkdir(parents=True)
            existing.write_text("old\n", encoding="utf-8")
            before = core._snapshot_filesystem(root)
            backup_root, backups = core._backup_existing_filesystem_contents([root])
            existing.write_text("new\n", encoding="utf-8")
            (root / "target" / "partial.tsv").write_text("partial\n", encoding="utf-8")
            core._remove_new_filesystem_entries(root, before)
            core._restore_existing_filesystem_contents(backup_root, backups)
            self.assertEqual(existing.read_text(encoding="utf-8"), "old\n")
            self.assertFalse((root / "target" / "partial.tsv").exists())

    def test_external_compressed_input_is_not_deleted(self):
        with TemporaryDirectory() as directory:
            import gzip

            root = Path(directory)
            source = root / "user_genome.fna.gz"
            destination = root / "managed" / "genome.fna"
            source.parent.mkdir(parents=True, exist_ok=True)
            with gzip.open(source, "wb") as handle:
                handle.write(b">chr1\nACGT\n")
            _decompress_once(source, destination)
            self.assertTrue(source.is_file())
            self.assertEqual(destination.read_text(encoding="utf-8"), ">chr1\nACGT\n")

    def test_annotation_database_validator_rejects_interrupted_staging_database(self):
        with TemporaryDirectory() as directory:
            database = Path(directory) / "annotation.sqlite"
            with sqlite3.connect(database) as connection:
                connection.execute("CREATE TABLE transcripts(id INTEGER PRIMARY KEY)")
                connection.execute("CREATE TABLE aliases(alias TEXT, transcript_id INTEGER)")
                connection.execute("CREATE TABLE exons(transcript_id INTEGER, seqid TEXT, start INTEGER, end INTEGER)")
                connection.execute("CREATE TABLE pending_exons(parent TEXT)")
            self.assertFalse(annotation_database_is_valid(database))

    def test_legacy_premrna_assets_are_migrated_to_the_assembly_directory(self):
        with TemporaryDirectory() as directory:
            root = Path(directory) / "UCB_Xtro_10.0"
            legacy = root / "premrna"
            legacy.mkdir(parents=True)
            (legacy / "genome.fna").write_text(">chr1\nACGT\n", encoding="utf-8")
            (legacy / "annotation.gff").write_text("##gff-version 3\n", encoding="utf-8")
            metadata = {
                "display_name": "Xenopus tropicalis",
                "assembly": "UCB_Xtro_10.0",
                "assembly_accession": "GCF_TEST.1",
                "reference_data_directory": str(legacy),
            }
            args = SimpleNamespace(email=None, api_key=None)
            from hcrprobeforge.premrna import _download_sources

            genome, annotation, _source = _download_sources(None, metadata, "xtr", args)
            self.assertEqual(genome, root / "genome.fna")
            self.assertEqual(annotation, root / "annotation.gff")
            self.assertTrue(genome.is_file())
            self.assertTrue(annotation.is_file())
            self.assertFalse(legacy.exists())

    def test_prepare_target_is_full_genomic_transcript_strand_aware_and_cached(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            genome = root / "genome.fna"
            genome.write_text(">chr1\nACGTACGTACGTACGTACGTAC\n", encoding="utf-8")
            annotation = root / "annotation.gff3"
            annotation.write_text(
                "##gff-version 3\n"
                "chr1\tncbi\tmRNA\t1\t22\t.\t-\t.\tID=rna1;transcript_id=NM_123.1;gene=sox9\n"
                "chr1\tncbi\texon\t1\t3\t.\t-\t.\tParent=rna1\n"
                "chr1\tncbi\texon\t10\t12\t.\t-\t.\tParent=rna1\n"
                "chr1\tncbi\texon\t20\t22\t.\t-\t.\tParent=rna1\n",
                encoding="utf-8",
            )
            source = {
                "reference_data_directory": str(root),
                "assembly_accession": "GCF_TEST.1",
            }
            args = SimpleNamespace(species="xla", index=None, premrna_introns="2")
            metadata = {
                "status": "ready",
                "species": "xla",
                "assembly": "test",
                "assembly_accession": "GCF_TEST.1",
                "reference_data_directory": str(root),
            }
            first = {
                "accession": "NM_123.1",
                "gene_symbol": "sox9",
                "length": 9,
                "title": "sox9",
            }
            second = dict(first)
            with (
                patch("hcrprobeforge.premrna.resolve_reference_metadata", return_value=metadata),
                patch("hcrprobeforge.premrna._download_sources", return_value=(genome, annotation, source)),
                patch("hcrprobeforge.references.get_species_preset", return_value=references.get_species_preset("xla")),
            ):
                prepared = prepare_target(None, first, args)
                with patch("hcrprobeforge.premrna._open_fasta", side_effect=AssertionError("cache should avoid FASTA extraction")):
                    cached = prepare_target(None, second, args)

        self.assertEqual(prepared["source"], "ncbi_premrna")
        self.assertEqual(prepared["premrna_target"]["target_representation"], "full_genomic_transcript")
        self.assertEqual(prepared["length"], 22)
        self.assertEqual(prepared["premrna_target"]["selected_introns"], [2])
        self.assertNotIn("N", prepared["_sequence"])
        regions = prepared["transcript_annotation"]["regions"]
        self.assertEqual([region["kind"] for region in regions], ["exon", "intron", "exon", "intron", "exon"])
        self.assertEqual([region["eligible"] for region in regions if region["kind"] == "intron"], [False, True])
        self.assertEqual(cached["_sequence"], prepared["_sequence"])

    def test_prepare_target_whole_mode_allows_exons_and_uses_separate_cache(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            genome = root / "genome.fna"
            genome.write_text(">chr1\nACGTACGTACGTACGTACGTAC\n", encoding="utf-8")
            annotation = root / "annotation.gff3"
            annotation.write_text(
                "chr1\tncbi\tmRNA\t1\t22\t.\t+\t.\tID=rna1;transcript_id=NM_123.1;gene=sox9\n"
                "chr1\tncbi\texon\t1\t3\t.\t+\t.\tParent=rna1\n"
                "chr1\tncbi\texon\t10\t12\t.\t+\t.\tParent=rna1\n"
                "chr1\tncbi\texon\t20\t22\t.\t+\t.\tParent=rna1\n",
                encoding="utf-8",
            )
            metadata = {
                "status": "ready",
                "species": "xtr",
                "display_name": "Xenopus tropicalis",
                "scientific_name": "Xenopus tropicalis",
                "assembly": "test",
                "assembly_accession": "GCF_TEST.1",
                "reference_data_directory": str(root),
            }
            source = {"reference_data_directory": str(root), "assembly": "test", "assembly_accession": "GCF_TEST.1"}
            args = SimpleNamespace(species="xtr", index=None, premrna_introns=None, target_type="pre-mrna-whole")
            record = {"accession": "NM_123.1", "gene_symbol": "sox9", "length": 9, "title": "sox9"}
            with (
                patch("hcrprobeforge.premrna.resolve_reference_metadata", return_value=metadata),
                patch("hcrprobeforge.premrna._download_sources", return_value=(genome, annotation, source)),
                patch("hcrprobeforge.references.get_species_preset", return_value=references.get_species_preset("xtr")),
                patch.dict("os.environ", {"HCRPROBEFORGE_CACHE_DIR": str(root / "cache")}),
            ):
                prepared = prepare_target(None, record, args)

        self.assertEqual(prepared["premrna_target"]["selection_mode"], "whole")
        self.assertEqual(prepared["premrna_target"]["mode"], "pre-mRNA-whole")
        self.assertEqual(prepared["length"], 22)
        self.assertEqual([region["kind"] for region in prepared["transcript_annotation"]["regions"]], ["exon", "intron", "exon", "intron", "exon"])
        self.assertTrue(all(region.get("eligible") for region in prepared["transcript_annotation"]["regions"] if region["kind"] == "intron"))

    def test_prepare_target_accepts_ready_local_reference_without_species_preset(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            genome = root / "genome.fna"
            genome.write_text(">chr1\nACGTACGTACGTACGTACGTAC\n", encoding="utf-8")
            annotation = root / "annotation.gff3"
            annotation.write_text(
                "chr1\tncbi\tmRNA\t1\t22\t.\t+\t.\tID=rna1;transcript_id=NM_123.1;gene=sox9\n"
                "chr1\tncbi\texon\t1\t3\t.\t+\t.\tParent=rna1\n"
                "chr1\tncbi\texon\t10\t12\t.\t+\t.\tParent=rna1\n",
                encoding="utf-8",
            )
            metadata = {
                "status": "ready",
                "species": "local_xtr",
                "display_name": "Local Xenopus reference",
                "assembly": "LocalAssembly",
                "assembly_accession": None,
                "reference_data_directory": str(root),
            }
            source = {"reference_data_directory": str(root), "assembly": "LocalAssembly"}
            args = SimpleNamespace(
                species="local_xtr",
                index=None,
                premrna_introns=None,
                target_type="pre-mrna",
            )
            record = {"accession": "NM_123.1", "gene_symbol": "sox9", "length": 6}
            with (
                patch("hcrprobeforge.premrna.resolve_reference_metadata", return_value=metadata),
                patch("hcrprobeforge.premrna._download_sources", return_value=(genome, annotation, source)),
                patch("hcrprobeforge.references.get_species_preset", return_value=None),
                patch.dict("os.environ", {"HCRPROBEFORGE_CACHE_DIR": str(root / "cache")}),
            ):
                prepared = prepare_target(None, record, args)

        self.assertEqual(prepared["source"], "ncbi_premrna")
        self.assertEqual(prepared["premrna_target"]["selection_mode"], "intronic")

    def test_whole_premrna_filter_keeps_exon_intron_and_boundary_candidates(self):
        annotation = {
            "selection_mode": "whole",
            "regions": [
                {"kind": "exon", "exon_number": 1, "start": 1, "end": 20},
                {"kind": "intron", "intron_number": 1, "start": 21, "end": 50, "genomic_start": 21, "genomic_end": 50},
                {"kind": "exon", "exon_number": 2, "start": 51, "end": 70},
            ],
        }
        rows = [
            {"candidate_id": "exon", "start": 2, "end": 12, "length": 11},
            {"candidate_id": "intron", "start": 25, "end": 35, "length": 11},
            {"candidate_id": "boundary", "start": 16, "end": 26, "length": 11},
        ]
        filtered = core.filter_premrna_candidates(
            rows,
            {"length": 70, "source": "ncbi_premrna", "premrna_target": annotation, "transcript_annotation": annotation},
        )
        self.assertEqual([row["transcript_region"] for row in filtered], ["exon_1", "intron_1", "intron_boundary"])

    def test_premrna_tier_qc_uses_only_eligible_idt_rows(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            probes = root / "probes.tsv"
            idt = root / "native_IDT.tsv"
            log = root / "design.log"
            params = root / "run_parameters.json"
            probes.write_text(
                "start\tlength\tGC\tTm\tdTm\tGibbsFE\tP1\tP2\n"
                "1\t10\t50\t60\t1\t-60\tACGTACGTAC\tTGCATGCATG\n"
                "15\t10\t50\t60\t1\t-60\tACGTACGTAA\tTTGCATGCAT\n",
                encoding="utf-8",
            )
            idt.write_text("name,sequence\nraw_P1,ACGT\nraw_P2,TGCA\n", encoding="utf-8")
            annotation = {
                "selection_mode": "intronic",
                "regions": [{"kind": "intron", "intron_number": 1, "start": 1, "end": 10, "eligible": True}],
            }
            args = SimpleNamespace(channel="B1")
            record = {
                "source": "ncbi_premrna",
                "accession": "NM_123.1",
                "length": 30,
                "transcript_annotation": annotation,
                "premrna_target": annotation,
            }
            run_metadata = {"outputs": {"idt_tsv": str(idt)}}
            with (
                patch.object(core, "run_design_probes", return_value=(probes, idt, log, params, run_metadata)),
                patch.object(core, "run_oligo_structure_qc", return_value={"status": "completed"}) as qc,
                patch.object(core, "clean_oligo_seq", side_effect=lambda sequence: sequence),
            ):
                result = core.run_design_tier(
                    fasta=root / "target.fa",
                    run_dir=root,
                    target="target",
                    args=args,
                    record=record,
                    sequence_length=30,
                )
                self.assertEqual(len(result["rows"]), 1)
                self.assertEqual(result["qc_input"].name, "target_eligible_IDT_order.csv")
                self.assertEqual(qc.call_args.args[0].name, "target_eligible_IDT_order.csv")
                self.assertEqual(result["eligible_idt_csv"].read_text(encoding="utf-8").count("_P1"), 1)

    def test_prepare_target_explains_inadequate_annotation(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            genome = root / "genome.fna"
            genome.write_text(">chr1\nACGTACGTACGT\n", encoding="utf-8")
            annotation = root / "annotation.gff3"
            annotation.write_text(
                "##gff-version 3\n"
                "chr1\tncbi\tmRNA\t1\t12\t.\t+\t.\tID=rna1;transcript_id=NM_123.1;gene=sox9\n"
                "chr1\tncbi\texon\t1\t12\t.\t+\t.\tParent=rna1\n",
                encoding="utf-8",
            )
            metadata = {
                "status": "ready",
                "species": "xtr",
                "display_name": "Xenopus tropicalis",
                "scientific_name": "Xenopus tropicalis",
                "assembly": "UCB_Xtro_10.0",
                "assembly_accession": "GCF_TEST.1",
                "reference_data_directory": str(root),
            }
            args = SimpleNamespace(species="xtr", index=None, premrna_introns=None)
            record = {"accession": "NM_123.1", "gene_symbol": "sox9"}
            source = {
                "reference_data_directory": str(root),
                "assembly_accession": "GCF_TEST.1",
            }
            with (
                patch("hcrprobeforge.premrna.resolve_reference_metadata", return_value=metadata),
                patch("hcrprobeforge.premrna._download_sources", return_value=(genome, annotation, source)),
                patch("hcrprobeforge.references.get_species_preset", return_value=references.get_species_preset("xtr")),
            ):
                with self.assertRaisesRegex(RuntimeError, "Pre-mRNA design cannot proceed") as raised:
                    prepare_target(None, record, args)

        message = str(raised.exception)
        self.assertIn("Xenopus tropicalis", message)
        self.assertIn("UCB_Xtro_10.0 (GCF_TEST.1)", message)
        self.assertIn("at least two linked exons", message)
        self.assertIn("matching local GFF3", message)
        self.assertIn("does not contain a usable multi-exon model", message)

    def test_missing_annotation_url_explains_annotation_is_unavailable(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            genome = root / "genome.fna"
            genome.write_text(">chr1\nACGT\n", encoding="utf-8")
            metadata = {
                "display_name": "Xenopus tropicalis",
                "assembly": "UCB_Xtro_10.0",
                "assembly_accession": "",
                "reference_data_directory": str(root),
            }
            args = SimpleNamespace(email=None, api_key=None)
            with self.assertRaisesRegex(RuntimeError, "no matching genomic annotation is available") as raised:
                from hcrprobeforge.premrna import _download_sources

                _download_sources(None, metadata, "xtr", args)

        message = str(raised.exception)
        self.assertIn("Xenopus tropicalis", message)
        self.assertIn("UCB_Xtro_10.0", message)
        self.assertIn("matching local GFF3", message)
