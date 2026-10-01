from __future__ import annotations

import json
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

import hcrprobeforge.references as references


class SpeciesDeletionTests(unittest.TestCase):
    def _metadata(self, root: Path, species: str, prefix: Path) -> Path:
        path = root / "refs" / species / "assembly_v1" / "reference.json"
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(
            json.dumps(
                {
                    "status": "ready",
                    "species": species,
                    "index_species_alias": species,
                    "index_prefix": str(prefix),
                    "assembly": "assembly_v1",
                }
            ),
            encoding="utf-8",
        )
        return path

    def _index_files(self, prefix: Path) -> list[Path]:
        files = [prefix.with_name(prefix.name + suffix) for suffix in references.BOWTIE2_INDEX_SUFFIXES]
        prefix.parent.mkdir(parents=True, exist_ok=True)
        for path in files:
            path.write_bytes(b"index")
        return files

    def _config(self, root: Path, alias: str, prefix: Path) -> Path:
        path = root / "design" / "HCRconfig.yaml"
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(
            "references:\n"
            f"  {alias}: {prefix}\n"
            "  built_in: /some/other/index\n",
            encoding="utf-8",
        )
        return path

    def test_custom_species_deletion_removes_private_data_and_exact_index_files(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "custom_genome.fna"
            source.write_text(">contig\nACGT\n", encoding="utf-8")
            prefix = root / "design" / "indices" / "octopus_v1" / "octopus_v1"
            with (
                patch.object(references, "species_data_root", return_value=root / "refs"),
                patch.object(references, "hcrprobedesign_data_root", return_value=root / "design"),
            ):
                references.register_custom_species(
                    "octopus_v1",
                    "Octopus vulgaris",
                    "Octopus vulgaris",
                    assembly_accession="GCF_000000001.1",
                    genome_fasta=source,
                )
                metadata = self._metadata(root, "octopus_v1", prefix)
                index_files = self._index_files(prefix)
                plan = references.custom_species_deletion_plan("octopus_v1")

                self.assertEqual(len(plan["index_files"]), len(references.BOWTIE2_INDEX_SUFFIXES))
                self.assertEqual(plan["preserved_index_files"], [])
                self.assertTrue(plan["saved_fasta"])
                references.delete_custom_species("octopus_v1")

                self.assertIsNone(references.get_species_preset("octopus_v1"))
                self.assertFalse(metadata.exists())
                self.assertTrue(all(not path.exists() for path in index_files))
                self.assertFalse(list((root / "refs" / "custom_species").rglob("*")))

    def test_custom_species_deletion_unregisters_alias_and_removes_empty_index_folder(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            prefix = root / "design" / "indices" / "c_elegans" / "c_elegans"
            config = self._config(root, "c_elegans", prefix)
            with (
                patch.object(references, "species_data_root", return_value=root / "refs"),
                patch.object(references, "hcrprobedesign_data_root", return_value=root / "design"),
            ):
                references.register_custom_species("c_elegans", "C. elegans", "Caenorhabditis elegans")
                self._metadata(root, "c_elegans", prefix)
                self._index_files(prefix)
                plan = references.custom_species_deletion_plan("c_elegans")

                self.assertEqual(plan["config_aliases_to_remove"], ["c_elegans"])
                self.assertEqual(plan["config_files"], [str(config)])
                references.delete_custom_species("c_elegans")

                config_text = config.read_text(encoding="utf-8")
                self.assertNotIn("c_elegans", config_text)
                self.assertIn("built_in", config_text)
                self.assertFalse(prefix.parent.exists())
                self.assertTrue(config.with_name("HCRconfig.yaml.hcrprobeforge.bak").is_file())

    def test_stale_custom_registration_is_repaired_before_rebuild(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            prefix = root / "design" / "indices" / "c_elegans" / "c_elegans"
            self._config(root, "c_elegans", prefix)
            prefix.parent.mkdir(parents=True, exist_ok=True)
            (prefix.parent / "c_elegans.1.bt2.tmp").write_bytes(b"partial")
            with (
                patch.object(references, "species_data_root", return_value=root / "refs"),
                patch.object(references, "hcrprobedesign_data_root", return_value=root / "design"),
            ):
                references.register_custom_species("c_elegans", "C. elegans", "Caenorhabditis elegans")
                self.assertTrue(
                    references._recover_stale_custom_registration(
                        "c_elegans",
                        custom_species_key="c_elegans",
                    )
                )
                self.assertNotIn(
                    "c_elegans",
                    (root / "design" / "HCRconfig.yaml").read_text(encoding="utf-8"),
                )
                self.assertFalse(prefix.parent.exists())

    def test_builtin_preset_cannot_delete_or_touch_its_index(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            prefix = root / "design" / "indices" / "xtr" / "xtr"
            index_files = self._index_files(prefix)
            with (
                patch.object(references, "species_data_root", return_value=root / "refs"),
                patch.object(references, "hcrprobedesign_data_root", return_value=root / "design"),
            ):
                with self.assertRaises(ValueError):
                    references.delete_custom_species("xtr")
            self.assertTrue(all(path.exists() for path in index_files))

    def test_builtin_index_removal_preserves_preset_and_reference_assets(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            prefix = root / "design" / "indices" / "xtr" / "xtr"
            index_files = self._index_files(prefix)
            reference_dir = root / "refs" / "xtr" / "UCB_Xtro_10.0"
            reference_dir.mkdir(parents=True, exist_ok=True)
            metadata = reference_dir / "reference.json"
            metadata.write_text(
                json.dumps(
                    {
                        "status": "ready",
                        "species": "xtr",
                        "index_species_alias": "xtr",
                        "index_prefix": str(prefix),
                        "assembly": "UCB_Xtro_10.0",
                    }
                ),
                encoding="utf-8",
            )
            genome = reference_dir / "genome.fna"
            annotation = reference_dir / "annotation.gff"
            database = reference_dir / "annotation.sqlite"
            for path in (genome, annotation, database):
                path.write_text("retained", encoding="utf-8")
            with (
                patch.object(references, "species_data_root", return_value=root / "refs"),
                patch.object(references, "hcrprobedesign_data_root", return_value=root / "design"),
            ):
                plan = references.builtin_species_index_deletion_plan("xtr")
                self.assertTrue(plan["builtin"])
                self.assertTrue(plan["preset_protected"])
                self.assertEqual(len(plan["index_files"]), len(index_files))
                references.delete_builtin_species_index("xtr")
                self.assertTrue(all(not path.exists() for path in index_files))
                self.assertTrue(metadata.exists())
                self.assertTrue(genome.exists())
                self.assertTrue(annotation.exists())
                self.assertTrue(database.exists())
                self.assertIsNotNone(references.get_species_preset("xtr"))

    def test_shared_index_is_retained_when_one_custom_species_is_removed(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            prefix = root / "design" / "indices" / "shared" / "shared"
            with (
                patch.object(references, "species_data_root", return_value=root / "refs"),
                patch.object(references, "hcrprobedesign_data_root", return_value=root / "design"),
            ):
                references.register_custom_species("species_a", "Species A", "Species A")
                references.register_custom_species("species_b", "Species B", "Species B")
                self._metadata(root, "species_a", prefix)
                self._metadata(root, "species_b", prefix)
                index_files = self._index_files(prefix)

                plan = references.custom_species_deletion_plan("species_a")
                self.assertEqual(plan["index_files"], [])
                self.assertEqual(len(plan["preserved_index_files"]), len(index_files))
                references.delete_custom_species("species_a")

                self.assertTrue(all(path.exists() for path in index_files))
                self.assertIsNone(references.get_species_preset("species_a"))
                self.assertIsNotNone(references.get_species_preset("species_b"))


if __name__ == "__main__":
    unittest.main()
