"""NCBI reference retrieval and HCRProbeDesign index management.

This module owns reference-data operations only.  It does not generate probe
candidates or select probe pairs.  Genome FASTA files are retrieved from the
NCBI Assembly record, then the installed ``buildGenomeIndex`` command creates
and registers the Bowtie2 index in HCRProbeDesign's normal data directory.
"""

from __future__ import annotations

import argparse
import gzip
import json
import os
import re
import shutil
import subprocess
import tempfile
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Callable
from urllib.parse import quote, unquote, urljoin, urlsplit, urlunsplit

from . import core


NCBI_EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
ASSEMBLY_ACCESSION_RE = re.compile(r"^(?:GCF|GCA)_[0-9]+(?:\.[0-9]+)?$", re.IGNORECASE)
ProgressCallback = Callable[[dict[str, object]], None]

# A Bowtie2 index is a small, fixed family of files.  Keep this list explicit
# so custom-species cleanup can remove only the files belonging to the exact
# registered prefix, never an entire HCRProbeDesign data directory.
BOWTIE2_INDEX_SUFFIXES: tuple[str, ...] = (
    ".1.bt2",
    ".2.bt2",
    ".3.bt2",
    ".4.bt2",
    ".rev.1.bt2",
    ".rev.2.bt2",
    ".1.bt2l",
    ".2.bt2l",
    ".3.bt2l",
    ".4.bt2l",
    ".rev.1.bt2l",
    ".rev.2.bt2l",
)
BOWTIE2_TEMP_SUFFIXES: tuple[str, ...] = tuple(
    f"{suffix}.tmp" for suffix in BOWTIE2_INDEX_SUFFIXES
)


class _AssemblyDownloadError(RuntimeError):
    """Download failure that retains the HTTP status for safe fallback logic."""

    def __init__(self, message: str, *, status_code: int | None = None) -> None:
        super().__init__(message)
        self.status_code = status_code


@dataclass(frozen=True)
class SpeciesPreset:
    """Approved organism/assembly pairing used by the automatic index flow."""

    key: str
    display_name: str
    scientific_name: str
    assembly_name: str
    assembly_accession: str
    aliases: tuple[str, ...] = ()
    local_genome_fasta: str | None = None


SPECIES_PRESETS: dict[str, SpeciesPreset] = {
    "xtr": SpeciesPreset(
        key="xtr",
        display_name="Xenopus tropicalis",
        scientific_name="Xenopus tropicalis",
        assembly_name="UCB_Xtro_10.0",
        assembly_accession="GCF_000004195.4",
        aliases=("xenopus_tropicalis", "xenopus-tropicalis"),
    ),
    "xla": SpeciesPreset(
        key="xla",
        display_name="Xenopus laevis",
        scientific_name="Xenopus laevis",
        assembly_name="Xenopus_laevis_v10.1",
        assembly_accession="GCF_017654675.1",
        aliases=("xenopus_laevis", "xenopus-laevis"),
    ),
    "zebrafish": SpeciesPreset(
        key="zebrafish",
        display_name="Zebrafish",
        scientific_name="Danio rerio",
        assembly_name="GRCz12ab",
        assembly_accession="GCF_052040795.1",
        aliases=("danio_rerio", "danio-rerio"),
    ),
    "mouse": SpeciesPreset(
        key="mouse",
        display_name="Mouse",
        scientific_name="Mus musculus",
        assembly_name="GRCm39",
        assembly_accession="GCF_000001635.27",
        aliases=("mus_musculus", "mus-musculus"),
    ),
    "chicken": SpeciesPreset(
        key="chicken",
        display_name="Chicken",
        scientific_name="Gallus gallus",
        assembly_name="bGalGal1.mat.broiler.GRCg7b",
        assembly_accession="GCF_016699485.2",
        aliases=("gallus_gallus", "gallus-gallus"),
    ),
    "human": SpeciesPreset(
        key="human",
        display_name="Human",
        scientific_name="Homo sapiens",
        assembly_name="GRCh38.p14",
        assembly_accession="GCF_000001405.40",
        aliases=("homo_sapiens", "homo-sapiens", "hg38", "grch38"),
    ),
}

# User-created HCRProbeDesign aliases remain accepted through the registered
# reference metadata path. Built-in aliases are intentionally limited to the
# identifiers declared by the current preset table.
COMPATIBLE_ORGANISM_ALIASES: dict[str, str] = {}

# These are directory names, not biological labels or HCRProbeDesign aliases.
# Keep the short internal workflow names accepted by the Python API while
# writing clearer, pluralized folders for users.
WORKFLOW_DIRECTORIES: dict[str, str] = {
    "design": "individual_designs",
    "individual_designs": "individual_designs",
    "manifest": "manifests",
    "manifests": "manifests",
    "plot": "plot",
    "qc": "qc",
    "index": "index",
}


def canonical_species(value: str | None) -> str:
    """Normalize a species alias without changing its biological label."""
    return str(value or "").strip().casefold().replace(" ", "_")


REMOVED_BUILTIN_ALIASES = frozenset({"xtr10"})


def species_aliases(value: str | None) -> tuple[str, ...]:
    """Return all equivalent identifiers for a preset or registered alias."""
    requested = canonical_species(value)
    if not requested:
        return ()
    preset = get_species_preset(requested)
    if preset is None:
        return (requested,)
    values = {
        canonical_species(item)
        for item in (preset.key, preset.display_name, preset.scientific_name, *preset.aliases)
        if item
    }
    values.add(requested)
    return tuple(sorted(values))


def normalize_assembly_accession(value: str | None) -> str:
    """Return a canonical, validated NCBI assembly accession."""
    accession = str(value or "").strip().upper()
    if not ASSEMBLY_ACCESSION_RE.fullmatch(accession):
        raise ValueError(
            "assembly_accession must look like GCF_<ASSEMBLY_ID>.<VERSION> or GCA_<ASSEMBLY_ID>.<VERSION>"
        )
    return accession


def _normalise_assembly_component(value: str | None) -> str:
    """Map an assembly label to NCBI's filename-safe underscore convention."""
    text = unquote(str(value or "")).strip()
    text = re.sub(r"\s+", "_", text)
    text = re.sub(r"[^A-Za-z0-9._-]+", "_", text)
    text = re.sub(r"_+", "_", text)
    return text.strip("._-")


def custom_species_registry_path() -> Path:
    """Return the user-level registry for webapp/CLI-created species."""
    return species_data_root() / "species_presets.json"


def _custom_species_records() -> list[dict[str, Any]]:
    path = custom_species_registry_path()
    if not path.is_file():
        return []
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError, TypeError):
        return []
    records = payload.get("species", []) if isinstance(payload, dict) else payload
    if not isinstance(records, list):
        return []
    return [record for record in records if isinstance(record, dict)]


def _custom_preset(record: dict[str, Any]) -> SpeciesPreset | None:
    key = canonical_species(str(record.get("key") or ""))
    display_name = str(record.get("display_name") or "").strip()
    scientific_name = str(record.get("scientific_name") or "").strip()
    if not key or not display_name or not scientific_name:
        return None
    aliases = tuple(
        canonical_species(str(alias))
        for alias in (record.get("aliases") or [])
        if str(alias).strip()
    )
    return SpeciesPreset(
        key=key,
        display_name=display_name,
        scientific_name=scientific_name,
        assembly_name=str(record.get("assembly_name") or "custom local assembly").strip(),
        assembly_accession=str(record.get("assembly_accession") or "").strip(),
        aliases=aliases,
        local_genome_fasta=str(record.get("local_genome_fasta") or "").strip() or None,
    )


def custom_species_presets() -> tuple[SpeciesPreset, ...]:
    """Return valid user-created presets in their saved order."""
    presets: list[SpeciesPreset] = []
    builtin_keys = {preset.key for preset in SPECIES_PRESETS.values()}
    for record in _custom_species_records():
        preset = _custom_preset(record)
        if (
            preset is not None
            and preset.key not in builtin_keys
            and preset.key not in REMOVED_BUILTIN_ALIASES
            and preset.key not in {item.key for item in presets}
        ):
            presets.append(preset)
    return tuple(presets)


def register_custom_species(
    key: str,
    display_name: str,
    scientific_name: str,
    *,
    assembly_accession: str | None = None,
    assembly_name: str | None = None,
    genome_fasta: Path | None = None,
) -> SpeciesPreset:
    """Create one persistent user species preset.

    The registry stores descriptive metadata and, when supplied, a private
    copy of the user's genome FASTA. Transcript lookup remains NCBI-backed
    through ``scientific_name``; the FASTA is used only by the index workflow.
    """
    normalized_key = canonical_species(key)
    if not normalized_key or normalized_key in {".", ".."} or core.safe_name(normalized_key) != normalized_key:
        raise ValueError("Species alias may contain only letters, numbers, dots, underscores, and hyphens")
    if normalized_key == "__add_new_species__":
        raise ValueError("Choose a different species alias")
    display = str(display_name or "").strip()
    organism = str(scientific_name or "").strip()
    if not display:
        raise ValueError("Display name is required")
    if not organism:
        raise ValueError("NCBI organism is required")
    if any(any(character in value for character in "\r\n\t") for value in (display, organism)):
        raise ValueError("Display name and NCBI organism cannot contain line breaks or tabs")
    if get_species_preset(normalized_key) is not None:
        raise ValueError(f"Species alias {normalized_key!r} already exists")
    accession = ""
    if assembly_accession:
        accession = normalize_assembly_accession(assembly_accession)
    name = str(assembly_name or "").strip() or "custom local assembly"
    if any(character in name for character in "\r\n\t"):
        raise ValueError("Assembly name cannot contain line breaks or tabs")
    source = Path(genome_fasta).expanduser().resolve() if genome_fasta is not None else None
    if source is not None and not source.is_file():
        raise ValueError(f"Genome FASTA not found: {source}")
    if not accession and source is None:
        # A preset can be created before its reference is available; the
        # index tab/CLI can receive the assembly or FASTA later.
        name = name or "custom assembly"

    records = _custom_species_records()
    species_root = species_data_root() / "custom_species" / normalized_key
    copied_fasta: Path | None = None
    temporary_fasta: Path | None = None
    try:
        if source is not None:
            species_root.mkdir(parents=True, exist_ok=True)
            copied_fasta = species_root / (core.safe_name(source.name) or "genome.fasta")
            temporary_fasta = copied_fasta.with_name(copied_fasta.name + f".tmp.{os.getpid()}")
            shutil.copyfile(source, temporary_fasta)
            os.replace(temporary_fasta, copied_fasta)
            temporary_fasta = None
        records.append(
            {
                "key": normalized_key,
                "display_name": display,
                "scientific_name": organism,
                "assembly_name": name,
                "assembly_accession": accession,
                "aliases": [],
                "local_genome_fasta": str(copied_fasta) if copied_fasta else None,
                "created_at_utc": core.utc_now_iso(),
            }
        )
        registry = custom_species_registry_path()
        registry.parent.mkdir(parents=True, exist_ok=True)
        temporary = registry.with_name(registry.name + f".tmp.{os.getpid()}")
        temporary.write_text(json.dumps({"schema_version": 1, "species": records}, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        os.replace(temporary, registry)
    except Exception:
        if temporary_fasta is not None:
            temporary_fasta.unlink(missing_ok=True)
        if copied_fasta is not None:
            copied_fasta.unlink(missing_ok=True)
        raise
    return _custom_preset(records[-1]) or SpeciesPreset(
        normalized_key, display, organism, name, accession, (), str(copied_fasta) if copied_fasta else None
    )


def _normalise_index_prefix(raw: str | Path | None, *, data_root: Path) -> Path | None:
    """Return a resolved Bowtie2 prefix, normalising a suffix if necessary."""
    text = str(raw or "").strip()
    if not text:
        return None
    candidate = Path(text).expanduser()
    if not candidate.is_absolute():
        candidate = data_root / candidate
    candidate = candidate.resolve()
    for suffix in sorted(BOWTIE2_INDEX_SUFFIXES, key=len, reverse=True):
        if candidate.name.endswith(suffix):
            candidate = candidate.with_name(candidate.name[: -len(suffix)])
            break
    return candidate


def _index_files_for_prefix(prefix: Path) -> list[Path]:
    """Return existing Bowtie2 files for one exact prefix."""
    return [
        prefix.with_name(prefix.name + suffix)
        for suffix in BOWTIE2_INDEX_SUFFIXES
        if prefix.with_name(prefix.name + suffix).is_file()
    ]


def _index_cleanup_files_for_prefix(prefix: Path) -> list[Path]:
    """Return final and interrupted-build files for one exact prefix."""
    return [
        prefix.with_name(prefix.name + suffix)
        for suffix in BOWTIE2_INDEX_SUFFIXES + BOWTIE2_TEMP_SUFFIXES
        if prefix.with_name(prefix.name + suffix).is_file()
    ]


def _build_rollback_backup(
    alias: str,
    metadata_path: Path,
    reference_directory: Path,
) -> tuple[Path, dict[Path, Path]]:
    """Back up managed files that an index build can replace or remove.

    Filesystem-entry snapshots are sufficient for a new build, but they cannot
    restore the contents of an existing genome, metadata, SQLite, Bowtie2, or
    configuration file after a forced rebuild fails halfway through. Keep a
    private byte-for-byte backup of those bounded managed targets.
    """
    backup_root = Path(tempfile.mkdtemp(prefix="hcrprobeforge-build-rollback-"))
    paths: list[Path] = [
        metadata_path,
        reference_directory / "genome.fna",
        reference_directory / "genome.fna.gz",
        reference_directory / "annotation.gff",
        reference_directory / "annotation.gff.gz",
        reference_directory / "annotation.sqlite",
        reference_directory / "premrna" / "genome.fna",
        reference_directory / "premrna" / "genome.fna.gz",
        reference_directory / "premrna" / "annotation.gff",
        reference_directory / "premrna" / "annotation.gff.gz",
        reference_directory / "premrna" / "annotation.sqlite",
        *_existing_hcrprobedesign_configs(),
    ]
    try:
        paths.extend(_index_cleanup_files_for_prefix(registered_index_prefix(alias)))
    except (OSError, ValueError):
        pass
    backups: dict[Path, Path] = {}
    for index, path in enumerate(dict.fromkeys(item.resolve() for item in paths)):
        if not path.is_file():
            continue
        backup = backup_root / f"{index:04d}_{core.safe_name(path.name) or 'file'}"
        shutil.copy2(path, backup)
        backups[path] = backup
    return backup_root, backups


def _restore_build_rollback_backup(
    backup_root: Path,
    backups: dict[Path, Path],
) -> None:
    """Restore backed-up managed files and remove the private backup tree."""
    try:
        for original, backup in backups.items():
            original.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(backup, original)
    finally:
        shutil.rmtree(backup_root, ignore_errors=True)


def _discard_build_rollback_backup(backup_root: Path) -> None:
    shutil.rmtree(backup_root, ignore_errors=True)


def _index_family_is_complete(prefix: Path) -> bool:
    """Return whether either the small or large Bowtie2 family is complete."""
    present = {path.name for path in _index_files_for_prefix(prefix)}
    small = {prefix.name + suffix for suffix in BOWTIE2_INDEX_SUFFIXES[:6]}
    large = {prefix.name + suffix for suffix in BOWTIE2_INDEX_SUFFIXES[6:]}
    return small.issubset(present) or large.issubset(present)


def _hcrprobedesign_config_candidates() -> tuple[Path, ...]:
    """Return known HCRProbeDesign config paths in preference order."""
    data_root = hcrprobedesign_data_root()
    return (
        data_root / "HCRconfig.yaml",
        data_root / "hcrconfig.yaml",
        data_root / "config.yaml",
    )


def _existing_hcrprobedesign_configs() -> list[Path]:
    paths: list[Path] = []
    for path in _hcrprobedesign_config_candidates():
        if path.is_file() and path not in paths:
            paths.append(path)
    return paths


def _yaml_module() -> Any:
    try:
        import yaml
    except ImportError as exc:  # pragma: no cover - PyYAML is a declared dependency
        raise RuntimeError(
            "HCRProbeDesign config cleanup requires PyYAML, which is a package dependency."
        ) from exc
    return yaml


def _load_hcrprobedesign_config(path: Path) -> dict[str, Any]:
    yaml = _yaml_module()
    try:
        payload = yaml.safe_load(path.read_text(encoding="utf-8"))
    except (OSError, ValueError, yaml.YAMLError) as exc:
        raise RuntimeError(f"Could not read HCRProbeDesign config {path}: {exc}") from exc
    if payload is None:
        return {}
    if not isinstance(payload, dict):
        raise RuntimeError(f"HCRProbeDesign config {path} must contain a YAML mapping")
    return payload


def _config_strings(value: Any) -> list[str]:
    if isinstance(value, str):
        return [value]
    if isinstance(value, dict):
        return [item for child in value.values() for item in _config_strings(child)]
    if isinstance(value, (list, tuple)):
        return [item for child in value for item in _config_strings(child)]
    return []


def _config_entry_matches(
    value: Any,
    *,
    expected_prefixes: set[Path],
    data_root: Path,
) -> bool:
    if not expected_prefixes:
        return True
    for raw in _config_strings(value):
        candidate = _normalise_index_prefix(raw, data_root=data_root)
        if candidate in expected_prefixes:
            return True
    return False


def _config_contains_alias(
    value: Any,
    alias: str,
    *,
    expected_prefixes: set[Path],
    data_root: Path,
) -> bool:
    normalized_alias = canonical_species(alias)
    if isinstance(value, dict):
        for key, child in value.items():
            if canonical_species(str(key)) == normalized_alias and _config_entry_matches(
                child, expected_prefixes=expected_prefixes, data_root=data_root
            ):
                return True
            if canonical_species(str(key)) in {"species", "alias", "name", "key"}:
                if canonical_species(str(child)) == normalized_alias and _config_entry_matches(
                    value, expected_prefixes=expected_prefixes, data_root=data_root
                ):
                    return True
            if _config_contains_alias(
                child,
                normalized_alias,
                expected_prefixes=expected_prefixes,
                data_root=data_root,
            ):
                return True
    elif isinstance(value, (list, tuple)):
        for child in value:
            if _config_contains_alias(
                child,
                normalized_alias,
                expected_prefixes=expected_prefixes,
                data_root=data_root,
            ):
                return True
    return False


def _remove_alias_from_config(
    value: Any,
    alias: str,
    *,
    expected_prefixes: set[Path],
    data_root: Path,
) -> int:
    """Remove matching mapping/list entries and return the number removed."""
    normalized_alias = canonical_species(alias)
    removed = 0
    if isinstance(value, dict):
        for key in list(value):
            child = value[key]
            if canonical_species(str(key)) == normalized_alias and _config_entry_matches(
                child, expected_prefixes=expected_prefixes, data_root=data_root
            ):
                del value[key]
                removed += 1
                continue
            if isinstance(child, dict):
                named_alias = next(
                    (
                        child_name
                        for child_name in ("species", "alias", "name", "key")
                        if canonical_species(str(child.get(child_name) or "")) == normalized_alias
                    ),
                    None,
                )
                if named_alias and _config_entry_matches(
                    child, expected_prefixes=expected_prefixes, data_root=data_root
                ):
                    del value[key]
                    removed += 1
                    continue
            removed += _remove_alias_from_config(
                child,
                normalized_alias,
                expected_prefixes=expected_prefixes,
                data_root=data_root,
            )
    elif isinstance(value, list):
        kept: list[Any] = []
        for child in value:
            if isinstance(child, dict):
                named_alias = next(
                    (
                        child_name
                        for child_name in ("species", "alias", "name", "key")
                        if canonical_species(str(child.get(child_name) or "")) == normalized_alias
                    ),
                    None,
                )
                if named_alias and _config_entry_matches(
                    child, expected_prefixes=expected_prefixes, data_root=data_root
                ):
                    removed += 1
                    continue
            removed += _remove_alias_from_config(
                child,
                normalized_alias,
                expected_prefixes=expected_prefixes,
                data_root=data_root,
            )
            kept.append(child)
        value[:] = kept
    return removed


def _config_paths_containing_alias(
    alias: str,
    *,
    expected_prefixes: set[Path],
) -> list[Path]:
    data_root = hcrprobedesign_data_root().resolve()
    paths: list[Path] = []
    for path in _existing_hcrprobedesign_configs():
        payload = _load_hcrprobedesign_config(path)
        if _config_contains_alias(
            payload,
            alias,
            expected_prefixes=expected_prefixes,
            data_root=data_root,
        ):
            paths.append(path)
    return paths


def _unregister_index_alias(
    alias: str,
    *,
    expected_prefixes: list[Path] | tuple[Path, ...] = (),
    allow_builtin: bool = False,
) -> dict[str, Any]:
    """Remove one verified alias from HCRProbeDesign config atomically.

    This deliberately does not call ``buildGenomeIndex --force``. It removes
    only a verified alias, preserving unrelated registrations. Built-in aliases
    are accepted only by the explicit index-removal workflow; that workflow
    never removes the built-in preset or its reference assets.
    """
    normalized_alias = canonical_species(alias)
    if not normalized_alias:
        raise ValueError("Index alias cannot be empty")
    preset = get_species_preset(normalized_alias)
    if preset is not None and preset.key in SPECIES_PRESETS and not allow_builtin:
        raise ValueError("Built-in HCRProbeDesign aliases cannot be unregistered")
    data_root = hcrprobedesign_data_root().resolve()
    normalized_prefixes = {
        prefix.resolve()
        for prefix in expected_prefixes
        if _within(prefix, data_root)
    }
    changed_configs: list[str] = []
    backup_paths: list[str] = []
    removed_entries = 0
    yaml = _yaml_module()
    for path in _existing_hcrprobedesign_configs():
        payload = _load_hcrprobedesign_config(path)
        if not _config_contains_alias(
            payload,
            normalized_alias,
            expected_prefixes=normalized_prefixes,
            data_root=data_root,
        ):
            continue
        removed = _remove_alias_from_config(
            payload,
            normalized_alias,
            expected_prefixes=normalized_prefixes,
            data_root=data_root,
        )
        if not removed:
            continue
        backup = path.with_name(path.name + ".hcrprobeforge.bak")
        if not backup.exists():
            shutil.copy2(path, backup)
        temporary = path.with_name(path.name + f".tmp.{os.getpid()}")
        try:
            temporary.write_text(
                yaml.safe_dump(payload, sort_keys=False, default_flow_style=False),
                encoding="utf-8",
            )
            os.replace(temporary, path)
        except Exception:
            temporary.unlink(missing_ok=True)
            raise RuntimeError(f"Could not update HCRProbeDesign config {path}") from None
        changed_configs.append(str(path))
        backup_paths.append(str(backup))
        removed_entries += removed
    return {
        "alias": normalized_alias,
        "config_files": changed_configs,
        "backup_files": backup_paths,
        "removed_entries": removed_entries,
    }


def unregister_custom_index_alias(
    alias: str,
    *,
    expected_prefixes: list[Path] | tuple[Path, ...] = (),
) -> dict[str, Any]:
    """Remove one custom alias while protecting all built-in aliases."""
    return _unregister_index_alias(alias, expected_prefixes=expected_prefixes)


def _within(path: Path, root: Path) -> bool:
    """Return whether ``path`` is contained by ``root`` after resolution."""
    try:
        path.resolve().relative_to(root.resolve())
    except ValueError:
        return False
    return True


def _reference_metadata_paths() -> list[Path]:
    root = species_data_root()
    if not root.exists():
        return []
    return sorted(root.glob("*/*/reference.json"))


def _prefix_referenced_by_other_species(prefix: Path, species_key: str) -> bool:
    """Return whether another species metadata record uses this prefix."""
    data_root = hcrprobedesign_data_root().resolve()
    excluded_root = (species_data_root() / species_key).resolve()
    for metadata_path in _reference_metadata_paths():
        if _within(metadata_path, excluded_root):
            continue
        metadata = _read_reference_metadata(metadata_path)
        if metadata is None:
            continue
        fallback_alias = str(
            metadata.get("species") or metadata.get("index_species_alias") or "reference"
        )
        other_prefix = _metadata_index_prefix(
            metadata,
            fallback_alias=fallback_alias,
            data_root=data_root,
        )
        if other_prefix == prefix:
            return True
    return False


def _remove_index_files_and_empty_directory(prefix: Path) -> dict[str, Any]:
    """Remove exact index/temp files and then only their empty alias folder."""
    cleanup_files = _index_cleanup_files_for_prefix(prefix)
    errors: list[str] = []
    for path in cleanup_files:
        try:
            path.unlink(missing_ok=True)
        except OSError as exc:
            errors.append(f"{path}: {exc}")
    removed_directory = ""
    directory = prefix.parent
    if not errors and directory.is_dir():
        try:
            if not any(directory.iterdir()):
                directory.rmdir()
                removed_directory = str(directory)
        except OSError as exc:
            errors.append(f"{directory}: {exc}")
    if errors:
        raise OSError("; ".join(errors))
    return {
        "files": [str(path) for path in cleanup_files],
        "directory": removed_directory,
    }


def _recover_stale_custom_registration(
    alias: str,
    *,
    custom_species_key: str,
) -> bool:
    """Repair an old custom alias left in config without a complete index."""
    normalized_alias = canonical_species(alias)
    preset = get_species_preset(custom_species_key)
    if preset is None or preset.key in SPECIES_PRESETS:
        return False
    data_root = hcrprobedesign_data_root().resolve()
    prefix = _normalise_index_prefix(
        registered_index_prefix(normalized_alias), data_root=data_root
    )
    if prefix is None or not _within(prefix, data_root):
        return False
    if _index_family_is_complete(prefix):
        return False
    if _prefix_referenced_by_other_species(prefix, preset.key):
        return False
    expected = {prefix}
    config_paths = _config_paths_containing_alias(
        normalized_alias,
        expected_prefixes=expected,
    )
    if not config_paths:
        return False
    unregister_custom_index_alias(
        normalized_alias,
        expected_prefixes=(prefix,),
    )
    _remove_index_files_and_empty_directory(prefix)
    return True


def _metadata_index_prefix(
    metadata: dict[str, Any],
    *,
    fallback_alias: str,
    data_root: Path,
) -> Path | None:
    raw_prefix = metadata.get("index_prefix")
    if raw_prefix:
        return _normalise_index_prefix(raw_prefix, data_root=data_root)
    alias = str(metadata.get("index_species_alias") or metadata.get("species") or fallback_alias)
    return _normalise_index_prefix(registered_index_prefix(alias), data_root=data_root)


def _read_reference_metadata(path: Path) -> dict[str, Any] | None:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError, TypeError):
        return None
    return payload if isinstance(payload, dict) else None


def _custom_species_record(key: str) -> tuple[str, dict[str, Any]]:
    normalized_key = canonical_species(key)
    if not normalized_key or normalized_key in SPECIES_PRESETS:
        raise ValueError("Only user-created species presets can be removed")
    for record in _custom_species_records():
        if canonical_species(str(record.get("key") or "")) == normalized_key:
            return normalized_key, record
    raise ValueError(f"Custom species {normalized_key!r} was not found")


def custom_species_deletion_plan(key: str) -> dict[str, Any]:
    """Describe the private data and exact index files removed for a preset.

    The plan is read-only and is used by both the confirmation dialog and the
    deletion operation.  Index files are eligible only when their prefix is
    inside HCRProbeDesign's data directory and no other installed reference
    metadata points to that same prefix.
    """
    normalized_key, removed = _custom_species_record(key)
    preset = _custom_preset(removed)
    display_name = preset.display_name if preset is not None else normalized_key
    reference_root = species_data_root()
    data_root = hcrprobedesign_data_root().resolve()
    metadata_root = reference_root / normalized_key
    metadata_paths = sorted(metadata_root.glob("*/reference.json")) if metadata_root.is_dir() else []
    target_metadata = set(metadata_paths)

    prefixes: list[Path] = []
    prefix_aliases: dict[Path, set[str]] = {}
    unsafe_prefixes: list[str] = []
    unsafe_aliases: set[str] = set()
    for metadata_path in metadata_paths:
        metadata = _read_reference_metadata(metadata_path)
        if metadata is None:
            continue
        alias = core.safe_name(
            str(metadata.get("index_species_alias") or metadata.get("species") or normalized_key)
        ) or normalized_key
        prefix = _metadata_index_prefix(metadata, fallback_alias=normalized_key, data_root=data_root)
        if prefix is None:
            continue
        if not _within(prefix, data_root):
            unsafe_prefixes.append(str(prefix))
            unsafe_aliases.add(alias)
            continue
        prefix_aliases.setdefault(prefix, set()).add(alias)
        if prefix not in prefixes:
            prefixes.append(prefix)

    # A custom species may have a registered index but no package metadata,
    # for example after an interrupted or manually completed build.
    fallback_prefix = _normalise_index_prefix(
        registered_index_prefix(normalized_key), data_root=data_root
    )
    if fallback_prefix is not None:
        if _within(fallback_prefix, data_root):
            prefix_aliases.setdefault(fallback_prefix, set()).add(normalized_key)
            if fallback_prefix not in prefixes:
                prefixes.append(fallback_prefix)
        else:
            unsafe_prefixes.append(str(fallback_prefix))
            unsafe_aliases.add(normalized_key)

    other_prefix_users: dict[Path, list[str]] = {}
    for metadata_path in _reference_metadata_paths():
        if metadata_path in target_metadata:
            continue
        metadata = _read_reference_metadata(metadata_path)
        if metadata is None:
            continue
        fallback_alias = str(metadata.get("species") or metadata.get("index_species_alias") or "reference")
        prefix = _metadata_index_prefix(metadata, fallback_alias=fallback_alias, data_root=data_root)
        if prefix is not None:
            other_prefix_users.setdefault(prefix, []).append(str(metadata_path))

    protected_aliases: set[str] = set(unsafe_aliases)
    for prefix, users in other_prefix_users.items():
        if users:
            protected_aliases.update(prefix_aliases.get(prefix, set()))
    index_files: list[Path] = []
    preserved_index_files: list[Path] = []
    shared_index_reasons: list[str] = []
    removable_prefixes: list[Path] = []
    for prefix in prefixes:
        files = _index_cleanup_files_for_prefix(prefix)
        users = other_prefix_users.get(prefix, [])
        aliases = prefix_aliases.get(prefix, set())
        if users or aliases.intersection(protected_aliases):
            preserved_index_files.extend(files)
            if users:
                shared_index_reasons.append(
                    f"{prefix} is also referenced by {', '.join(users)}"
                )
        else:
            index_files.extend(files)
            removable_prefixes.append(prefix)

    removable_aliases = {
        alias
        for prefix in removable_prefixes
        for alias in prefix_aliases.get(prefix, set())
        if alias not in protected_aliases
    }
    config_files: list[str] = []
    for alias in sorted(removable_aliases):
        expected = {
            prefix
            for prefix in removable_prefixes
            if alias in prefix_aliases.get(prefix, set())
        }
        config_files.extend(
            str(path)
            for path in _config_paths_containing_alias(alias, expected_prefixes=expected)
        )
    cleanup_names = {
        prefix.name + suffix
        for prefix in removable_prefixes
        for suffix in BOWTIE2_INDEX_SUFFIXES + BOWTIE2_TEMP_SUFFIXES
    }
    index_directories = sorted(
        {
            str(prefix.parent)
            for prefix in removable_prefixes
            if prefix.parent.is_dir()
            and all(child.name in cleanup_names for child in prefix.parent.iterdir())
        }
    )

    saved = Path(str(removed.get("local_genome_fasta") or "")).expanduser()
    expected_root = (reference_root / "custom_species" / normalized_key).resolve()
    saved_fasta = str(saved.resolve()) if saved.is_file() and _within(saved, expected_root) else ""
    return {
        "species": normalized_key,
        "display_name": display_name,
        "registry_file": str(custom_species_registry_path()),
        "saved_fasta": saved_fasta,
        "reference_metadata": [str(path) for path in metadata_paths],
        "index_root": str(data_root),
        "index_prefixes": [str(prefix) for prefix in prefixes],
        "index_aliases": sorted(
            {alias for prefix in prefixes for alias in prefix_aliases.get(prefix, set())}
        ),
        "config_aliases_to_remove": sorted(removable_aliases),
        "config_alias_prefixes": {
            alias: [
                str(prefix)
                for prefix in removable_prefixes
                if alias in prefix_aliases.get(prefix, set())
            ]
            for alias in sorted(removable_aliases)
        },
        "config_files": sorted(set(config_files)),
        "index_files": [str(path) for path in index_files],
        "preserved_index_files": [str(path) for path in preserved_index_files],
        "removable_index_prefixes": [str(prefix) for prefix in removable_prefixes],
        "index_directories": index_directories,
        "shared_index_reasons": shared_index_reasons,
        "unsafe_index_prefixes": sorted(set(unsafe_prefixes)),
    }


def builtin_species_index_deletion_plan(key: str) -> dict[str, Any]:
    """Describe removal of only a built-in species' installed index files.

    Built-in presets and their downloaded genome/annotation assets are retained.
    Only unshared Bowtie2 index families, their empty alias folders, and the
    matching HCRProbeDesign registrations are eligible for removal.
    """
    preset = get_species_preset(key)
    if preset is None or preset.key not in SPECIES_PRESETS:
        raise ValueError("Only built-in species indexes can be removed this way")
    normalized_key = preset.key
    aliases = set(species_aliases(normalized_key))
    aliases.add(normalized_key)
    data_root = hcrprobedesign_data_root().resolve()
    target_metadata: set[Path] = set()
    target_rows: list[tuple[Path, dict[str, Any]]] = []
    for metadata_path in _reference_metadata_paths():
        metadata = _read_reference_metadata(metadata_path)
        if metadata is None:
            continue
        row_aliases = {
            canonical_species(metadata.get(name))
            for name in ("index_species_alias", "species", "display_name", "scientific_name")
            if metadata.get(name)
        }
        if aliases.intersection(row_aliases):
            target_metadata.add(metadata_path)
            target_rows.append((metadata_path, metadata))

    prefixes: list[Path] = []
    prefix_aliases: dict[Path, set[str]] = {}
    unsafe_prefixes: list[str] = []
    for metadata_path, metadata in target_rows:
        alias = core.safe_name(
            str(metadata.get("index_species_alias") or metadata.get("species") or normalized_key)
        ) or normalized_key
        prefix = _metadata_index_prefix(metadata, fallback_alias=alias, data_root=data_root)
        if prefix is None:
            continue
        if not _within(prefix, data_root):
            unsafe_prefixes.append(str(prefix))
            continue
        prefix_aliases.setdefault(prefix, set()).add(alias)
        if prefix not in prefixes:
            prefixes.append(prefix)

    # Include a registered index that predates HCRProbeForge metadata while
    # respecting the currently supported alias set.
    for alias in sorted(aliases):
        prefix = _normalise_index_prefix(registered_index_prefix(alias), data_root=data_root)
        if prefix is None:
            continue
        if not _within(prefix, data_root):
            unsafe_prefixes.append(str(prefix))
            continue
        prefix_aliases.setdefault(prefix, set()).add(alias)
        if prefix not in prefixes:
            prefixes.append(prefix)

    other_prefix_users: dict[Path, list[str]] = {}
    for metadata_path in _reference_metadata_paths():
        if metadata_path in target_metadata:
            continue
        metadata = _read_reference_metadata(metadata_path)
        if metadata is None:
            continue
        fallback_alias = str(metadata.get("species") or metadata.get("index_species_alias") or "reference")
        prefix = _metadata_index_prefix(metadata, fallback_alias=fallback_alias, data_root=data_root)
        if prefix is not None:
            other_prefix_users.setdefault(prefix, []).append(str(metadata_path))

    index_files: list[Path] = []
    preserved_index_files: list[Path] = []
    shared_index_reasons: list[str] = []
    removable_prefixes: list[Path] = []
    for prefix in prefixes:
        files = _index_cleanup_files_for_prefix(prefix)
        users = other_prefix_users.get(prefix, [])
        if users:
            preserved_index_files.extend(files)
            shared_index_reasons.append(f"{prefix} is also referenced by {', '.join(users)}")
        else:
            index_files.extend(files)
            removable_prefixes.append(prefix)

    removable_aliases = sorted(
        {
            alias
            for prefix in removable_prefixes
            for alias in prefix_aliases.get(prefix, set())
        }
    )
    config_files: list[str] = []
    config_alias_prefixes: dict[str, list[str]] = {}
    for alias in removable_aliases:
        expected = {
            prefix for prefix in removable_prefixes if alias in prefix_aliases.get(prefix, set())
        }
        config_alias_prefixes[alias] = [str(prefix) for prefix in sorted(expected)]
        config_files.extend(
            str(path) for path in _config_paths_containing_alias(alias, expected_prefixes=expected)
        )
    cleanup_names = {
        prefix.name + suffix
        for prefix in removable_prefixes
        for suffix in BOWTIE2_INDEX_SUFFIXES + BOWTIE2_TEMP_SUFFIXES
    }
    index_directories = sorted(
        {
            str(prefix.parent)
            for prefix in removable_prefixes
            if prefix.parent.is_dir()
            and all(child.name in cleanup_names for child in prefix.parent.iterdir())
        }
    )
    preserved_assets = sorted({str(path.parent) for path in target_metadata})
    return {
        "species": normalized_key,
        "display_name": preset.display_name,
        "builtin": True,
        "preset_protected": True,
        "reference_assets_preserved": preserved_assets,
        "reference_metadata": [str(path) for path in sorted(target_metadata)],
        "index_root": str(data_root),
        "index_prefixes": [str(prefix) for prefix in prefixes],
        "index_aliases": sorted({alias for prefix in prefixes for alias in prefix_aliases.get(prefix, set())}),
        "config_aliases_to_remove": removable_aliases,
        "config_alias_prefixes": config_alias_prefixes,
        "config_files": sorted(set(config_files)),
        "index_files": [str(path) for path in index_files],
        "preserved_index_files": [str(path) for path in preserved_index_files],
        "removable_index_prefixes": [str(prefix) for prefix in removable_prefixes],
        "index_directories": index_directories,
        "shared_index_reasons": shared_index_reasons,
        "unsafe_index_prefixes": sorted(set(unsafe_prefixes)),
    }


def species_deletion_plan(key: str) -> dict[str, Any]:
    """Return the appropriate protected cleanup plan for built-in or custom species."""
    preset = get_species_preset(key)
    if preset is not None and preset.key in SPECIES_PRESETS:
        return builtin_species_index_deletion_plan(key)
    return custom_species_deletion_plan(key)


def delete_builtin_species_index(key: str) -> None:
    """Remove a built-in species' unshared index files, preserving its preset/data."""
    plan = builtin_species_index_deletion_plan(key)
    for alias in plan["config_aliases_to_remove"]:
        _unregister_index_alias(
            alias,
            expected_prefixes=tuple(Path(prefix) for prefix in plan["config_alias_prefixes"].get(alias, [])),
            allow_builtin=True,
        )
    errors: list[str] = []
    for raw_prefix in plan["removable_index_prefixes"]:
        try:
            _remove_index_files_and_empty_directory(Path(raw_prefix))
        except OSError as exc:
            errors.append(str(exc))
    if errors:
        raise OSError(
            f"Built-in species {plan['display_name']!r} index cleanup failed: " + "; ".join(errors)
        )


def delete_species(key: str) -> None:
    """Apply the protected deletion operation selected by :func:`species_deletion_plan`."""
    preset = get_species_preset(key)
    if preset is not None and preset.key in SPECIES_PRESETS:
        delete_builtin_species_index(key)
    else:
        delete_custom_species(key)


def delete_custom_species(key: str) -> None:
    """Remove a custom preset and its exact private/index data.

    Built-in presets cannot reach this path.  For a custom preset, the saved
    FASTA, HCRProbeForge reference metadata, and unshared Bowtie2 files for
    its registered prefixes are removed.  The matching custom aliases are
    also unregistered from HCRProbeDesign, while built-in aliases and any
    index shared with another installed reference are retained.
    """
    normalized_key, removed = _custom_species_record(key)
    plan = custom_species_deletion_plan(normalized_key)

    # Remove stale HCRProbeDesign registrations before deleting files.  This
    # is what allows the same custom alias to be built again without --force.
    for alias in plan["config_aliases_to_remove"]:
        unregister_custom_index_alias(
            alias,
            expected_prefixes=tuple(
                Path(prefix)
                for prefix in plan["config_alias_prefixes"].get(alias, [])
            ),
        )

    index_cleanup_errors: list[str] = []
    for raw_prefix in plan["removable_index_prefixes"]:
        try:
            _remove_index_files_and_empty_directory(Path(raw_prefix))
        except OSError as exc:
            index_cleanup_errors.append(str(exc))
    if index_cleanup_errors:
        raise OSError(
            f"Custom species {normalized_key!r} was not fully removed because index cleanup failed: "
            + "; ".join(index_cleanup_errors)
        )

    records = _custom_species_records()
    remaining = [
        record
        for record in records
        if canonical_species(str(record.get("key") or "")) != normalized_key
    ]
    registry = custom_species_registry_path()
    registry.parent.mkdir(parents=True, exist_ok=True)
    temporary = registry.with_name(registry.name + f".tmp.{os.getpid()}")
    temporary.write_text(
        json.dumps({"schema_version": 1, "species": remaining}, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    os.replace(temporary, registry)

    saved = Path(str(removed.get("local_genome_fasta") or "")).expanduser()
    expected_root = (species_data_root() / "custom_species" / normalized_key).resolve()
    try:
        if saved.is_file() and _within(saved, expected_root):
            saved.unlink()
        if expected_root.is_dir() and not any(expected_root.iterdir()):
            expected_root.rmdir()
    except OSError:
        # The registry is authoritative; stale private data is harmless.
        pass

    metadata_root = species_data_root() / normalized_key
    if metadata_root.is_dir():
        for metadata_file in metadata_root.glob("*/reference.json"):
            try:
                metadata_file.unlink()
            except OSError:
                pass
        for assembly_dir in sorted(
            (item for item in metadata_root.iterdir() if item.is_dir()), reverse=True
        ):
            try:
                assembly_dir.rmdir()
            except OSError:
                pass
        try:
            metadata_root.rmdir()
        except OSError:
            pass



def get_species_preset(value: str | None) -> SpeciesPreset | None:
    key = canonical_species(value)
    for preset in SPECIES_PRESETS.values():
        if key == preset.key or key in {canonical_species(alias) for alias in preset.aliases}:
            return preset
    for preset in custom_species_presets():
        if key == preset.key or key in {canonical_species(alias) for alias in preset.aliases}:
            return preset
    return None


def supported_species() -> tuple[SpeciesPreset, ...]:
    return tuple(SPECIES_PRESETS.values()) + custom_species_presets()


def species_data_root() -> Path:
    override = os.getenv("HCRPROBEFORGE_REFERENCE_DIR")
    if override:
        return Path(override).expanduser()
    return Path.home() / ".hcrprobeforge" / "references"


def species_output_root(outdir: Path, species: str) -> Path:
    """Return the species-scoped results directory.

    A caller may already provide a species-scoped path (batch runners do this),
    so this helper is idempotent for the canonical species aliases.
    """
    base = Path(outdir).expanduser().resolve()
    component = species_component(species)
    if base.name == component:
        return base
    return base / component


def species_component(species: str) -> str:
    """Return the readable filesystem component used for one species.

    The HCRProbeDesign alias and the user-facing species directory are kept
    deliberately separate. The display name is used for output/cache folders;
    aliases remain internal identifiers for indexes and CLI resolution.
    """
    preset = get_species_preset(species)
    if preset is not None:
        return core.safe_name(preset.display_name) or core.safe_name(preset.scientific_name) or "reference"
    return core.safe_name(canonical_species(species)) or "reference"


def results_root(outdir: Path) -> Path:
    """Return the project-level HCRProbeForge results root."""
    return Path(outdir).expanduser().resolve()


def workflow_directory(workflow: str) -> str:
    """Return the user-facing directory for an internal workflow name."""
    value = canonical_species(workflow)
    return core.safe_name(WORKFLOW_DIRECTORIES.get(value, value)) or "workflow"


def species_run_root(outdir: Path, species: str, workflow: str = "design") -> Path:
    """Return ``<project>/runs/<species>/<workflow>``."""
    workflow_component = workflow_directory(workflow)
    return results_root(outdir) / "runs" / species_component(species) / workflow_component


def species_cache_root(outdir: Path, species: str, workflow: str = "design", run_name: str | None = None) -> Path:
    """Return the shared species/workflow cache directory.

    Caches are intentionally separate from visible run outputs, but all cache
    types now share one predictable project-level root.
    """
    root = results_root(outdir) / "cache" / species_component(species) / workflow_directory(workflow)
    if run_name:
        root /= core.safe_name(run_name) or "run"
    return root


def index_metadata_root(outdir: Path) -> Path:
    """Return the project-level folder for webapp index metadata records."""
    return results_root(outdir).parent / "index_metadata"


def hcrprobedesign_data_root() -> Path:
    """Return HCRProbeDesign's persistent data directory."""
    override = os.getenv("HCRPROBEDESIGN_DATA_DIR")
    if override:
        return Path(override).expanduser().resolve()
    return (Path.home() / ".hcrprobedesign").resolve()


def registered_index_prefix(alias: str) -> Path:
    """Find the registered Bowtie2 prefix for an HCRProbeDesign alias.

    Recent HCRProbeDesign releases record the prefix in ``HCRconfig.yaml``.
    The deterministic fallback matches the documented default layout and is
    still useful immediately after a successful build when a config parser is
    not available.
    """
    requested_alias = core.safe_name(alias) or "reference"
    equivalent_aliases = {
        core.safe_name(value) or requested_alias
        for value in species_aliases(requested_alias)
    }
    equivalent_aliases.add(requested_alias)
    data_root = hcrprobedesign_data_root()
    config_candidates = _hcrprobedesign_config_candidates()
    try:
        import yaml
    except ImportError:  # pragma: no cover - PyYAML is a declared dependency
        yaml = None
    if yaml is not None:
        for config_path in config_candidates:
            if not config_path.is_file():
                continue
            try:
                payload = yaml.safe_load(config_path.read_text(encoding="utf-8"))
            except (OSError, ValueError, yaml.YAMLError):
                continue

            def walk(value: Any) -> list[str]:
                found: list[str] = []
                if isinstance(value, dict):
                    for key, item in value.items():
                        if canonical_species(str(key)) in equivalent_aliases:
                            found.extend(collect_strings(item))
                        found.extend(walk(item))
                elif isinstance(value, (list, tuple)):
                    for item in value:
                        found.extend(walk(item))
                elif isinstance(value, str):
                    # Match an alias as a path/name component instead of as a
                    # substring, so short aliases cannot discover unrelated
                    # registered names.
                    if any(
                        re.search(
                            rf"(?<![A-Za-z0-9]){re.escape(item.casefold())}(?![A-Za-z0-9])",
                            value.casefold(),
                        )
                        for item in equivalent_aliases
                    ):
                        found.append(value)
                return found

            def collect_strings(value: Any) -> list[str]:
                if isinstance(value, str):
                    return [value]
                if isinstance(value, dict):
                    return [item for child in value.values() for item in collect_strings(child)]
                if isinstance(value, (list, tuple)):
                    return [item for child in value for item in collect_strings(child)]
                return []

            for raw in walk(payload):
                candidate = Path(raw).expanduser()
                if candidate.suffix in {".bt2", ".bt2l"}:
                    candidate = candidate.with_name(re.sub(r"\.(?:[1-6])\.bt2l?$", "", candidate.name))
                if candidate.is_absolute():
                    return candidate.resolve()
                relative = (data_root / candidate).resolve()
                if relative.is_file() or any(
                    relative.with_name(relative.name + suffix).is_file()
                    for suffix in (".1.bt2", ".1.bt2l")
                ):
                    return relative
    # Prefer an existing registered layout before choosing the canonical
    # layout. This allows an index created by an earlier layout to remain
    # usable when its exact current alias is present in the configuration.
    for candidate_alias in sorted(equivalent_aliases, key=lambda item: (item != requested_alias, item)):
        candidate = (data_root / "indices" / candidate_alias / candidate_alias).resolve()
        if _index_files_for_prefix(candidate):
            return candidate
    return (data_root / "indices" / requested_alias / requested_alias).resolve()


def registered_index_is_ready(
    alias: str,
    *,
    metadata: dict[str, Any] | None = None,
) -> bool:
    """Return whether a registered alias appears to have Bowtie2 files."""
    # Callers that have already enumerated the reference registry (notably the
    # web setup page) can provide the matching row. Avoiding a second full
    # metadata-directory scan keeps returning from a completed index build
    # responsive without changing the readiness test itself.
    if metadata is None:
        metadata = find_installed_reference(alias)
    if metadata is not None and metadata.get("status") == "ready":
        prefix = _normalise_index_prefix(
            metadata.get("index_prefix") or registered_index_prefix(alias),
            data_root=hcrprobedesign_data_root().resolve(),
        ) or registered_index_prefix(alias)
    else:
        prefix = registered_index_prefix(alias)
    return bool(_index_files_for_prefix(prefix))


def _value(record: dict[str, Any], *names: str) -> Any:
    lowered = {str(key).casefold(): value for key, value in record.items()}
    for name in names:
        if name.casefold() in lowered:
            return lowered[name.casefold()]
    return None


def _base_params(email: str | None, api_key: str | None) -> dict[str, str]:
    params = {"tool": "hcrprobeforge", "retmode": "json"}
    if email:
        params["email"] = email
    if api_key:
        params["api_key"] = api_key
    return params


def _ncbi_json(
    session: Any,
    endpoint: str,
    params: dict[str, str],
    *,
    attempts: int = 3,
) -> dict[str, Any]:
    """Request one NCBI JSON response with bounded transient-error retries."""
    last_error: Exception | None = None
    for attempt in range(attempts):
        try:
            response = session.get(endpoint, params=params, timeout=60)
            status = int(getattr(response, "status_code", 200) or 200)
            if status == 429 or status >= 500:
                raise RuntimeError(f"NCBI returned HTTP {status}.")
            response.raise_for_status()
            payload = response.json()
            if not isinstance(payload, dict):
                raise RuntimeError("NCBI returned a non-object JSON response.")
            return {str(key): value for key, value in payload.items()}
        except Exception as exc:  # requests exceptions vary by installed requests version
            last_error = exc
            if attempt + 1 >= attempts:
                break
            core.raise_if_cancelled()
            time.sleep(1.5 * (attempt + 1))
    raise RuntimeError(f"NCBI request failed after {attempts} attempts: {last_error}") from last_error


def _assembly_uids(
    accession: str,
    *,
    email: str | None,
    api_key: str | None,
    session: Any,
) -> list[str]:
    """Resolve an Assembly accession to candidate Entrez UIDs.

    ESearch normally returns one UID for an exact accession.  Keeping the
    bounded list lets :func:`_assembly_record` verify the accession in the
    returned ESummary instead of trusting the first search hit if NCBI ever
    returns a stale, broad, or otherwise unexpected result.
    """
    accession = normalize_assembly_accession(accession)
    found: list[str] = []
    for term in (f'"{accession}"[ASAC]', accession):
        params = _base_params(email, api_key) | {
            "db": "assembly",
            "term": term,
            "retmode": "json",
            "retmax": "20",
        }
        payload = _ncbi_json(session, f"{NCBI_EUTILS}/esearch.fcgi", params)
        result = payload.get("esearchresult")
        if not isinstance(result, dict):
            continue
        identifiers = result.get("idlist") or result.get("IdList") or []
        for identifier in identifiers:
            value = str(identifier).strip()
            if value and value not in found:
                found.append(value)
    return found


def _assembly_uid(
    accession: str,
    *,
    email: str | None,
    api_key: str | None,
    session: Any,
) -> str:
    """Resolve an Assembly accession to the first candidate Entrez UID."""
    return next(iter(_assembly_uids(accession, email=email, api_key=api_key, session=session)), "")


def _assembly_record(
    accession: str,
    *,
    email: str | None = None,
    api_key: str | None = None,
    session: Any = None,
) -> dict[str, Any]:
    """Fetch and validate one NCBI Assembly record.

    NCBI's Assembly ESummary endpoint is UID-based.  Resolve the accession
    with ESearch first, then request the matching document summary.
    """
    if session is None:
        if core.requests is None:
            raise RuntimeError("NCBI reference retrieval requires requests.")
        session = core.requests.Session()
        session.headers.update({"User-Agent": "hcrprobeforge"})
    accession = normalize_assembly_accession(accession)
    uids = _assembly_uids(
        accession,
        email=email,
        api_key=api_key,
        session=session,
    )
    if not uids:
        raise RuntimeError(f"NCBI Assembly accession {accession} was not found.")
    mismatches: list[str] = []
    for uid in uids:
        params = _base_params(email, api_key) | {
            "db": "assembly",
            "id": uid,
            "version": "2.0",
        }
        payload = _ncbi_json(session, f"{NCBI_EUTILS}/esummary.fcgi", params)
        result = payload.get("result") if isinstance(payload, dict) else None
        if not isinstance(result, dict):
            continue
        summary_uids = result.get("uids") or []
        summary_uid = str(summary_uids[0]) if summary_uids else str(uid)
        raw = result.get(summary_uid) if summary_uid else None
        if not isinstance(raw, dict):
            raw = next(
                (value for key, value in result.items() if str(key).isdigit() and isinstance(value, dict)),
                None,
            )
        if not isinstance(raw, dict):
            continue
        record = {str(key): value for key, value in raw.items()}
        record_accession = str(_value(record, "assemblyaccession", "AssemblyAccession") or "").strip().upper()
        if not record_accession:
            mismatches.append(f"UID {uid} had no assembly accession in its summary")
            continue
        if record_accession != accession:
            mismatches.append(f"UID {uid} returned {record_accession}")
            continue
        record["uid"] = summary_uid
        record["requested_accession"] = accession
        return record
    detail = "; ".join(mismatches) if mismatches else "no matching ESummary record"
    raise RuntimeError(
        f"NCBI did not return an Assembly record matching {accession}: {detail}. "
        "Check the accession or provide a local genome FASTA."
    )


def _assembly_ftp_url(record: dict[str, Any]) -> str:
    """Return a normalized HTTPS directory containing the genomic FASTA.

    NCBI's ESummary field is the authoritative directory location.  The
    returned URL is normalized only as a URL; the actual FASTA filename is
    derived from the directory basename below rather than from the display
    assembly name, whose spaces and punctuation do not always match NCBI's
    filename convention.
    """
    for name in (
        "FtpPath_RefSeq",
        "ftp_path_refseq",
        "FtpPath_GenBank",
        "ftp_path_genbank",
    ):
        value = str(_value(record, name) or "").strip()
        if not value or value.casefold() in {"na", "n/a", "none", "null"}:
            continue
        if value.startswith("ftp://"):
            value = "https://" + value[len("ftp://") :]
        parsed = urlsplit(value)
        if parsed.scheme not in {"http", "https"} or not parsed.netloc:
            continue
        path = parsed.path.rstrip("/")
        if not path or path == "/":
            continue
        # Decode and re-encode once so spaces or already escaped characters
        # cannot create a malformed request URL.  Underscores remain literal.
        encoded_path = quote(unquote(path), safe="/._-~")
        return urlunsplit(("https", parsed.netloc, encoded_path, "", ""))
    raise RuntimeError("The NCBI Assembly record did not provide a downloadable FTP path.")


def _assembly_download_urls(
    record: dict[str, Any],
    accession: str,
    assembly_name: str,
) -> list[str]:
    """Return deterministic FASTA URL candidates for one Assembly record.

    NCBI documents these files as ``*_genomic.fna.gz``.  The directory
    basename is preferred because it is the server's canonical, already
    normalized assembly label.  Display-name and whitespace-normalized
    fallbacks cover older or unusual records whose summary label differs from
    the FTP directory basename.
    """
    directory = _assembly_ftp_url(record).rstrip("/") + "/"
    parsed = urlsplit(directory)
    directory_name = unquote(parsed.path.rstrip("/").rsplit("/", 1)[-1])
    if not directory_name:
        raise RuntimeError("NCBI returned an FTP path without an assembly directory name.")
    actual_accession = str(
        _value(record, "assemblyaccession", "AssemblyAccession") or accession
    ).strip().upper()
    components: list[str] = []

    def add_component(value: str) -> None:
        value = str(value or "").strip()
        if value and value not in components:
            components.append(value)

    add_component(_normalise_assembly_component(directory_name))
    add_component(directory_name)
    normalized_name = _normalise_assembly_component(assembly_name)
    if normalized_name:
        add_component(f"{actual_accession}_{normalized_name}")
    if str(assembly_name or "").strip():
        add_component(f"{actual_accession}_{assembly_name}")
    return [
        urljoin(directory, quote(f"{component}_genomic.fna.gz", safe="._-~"))
        for component in components
    ]


def _assembly_annotation_download_urls(
    record: dict[str, Any],
    accession: str,
    assembly_name: str,
) -> list[str]:
    """Return deterministic GFF3/GFF download candidates for an assembly."""
    directory = _assembly_ftp_url(record).rstrip("/") + "/"
    parsed = urlsplit(directory)
    directory_name = unquote(parsed.path.rstrip("/").rsplit("/", 1)[-1])
    actual_accession = str(
        _value(record, "assemblyaccession", "AssemblyAccession") or accession
    ).strip().upper()
    components: list[str] = []

    def add(value: str) -> None:
        normalized = _normalise_assembly_component(value)
        if normalized and normalized not in components:
            components.append(normalized)

    add(directory_name)
    add(f"{actual_accession}_{assembly_name}")
    suffixes = ("_genomic.gff.gz", "_genomic.gff3.gz")
    return [
        urljoin(directory, quote(f"{component}{suffix}", safe="._-~"))
        for component in components
        for suffix in suffixes
    ]


def _assembly_annotation_download_urls_from_genome_url(genome_url: str) -> list[str]:
    """Derive annotation candidates from a validated genomic FASTA URL."""
    parsed = urlsplit(str(genome_url or ""))
    if not parsed.path:
        return []
    filename = unquote(parsed.path.rsplit("/", 1)[-1])
    stem = re.sub(r"_genomic\.fna(?:\.gz)?$", "", filename, flags=re.IGNORECASE)
    if not stem:
        return []
    directory = str(genome_url).rsplit("/", 1)[0] + "/"
    return [
        urljoin(directory, quote(stem + suffix, safe="._-~"))
        for suffix in ("_genomic.gff.gz", "_genomic.gff3.gz")
    ]


def _discover_assembly_download_urls(record: dict[str, Any], *, session: Any = None) -> list[str]:
    """Discover genomic FASTA links from an NCBI directory listing.

    This is a last-resort fallback after the deterministic names return 404.
    It protects future assemblies whose directory contains a valid genomic
    FASTA but whose naming differs from the usual accession/assembly pattern.
    """
    directory = _assembly_ftp_url(record).rstrip("/") + "/"
    if session is None:
        if core.requests is None:
            return []
        session = core.requests.Session()
        session.headers.update({"User-Agent": "hcrprobeforge"})
    try:
        response = session.get(directory, timeout=60)
        status = int(getattr(response, "status_code", 200) or 200)
        if status != 200:
            return []
        names = re.findall(
            r"(?:href=[\"']?)?([^\"'<>\s]+_genomic\.fna\.gz)(?:[\"']|\s|$)",
            response.text,
            flags=re.IGNORECASE,
        )
    except Exception:
        return []
    urls: list[str] = []
    for name in names:
        filename = unquote(name.rsplit("/", 1)[-1])
        url = urljoin(directory, quote(filename, safe="._-~"))
        if url not in urls:
            urls.append(url)
    return urls


def _assembly_metadata(record: dict[str, Any], preset: SpeciesPreset | None) -> dict[str, Any]:
    accession = str(
        _value(record, "assemblyaccession", "AssemblyAccession")
        or record.get("requested_accession")
        or ""
    )
    return {
        "schema_version": 1,
        "source": "NCBI Assembly ESummary",
        "species": preset.key if preset else None,
        "display_name": preset.display_name if preset else None,
        "scientific_name": preset.scientific_name if preset else str(_value(record, "speciesname", "organism") or ""),
        "assembly": str(_value(record, "assemblyname", "AssemblyName") or ""),
        "assembly_accession": accession,
        "ncbi_uid": record.get("uid"),
        "refseq_category": _value(record, "refseq_category", "RefSeq_category"),
        "assembly_status": _value(record, "assemblystatus", "AssemblyStatus"),
        "ftp_directory": _assembly_ftp_url(record),
        "record": core.json_ready(record),
        "retrieved_at_utc": core.utc_now_iso(),
    }


def _download(
    url: str,
    destination: Path,
    *,
    progress_callback: ProgressCallback | None = None,
    session: Any = None,
) -> None:
    if session is None:
        if core.requests is None:
            raise RuntimeError("Reference retrieval requires requests.")
        session = core.requests.Session()
        session.headers.update({"User-Agent": "hcrprobeforge"})
    destination.parent.mkdir(parents=True, exist_ok=True)
    last_error: Exception | None = None
    for attempt in range(3):
        temporary = destination.with_name(destination.name + f".part.{os.getpid()}")
        try:
            response = session.get(url, stream=True, timeout=120)
            status = int(getattr(response, "status_code", 200) or 200)
            if status in {404, 410}:
                raise _AssemblyDownloadError(
                    f"NCBI did not provide the requested genome file (HTTP {status}): {url}",
                    status_code=status,
                )
            if status == 429 or status >= 500:
                raise _AssemblyDownloadError(
                    f"NCBI temporarily rejected the genome download (HTTP {status}): {url}",
                    status_code=status,
                )
            response.raise_for_status()
            total = int(response.headers.get("Content-Length") or 0)
            received = 0
            with temporary.open("wb") as handle:
                for chunk in response.iter_content(chunk_size=1024 * 1024):
                    core.raise_if_cancelled()
                    if not chunk:
                        continue
                    handle.write(chunk)
                    received += len(chunk)
                    if progress_callback is not None:
                        progress_callback(
                            {
                                "phase": "download",
                                "completed": received,
                                "total": total or None,
                                "message": f"Downloaded {received / 1024 / 1024:.1f} MB"
                                + (f" of {total / 1024 / 1024:.1f} MB" if total else ""),
                            }
                        )
            if not received:
                raise _AssemblyDownloadError(
                    f"NCBI returned an empty genome file: {url}",
                    status_code=status,
                )
            if total and received != total:
                raise _AssemblyDownloadError(
                    f"Genome download was incomplete ({received} of {total} bytes): {url}",
                    status_code=status,
                )
            os.replace(temporary, destination)
            return
        except _AssemblyDownloadError as exc:
            last_error = exc
            temporary.unlink(missing_ok=True)
            if exc.status_code in {404, 410}:
                raise
            if attempt + 1 >= 3:
                raise RuntimeError(
                    f"Genome download failed after 3 attempts: {url}: {exc}"
                ) from exc
        except Exception as exc:
            last_error = exc
            temporary.unlink(missing_ok=True)
            if attempt + 1 >= 3:
                raise RuntimeError(
                    f"Genome download failed after 3 attempts: {url}: {exc}"
                ) from exc
        core.raise_if_cancelled()
        time.sleep(1.5 * (attempt + 1))
    raise RuntimeError(f"Genome download failed: {url}: {last_error}") from last_error


def _decompress_genome(source: Path, destination: Path) -> Path:
    if source.suffix != ".gz":
        return source
    with gzip.open(source, "rb") as input_handle, destination.open("wb") as output_handle:
        shutil.copyfileobj(input_handle, output_handle, length=1024 * 1024)
    return destination


def _run_build_genome_index(
    *,
    species: str,
    fasta: Path,
    threads: int,
    force: bool,
    progress_callback: ProgressCallback | None,
    custom_species_key: str | None = None,
) -> tuple[int, str, str]:
    executable = shutil.which("buildGenomeIndex")
    if not executable:
        raise RuntimeError(
            "buildGenomeIndex was not found. Install HCRProbeDesign in the active environment."
        )
    if custom_species_key and not force:
        recovered = _recover_stale_custom_registration(
            species,
            custom_species_key=custom_species_key,
        )
        if recovered and progress_callback:
            progress_callback(
                {
                    "phase": "build",
                    "completed": 0,
                    "total": None,
                    "message": "Removed a stale custom index registration before rebuilding",
                }
            )
    command = [executable, "--species", species, "--fasta", str(fasta), "--threads", str(threads)]
    if force:
        command.append("--force")
    if progress_callback:
        progress_callback({"phase": "build", "completed": 0, "total": None, "message": "Building and registering the Bowtie2 index"})
    process = subprocess.Popen(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        start_new_session=(os.name != "nt"),
    )
    try:
        while process.poll() is None:
            core.raise_if_cancelled()
            time.sleep(0.5)
    except core.RunCancelled:
        core._stop_design_process(process, force=True)
        process.communicate()
        raise
    stdout, stderr = process.communicate()
    core.raise_if_cancelled()
    if process.returncode:
        raise RuntimeError(
            f"buildGenomeIndex failed with exit code {process.returncode}: {stderr.strip()}"
        )
    if progress_callback:
        progress_callback({"phase": "build", "completed": 1, "total": 1, "message": "Index built and registered"})
    return int(process.returncode or 0), stdout, stderr


def _metadata_path(preset: SpeciesPreset, assembly: str) -> Path:
    return species_data_root() / preset.key / core.safe_name(assembly) / "reference.json"


def _write_metadata(path: Path, metadata: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + f".tmp.{os.getpid()}")
    try:
        temporary.write_text(
            json.dumps(core.json_ready(metadata), indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def fetch_and_build_index(
    species: str,
    *,
    assembly_accession: str | None = None,
    assembly_name: str | None = None,
    index_alias: str | None = None,
    threads: int = 4,
    force: bool = False,
    build_annotation_database: bool = False,
    email: str | None = None,
    api_key: str | None = None,
    progress_callback: ProgressCallback | None = None,
) -> dict[str, Any]:
    """Download an NCBI assembly and register its Bowtie2 index.

    ``assembly_accession`` is optional for the built-in reference.  Supplying
    it makes a separate assembly entry, which is useful when a new build is
    released or when a project must reproduce an older build.  Unless an
    explicit ``index_alias`` is supplied, alternate builds receive a unique
    alias derived from the preset and assembly name.
    """
    preset = get_species_preset(species)
    if preset is None:
        supported = ", ".join(preset.key for preset in supported_species())
        raise ValueError(f"Automatic index retrieval supports: {supported}. Use index build for a local FASTA.")
    if threads < 1:
        raise ValueError("threads must be at least 1")
    if core.requests is None:
        raise RuntimeError("Automatic reference retrieval requires requests.")
    session = core.requests.Session()
    session.headers.update({"User-Agent": "hcrprobeforge"})
    raw_requested_accession = assembly_accession or preset.assembly_accession
    if not raw_requested_accession:
        raise ValueError(
            f"Species {preset.display_name} has no saved NCBI assembly accession. "
            "Provide --assembly-accession or build the index from a local genome FASTA."
        )
    requested_accession = normalize_assembly_accession(raw_requested_accession)
    try:
        record = _assembly_record(
            requested_accession,
            email=email,
            api_key=api_key,
            session=session,
        )
    except RuntimeError as exc:
        raise RuntimeError(
            f"Could not retrieve and validate NCBI Assembly {requested_accession}. "
            "Verify the accession in NCBI Assembly or download its genomic FASTA "
            "and use the local index-build workflow. Details: "
            f"{exc}"
        ) from exc
    metadata = _assembly_metadata(record, preset)
    assembly = str(metadata["assembly"] or assembly_name or preset.assembly_name)
    if assembly_name and _normalise_assembly_component(assembly_name).casefold() != _normalise_assembly_component(assembly).casefold():
        raise ValueError(
            f"NCBI accession {requested_accession} is named {assembly!r}, not {assembly_name!r}."
        )
    preset_accession = (
        normalize_assembly_accession(preset.assembly_accession)
        if preset.assembly_accession
        else ""
    )
    default_request = requested_accession == preset_accession and not index_alias
    alias = core.safe_name(index_alias or (preset.key if default_request else f"{preset.key}_{assembly}"))
    if not alias:
        raise ValueError("index alias cannot be empty")
    metadata_path = _metadata_path(preset, assembly)
    if metadata_path.exists() and not force:
        try:
            existing = json.loads(metadata_path.read_text(encoding="utf-8"))
            if existing.get("status") == "ready" and (
                not index_alias
                or core.safe_name(str(existing.get("index_species_alias") or "")) == alias
            ):
                existing_database = Path(
                    str(existing.get("annotation_database_path") or metadata_path.parent / "annotation.sqlite")
                )
                from . import premrna

                # A stale or interrupted SQLite file is not a usable database.
                # Treat it exactly like a missing database so the explicit
                # index-build opt-in repairs it instead of returning a false
                # ready state.
                if build_annotation_database and not premrna.annotation_database_is_valid(existing_database):
                    existing["metadata_path"] = str(metadata_path)
                    existing["reference_data_directory"] = str(metadata_path.parent)
                    source_args = argparse.Namespace(email=email, api_key=api_key)
                    _genome, annotation, _source = premrna._download_sources(
                        session,
                        existing,
                        preset.key,
                        source_args,
                    )
                    database = premrna.annotation_index_path(annotation)
                    if not premrna.annotation_database_is_valid(database):
                        premrna.build_annotation_index(
                            annotation,
                            database,
                            progress_callback=progress_callback,
                        )
                    existing.update(
                        {
                            "annotation_gff_path": str(annotation),
                            "annotation_database_path": str(database),
                            "annotation_database_status": "ready",
                            "annotation_database_requested": True,
                            "metadata_path": str(metadata_path),
                        }
                    )
                    _write_metadata(metadata_path, existing)
                elif build_annotation_database:
                    existing.update(
                        {
                            "annotation_database_path": str(existing_database),
                            "annotation_database_status": "ready",
                            "annotation_database_requested": True,
                            "metadata_path": str(metadata_path),
                        }
                    )
                    _write_metadata(metadata_path, existing)
                return existing | {
                    "reused": True,
                    "metadata_path": str(metadata_path),
                    "index_prefix": str(existing.get("index_prefix") or registered_index_prefix(alias)),
                }
        except (OSError, ValueError, TypeError):
            pass
    previous_metadata = metadata_path.read_bytes() if metadata_path.is_file() else None
    reference_root = species_data_root()
    design_data_root = hcrprobedesign_data_root()
    reference_snapshot = core._snapshot_filesystem(reference_root)
    design_snapshot = core._snapshot_filesystem(design_data_root)
    rollback_root, rollback_files = _build_rollback_backup(
        alias,
        metadata_path,
        metadata_path.parent,
    )
    try:
        with tempfile.TemporaryDirectory(prefix=f"hcrprobeforge-{preset.key}-") as temporary:
            temporary_dir = Path(temporary)
            genome_url = ""
            genome_fasta: Path | None = None
            download_errors: list[str] = []
            candidate_urls = _assembly_download_urls(record, requested_accession, assembly)
            for candidate_index, candidate_url in enumerate(candidate_urls, start=1):
                compressed = temporary_dir / f"assembly_candidate_{candidate_index}.fna.gz"
                decompressed = temporary_dir / f"assembly_candidate_{candidate_index}.fna"
                try:
                    if progress_callback:
                        progress_callback({"phase": "download", "message": f"Downloading genome FASTA candidate {candidate_index}…"})
                    _download(candidate_url, compressed, progress_callback=progress_callback, session=session)
                    candidate_fasta = _decompress_genome(compressed, decompressed)
                    if candidate_fasta.stat().st_size <= 0:
                        raise RuntimeError("NCBI returned an empty decompressed genome FASTA")
                    genome_url = candidate_url
                    genome_fasta = candidate_fasta
                    break
                except _AssemblyDownloadError as exc:
                    download_errors.append(str(exc))
                    compressed.unlink(missing_ok=True)
                    decompressed.unlink(missing_ok=True)
                    if exc.status_code not in {404, 410}:
                        raise
                except (OSError, EOFError, gzip.BadGzipFile, RuntimeError) as exc:
                    download_errors.append(f"{candidate_url}: {exc}")
                    compressed.unlink(missing_ok=True)
                    decompressed.unlink(missing_ok=True)
            if genome_fasta is None:
                discovered_urls = _discover_assembly_download_urls(record, session=session)
                for candidate_index, candidate_url in enumerate(discovered_urls, start=len(candidate_urls) + 1):
                    compressed = temporary_dir / f"assembly_candidate_{candidate_index}.fna.gz"
                    decompressed = temporary_dir / f"assembly_candidate_{candidate_index}.fna"
                    try:
                        _download(candidate_url, compressed, progress_callback=progress_callback, session=session)
                        candidate_fasta = _decompress_genome(compressed, decompressed)
                        if candidate_fasta.stat().st_size <= 0:
                            raise RuntimeError("NCBI returned an empty decompressed genome FASTA")
                        genome_url = candidate_url
                        genome_fasta = candidate_fasta
                        break
                    except _AssemblyDownloadError as exc:
                        download_errors.append(str(exc))
                        compressed.unlink(missing_ok=True)
                        decompressed.unlink(missing_ok=True)
                        if exc.status_code not in {404, 410}:
                            raise
                    except (OSError, EOFError, gzip.BadGzipFile, RuntimeError) as exc:
                        download_errors.append(f"{candidate_url}: {exc}")
                        compressed.unlink(missing_ok=True)
                        decompressed.unlink(missing_ok=True)
            if genome_fasta is None:
                detail = "\n".join(f"- {item}" for item in download_errors[-8:])
                raise RuntimeError(
                    f"NCBI provided Assembly {requested_accession}, but its genomic FASTA could not be located. "
                    "The FTP directory was valid but none of the expected genomic FASTA names were available. "
                    "Verify the accession at NCBI or use the Build index workflow with a downloaded FASTA."
                    + (f"\nTried:\n{detail}" if detail else "")
                )

            # Reference assets are persistent assembly data, not request
            # staging files. Keep only the uncompressed forms: pysam/Bowtie2
            # use them directly and retaining duplicate .gz files wastes disk.
            reference_directory = metadata_path.parent
            reference_directory.mkdir(parents=True, exist_ok=True)
            persistent_genome = reference_directory / "genome.fna"
            shutil.copyfile(genome_fasta, persistent_genome)
            genome_fasta = persistent_genome
            annotation_url = ""
            annotation_path = reference_directory / "annotation.gff"
            annotation_errors: list[str] = []
            annotation_candidates = _assembly_annotation_download_urls(record, requested_accession, assembly)
            for candidate_index, candidate_url in enumerate(annotation_candidates, start=1):
                compressed = temporary_dir / f"annotation_candidate_{candidate_index}.gff.gz"
                decompressed = temporary_dir / f"annotation_candidate_{candidate_index}.gff"
                try:
                    if progress_callback:
                        progress_callback({"phase": "download", "message": "Downloading genomic annotation…"})
                    _download(candidate_url, compressed, progress_callback=progress_callback, session=session)
                    candidate_annotation = _decompress_genome(compressed, decompressed)
                    if candidate_annotation.stat().st_size <= 0:
                        raise RuntimeError("NCBI returned an empty genomic annotation")
                    shutil.copyfile(candidate_annotation, annotation_path)
                    annotation_url = candidate_url
                    break
                except _AssemblyDownloadError as exc:
                    annotation_errors.append(str(exc))
                    compressed.unlink(missing_ok=True)
                    decompressed.unlink(missing_ok=True)
                    if exc.status_code not in {404, 410}:
                        raise
                except (OSError, EOFError, gzip.BadGzipFile, RuntimeError) as exc:
                    annotation_errors.append(f"{candidate_url}: {exc}")
                    compressed.unlink(missing_ok=True)
                    decompressed.unlink(missing_ok=True)
            annotation_available = annotation_path.is_file()
            if not annotation_available:
                # A genome index remains valid for the mature-transcript
                # workflow even when an assembly has no downloadable GFF3 at
                # this moment. Preserve the diagnostic details so the first
                # Pre-mRNA run can explain/retry the annotation fetch.
                annotation_path = None

            metadata.update(
                {
                    "genome_url": genome_url,
                    "annotation_url": annotation_url,
                    "genome_fasta_path": str(genome_fasta),
                    "annotation_gff_path": str(annotation_path) if annotation_path else None,
                    "reference_data_directory": str(reference_directory),
                    "genome_fasta_sha256": core.sha256_file(genome_fasta),
                    "genome_fasta_bytes": genome_fasta.stat().st_size,
                    "annotation_gff_bytes": annotation_path.stat().st_size if annotation_path else None,
                    "annotation_status": "ready" if annotation_available else "unavailable_at_index_build",
                    "annotation_database_path": None,
                    "annotation_database_status": "deferred" if annotation_available else "unavailable",
                    "annotation_database_requested": bool(build_annotation_database),
                    "annotation_download_errors": annotation_errors[-8:],
                    "index_species_alias": alias,
                    "requested_assembly_accession": requested_accession,
                    "status": "building",
                    "metadata_path": str(metadata_path),
                }
            )
            _write_metadata(metadata_path, metadata)
            from . import premrna

            if annotation_path is not None and build_annotation_database:
                premrna.build_annotation_index(annotation_path, progress_callback=progress_callback)
                metadata["annotation_database_path"] = str(premrna.annotation_index_path(annotation_path))
                metadata["annotation_database_status"] = "ready"
                _write_metadata(metadata_path, metadata)
            _run_build_genome_index(
                species=alias,
                fasta=genome_fasta,
                threads=threads,
                force=force,
                custom_species_key=(preset.key if preset.key not in SPECIES_PRESETS else None),
                progress_callback=progress_callback,
            )
            if annotation_path is not None and not build_annotation_database:
                # A rebuilt annotation must never be paired with a stale
                # SQLite lookup. The next Pre-mRNA run will recreate it.
                premrna.annotation_index_path(annotation_path).unlink(missing_ok=True)
        metadata["status"] = "ready"
        metadata["completed_at_utc"] = core.utc_now_iso()
        metadata["metadata_path"] = str(metadata_path)
        metadata["index_prefix"] = str(registered_index_prefix(alias))
        metadata["reused"] = False
        _write_metadata(metadata_path, metadata)
        _discard_build_rollback_backup(rollback_root)
        return metadata
    except Exception:
        # A failed build must not leave an apparently usable registry entry or
        # a newly created partial HCRProbeDesign data tree behind. Existing
        # references are restored byte-for-byte; files belonging to unrelated
        # aliases are never touched.
        core._remove_new_filesystem_entries(reference_root, reference_snapshot)
        core._remove_new_filesystem_entries(design_data_root, design_snapshot)
        _restore_build_rollback_backup(rollback_root, rollback_files)
        if previous_metadata is not None:
            metadata_path.parent.mkdir(parents=True, exist_ok=True)
            metadata_path.write_bytes(previous_metadata)
        raise


def build_local_index(
    species: str,
    fasta: Path,
    *,
    threads: int = 4,
    force: bool = False,
    assembly: str | None = None,
    annotation: Path | None = None,
    build_annotation_database: bool = False,
    progress_callback: ProgressCallback | None = None,
) -> dict[str, Any]:
    """Build/register an index from local FASTA and optional matching GFF3."""
    fasta = Path(fasta).expanduser().resolve()
    if not fasta.is_file():
        raise ValueError(f"genome FASTA not found: {fasta}")
    if threads < 1:
        raise ValueError("threads must be at least 1")
    alias = canonical_species(species)
    if not alias:
        raise ValueError("species alias cannot be empty")
    preset = get_species_preset(alias) or SpeciesPreset(alias, alias, "", assembly or alias, "")
    assembly_name = assembly or preset.assembly_name
    metadata_path = _metadata_path(preset, assembly_name)
    reference_directory = metadata_path.parent
    persistent_genome = reference_directory / "genome.fna"
    persistent_annotation: Path | None = None
    if annotation is not None:
        annotation = Path(annotation).expanduser().resolve()
        if not annotation.is_file():
            raise ValueError(f"genomic annotation not found: {annotation}")
        persistent_annotation = reference_directory / "annotation.gff"
    if build_annotation_database and persistent_annotation is None:
        raise ValueError(
            "Building the annotation database requires a matching local GFF3 annotation. "
            "Provide one, or leave the option off and let the first Pre-mRNA run retrieve/use annotation when available."
        )
    metadata = {
        "schema_version": 2,
        "status": "building",
        "source": "local FASTA and GFF3" if persistent_annotation else "local FASTA",
        "species": alias,
        "display_name": preset.display_name,
        "scientific_name": preset.scientific_name,
        "assembly": assembly_name,
        "assembly_accession": None,
        "local_fasta": str(fasta),
        "genome_fasta_path": str(persistent_genome),
        "genome_fasta_sha256": None,
        "genome_fasta_bytes": None,
        "annotation_gff_path": str(persistent_annotation) if persistent_annotation else None,
        "annotation_gff_bytes": None,
        "annotation_database_path": None,
        "annotation_database_status": "deferred" if persistent_annotation else "unavailable",
        "annotation_database_requested": bool(build_annotation_database),
        "reference_data_directory": str(reference_directory),
        "index_species_alias": alias,
        "retrieved_at_utc": core.utc_now_iso(),
    }
    previous_metadata = metadata_path.read_bytes() if metadata_path.is_file() else None
    reference_root = species_data_root()
    design_data_root = hcrprobedesign_data_root()
    reference_snapshot = core._snapshot_filesystem(reference_root)
    design_snapshot = core._snapshot_filesystem(design_data_root)
    rollback_root, rollback_files = _build_rollback_backup(
        alias,
        metadata_path,
        reference_directory,
    )
    staged_genome = reference_directory / f".genome.fna.tmp.{os.getpid()}"
    staged_annotation = reference_directory / f".annotation.gff.tmp.{os.getpid()}" if annotation is not None else None
    staged_database = reference_directory / f".annotation.sqlite.tmp.{os.getpid()}" if annotation is not None and build_annotation_database else None
    try:
        reference_directory.mkdir(parents=True, exist_ok=True)
        if fasta.suffix.casefold() == ".gz":
            _decompress_genome(fasta, staged_genome)
        else:
            shutil.copyfile(fasta, staged_genome)
        if annotation is not None and staged_annotation is not None:
            if annotation.suffix.casefold() == ".gz":
                _decompress_genome(annotation, staged_annotation)
            else:
                shutil.copyfile(annotation, staged_annotation)
        metadata["genome_fasta_sha256"] = core.sha256_file(staged_genome)
        metadata["genome_fasta_bytes"] = staged_genome.stat().st_size
        if staged_annotation is not None:
            metadata["annotation_gff_bytes"] = staged_annotation.stat().st_size
        _write_metadata(metadata_path, metadata)
        if persistent_annotation is not None and build_annotation_database:
            from . import premrna

            # Prepare the optional annotation lookup before invoking the
            # external genome-index builder. This keeps the web progress phases
            # truthful and avoids an expensive Bowtie2 build when the supplied
            # annotation is malformed or unusable.
            premrna.build_annotation_index(staged_annotation, staged_database, progress_callback=progress_callback)
        _run_build_genome_index(
            species=alias,
            fasta=staged_genome,
            threads=threads,
            force=force,
            custom_species_key=(preset.key if preset.key not in SPECIES_PRESETS else None),
            progress_callback=progress_callback,
        )
        os.replace(staged_genome, persistent_genome)
        if staged_annotation is not None and persistent_annotation is not None:
            os.replace(staged_annotation, persistent_annotation)
        if staged_database is not None:
            os.replace(staged_database, persistent_annotation.with_name("annotation.sqlite"))
            metadata["annotation_database_path"] = str(persistent_annotation.with_name("annotation.sqlite"))
            metadata["annotation_database_status"] = "ready"
        elif not build_annotation_database:
            # Do not leave a lookup built from an older local annotation or
            # genome beside the newly registered reference.
            persistent_annotation_path = persistent_annotation or reference_directory / "annotation.gff"
            persistent_annotation_path.with_name("annotation.sqlite").unlink(missing_ok=True)
            if persistent_annotation is None:
                persistent_annotation_path.unlink(missing_ok=True)
        metadata["status"] = "ready"
        metadata["completed_at_utc"] = core.utc_now_iso()
        metadata["metadata_path"] = str(metadata_path)
        metadata["index_prefix"] = str(registered_index_prefix(alias))
        _write_metadata(metadata_path, metadata)
        _discard_build_rollback_backup(rollback_root)
        return metadata
    except Exception:
        core._remove_new_filesystem_entries(reference_root, reference_snapshot)
        core._remove_new_filesystem_entries(design_data_root, design_snapshot)
        _restore_build_rollback_backup(rollback_root, rollback_files)
        if previous_metadata is not None:
            metadata_path.parent.mkdir(parents=True, exist_ok=True)
            metadata_path.write_bytes(previous_metadata)
        raise
    finally:
        staged_genome.unlink(missing_ok=True)
        if staged_annotation is not None:
            staged_annotation.unlink(missing_ok=True)
        if staged_database is not None:
            staged_database.unlink(missing_ok=True)


def list_installed_references() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for path in _reference_metadata_paths():
        data = _read_reference_metadata(path)
        if data is not None:
            rows.append(data)
    return rows


def find_installed_reference(species: str | None) -> dict[str, Any] | None:
    """Return ready metadata for a registered genome alias, if available."""
    requested = canonical_species(species)
    if not requested or requested in REMOVED_BUILTIN_ALIASES:
        return None
    # A custom preset is selected by its stable key/display name, while the
    # registered reference may be keyed by an assembly alias. Treat all of
    # the preset's identifiers as equivalent when resolving its metadata.
    requested_keys = {requested}
    preset = get_species_preset(species)
    if preset is not None:
        requested_keys.update(
            canonical_species(value)
            for value in (
                preset.key,
                preset.display_name,
                preset.scientific_name,
                *preset.aliases,
            )
            if value
        )
    for row in list_installed_references():
        aliases = {
            canonical_species(row.get("index_species_alias")),
            canonical_species(row.get("species")),
            canonical_species(row.get("display_name")),
            canonical_species(row.get("scientific_name")),
        }
        if requested_keys.intersection(aliases) and row.get("status") == "ready":
            return row
    return None


def _progress_print(update: dict[str, object]) -> None:
    message = str(update.get("message") or update.get("phase") or "working")
    print(message)


def species_main(argv: list[str] | None = None) -> int:
    """Create and inspect persistent user-defined species presets."""
    parser = argparse.ArgumentParser(
        prog="hcrprobeforge species",
        description="Create species presets for organisms not included in the built-in list.",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)
    list_parser = subparsers.add_parser("list", help="List built-in and user-created species presets")
    list_parser.set_defaults(_command="list")
    add_parser = subparsers.add_parser("add", help="Save a user species preset")
    add_parser.add_argument("--alias", required=True, help="Stable filesystem/index alias, e.g. octopus_v1")
    add_parser.add_argument("--display-name", required=True, help="Name shown in the webapp")
    add_parser.add_argument("--organism", required=True, help="Exact NCBI organism name used for transcript lookup")
    add_parser.add_argument("--assembly-accession", help="Optional NCBI GCF_/GCA_ assembly accession")
    add_parser.add_argument("--assembly-name", help="Optional assembly name")
    add_parser.add_argument("--genome-fasta", type=Path, help="Optional local genome FASTA to keep for the index workflow")
    add_parser.set_defaults(_command="add")
    args = parser.parse_args(argv)
    try:
        if args._command == "list":
            for preset in supported_species():
                source = "built-in" if preset.key in SPECIES_PRESETS else "user"
                print(f"{preset.key}: {preset.display_name} [{preset.scientific_name}] ({source})")
            return 0
        preset = register_custom_species(
            args.alias,
            args.display_name,
            args.organism,
            assembly_accession=args.assembly_accession,
            assembly_name=args.assembly_name,
            genome_fasta=args.genome_fasta,
        )
        print(json.dumps(asdict(preset), indent=2, sort_keys=True))
        return 0
    except (OSError, ValueError) as exc:
        print(f"ERROR: {exc}", file=os.sys.stderr)
        return 1


def index_main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="hcrprobeforge index",
        description="Download NCBI genome references and build/register Bowtie2 indexes.",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)
    list_parser = subparsers.add_parser("list", help="List HCRProbeForge-managed references")
    list_parser.set_defaults(_command="list")
    for name, help_text in (("fetch", "Download an approved NCBI assembly and build its index"), ("build", "Build/register an index from a local genome FASTA")):
        child = subparsers.add_parser(name, help=help_text)
        child.add_argument(
            "--species",
            required=True,
            help="Built-in or user-created species alias; registered custom aliases are accepted",
        )
        child.add_argument("--threads", type=int, default=4)
        child.add_argument("--force", action="store_true", help="Rebuild and replace the registered entry")
        child.add_argument(
            "--with-annotation-db",
            action="store_true",
            help="Also build the optional genomic annotation SQLite database for future Pre-mRNA/intron designs",
        )
        child.add_argument("--email", default=os.getenv("NCBI_EMAIL"))
        child.add_argument("--api-key", default=os.getenv("NCBI_API_KEY"))
        if name == "fetch":
            child.add_argument("--assembly-accession", help="Optional NCBI assembly accession; default is the approved reference")
            child.add_argument("--assembly-name", help="Optional expected NCBI assembly name; prevents accidental mismatches")
            child.add_argument("--index-alias", help="Optional HCRProbeDesign alias for this assembly")
            child.set_defaults(_command="fetch")
        else:
            child.add_argument("--fasta", type=Path, required=True)
            child.add_argument(
                "--annotation",
                type=Path,
                help="Optional matching genomic GFF3 or compressed GFF3 for Pre-mRNA design",
            )
            child.add_argument("--assembly", help="Assembly display name for a local FASTA")
            child.set_defaults(_command="build")
    args = parser.parse_args(argv)
    try:
        if args._command == "list":
            rows = list_installed_references()
            if not rows:
                print("No HCRProbeForge-managed references are installed.")
            for row in rows:
                alias = row.get("index_species_alias") or row.get("species")
                organism = row.get("scientific_name") or row.get("display_name") or ""
                print(
                    f"{alias}: {organism}; {row.get('assembly')} "
                    f"({row.get('assembly_accession') or 'local'}) - {row.get('status')}"
                )
            return 0
        callback = _progress_print
        if args._command == "fetch":
            result = fetch_and_build_index(
                args.species,
                assembly_accession=args.assembly_accession,
                assembly_name=args.assembly_name,
                index_alias=args.index_alias,
                threads=args.threads,
                force=args.force,
                build_annotation_database=args.with_annotation_db,
                email=args.email,
                api_key=args.api_key,
                progress_callback=callback,
            )
        else:
            result = build_local_index(
                args.species,
                args.fasta,
                threads=args.threads,
                force=args.force,
                assembly=args.assembly,
                annotation=args.annotation,
                build_annotation_database=args.with_annotation_db,
                progress_callback=callback,
            )
        print(json.dumps(core.json_ready(result), indent=2, sort_keys=True))
        return 0
    except (OSError, RuntimeError, ValueError) as exc:
        print(f"ERROR: {exc}", file=os.sys.stderr)
        return 1


def main(argv: list[str] | None = None) -> int:
    return index_main(argv)
