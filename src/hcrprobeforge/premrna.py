"""Genomic Pre-mRNA target preparation.

Pre-mRNA mode represents the selected transcript as its full genomic
transcript span (exons and introns in transcript 5'-to-3' order). Candidate
selection is subsequently restricted to complete intronic intervals. The
ordinary mature-transcript design path is not used to make this decision.

NCBI GFF3 files can be hundreds of megabytes. A compact SQLite lookup is
therefore created once while a reference is prepared. A design then performs
an indexed accession query, uses a random-access FASTA reader, and caches the
assembled target.
"""

from __future__ import annotations

import gzip
import json
import os
import re
import shutil
import sqlite3
from pathlib import Path
from typing import Any
from urllib.parse import unquote


_ACCESSION_RE = re.compile(
    r"(?:^|[^A-Za-z])((?:NM|NR|XM|XR)_\d+(?:\.\d+)?)(?:$|[^0-9])",
    re.IGNORECASE,
)
_REFSEQ_RE = re.compile(
    r"(?:RefSeq:|refseq:)((?:NM|NR|XM|XR)_\d+(?:\.\d+)?)",
    re.IGNORECASE,
)
_TRANSCRIPT_FEATURES = {
    "transcript",
    "mrna",
    "lnc_rna",
    "ncrna",
    "primary_transcript",
}

ANNOTATION_SCHEMA_VERSION = "3"
INTRONIC_TARGET_TYPE = "pre-mrna"
WHOLE_TARGET_TYPE = "pre-mrna-whole"


def is_premrna_target_type(value: str | None) -> bool:
    """Return whether a target type uses the genomic Pre-mRNA reference."""
    return str(value or "").strip().casefold() in {
        INTRONIC_TARGET_TYPE,
        WHOLE_TARGET_TYPE,
    }


def is_whole_premrna_target_type(value: str | None) -> bool:
    """Return whether a Pre-mRNA design includes exons and introns as eligible."""
    return str(value or "").strip().casefold() == WHOLE_TARGET_TYPE


def premrna_selection_mode(value: str | None) -> str:
    """Return the stable metadata/cache name for a Pre-mRNA selection mode."""
    return "whole" if is_whole_premrna_target_type(value) else "intronic"


def _core():
    from . import core

    return core


def _references():
    from . import references

    return references


def parse_intron_selection(value: str | None) -> tuple[int, ...] | None:
    """Parse a 1-based comma-separated intron list; blank means all introns."""
    if value is None or not str(value).strip():
        return None
    selected: list[int] = []
    for token in str(value).split(","):
        token = token.strip()
        if not token:
            continue
        if not token.isdigit() or int(token) < 1:
            raise ValueError(
                "--premrna-introns must contain comma-separated positive intron "
                "numbers, for example 1,3,5."
            )
        number = int(token)
        if number not in selected:
            selected.append(number)
    return tuple(selected) or None


def _attributes(text: str) -> dict[str, str]:
    result: dict[str, str] = {}
    for field in text.rstrip(";\n").split(";"):
        if not field:
            continue
        if "=" in field:
            key, value = field.split("=", 1)
        elif " " in field:
            key, value = field.split(" ", 1)
            value = value.strip('"')
        else:
            continue
        result[key.strip()] = unquote(value.strip())
    return result


def _attribute_values(attributes: dict[str, str], *keys: str) -> list[str]:
    wanted = {key.casefold() for key in keys}
    values: list[str] = []
    for key, value in attributes.items():
        if key.casefold() not in wanted:
            continue
        for item in str(value).split(","):
            item = item.strip()
            if item and item not in values:
                values.append(item)
    return values


def _refseq_aliases(value: str) -> list[str]:
    aliases = [value]
    aliases.extend(match.group(1) for match in _REFSEQ_RE.finditer(value))
    aliases.extend(match.group(1) for match in _ACCESSION_RE.finditer(value))
    return aliases


def _normalise_alias(value: str) -> str:
    return str(value).strip().casefold()


def _accession_matches(value: str, accession: str) -> bool:
    left = _normalise_alias(value)
    right = _normalise_alias(accession)
    return bool(left and right and left == right)


def _accession_compatible(value: str, accession: str) -> bool:
    """Match an exact version, or an explicitly unversioned annotation alias.

    A versioned annotation must never silently satisfy a different requested
    transcript version. Some GFF3 files omit the version entirely, so an
    unversioned annotation alias remains a deliberate compatibility fallback.
    """
    left = _normalise_alias(value)
    right = _normalise_alias(accession)
    if not left or not right:
        return False
    if left == right:
        return True
    return "." not in left and left == right.split(".", 1)[0]


def _gff_transcript_model(
    seqid: str,
    start: int,
    end: int,
    strand: str,
    attributes: dict[str, str],
) -> dict[str, Any] | None:
    values: list[str] = []
    values.extend(_attribute_values(attributes, "ID", "transcript_id", "Name", "gene_id"))
    for value in list(values):
        values.extend(_refseq_aliases(value))
    values = list(dict.fromkeys(values))
    if not values:
        return None
    gene_id_values = _attribute_values(attributes, "gene_id")
    if not gene_id_values:
        gene_id_values = [value for value in _attribute_values(attributes, "Dbxref") if "GeneID:" in value]
    gene_symbol_values = _attribute_values(attributes, "gene", "gene_name", "locus_tag")
    return {
        "seqid": seqid,
        "strand": strand if strand in {"+", "-"} else "+",
        "start": start,
        "end": end,
        "aliases": values,
        "gene_id": gene_id_values[0] if gene_id_values else None,
        "gene_symbol": gene_symbol_values[0] if gene_symbol_values else None,
        "exons": [],
    }


def _target_model_match(model: dict[str, Any], record: dict[str, Any]) -> tuple[bool, bool]:
    accession = str(record.get("accession") or "")
    exact = bool(accession) and any(
        _accession_compatible(alias, accession) for alias in model.get("aliases", [])
    )
    if exact:
        return True, False
    gene_id = str(record.get("gene_id") or "").casefold()
    model_gene_id = str(model.get("gene_id") or "").casefold()
    gene_symbol = str(record.get("resolved_gene_symbol") or record.get("gene_symbol") or "").casefold()
    model_gene_symbol = str(model.get("gene_symbol") or "").casefold()
    fallback = bool(
        (gene_id and gene_id == model_gene_id)
        or (gene_symbol and gene_symbol == model_gene_symbol)
    )
    return False, fallback


def _normalise_model(model: dict[str, Any]) -> dict[str, Any] | None:
    exons = sorted(
        {
            (int(exon["start"]), int(exon["end"]))
            for exon in model.get("exons", [])
            if int(exon["start"]) <= int(exon["end"])
        }
    )
    if not exons:
        return None
    model["exons"] = [{"start": start, "end": end} for start, end in exons]
    model["start"] = min(int(model["start"]), exons[0][0])
    model["end"] = max(int(model["end"]), exons[-1][1])
    return model


def parse_gff3_models(path: Path) -> list[dict[str, Any]]:
    """Parse all transcript/exon relationships (legacy/general utility)."""
    transcripts: list[dict[str, Any]] = []
    aliases: dict[str, dict[str, Any]] = {}
    pending_exons: list[tuple[list[str], dict[str, Any]]] = []
    opener = gzip.open if path.suffix.casefold() == ".gz" else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) != 9:
                continue
            seqid, _source, feature, start, end, _score, strand, _phase, attribute_text = fields
            try:
                start_i, end_i = int(start), int(end)
            except ValueError:
                continue
            attributes = _attributes(attribute_text)
            if feature.casefold() in _TRANSCRIPT_FEATURES:
                model = _gff_transcript_model(seqid, start_i, end_i, strand, attributes)
                if model is None:
                    continue
                transcripts.append(model)
                for alias in model["aliases"]:
                    aliases[_normalise_alias(alias)] = model
            elif feature.casefold() == "exon":
                parents = _attribute_values(attributes, "Parent", "transcript_id", "ID")
                if parents:
                    pending_exons.append((parents, {"seqid": seqid, "start": start_i, "end": end_i}))

    for parents, exon in pending_exons:
        matched: dict[int, dict[str, Any]] = {}
        for parent in parents:
            normalized = _normalise_alias(parent)
            model = aliases.get(normalized) or aliases.get(normalized.split(".", 1)[0])
            if model is not None:
                matched[id(model)] = model
        for model in matched.values():
            if model["seqid"] == exon["seqid"]:
                model["exons"].append({"start": exon["start"], "end": exon["end"]})
    return [usable for model in transcripts if (usable := _normalise_model(model)) is not None]


def annotation_index_path(gff_path: Path) -> Path:
    """Return the compact lookup path stored beside one assembly annotation."""
    return gff_path.with_name("annotation.sqlite")


def annotation_database_is_valid(database_path: Path, *, deep: bool = False) -> bool:
    """Validate a completed annotation lookup before it is reused.

    Existence alone is not sufficient: an interrupted older build can leave a
    readable SQLite file containing staging tables but no metadata.  Requiring
    the final schema marker and indexes lets callers rebuild those files
    atomically without affecting mature-transcript designs.  The normal check
    uses SQLite's fast structural ``quick_check`` so a webapp status badge or
    ordinary Pre-mRNA run does not scan every row repeatedly.  Callers doing a
    release/diagnostic audit can request the slower full ``integrity_check``.
    """
    database_path = Path(database_path)
    if not database_path.is_file():
        return False
    required_tables = {"transcripts", "aliases", "exons", "metadata"}
    required_indexes = {"aliases_alias", "exons_transcript", "transcripts_symbol"}
    connection: sqlite3.Connection | None = None
    try:
        with sqlite3.connect(database_path) as connection:
            pragma = "integrity_check" if deep else "quick_check"
            if connection.execute(f"PRAGMA {pragma}").fetchone()[0] != "ok":
                return False
            tables = {
                str(row[0])
                for row in connection.execute(
                    "SELECT name FROM sqlite_master WHERE type='table'"
                )
            }
            indexes = {
                str(row[0])
                for row in connection.execute(
                    "SELECT name FROM sqlite_master WHERE type='index'"
                )
            }
            if not required_tables.issubset(tables) or not required_indexes.issubset(indexes):
                return False
            schema = connection.execute(
                "SELECT value FROM metadata WHERE key='schema_version'"
            ).fetchone()
            if schema is None or str(schema[0]) != ANNOTATION_SCHEMA_VERSION:
                return False
            if connection.execute("SELECT COUNT(*) FROM transcripts").fetchone()[0] < 1:
                return False
            # A database with no exon rows cannot support either Pre-mRNA mode.
            if connection.execute("SELECT COUNT(*) FROM exons").fetchone()[0] < 1:
                return False
    except (OSError, sqlite3.Error, TypeError, ValueError):
        return False
    return True


def annotation_database_is_structurally_ready(database_path: Path) -> bool:
    """Perform a fast, read-only readiness check for a webapp status badge.

    The full validator intentionally runs SQLite ``quick_check`` and is used
    before a database is consumed by a Pre-mRNA design.  Running that scan for
    every species while rendering the setup page makes returning from a result
    page unnecessarily slow on machines with several large references.  The
    webapp only needs to identify a completed final-format database here; the
    design path still performs the full validator before reuse.
    """
    database_path = Path(database_path)
    if not database_path.is_file():
        return False
    required_tables = {"transcripts", "aliases", "exons", "metadata"}
    required_indexes = {"aliases_alias", "exons_transcript", "transcripts_symbol"}
    try:
        with sqlite3.connect(database_path) as connection:
            tables = {
                str(row[0])
                for row in connection.execute(
                    "SELECT name FROM sqlite_master WHERE type='table'"
                )
            }
            indexes = {
                str(row[0])
                for row in connection.execute(
                    "SELECT name FROM sqlite_master WHERE type='index'"
                )
            }
            if not required_tables.issubset(tables) or not required_indexes.issubset(indexes):
                return False
            schema = connection.execute(
                "SELECT value FROM metadata WHERE key='schema_version'"
            ).fetchone()
            if schema is None or str(schema[0]) != ANNOTATION_SCHEMA_VERSION:
                return False
            return bool(
                connection.execute("SELECT 1 FROM transcripts LIMIT 1").fetchone()
                and connection.execute("SELECT 1 FROM exons LIMIT 1").fetchone()
            )
    except (OSError, sqlite3.Error, TypeError, ValueError):
        return False


def _gff_opener(path: Path):
    return gzip.open if path.suffix.casefold() == ".gz" else open


def _remove_annotation_sqlite_temporary(path: Path) -> None:
    """Remove an interrupted temporary SQLite database and sidecars."""
    for suffix in ("", "-journal", "-wal", "-shm"):
        path.with_name(path.name + suffix).unlink(missing_ok=True)


def _annotation_progress(
    progress_callback: Any,
    *,
    completed: int,
    total: int | None,
    line_count: int,
    message: str,
) -> None:
    if progress_callback is None:
        return
    progress_callback(
        {
            "phase": "annotation",
            "completed": completed,
            "total": total,
            "records_scanned": line_count,
            "message": message,
        }
    )


def build_annotation_index(
    gff_path: Path,
    database_path: Path | None = None,
    *,
    progress_callback: Any = None,
) -> Path:
    """Build an atomic, accession-indexed SQLite transcript/exon lookup.

    The input is streamed once. Transcript IDs and aliases are resolved in
    memory, while rows are written to SQLite in batches. This preserves the
    exon-before-transcript handling required by GFF3 while avoiding a SQL
    lookup and transaction-visible insert for every individual record.
    """
    gff_path = Path(gff_path)
    database_path = Path(database_path or annotation_index_path(gff_path))
    database_path.parent.mkdir(parents=True, exist_ok=True)
    temporary = database_path.with_name(database_path.name + f".tmp.{os.getpid()}")
    _remove_annotation_sqlite_temporary(temporary)
    # Report physical bytes read for both plain and gzip-compressed GFF3.  A
    # compressed annotation cannot expose a useful uncompressed byte total,
    # but its underlying file position still gives an honest download/read
    # progress signal instead of leaving the webapp at one fixed percentage.
    file_total = gff_path.stat().st_size
    batch_size = 10_000
    transcript_count = 0
    exon_count = 0
    line_count = 0
    try:
        _annotation_progress(
            progress_callback,
            completed=0,
            total=file_total,
            line_count=0,
            message="Reading genomic annotation…",
        )
        connection = sqlite3.connect(temporary)
        with connection:
            connection.executescript(
                """
                PRAGMA journal_mode=OFF;
                PRAGMA synchronous=OFF;
                PRAGMA temp_store=MEMORY;
                CREATE TABLE transcripts (
                    id INTEGER PRIMARY KEY,
                    seqid TEXT NOT NULL,
                    start INTEGER NOT NULL,
                    end INTEGER NOT NULL,
                    strand TEXT NOT NULL,
                    gene_id TEXT,
                    gene_symbol TEXT,
                    aliases_json TEXT NOT NULL
                );
                CREATE TABLE aliases(alias TEXT NOT NULL, transcript_id INTEGER NOT NULL);
                CREATE TABLE exons(transcript_id INTEGER NOT NULL, seqid TEXT NOT NULL, start INTEGER NOT NULL, end INTEGER NOT NULL);
                CREATE TABLE metadata(key TEXT PRIMARY KEY, value TEXT NOT NULL);
                """
            )
            connection.execute(
                "INSERT INTO metadata(key,value) VALUES(?,?)",
                ("schema_version", ANNOTATION_SCHEMA_VERSION),
            )
            transcript_rows: list[tuple[Any, ...]] = []
            alias_rows: list[tuple[str, int]] = []
            exon_rows: list[tuple[int, str, int, int]] = []
            pending_exons: list[tuple[list[str], str, int, int]] = []
            alias_to_transcripts: dict[str, set[int]] = {}
            next_transcript_id = 1

            def flush_rows() -> None:
                if transcript_rows:
                    connection.executemany(
                        "INSERT INTO transcripts(id,seqid,start,end,strand,gene_id,gene_symbol,aliases_json) VALUES(?,?,?,?,?,?,?,?)",
                        transcript_rows,
                    )
                    transcript_rows.clear()
                if alias_rows:
                    connection.executemany(
                        "INSERT INTO aliases(alias,transcript_id) VALUES(?,?)",
                        alias_rows,
                    )
                    alias_rows.clear()
                if exon_rows:
                    connection.executemany(
                        "INSERT INTO exons(transcript_id,seqid,start,end) VALUES(?,?,?,?)",
                        exon_rows,
                    )
                    exon_rows.clear()

            with _gff_opener(gff_path)(gff_path, "rt", encoding="utf-8", errors="replace") as handle:
                for raw in handle:
                    line_count += 1
                    if line_count % 50_000 == 0:
                        _core().raise_if_cancelled()
                        completed = line_count
                        total = None
                        try:
                            if gff_path.suffix.casefold() == ".gz":
                                raw_handle = getattr(getattr(handle, "buffer", None), "fileobj", None)
                                completed = int(raw_handle.tell()) if raw_handle is not None else int(handle.tell())
                            else:
                                completed = int(handle.tell())
                            completed = min(file_total, max(0, completed))
                            total = file_total
                        except (OSError, ValueError, AttributeError):
                            completed = line_count
                        suffix = f" ({completed / total:.0%} of file)" if total else ""
                        _annotation_progress(
                            progress_callback,
                            completed=completed,
                            total=total,
                            line_count=line_count,
                            message=f"Indexing genomic annotation{suffix}; {line_count:,} records scanned…",
                        )
                    if not raw.strip() or raw.startswith("#"):
                        continue
                    fields = raw.rstrip("\r\n").split("\t")
                    if len(fields) != 9:
                        continue
                    feature_lower = fields[2].casefold()
                    # Most genomic GFF3 rows are not needed for the transcript
                    # lookup. Avoid parsing their attributes; on large NCBI
                    # annotations this is a substantial CPU and allocation
                    # reduction without changing the retained feature set.
                    if feature_lower not in _TRANSCRIPT_FEATURES and feature_lower != "exon":
                        continue
                    seqid, _source, feature, start, end, _score, strand, _phase, attribute_text = fields
                    try:
                        start_i, end_i = int(start), int(end)
                    except ValueError:
                        continue
                    attributes = _attributes(attribute_text)
                    if feature_lower in _TRANSCRIPT_FEATURES:
                        model = _gff_transcript_model(seqid, start_i, end_i, strand, attributes)
                        if model is None:
                            continue
                        transcript_id = next_transcript_id
                        next_transcript_id += 1
                        transcript_rows.append(
                            (
                                transcript_id,
                                model["seqid"], model["start"], model["end"], model["strand"],
                                model.get("gene_id"), model.get("gene_symbol"),
                                json.dumps(model["aliases"], ensure_ascii=False),
                            )
                        )
                        for alias in model["aliases"]:
                            normalized = _normalise_alias(alias)
                            alias_to_transcripts.setdefault(normalized, set()).add(transcript_id)
                            alias_rows.append((normalized, transcript_id))
                            base_alias = normalized.split(".", 1)[0]
                            if base_alias != normalized:
                                alias_to_transcripts.setdefault(base_alias, set()).add(transcript_id)
                                alias_rows.append((base_alias, transcript_id))
                        transcript_count += 1
                    elif feature_lower == "exon":
                        parents = _attribute_values(attributes, "Parent", "transcript_id", "ID")
                        matched: set[int] = set()
                        for parent in parents:
                            normalized_parent = _normalise_alias(parent)
                            matched.update(alias_to_transcripts.get(normalized_parent, set()))
                            matched.update(alias_to_transcripts.get(normalized_parent.split(".", 1)[0], set()))
                        if matched:
                            for transcript_id in matched:
                                exon_rows.append((transcript_id, seqid, start_i, end_i))
                                exon_count += 1
                        else:
                            pending_exons.append((parents, seqid, start_i, end_i))
                    if len(transcript_rows) >= batch_size or len(alias_rows) >= batch_size or len(exon_rows) >= batch_size:
                        flush_rows()

            for parents, seqid, start_i, end_i in pending_exons:
                matched: set[int] = set()
                for parent in parents:
                    normalized_parent = _normalise_alias(parent)
                    matched.update(alias_to_transcripts.get(normalized_parent, set()))
                    matched.update(alias_to_transcripts.get(normalized_parent.split(".", 1)[0], set()))
                for transcript_id in matched:
                    exon_rows.append((transcript_id, seqid, start_i, end_i))
                    exon_count += 1
                if len(exon_rows) >= batch_size:
                    flush_rows()
            flush_rows()
            exon_count = int(connection.execute("SELECT COUNT(*) FROM exons").fetchone()[0])
            connection.executescript(
                "CREATE INDEX aliases_alias ON aliases(alias);"
                "CREATE INDEX exons_transcript ON exons(transcript_id);"
                "CREATE INDEX transcripts_symbol ON transcripts(gene_symbol);"
            )
            connection.commit()
        # The context manager commits/rolls back but does not close the
        # connection. Close it explicitly before replacing the database so
        # readers can reopen the completed file immediately. In particular,
        # never leave a persistent annotation database in EXCLUSIVE mode.
        connection.close()
        connection = None
        os.replace(temporary, database_path)
    except Exception:
        if connection is not None:
            connection.close()
        _remove_annotation_sqlite_temporary(temporary)
        raise
    _remove_annotation_sqlite_temporary(temporary)
    if progress_callback is not None:
        progress_callback(
            {
                "phase": "annotation",
                "completed": 1,
                "total": 1,
                "records_scanned": line_count,
                "message": f"Annotation lookup ready ({transcript_count:,} transcripts, {exon_count:,} exons)",
            }
        )
    return database_path


def _model_from_row(connection: sqlite3.Connection, row: sqlite3.Row) -> dict[str, Any] | None:
    model = {
        "seqid": row["seqid"],
        "strand": row["strand"],
        "start": int(row["start"]),
        "end": int(row["end"]),
        "aliases": json.loads(row["aliases_json"]),
        "gene_id": row["gene_id"],
        "gene_symbol": row["gene_symbol"],
        "exons": [],
    }
    for exon in connection.execute(
        "SELECT seqid,start,end FROM exons WHERE transcript_id=? ORDER BY start,end",
        (int(row["id"]),),
    ):
        if exon["seqid"] == model["seqid"]:
            model["exons"].append({"start": int(exon["start"]), "end": int(exon["end"])})
    return _normalise_model(model)


def _sqlite_models(database_path: Path, record: dict[str, Any]) -> list[dict[str, Any]]:
    if not database_path.is_file():
        return []
    connection = sqlite3.connect(database_path)
    connection.row_factory = sqlite3.Row
    try:
        accession = str(record.get("accession") or "")
        normalized = _normalise_alias(accession)
        exact_ids: set[int] = set()
        for alias in (normalized, normalized.split(".", 1)[0]):
            exact_ids.update(
                int(row["transcript_id"])
                for row in connection.execute("SELECT transcript_id FROM aliases WHERE alias=?", (alias,))
            )
        rows: list[sqlite3.Row] = []
        if exact_ids:
            placeholders = ",".join("?" for _ in exact_ids)
            rows = list(connection.execute(
                f"SELECT * FROM transcripts WHERE id IN ({placeholders})", tuple(exact_ids)
            ))
            compatible_models = [
                model
                for row in rows
                if (model := _model_from_row(connection, row)) is not None
                and any(
                    _accession_compatible(alias, accession)
                    for alias in model.get("aliases", [])
                )
            ]
            if compatible_models:
                return compatible_models
            rows = []
        if not rows:
            gene_id = str(record.get("gene_id") or "")
            gene_symbol = str(record.get("resolved_gene_symbol") or record.get("gene_symbol") or "")
            if gene_id:
                rows.extend(connection.execute(
                    "SELECT * FROM transcripts WHERE lower(gene_id)=lower(?)", (gene_id,)
                ))
            if not rows and gene_symbol:
                rows.extend(connection.execute(
                    "SELECT * FROM transcripts WHERE lower(gene_symbol)=lower(?)", (gene_symbol,)
                ))
        return [model for row in rows if (model := _model_from_row(connection, row)) is not None]
    finally:
        connection.close()


def parse_gff3_model_for_record(
    path: Path,
    record: dict[str, Any],
    *,
    database_path: Path | None = None,
    progress_callback: Any = None,
) -> list[dict[str, Any]]:
    """Return only models matching one accession/gene, using the compact index."""
    database_path = Path(database_path or annotation_index_path(path))
    try:
        if annotation_database_is_valid(database_path):
            models = _sqlite_models(database_path, record)
            if models:
                return models
        build_annotation_index(path, database_path, progress_callback=progress_callback)
        return _sqlite_models(database_path, record)
    except (OSError, sqlite3.Error, ValueError, json.JSONDecodeError):
        pass

    # Small/local fallback when an SQLite database cannot be created.
    exact_models: list[dict[str, Any]] = []
    fallback_models: list[dict[str, Any]] = []
    aliases: dict[str, list[dict[str, Any]]] = {}
    pending_exons: list[tuple[list[str], dict[str, Any]]] = []
    line_count = 0
    with _gff_opener(path)(path, "rt", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line_count += 1
            if line_count % 100_000 == 0:
                _core().raise_if_cancelled()
                if progress_callback is not None:
                    progress_callback({
                        "phase": "annotation",
                        "completed": line_count,
                        "total": None,
                        "message": "Scanning the genomic annotation for the selected transcript…",
                    })
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) != 9:
                continue
            seqid, _source, feature, start, end, _score, strand, _phase, attribute_text = fields
            try:
                start_i, end_i = int(start), int(end)
            except ValueError:
                continue
            attributes = _attributes(attribute_text)
            if feature.casefold() == "exon":
                parents = _attribute_values(attributes, "Parent", "transcript_id", "ID")
                if parents:
                    pending_exons.append((parents, {"seqid": seqid, "start": start_i, "end": end_i}))
                continue
            if feature.casefold() not in _TRANSCRIPT_FEATURES:
                continue
            model = _gff_transcript_model(seqid, start_i, end_i, strand, attributes)
            if model is None:
                continue
            exact, fallback = _target_model_match(model, record)
            if exact:
                exact_models.append(model)
            elif fallback:
                fallback_models.append(model)
            for alias in model["aliases"]:
                aliases.setdefault(_normalise_alias(alias), []).append(model)
                aliases.setdefault(_normalise_alias(alias).split(".", 1)[0], []).append(model)
    selected = exact_models or fallback_models
    selected_ids = {id(model) for model in selected}
    for parents, exon in pending_exons:
        matched: dict[int, dict[str, Any]] = {}
        for parent in parents:
            for model in aliases.get(_normalise_alias(parent), []):
                if id(model) in selected_ids:
                    matched[id(model)] = model
            for model in aliases.get(_normalise_alias(parent).split(".", 1)[0], []):
                if id(model) in selected_ids:
                    matched[id(model)] = model
        for model in matched.values():
            if model["seqid"] == exon["seqid"]:
                model["exons"].append({"start": exon["start"], "end": exon["end"]})
    return [usable for model in selected if (usable := _normalise_model(model)) is not None]


def _model_cache_path(gff_path: Path, record: dict[str, Any]) -> Path:
    """Return the project-visible cache path for one transcript model."""
    accession = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(record.get("accession") or "transcript"))
    assembly = _core().safe_name(gff_path.parent.name) or "reference"
    root = _core().hcrprobeforge_cache_root() / "premrna" / "models" / assembly
    return root / f"transcript_{accession}_model.json"


def _legacy_model_cache_path(gff_path: Path, record: dict[str, Any]) -> Path:
    accession = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(record.get("accession") or "transcript"))
    return gff_path.with_name(f"transcript_{accession}_model.json")


def _load_or_parse_target_models(
    gff_path: Path,
    record: dict[str, Any],
    *,
    database_path: Path | None = None,
    progress_callback: Any = None,
) -> list[dict[str, Any]]:
    cache_path = _model_cache_path(gff_path, record)
    legacy_cache_path = _legacy_model_cache_path(gff_path, record)
    try:
        stat = gff_path.stat()
    except OSError:
        stat = None
    for candidate_path in (cache_path, legacy_cache_path):
        try:
            if stat is None:
                break
            raw_cache = candidate_path.read_text(encoding="utf-8")
            cached = json.loads(raw_cache)
            if (
                cached.get("accession") == record.get("accession")
                and int(cached.get("gff3_size", -1)) == int(stat.st_size)
                and int(cached.get("gff3_mtime_ns", -1)) == int(stat.st_mtime_ns)
                and isinstance(cached.get("models"), list)
                and cached["models"]
            ):
                models = cached["models"]
                if candidate_path != cache_path or "\n" not in raw_cache:
                    try:
                        cache_path.parent.mkdir(parents=True, exist_ok=True)
                        cache_path.write_text(
                            json.dumps(cached, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
                            encoding="utf-8",
                        )
                        candidate_path.unlink(missing_ok=True)
                    except OSError:
                        pass
                legacy_cache_path.unlink(missing_ok=True)
                return models
        except (OSError, TypeError, ValueError, json.JSONDecodeError):
            continue
    models = parse_gff3_model_for_record(
        gff_path,
        record,
        database_path=database_path,
        progress_callback=progress_callback,
    )
    if models:
        temporary: Path | None = None
        try:
            stat = gff_path.stat()
            cache_path.parent.mkdir(parents=True, exist_ok=True)
            temporary = cache_path.with_name(cache_path.name + f".tmp.{os.getpid()}")
            temporary.write_text(
                json.dumps({
                    "schema_version": 2,
                    "accession": record.get("accession"),
                    "gff3_size": stat.st_size,
                    "gff3_mtime_ns": stat.st_mtime_ns,
                    "models": models,
                }, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
                encoding="utf-8",
            )
            os.replace(temporary, cache_path)
            if legacy_cache_path != cache_path:
                legacy_cache_path.unlink(missing_ok=True)
        except OSError:
            if temporary is not None:
                temporary.unlink(missing_ok=True)
        legacy_cache_path.unlink(missing_ok=True)
    return models


def _open_fasta(path: Path):
    """Return a random-access FASTA reader, creating its .fai when needed."""
    try:
        import pysam  # type: ignore

        index = Path(str(path) + ".fai")
        if not index.exists():
            _core().notify_progress("fasta", "Creating a random-access FASTA index…")
            pysam.faidx(str(path))
        handle = pysam.FastaFile(str(path))

        class _PysamAccessor:
            def __enter__(self):
                return self

            def __exit__(self, *_args):
                handle.close()

            def fetch(self, seqid: str, start: int, end: int) -> str:
                try:
                    return str(handle.fetch(seqid, start - 1, end)).upper()
                except (KeyError, ValueError) as exc:
                    raise RuntimeError(
                        f"Genome FASTA does not contain contig {seqid!r} referenced by the annotation."
                    ) from exc

        return _PysamAccessor()
    except Exception:
        pass

    try:
        from Bio import SeqIO  # type: ignore

        indexed = SeqIO.index(str(path), "fasta")

        class _BiopythonAccessor:
            def __enter__(self):
                return self

            def __exit__(self, *_args):
                indexed.close()

            def fetch(self, seqid: str, start: int, end: int) -> str:
                try:
                    return str(indexed[seqid].seq[start - 1 : end]).upper()
                except KeyError as exc:
                    raise RuntimeError(
                        f"Genome FASTA does not contain contig {seqid!r} referenced by the annotation."
                    ) from exc

        return _BiopythonAccessor()
    except Exception:
        records: dict[str, str] = {}
        current: str | None = None
        chunks: list[str] = []
        with path.open(encoding="utf-8", errors="replace") as handle:
            for raw in handle:
                line = raw.strip()
                if line.startswith(">"):
                    if current is not None:
                        records[current] = "".join(chunks).upper()
                    current = line[1:].split()[0]
                    chunks = []
                elif current is not None:
                    chunks.append(line)
            if current is not None:
                records[current] = "".join(chunks).upper()

        class _MemoryAccessor:
            def __enter__(self):
                return self

            def __exit__(self, *_args):
                return None

            def fetch(self, seqid: str, start: int, end: int) -> str:
                try:
                    return records[seqid][start - 1 : end]
                except KeyError as exc:
                    raise RuntimeError(
                        f"Genome FASTA does not contain contig {seqid!r} referenced by the annotation."
                    ) from exc

        return _MemoryAccessor()


def _reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1].upper()


def _select_model(models: list[dict[str, Any]], record: dict[str, Any]) -> dict[str, Any]:
    accession = str(record.get("accession") or "")
    exact = [
        model for model in models
        if any(_accession_compatible(alias, accession) for alias in model.get("aliases", []))
    ]
    if exact:
        exact.sort(key=lambda model: (-(len(model.get("exons") or [])), model["seqid"], model["start"]))
        return exact[0]
    gene_id = str(record.get("gene_id") or "").casefold()
    gene_symbol = str(record.get("resolved_gene_symbol") or record.get("gene_symbol") or "").casefold()
    fallback = [
        model for model in models
        if (gene_id and gene_id == str(model.get("gene_id") or "").casefold())
        or (gene_symbol and gene_symbol == str(model.get("gene_symbol") or "").casefold())
    ]
    if len(fallback) == 1:
        return fallback[0]
    raise RuntimeError(
        f"The exact RefSeq transcript {accession or 'selected transcript'} was not found in the matching genomic annotation. "
        "The assembly may have a different annotation release; rebuild the species index from the matching NCBI assembly or provide a mature transcript target."
    )


def _reference_directory(metadata: dict[str, Any], species: str) -> Path:
    raw = str(metadata.get("reference_data_directory") or "").strip()
    if raw:
        path = Path(raw).expanduser()
        # Releases before the persistent reference layout stored genome and
        # annotation assets below ``<assembly>/premrna``. Treat that folder
        # as a legacy cache, not as the canonical reference directory.
        return path.parent if path.name.casefold() == "premrna" else path
    metadata_path = str(metadata.get("metadata_path") or "").strip()
    if metadata_path:
        return Path(metadata_path).expanduser().resolve().parent
    references = _references()
    assembly = references.core.safe_name(str(metadata.get("assembly") or metadata.get("assembly_accession") or "assembly"))
    return references.species_data_root() / (references.core.safe_name(species) or "species") / (assembly or "assembly")


def _reference_label(metadata: dict[str, Any], species: str) -> str:
    species_label = str(
        metadata.get("display_name")
        or metadata.get("scientific_name")
        or species
        or "the selected species"
    ).strip()
    assembly = str(metadata.get("assembly") or metadata.get("assembly_accession") or "selected assembly").strip()
    accession = str(metadata.get("assembly_accession") or "").strip()
    suffix = f" ({accession})" if accession else ""
    return f"{species_label} · {assembly}{suffix}"


def _annotation_error(metadata: dict[str, Any], species: str, detail: str) -> RuntimeError:
    return RuntimeError(
        f"Pre-mRNA design cannot proceed for {_reference_label(metadata, species)}: {detail} "
        "Intronic Pre-mRNA design requires a genomic GFF3 annotation from the same "
        "assembly with the selected transcript and at least two linked exons; whole "
        "Pre-mRNA can use a valid single-exon model but still needs the matching "
        "annotation. Rebuild the matching NCBI assembly index, or provide a matching "
        "local GFF3 when building a local index. Mature-transcript design does not "
        "require this annotation."
    )


def _decompress_once(
    source: Path,
    destination: Path,
    *,
    remove_source: bool = False,
) -> Path:
    """Decompress *source* atomically without mutating external inputs.

    Only compressed files created by HCRProbeForge in a managed reference
    directory may be removed, and callers must opt into that cleanup.
    """
    temporary = destination.with_name(destination.name + f".tmp.{os.getpid()}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    try:
        with gzip.open(source, "rb") as input_handle, temporary.open("wb") as output_handle:
            shutil.copyfileobj(input_handle, output_handle, length=1024 * 1024)
        os.replace(temporary, destination)
        if remove_source:
            source.unlink(missing_ok=True)
    finally:
        temporary.unlink(missing_ok=True)
    return destination


def _canonical_asset_path(
    metadata: dict[str, Any],
    reference_dir: Path,
    metadata_key: str,
    filename: str,
) -> Path:
    """Reuse/migrate legacy Pre-mRNA assets into the canonical assembly dir."""
    destination = reference_dir / filename
    configured = str(metadata.get(metadata_key) or "").strip()
    candidates: list[Path] = []
    if configured:
        candidates.append(Path(configured).expanduser())
    candidates.extend(
        [
            destination,
            reference_dir / "premrna" / filename,
            reference_dir / "premrna" / f"{filename}.gz",
        ]
    )
    seen: set[str] = set()
    unique_candidates: list[Path] = []
    for candidate in candidates:
        key = str(candidate)
        if key not in seen:
            seen.add(key)
            unique_candidates.append(candidate)

    if destination.is_file():
        # If metadata still points to the old managed location, remove only
        # that stale duplicate after the canonical file has been found.
        if configured and Path(configured).expanduser() != destination:
            configured_path = Path(configured).expanduser()
            if configured_path.parent == reference_dir / "premrna":
                configured_path.unlink(missing_ok=True)
        return destination

    for candidate in unique_candidates:
        if not candidate.is_file() or candidate == destination:
            continue
        if candidate.suffix.casefold() == ".gz":
            _decompress_once(
                candidate,
                destination,
                remove_source=candidate.parent == reference_dir / "premrna",
            )
        else:
            temporary = destination.with_name(destination.name + f".tmp.{os.getpid()}")
            try:
                shutil.copyfile(candidate, temporary)
                os.replace(temporary, destination)
            finally:
                temporary.unlink(missing_ok=True)
            if candidate.parent == reference_dir / "premrna":
                candidate.unlink(missing_ok=True)
        return destination
    return destination


def _remove_legacy_reference_duplicates(reference_dir: Path) -> None:
    """Remove only old managed genome/annotation copies under ``premrna``."""
    legacy_dir = reference_dir / "premrna"
    if not legacy_dir.is_dir():
        return
    for filename in (
        "genome.fna",
        "genome.fna.gz",
        "annotation.gff",
        "annotation.gff.gz",
    ):
        (legacy_dir / filename).unlink(missing_ok=True)
    try:
        legacy_dir.rmdir()
    except OSError:
        # Leave unrelated legacy/cache files untouched.
        pass


def _forward_annotation_progress(update: dict[str, Any]) -> None:
    phase = str(update.get("phase") or "annotation")
    message = str(update.get("message") or "Indexing genomic annotation…")
    details = {
        key: value
        for key, value in update.items()
        if key not in {"phase", "message"}
    }
    _core().notify_progress(phase, message, **details)


def _download_sources(
    session: Any,
    metadata: dict[str, Any],
    species: str,
    args: Any,
) -> tuple[Path, Path, dict[str, Any]]:
    references = _references()
    def download_progress(update: dict[str, Any]) -> None:
        _core().notify_progress(
            str(update.get("phase") or "download"),
            str(update.get("message") or "Downloading genomic reference data…"),
            **{
                key: value
                for key, value in update.items()
                if key not in {"phase", "message"}
            },
        )

    reference_dir = _reference_directory(metadata, species)
    reference_dir.mkdir(parents=True, exist_ok=True)
    genome_path = _canonical_asset_path(
        metadata,
        reference_dir,
        "genome_fasta_path",
        "genome.fna",
    )
    gff_path = _canonical_asset_path(
        metadata,
        reference_dir,
        "annotation_gff_path",
        "annotation.gff",
    )
    _remove_legacy_reference_duplicates(reference_dir)

    local = str(metadata.get("local_fasta") or metadata.get("local_genome_fasta") or "").strip()
    if local and not genome_path.is_file():
        local_path = Path(local).expanduser()
        if local_path.is_file():
            if local_path.suffix.casefold() == ".gz":
                _decompress_once(local_path, genome_path)
            else:
                shutil.copyfile(local_path, genome_path)

    genome_url = str(metadata.get("genome_url") or "").strip()
    annotation_url = str(metadata.get("annotation_url") or "").strip()
    ftp_directory = str(metadata.get("ftp_directory") or "").strip()
    assembly_accession = str(metadata.get("assembly_accession") or "").strip()
    assembly = str(metadata.get("assembly") or "").strip()
    # Do not refresh NCBI metadata when both persistent assets are already
    # present. This keeps a first Pre-mRNA run offline and avoids needless
    # network work for references prepared by an older release.
    if (not genome_path.is_file() or not gff_path.is_file()) and assembly_accession:
        try:
            assembly_record = references._assembly_record(
                assembly_accession,
                email=getattr(args, "email", None),
                api_key=getattr(args, "api_key", None),
                session=session,
            )
            ftp_directory = ftp_directory or references._assembly_ftp_url(assembly_record)
            genome_candidates = references._assembly_download_urls(assembly_record, assembly_accession, assembly)
            genome_url = genome_url or (genome_candidates[0] if genome_candidates else "")
            annotation_candidates = references._assembly_annotation_download_urls(
                assembly_record, assembly_accession, assembly
            )
            annotation_url = annotation_url or (annotation_candidates[0] if annotation_candidates else "")
        except Exception as exc:
            if not genome_path.is_file():
                raise RuntimeError(
                    f"Could not refresh the genomic reference metadata for {species}: {exc}"
                ) from exc
            if not gff_path.is_file():
                raise _annotation_error(
                    metadata,
                    species,
                    f"the assembly metadata could not be refreshed while locating its annotation ({exc}).",
                ) from exc

    if not genome_path.is_file():
        if not genome_url:
            raise RuntimeError(
                f"Pre-mRNA design needs the genomic FASTA for {species}. Build the matching genome index first or provide a local genome FASTA."
            )
        compressed = reference_dir / "genome.fna.gz"
        _core().notify_progress("download", f"Downloading genome FASTA for {species}…")
        references._download(genome_url, compressed, progress_callback=download_progress, session=session)
        _decompress_once(compressed, genome_path, remove_source=True)

    if not gff_path.is_file():
        if not annotation_url:
            annotation_url = next(iter(references._assembly_annotation_download_urls_from_genome_url(genome_url)), "")
        if not annotation_url:
            raise _annotation_error(
                metadata,
                species,
                "no matching genomic annotation is available for this registered reference.",
            )
        compressed = reference_dir / "annotation.gff.gz"
        _core().notify_progress("download", f"Downloading genomic annotation for {species}…")
        try:
            references._download(annotation_url, compressed, progress_callback=download_progress, session=session)
            _decompress_once(compressed, gff_path, remove_source=True)
        except _core().RunCancelled:
            raise
        except Exception as exc:
            compressed.unlink(missing_ok=True)
            raise _annotation_error(
                metadata,
                species,
                f"the matching genomic annotation could not be downloaded ({exc}).",
            ) from exc

    metadata.update({
        "reference_data_directory": str(reference_dir),
        "genome_fasta_path": str(genome_path),
        "annotation_gff_path": str(gff_path),
        "genome_url": genome_url or metadata.get("genome_url"),
        "annotation_url": annotation_url or metadata.get("annotation_url"),
    })
    metadata_path = str(metadata.get("metadata_path") or "").strip()
    if metadata_path:
        try:
            references._write_metadata(Path(metadata_path), metadata)
        except OSError:
            pass
    return genome_path, gff_path, {
        "genome_fasta": str(genome_path),
        "gff3": str(gff_path),
        "genome_url": genome_url or None,
        "annotation_url": annotation_url or None,
        "assembly": assembly,
        "assembly_accession": assembly_accession or None,
        "reference_data_directory": str(reference_dir),
    }


def resolve_reference_metadata(references: Any, species: str) -> dict[str, Any] | None:
    """Resolve assembly metadata for any preset with a ready index."""
    metadata = references.find_installed_reference(species)
    if metadata is not None:
        return metadata
    preset = references.get_species_preset(species)
    if preset is None or not references.registered_index_is_ready(preset.key):
        return None
    metadata_path = (
        references.species_data_root()
        / preset.key
        / references.core.safe_name(preset.assembly_name)
        / "reference.json"
    )
    return {
        "schema_version": 2,
        "status": "ready",
        "source": "HCRProbeDesign registered index; HCRProbeForge species preset metadata",
        "metadata_fallback": True,
        "species": preset.key,
        "index_species_alias": preset.key,
        "display_name": preset.display_name,
        "scientific_name": preset.scientific_name,
        "assembly": preset.assembly_name,
        "assembly_accession": preset.assembly_accession,
        "index_prefix": str(references.registered_index_prefix(preset.key)),
        "metadata_path": str(metadata_path),
    }


def _target_cache_paths(
    reference_dir: Path,
    record: dict[str, Any],
    selected: list[int],
    *,
    target_type: str = INTRONIC_TARGET_TYPE,
    assembly: str | None = None,
) -> tuple[Path, Path]:
    """Return project-level FASTA/metadata paths for a prepared target."""
    accession = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(record.get("accession") or "transcript"))
    selection = "all" if not selected else "intron_" + "_".join(str(item) for item in selected)
    assembly_component = _core().safe_name(assembly or reference_dir.name) or "reference"
    mode = premrna_selection_mode(target_type)
    root = _core().hcrprobeforge_cache_root() / "premrna" / "targets" / assembly_component / accession
    stem = f"{mode}_{selection}"
    return root / f"{stem}.fa", root / f"{stem}.json"


def _legacy_target_cache_paths(
    reference_dir: Path,
    record: dict[str, Any],
    selected: list[int],
) -> tuple[Path, Path]:
    accession = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(record.get("accession") or "transcript"))
    selection = "all" if not selected else "intron_" + "_".join(str(item) for item in selected)
    root = reference_dir / "premrna" / "targets"
    return root / f"{accession}_{selection}.fa", root / f"{accession}_{selection}.json"


def _write_target_cache(
    sequence_path: Path,
    metadata_path: Path,
    sequence: str,
    target: dict[str, Any],
    source_info: dict[str, Any],
) -> None:
    sequence_path.parent.mkdir(parents=True, exist_ok=True)
    temporary = sequence_path.with_name(sequence_path.name + f".tmp.{os.getpid()}")
    metadata_temporary = metadata_path.with_name(metadata_path.name + f".tmp.{os.getpid()}")
    try:
        temporary.write_text(
            ">pre-mRNA\n" + "\n".join(
                sequence[index:index + 80] for index in range(0, len(sequence), 80)
            ) + "\n",
            encoding="utf-8",
        )
        os.replace(temporary, sequence_path)
        metadata_temporary.write_text(
            json.dumps({"schema_version": 2, "target": target, "source": source_info}, indent=2, sort_keys=True),
            encoding="utf-8",
        )
        os.replace(metadata_temporary, metadata_path)
    finally:
        temporary.unlink(missing_ok=True)
        metadata_temporary.unlink(missing_ok=True)


def prepare_target(session: Any, record: dict[str, Any], args: Any) -> dict[str, Any]:
    """Prepare a full genomic transcript for intronic or whole-target design."""
    species = str(getattr(args, "species", "")).casefold()
    target_type = str(getattr(args, "target_type", INTRONIC_TARGET_TYPE) or INTRONIC_TARGET_TYPE)
    if not is_premrna_target_type(target_type):
        raise ValueError(f"Unsupported Pre-mRNA target type: {target_type}")
    whole_target = is_whole_premrna_target_type(target_type)
    mature_transcript_length = (
        record.get("mature_transcript_length_nt")
        or record.get("mature_length_nt")
        or record.get("length")
    )
    references = _references()
    preset = references.get_species_preset(species)
    metadata = resolve_reference_metadata(references, species)
    if metadata is None:
        raise RuntimeError(
            f"Pre-mRNA design requires a ready registered genome index and matching annotation for {species or 'the selected species'}. Build the species index first."
        )
    if str(getattr(args, "index", "") or "").strip():
        raise ValueError(
            "Pre-mRNA design uses the registered species assembly so genomic annotation and specificity screening cannot disagree; omit --index."
        )
    reference_species = str(metadata.get("species") or species).strip() or species
    reference_display_name = str(
        metadata.get("display_name")
        or metadata.get("scientific_name")
        or (preset.display_name if preset is not None else species)
        or "the selected species"
    ).strip()
    _core().notify_progress(
        "reference",
        f"Checking the cached genome and annotation for {reference_display_name}…",
    )
    genome_path, gff_path, source_info = _download_sources(session, metadata, reference_species, args)
    database_path = annotation_index_path(gff_path)
    database_ready = annotation_database_is_valid(database_path)
    if not database_ready:
        _core().notify_progress(
            "annotation",
            (
                f"Preparing the searchable genomic annotation for {reference_display_name}…"
                if not database_path.exists()
                else f"Refreshing the incomplete genomic annotation database for {reference_display_name}…"
            ),
        )
        try:
            build_annotation_index(
                gff_path,
                database_path,
                progress_callback=_forward_annotation_progress,
            )
        except _core().RunCancelled:
            raise
        except Exception as exc:
            raise _annotation_error(
                metadata,
                species,
                f"the genomic annotation could not be indexed ({exc}).",
            ) from exc
        metadata["annotation_database_path"] = str(database_path)
        metadata["annotation_database_status"] = "ready"
        metadata_path = str(metadata.get("metadata_path") or "").strip()
        if metadata_path:
            try:
                references._write_metadata(Path(metadata_path), metadata)
            except OSError:
                pass
    else:
        metadata["annotation_database_path"] = str(database_path)
        metadata["annotation_database_status"] = "ready"
        metadata_path = str(metadata.get("metadata_path") or "").strip()
        if metadata_path:
            try:
                references._write_metadata(Path(metadata_path), metadata)
            except OSError:
                pass
    source_info["annotation_database_path"] = str(database_path)
    source_info["annotation_database_status"] = "ready"
    _core().notify_progress("annotation", "Looking up the selected transcript in the genomic annotation…")
    models = _load_or_parse_target_models(gff_path, record, database_path=database_path)
    if not models:
        requirement = "usable transcript model" if whole_target else "usable multi-exon model"
        raise _annotation_error(
            metadata,
            species,
            f"the annotation does not contain a {requirement} for RefSeq transcript {record.get('accession') or 'the selected transcript'}.",
        )
    try:
        model = _select_model(models, record)
    except RuntimeError as exc:
        raise _annotation_error(metadata, species, str(exc)) from exc
    genomic_exons = sorted(model["exons"], key=lambda item: (int(item["start"]), int(item["end"])))
    genomic_introns: list[dict[str, int]] = []
    for left, right in zip(genomic_exons, genomic_exons[1:]):
        start = int(left["end"]) + 1
        end = int(right["start"]) - 1
        if start <= end:
            genomic_introns.append({"start": start, "end": end})
    if model.get("strand") == "-":
        genomic_introns.reverse()
    if not genomic_introns and not whole_target:
        raise _annotation_error(
            metadata,
            species,
            f"the annotation model for transcript {record.get('accession') or 'the selected transcript'} does not contain a usable multi-exon model or resolvable genomic intron interval.",
        )
    requested = parse_intron_selection(getattr(args, "premrna_introns", None))
    if whole_target and requested is not None:
        raise ValueError("--premrna-introns applies only to the intronic Pre-mRNA target type.")
    available = set(range(1, len(genomic_introns) + 1))
    if whole_target:
        selected_numbers: list[int] = []
    elif requested is not None:
        missing = [number for number in requested if number not in available]
        if missing:
            raise ValueError(
                f"Requested intron number(s) {', '.join(map(str, missing))} are not present in {record.get('accession')}; the transcript has {len(genomic_introns)} intron(s)."
            )
        selected_numbers = sorted(requested)
    else:
        selected_numbers = list(range(1, len(genomic_introns) + 1))

    ordered_segments: list[tuple[str, int | None, int, int]] = []
    exon_order = genomic_exons if model.get("strand") != "-" else list(reversed(genomic_exons))
    for index, exon in enumerate(exon_order):
        ordered_segments.append(("exon", index + 1, int(exon["start"]), int(exon["end"])))
        if index < len(genomic_introns):
            intron = genomic_introns[index]
            ordered_segments.append(("intron", index + 1, int(intron["start"]), int(intron["end"])))

    reference_dir = Path(source_info["reference_data_directory"])
    target_fasta, target_json = _target_cache_paths(
        reference_dir,
        record,
        selected_numbers if requested is not None else [],
        target_type=target_type,
        assembly=str(source_info.get("assembly") or reference_dir.name),
    )
    legacy_target_fasta, legacy_target_json = _legacy_target_cache_paths(
        reference_dir, record, selected_numbers if requested is not None else []
    )
    source_fingerprint = {
        "genome_size": genome_path.stat().st_size,
        "genome_mtime_ns": genome_path.stat().st_mtime_ns,
        "annotation_size": gff_path.stat().st_size,
        "annotation_mtime_ns": gff_path.stat().st_mtime_ns,
        "assembly_accession": source_info.get("assembly_accession"),
    }
    cached_target: dict[str, Any] | None = None
    cached_fasta = target_fasta
    cached_json = target_json
    for candidate_fasta, candidate_json in (
        (target_fasta, target_json),
        (legacy_target_fasta, legacy_target_json),
    ):
        try:
            payload = json.loads(candidate_json.read_text(encoding="utf-8"))
            target_payload = payload.get("target") if isinstance(payload.get("target"), dict) else None
            mode_matches = (
                str(target_payload.get("selection_mode") or "intronic")
                == premrna_selection_mode(target_type)
            ) if target_payload is not None else False
            if (
                payload.get("source", {}).get("fingerprint") == source_fingerprint
                and candidate_fasta.is_file()
                and mode_matches
            ):
                cached_target = target_payload
                cached_fasta = candidate_fasta
                cached_json = candidate_json
                _core().notify_progress("premrna", "Reusing the cached full Pre-mRNA target sequence…")
                break
        except (OSError, ValueError, TypeError, json.JSONDecodeError):
            continue

    regions: list[dict[str, Any]] = []
    sequence_parts: list[str] = []
    cursor = 1
    if cached_target is not None:
        sequence = "".join(
            line.strip() for line in cached_fasta.read_text(encoding="utf-8").splitlines()
            if line and not line.startswith(">")
        ).upper()
        regions = list(cached_target.get("regions") or [])
        # Do not trust an interrupted or older cache merely because its files
        # exist. Require a complete coordinate backbone before reusing it.
        region_end = max((int(region.get("end", 0)) for region in regions), default=0)
        if (
            not sequence
            or len(sequence) != int(cached_target.get("sequence_length_nt", -1))
            or region_end != len(sequence)
            or not regions
        ):
            cached_target = None
            cached_fasta = target_fasta
            cached_json = target_json
            regions = []
            sequence = ""
            _core().notify_progress("fasta", "Refreshing the incomplete cached Pre-mRNA target sequence…")
    else:
        _core().notify_progress("fasta", "Extracting the selected transcript's full genomic sequence…")
        with _open_fasta(genome_path) as genome:
            for kind, number, start, end in ordered_segments:
                sequence_piece = genome.fetch(model["seqid"], start, end)
                if model.get("strand") == "-":
                    sequence_piece = _reverse_complement(sequence_piece)
                if not sequence_piece or re.search(r"[^ACGTN]", sequence_piece):
                    raise RuntimeError(f"The genomic sequence for {kind} {number or ''} contains unsupported bases or is empty.")
                segment_start = cursor
                sequence_parts.append(sequence_piece)
                cursor += len(sequence_piece)
                region: dict[str, Any] = {
                    "kind": kind,
                    "label": f"{kind} {number}",
                    "start": segment_start,
                    "end": cursor - 1,
                    "length_nt": len(sequence_piece),
                    "seqid": model["seqid"],
                    "genomic_start": start,
                    "genomic_end": end,
                    "strand": model["strand"],
                }
                if kind == "intron":
                    region["intron_number"] = int(number or 0)
                    region["eligible"] = whole_target or int(number or 0) in selected_numbers
                else:
                    region["exon_number"] = int(number or 0)
                regions.append(region)
        sequence = "".join(sequence_parts)

    for region in regions:
        if str(region.get("kind")) == "intron":
            region["eligible"] = whole_target or int(region.get("intron_number", 0)) in selected_numbers
    if not sequence or re.search(r"[^ACGTN]", sequence):
        raise RuntimeError("The assembled Pre-mRNA sequence is empty or contains unsupported bases.")

    selection_mode = premrna_selection_mode(target_type)
    target_note = (
        "The target contains the full genomic transcript span. Probe selection includes eligible exonic and intronic sequence; candidates are annotated by region."
        if whole_target
        else "The target contains the full genomic transcript span. Probe selection is restricted to complete intronic regions; exonic and intron-boundary candidates are excluded."
    )
    target = {
        "mode": "pre-mRNA-whole" if whole_target else "pre-mRNA-intronic",
        "selection_mode": selection_mode,
        "target_representation": "full_genomic_transcript",
        "distribution": "whole_target" if whole_target else ("even_by_intron" if requested is None else "selected_introns"),
        "selected_introns": selected_numbers,
        "available_introns": len(genomic_introns),
        "transcript_accession": record.get("accession"),
        "transcript_model": {
            "seqid": model["seqid"],
            "strand": model["strand"],
            "genomic_start": min(int(exon["start"]) for exon in genomic_exons),
            "genomic_end": max(int(exon["end"]) for exon in genomic_exons),
            "exon_count": len(genomic_exons),
        },
        "mature_transcript_length_nt": int(mature_transcript_length) if mature_transcript_length else None,
        "sequence_length_nt": len(sequence),
        "regions": regions,
        "source": source_info,
        "note": target_note,
    }
    source_info["fingerprint"] = source_fingerprint
    if cached_target is None or cached_fasta != target_fasta or cached_target != target:
        _write_target_cache(target_fasta, target_json, sequence, target, source_info)
    if cached_fasta != target_fasta:
        cached_fasta.unlink(missing_ok=True)
        cached_json.unlink(missing_ok=True)
    # Remove only the old managed cache locations. The canonical project cache
    # above is retained and remains reusable by later runs.
    if legacy_target_fasta != target_fasta:
        legacy_target_fasta.unlink(missing_ok=True)
        legacy_target_json.unlink(missing_ok=True)

    record["source"] = "ncbi_premrna"
    record["_sequence"] = sequence
    record["length"] = len(sequence)
    record["title"] = f"{record.get('title') or record.get('accession')} · Pre-mRNA"
    record["premrna_target"] = target
    record["transcript_annotation"] = {
        "status": "resolved_premrna",
        "mode": "pre-mRNA-whole" if whole_target else "pre-mRNA-intronic",
        "selection_mode": selection_mode,
        "target_representation": "full_genomic_transcript",
        "accession": record.get("accession"),
        "transcript_length_nt": len(sequence),
        "coordinate_system": "1-based inclusive genomic-transcript coordinates",
        "regions": regions,
        "source": source_info,
        "note": target["note"],
    }
    return record
