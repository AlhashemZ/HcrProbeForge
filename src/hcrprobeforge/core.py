#!/usr/bin/env python3
"""HCRProbeForge: design, QC, auto-curate, plot, and export HCR probe sets."""

from __future__ import annotations

import argparse
import copy
import csv
import hashlib
import itertools
import importlib.metadata
import json
import math
import os
import platform
import re
import signal
import shutil
import subprocess
import sys
import tempfile
import time
import xml.etree.ElementTree as ET
from contextlib import contextmanager
from contextvars import ContextVar
from datetime import datetime, timezone
from pathlib import Path
from bisect import bisect_right
from typing import Any, Callable, Iterator, Tuple

try:
    import requests
except ImportError:  # Keep --help, --version, plot-only, and QC-only importable.
    requests = None  # type: ignore[assignment]

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
__version__ = "1.3.7"
SCRIPT_BUILD = __version__


class RetryableNCBIError(RuntimeError):
    """NCBI/network failure that can be retried safely by batch commands."""


class RunCancelled(RuntimeError):
    """Raised when a user cancels a running design from the local webapp."""


_NCBI_LAST_REQUEST_MONOTONIC = 0.0
_JSON_MISSING = object()
_CANCEL_CHECK: ContextVar[Callable[[], bool] | None] = ContextVar("hcrprobeforge_cancel_check", default=None)
_PROGRESS_CALLBACK: ContextVar[Callable[[dict[str, object]], None] | None] = ContextVar(
    "hcrprobeforge_progress_callback",
    default=None,
)


@contextmanager
def cancellation_scope(check: Callable[[], bool] | None) -> Iterator[None]:
    """Install a cooperative cancellation callback for the current run.

    Normal CLI calls leave this unset, so ordinary scientific execution keeps
    the same behavior. The webapp uses it to stop external design processes
    and to avoid starting the next channel or manifest job after cancellation.
    """
    token = _CANCEL_CHECK.set(check)
    try:
        yield
    finally:
        _CANCEL_CHECK.reset(token)


@contextmanager
def progress_scope(callback: Callable[[dict[str, object]], None] | None) -> Iterator[None]:
    """Install an optional progress sink for the webapp or another caller."""
    token = _PROGRESS_CALLBACK.set(callback)
    try:
        yield
    finally:
        _PROGRESS_CALLBACK.reset(token)


def notify_progress(phase: str, message: str, **details: object) -> None:
    """Report a non-scientific pipeline phase without changing the CLI path."""
    callback = _PROGRESS_CALLBACK.get()
    if callback is None:
        return
    update: dict[str, object] = {"phase": phase, "message": message}
    update.update(details)
    try:
        callback(update)
    except RunCancelled:
        raise
    except Exception:
        # Progress presentation must never alter a scientific run.
        return


def raise_if_cancelled() -> None:
    check = _CANCEL_CHECK.get()
    if check is not None and check():
        raise RunCancelled("Run cancelled by user.")


def _sleep_with_cancellation(seconds: float) -> None:
    deadline = time.monotonic() + max(0.0, seconds)
    while True:
        raise_if_cancelled()
        remaining = deadline - time.monotonic()
        if remaining <= 0:
            return
        time.sleep(min(0.1, remaining))


def _stop_design_process(process: subprocess.Popen[str], *, force: bool = False) -> None:
    """Stop designProbes and its Bowtie2 descendants during webapp cancel."""
    if os.name != "nt":
        try:
            os.killpg(process.pid, signal.SIGKILL if force else signal.SIGTERM)
            return
        except (OSError, ProcessLookupError):
            pass
    if force and os.name == "nt":
        try:
            subprocess.run(
                ["taskkill", "/F", "/T", "/PID", str(process.pid)],
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
                check=False,
            )
        except OSError:
            pass
    try:
        (process.kill if force else process.terminate)()
    except ProcessLookupError:
        pass


def safe_name(text: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", text).strip("_")


def designprobes_target_name(target: str) -> str:
    """Return a conservative target name for the external design engine.

    Human-facing Pre-mRNA folders intentionally retain a comma-separated
    intron scope (for example ``premrna_introns1,2``). ``designProbes`` is an
    external program and uses ``--targetName`` in generated identifiers and
    paths, so pass it the same readable label with punctuation normalized.
    """
    return safe_name(target) or "target"





def resolve_species_arguments(args: argparse.Namespace, argv: list[str]) -> None:
    """Resolve the biological organism from the selected genome alias.

    Known aliases are canonicalized before they reach designProbes so
    transcript retrieval and genome masking cannot silently use different
    organisms. Existing locally registered aliases remain supported when the
    user supplies their NCBI organism explicitly.
    """
    from . import references

    supplied_organism = any(
        item == "--organism" or item.startswith("--organism=") for item in argv
    )
    preset = references.get_species_preset(args.species)
    if preset is not None:
        if supplied_organism and args.organism and str(args.organism).casefold() != preset.scientific_name.casefold():
            raise ValueError(
                f"--species {args.species} is paired with NCBI organism {preset.scientific_name!r}; "
                f"received --organism {args.organism!r}. Choose a matching species or use a custom alias."
            )
        args.species = preset.key
        args.organism = preset.scientific_name
        args.reference = {
            "species": preset.key,
            "display_name": preset.display_name,
            "scientific_name": preset.scientific_name,
            "assembly": preset.assembly_name,
            "assembly_accession": preset.assembly_accession,
            "source": "NCBI Assembly",
        }
        return
    compatible = references.COMPATIBLE_ORGANISM_ALIASES.get(str(args.species).casefold())
    if compatible and not args.organism:
        args.organism = compatible
    installed = references.find_installed_reference(args.species)
    if installed is not None:
        installed_organism = str(installed.get("scientific_name") or "").strip()
        if supplied_organism and installed_organism and args.organism:
            if str(args.organism).casefold() != installed_organism.casefold():
                raise ValueError(
                    f"--species {args.species} is registered for NCBI organism {installed_organism!r}; "
                    f"received --organism {args.organism!r}. Choose a matching reference."
                )
        if installed_organism:
            args.organism = installed_organism
        args.species = str(installed.get("index_species_alias") or args.species)
        args.reference = {
            "species": installed.get("species") or args.species,
            "index_species_alias": args.species,
            "display_name": installed.get("display_name") or args.species,
            "scientific_name": args.organism,
            "assembly": installed.get("assembly"),
            "assembly_accession": installed.get("assembly_accession"),
            "source": installed.get("source") or "registered HCRProbeDesign reference",
        }
        return
    if not args.organism:
        raise ValueError(
            f"Genome alias {args.species!r} is not one of the automatic species presets. "
            "Provide --organism for a locally registered custom alias."
        )
    args.reference = {
        "species": str(args.species),
        "display_name": str(args.species),
        "scientific_name": str(args.organism),
        "assembly": None,
        "assembly_accession": None,
        "source": "user-registered HCRProbeDesign reference",
    }


def _argv_option_value(argv: list[str], option: str) -> str | None:
    """Return one simple CLI option value without reparsing a failed command."""
    for index, item in enumerate(argv):
        if item == option and index + 1 < len(argv):
            return str(argv[index + 1])
        if item.startswith(option + "="):
            return item.split("=", 1)[1]
    return None


def design_process_failure_message(return_code: int, command: list[str]) -> str:
    """Turn a designProbes failure into an actionable, honest diagnostic.

    Failed runs are rolled back by :func:`main`, so telling a user to inspect a
    log that no longer exists is misleading. A missing registered index is a
    common and recoverable failure, while all other failures retain the exit
    code and explain that no failed-run log was kept.
    """
    species = _argv_option_value(command, "--species") or ""
    if "--no-genomemask" not in command and "--index" not in command and species:
        from . import references

        if not references.registered_index_is_ready(species):
            return missing_genome_index_message(species)
    return (
        f"designProbes failed with exit code {return_code}. "
        "The failed run was rolled back, so no log file was retained."
    )


def missing_genome_index_message(species: str) -> str:
    """Return the actionable error shared by CLI and webapp design paths."""
    from . import references

    preset = references.get_species_preset(species)
    display = preset.display_name if preset is not None else species
    alias = preset.key if preset is not None else species
    if preset is not None and preset.local_genome_fasta:
        source_hint = (
            f"run 'hcrprobeforge index build --species {alias} "
            f"--fasta {preset.local_genome_fasta}'"
        )
    elif preset is not None and preset.assembly_accession:
        source_hint = f"run 'hcrprobeforge index fetch --species {alias}'"
    else:
        source_hint = (
            "provide a versioned NCBI assembly accession or local genome FASTA "
            "in the Build a genome index workflow"
        )
    return (
        f"No Bowtie2 genome index is registered for {display} (index alias {alias}). "
        "Build it first with the Build a genome index workflow or "
        f"{source_hint}, then rerun the design. "
        "The failed run was rolled back, so no design output or log file was retained."
    )


def utc_now_iso() -> str:
    """Return current UTC time as an ISO-8601 string."""
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def json_ready(value: Any) -> Any:
    """Convert common Python objects into JSON-serializable values."""
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, dict):
        return {str(key): json_ready(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_ready(item) for item in value]
    return value


def write_json(path: Path, data: dict[str, Any]) -> None:
    # Keep biological labels such as 5′ UTR readable in the files users open
    # directly.  This is only a serialization change; the parsed values and
    # scientific pipeline are unchanged.
    path.write_text(
        json.dumps(json_ready(data), indent=2, sort_keys=True, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def package_version(package: str) -> str | None:
    try:
        return importlib.metadata.version(package)
    except importlib.metadata.PackageNotFoundError:
        return None


def executable_version(executable: str) -> str | None:
    """Return the first version line reported by an optional executable."""
    resolved = shutil.which(executable)
    if not resolved:
        return None
    try:
        completed = subprocess.run(
            [resolved, "--version"],
            check=False,
            capture_output=True,
            text=True,
            timeout=5,
        )
    except (OSError, subprocess.SubprocessError):
        return None
    output = (completed.stdout or completed.stderr or "").strip()
    return output.splitlines()[0].strip() if output else None


def sanitized_args(args: argparse.Namespace) -> dict[str, Any]:
    """Return CLI arguments without exposing secrets such as NCBI_API_KEY."""
    hidden = {"api_key"}
    data: dict[str, Any] = {}
    for key, value in vars(args).items():
        if key.startswith("_"):
            continue
        if key in hidden:
            data[key] = "<provided>" if value else None
        else:
            data[key] = json_ready(value)
    return data


def sanitized_argv(argv: list[str]) -> list[str]:
    """Return command-line arguments with secret values redacted."""
    cleaned: list[str] = []
    redact_next = False
    for item in argv:
        if redact_next:
            cleaned.append("<provided>")
            redact_next = False
            continue
        if item == "--api-key":
            cleaned.append(item)
            redact_next = True
        elif item.startswith("--api-key="):
            cleaned.append("--api-key=<provided>")
        else:
            cleaned.append(item)
    return cleaned


def _ncbi_min_request_interval(params: dict[str, Any]) -> float:
    """Return a conservative per-process E-utilities request interval.

    NCBI's standard guidance allows a higher request rate when an API key is
    supplied.  The small guard here is intentionally per process; it protects
    long sequential catalogue runs without changing any scientific behavior.
    """
    return 0.11 if params.get("api_key") else 0.34


def _pace_ncbi_request(params: dict[str, Any]) -> None:
    global _NCBI_LAST_REQUEST_MONOTONIC
    interval = _ncbi_min_request_interval(params)
    now = time.monotonic()
    wait = interval - (now - _NCBI_LAST_REQUEST_MONOTONIC)
    if wait > 0:
        _sleep_with_cancellation(wait)
    _NCBI_LAST_REQUEST_MONOTONIC = time.monotonic()


def ncbi_get(session: requests.Session, endpoint: str, params: dict[str, Any]) -> requests.Response:
    """GET one NCBI E-utilities endpoint with HTTP-level retry/backoff."""
    last_error: Exception | None = None
    for attempt in range(4):
        try:
            raise_if_cancelled()
            _pace_ncbi_request(params)
            response = session.get(f"{EUTILS}/{endpoint}", params=params, timeout=45)
            response.raise_for_status()
            return response
        except requests.RequestException as exc:
            last_error = exc
            if attempt == 3:
                break
            _sleep_with_cancellation(2**attempt)
    raise RetryableNCBIError(f"NCBI request failed after retries: {last_error}")


def _json_path_value(data: dict[str, Any], path: tuple[str, ...]) -> Any:
    value: Any = data
    for key in path:
        if not isinstance(value, dict) or key not in value:
            return _JSON_MISSING
        value = value[key]
    return value


def _ncbi_payload_error_message(data: dict[str, Any]) -> str | None:
    """Extract an explicit NCBI error string from common JSON response shapes."""
    for key in ("error", "ERROR"):
        value = data.get(key)
        if value:
            return str(value)
    esearch = data.get("esearchresult")
    if isinstance(esearch, dict):
        for key in ("error", "ERROR"):
            value = esearch.get(key)
            if value:
                return str(value)
    return None


def _decode_ncbi_json(response: requests.Response, endpoint: str) -> dict[str, Any]:
    """Decode an NCBI JSON response, repairing literal control characters."""
    try:
        data = response.json()
        if not isinstance(data, dict):
            raise ValueError(
                f"NCBI returned unexpected JSON type from {endpoint}: {type(data).__name__}."
            )
        return data
    except (requests.exceptions.JSONDecodeError, json.JSONDecodeError, ValueError) as strict_exc:
        repaired = re.sub(r"[\x00-\x1f\x7f]", " ", response.text)
        try:
            data = json.loads(repaired)
            if not isinstance(data, dict):
                raise ValueError(
                    f"NCBI returned unexpected JSON type from {endpoint}: {type(data).__name__}."
                )
        except (json.JSONDecodeError, ValueError) as repaired_exc:
            raise ValueError(
                f"could not decode JSON from {endpoint}: {repaired_exc}"
            ) from strict_exc
        print(
            f"WARNING: NCBI {endpoint} JSON contained an unescaped control "
            "character; the affected metadata text was safely normalized.",
            file=sys.stderr,
        )
        return data


def ncbi_json_get(
    session: requests.Session,
    endpoint: str,
    params: dict[str, Any],
    *,
    required_path: tuple[str, ...] | None = None,
) -> dict[str, Any]:
    """Return validated JSON from NCBI with semantic retry/backoff.

    In addition to malformed JSON, NCBI can occasionally return a syntactically
    valid service/error payload that lacks the fields expected by ESearch or
    ESummary.  Such a response must not be treated as a successful lookup.  The
    payload is therefore validated after control-character repair and retried
    when a required response path is absent.
    """
    last_error: Exception | None = None
    for attempt in range(4):
        try:
            response = ncbi_get(session, endpoint, params)
            data = _decode_ncbi_json(response, endpoint)
            reported_error = _ncbi_payload_error_message(data)
            if reported_error:
                raise ValueError(f"NCBI reported an error: {reported_error}")
            if required_path is not None:
                value = _json_path_value(data, required_path)
                if value is _JSON_MISSING:
                    raise ValueError(
                        "NCBI returned an incomplete JSON response; missing "
                        + ".".join(required_path)
                    )
            return data
        except (RetryableNCBIError, ValueError, json.JSONDecodeError) as exc:
            last_error = exc
            if attempt == 3:
                break
            delay = 2**attempt
            print(
                f"WARNING: NCBI {endpoint} response was incomplete or transient "
                f"({exc}); retrying in {delay}s ({attempt + 1}/4).",
                file=sys.stderr,
            )
            _sleep_with_cancellation(delay)

    raise RetryableNCBIError(
        f"NCBI {endpoint} did not return a complete valid response after retries: {last_error}"
    )

def base_params(args: argparse.Namespace) -> dict[str, str]:
    params = {"tool": "hcrprobeforge"}
    if args.email:
        params["email"] = args.email
    if args.api_key:
        params["api_key"] = args.api_key
    return params


def search_ids(
    session: requests.Session,
    db: str,
    term: str,
    args: argparse.Namespace,
    retmax: int = 50,
) -> list[str]:
    params = base_params(args) | {
        "db": db,
        "term": term,
        "retmode": "json",
        "retmax": str(retmax),
    }
    data = ncbi_json_get(
        session, "esearch.fcgi", params, required_path=("esearchresult", "idlist")
    )
    return list(data["esearchresult"]["idlist"])


def gene_summaries(
    session: requests.Session, ids: list[str], args: argparse.Namespace
) -> list[dict[str, Any]]:
    """Return compact NCBI Gene summaries for candidate Gene IDs."""
    if not ids:
        return []
    params = base_params(args) | {
        "db": "gene",
        "id": ",".join(ids),
        "retmode": "json",
    }
    data = ncbi_json_get(session, "esummary.fcgi", params, required_path=("result", "uids"))["result"]
    records: list[dict[str, Any]] = []
    for uid in data.get("uids", []):
        item = data.get(uid, {})
        if not isinstance(item, dict) or not item:
            raise RetryableNCBIError(f"NCBI ESummary result was incomplete for Gene UID {uid}.")
        organism = item.get("organism", {})
        if isinstance(organism, dict):
            scientific_name = organism.get("scientificname", "")
        else:
            scientific_name = str(organism or "")
        records.append(
            {
                "gene_id": str(uid),
                "symbol": str(item.get("name") or item.get("nomenclaturesymbol") or ""),
                "nomenclature_symbol": str(item.get("nomenclaturesymbol") or ""),
                "description": str(item.get("description") or ""),
                "organism": scientific_name,
                "status": str(item.get("status") or ""),
                "other_aliases": str(item.get("otheraliases") or ""),
            }
        )
    return records


def _prompt_number(title: str, options: list[str], default_index: int = 0) -> int | None:
    """Prompt for a 0-based option index when stdin is interactive."""
    if not sys.stdin.isatty():
        return None
    print(title)
    for index, option in enumerate(options, start=1):
        marker = " (default)" if index - 1 == default_index else ""
        print(f"  {index}. {option}{marker}")
    while True:
        response = input(f"Choose [default {default_index + 1}]: ").strip()
        if not response:
            return default_index
        try:
            chosen = int(response) - 1
        except ValueError:
            print("Enter a number from the list.")
            continue
        if 0 <= chosen < len(options):
            return chosen
        print("Enter a number from the list.")


def resolve_gene_record(
    session: requests.Session, gene: str, organism: str, args: argparse.Namespace
) -> tuple[str, dict[str, Any]]:
    """Resolve a gene symbol to one NCBI Gene record without requiring a rerun.

    Exact official-symbol matches are preferred. Remaining ties are ranked by
    active status and availability of linked RefSeq RNA records. An interactive
    terminal is offered an in-place choice; non-interactive execution uses a
    deterministic GeneID tie-break and records the ambiguity in metadata.
    """
    ids = search_ids(session, "gene", f'{gene}[sym] AND "{organism}"[orgn]', args)
    search_mode = "symbol"
    if not ids:
        ids = search_ids(session, "gene", f'"{gene}"[All Fields] AND "{organism}"[orgn]', args)
        search_mode = "all_fields"
    if not ids:
        raise RuntimeError(
            f"No NCBI Gene match was found for {gene!r} in {organism!r}. "
            "For species with limited annotation, retry with the gene's NCBI/RefSeq accession or provide a FASTA file."
        )

    records = gene_summaries(session, ids, args)
    if not records:
        raise RuntimeError(
            f"NCBI returned Gene IDs for {gene!r}, but no summaries could be read. "
            "Retry with the gene's NCBI/RefSeq accession or provide a FASTA file."
        )

    query_fold = gene.casefold()
    for record in records:
        symbols = {record.get("symbol", "").casefold(), record.get("nomenclature_symbol", "").casefold()}
        record["exact_official_symbol"] = query_fold in symbols
        record["is_discontinued"] = "discontinued" in record.get("status", "").casefold()
        record["linked_refseq_rna_count"] = None

    # Only query RefSeq links for the best symbol/status group. This avoids one
    # additional NCBI request for every broad all-fields match.
    best_symbol_status = min(
        (not bool(record["exact_official_symbol"]), bool(record["is_discontinued"]))
        for record in records
    )
    link_rank_candidates = [
        record
        for record in records
        if (not bool(record["exact_official_symbol"]), bool(record["is_discontinued"])) == best_symbol_status
    ]
    for record in link_rank_candidates:
        record["linked_refseq_rna_count"] = len(linked_refseq_rna_ids(session, record["gene_id"], args))

    ranked = sorted(
        records,
        key=lambda record: (
            not bool(record["exact_official_symbol"]),
            bool(record["is_discontinued"]),
            -int(record["linked_refseq_rna_count"] or 0),
            int(record["gene_id"]),
        ),
    )

    primary = lambda record: (
        not bool(record["exact_official_symbol"]),
        bool(record["is_discontinued"]),
        -int(record["linked_refseq_rna_count"] or 0),
    )
    tied = [record for record in ranked if primary(record) == primary(ranked[0])]
    ambiguity_note = None

    if args.transcript_policy == "require-accession" and len(records) != 1:
        candidates = ", ".join(f"{r['symbol'] or '?'} (GeneID {r['gene_id']})" for r in ranked)
        raise RuntimeError(
            f"Multiple NCBI Gene matches were found for {gene!r}: {candidates}. "
            "Supply --accession or use the default automatic transcript policy."
        )

    selected = ranked[0]
    if len(tied) > 1 and args.transcript_policy in {"auto", "interactive"}:
        choice = _prompt_number(
            f"Multiple plausible NCBI Gene records match {gene!r} in {organism!r}:",
            [
                f"{r['symbol'] or '?'}; GeneID {r['gene_id']}; {r['description']}; "
                f"RefSeq RNAs={r['linked_refseq_rna_count']}"
                for r in tied
            ],
        )
        if choice is not None:
            selected = tied[choice]
        else:
            ambiguity_note = (
                "Several Gene records remained tied after symbol/status/RefSeq ranking in non-interactive mode; "
                f"the lowest deterministic GeneID ({selected['gene_id']}) was selected."
            )
            print(f"WARNING: {ambiguity_note}", file=sys.stderr)
    elif len(tied) > 1:
        ambiguity_note = (
            "Several Gene records remained tied after symbol/status/RefSeq ranking; "
            f"the lowest deterministic GeneID ({selected['gene_id']}) was selected."
        )
        print(f"WARNING: {ambiguity_note}", file=sys.stderr)

    reason_parts = []
    if selected["exact_official_symbol"]:
        reason_parts.append("exact official-symbol match")
    else:
        reason_parts.append(f"best {search_mode.replace('_', ' ')} match")
    if not selected["is_discontinued"]:
        reason_parts.append("active Gene record")
    reason_parts.append(f"{selected['linked_refseq_rna_count']} linked RefSeq RNA record(s)")

    resolution = {
        "query_gene_symbol": gene,
        "organism": organism,
        "search_mode": search_mode,
        "candidate_gene_count": len(records),
        "selected_gene_id": selected["gene_id"],
        "selected_gene_symbol": selected.get("symbol"),
        "selected_gene_description": selected.get("description"),
        "selection_reason": "; ".join(reason_parts),
        "ambiguity_note": ambiguity_note,
        "gene_candidates": [
            {
                "gene_id": r["gene_id"],
                "symbol": r.get("symbol"),
                "description": r.get("description"),
                "status": r.get("status"),
                "exact_official_symbol": bool(r.get("exact_official_symbol")),
                "linked_refseq_rna_count": int(r.get("linked_refseq_rna_count") or 0),
            }
            for r in ranked
        ],
    }
    return selected["gene_id"], resolution


def linked_refseq_rna_ids(
    session: requests.Session, gene_id: str, args: argparse.Namespace
) -> list[str]:
    params = base_params(args) | {
        "dbfrom": "gene",
        "db": "nuccore",
        "id": gene_id,
        "linkname": "gene_nuccore_refseqrna",
        "retmode": "xml",
    }
    last_error: Exception | None = None
    root: ET.Element | None = None
    for attempt in range(4):
        try:
            root = ET.fromstring(ncbi_get(session, "elink.fcgi", params).text)
            break
        except (RetryableNCBIError, ET.ParseError) as exc:
            last_error = exc
            if attempt == 3:
                break
            delay = 2**attempt
            print(
                f"WARNING: NCBI elink.fcgi response was incomplete or transient ({exc}); "
                f"retrying in {delay}s ({attempt + 1}/4).",
                file=sys.stderr,
            )
            time.sleep(delay)
    if root is None:
        raise RetryableNCBIError(
            f"NCBI elink.fcgi did not return valid XML after retries: {last_error}"
        )
    ids: list[str] = []
    for block in root.findall(".//LinkSetDb"):
        name = block.findtext("LinkName", default="")
        if name == "gene_nuccore_refseqrna":
            ids.extend(node.text for node in block.findall("./Link/Id") if node.text)
    return ids


def nucleotide_summaries(
    session: requests.Session, ids: list[str], args: argparse.Namespace
) -> list[dict[str, Any]]:
    if not ids:
        return []
    params = base_params(args) | {
        "db": "nuccore",
        "id": ",".join(ids),
        "retmode": "json",
    }
    data = ncbi_json_get(session, "esummary.fcgi", params, required_path=("result", "uids"))["result"]
    records = []
    for uid in data.get("uids", []):
        item = data.get(uid, {})
        if not isinstance(item, dict) or not item:
            raise RetryableNCBIError(f"NCBI ESummary result was incomplete for nuccore UID {uid}.")
        accession = item.get("accessionversion") or item.get("caption")
        if accession:
            records.append(
                {
                    "uid": uid,
                    "accession": accession,
                    "title": item.get("title", ""),
                    "length": int(item.get("slen", 0)),
                }
            )
    return records


def refseq_select_accessions(
    session: requests.Session, gene: str, organism: str, args: argparse.Namespace
) -> tuple[set[str], dict[str, Any]]:
    """Best-effort lookup of RefSeq Select records for one gene.

    RefSeq Select is not currently available for every organism. A missing or
    unavailable Select record is therefore metadata, not a fatal lookup error.
    """
    query = f'"{gene}"[gene] AND "{organism}"[orgn] AND Refseq_select[filter]'
    try:
        ids = search_ids(session, "nuccore", query, args, retmax=50)
        records = nucleotide_summaries(session, ids, args)
    except (RuntimeError, ValueError, KeyError) as exc:
        return set(), {
            "status": "lookup_failed",
            "query": query,
            "accessions": [],
            "note": str(exc),
        }
    accessions = {str(record["accession"]) for record in records if record.get("accession")}
    return accessions, {
        "status": "found" if accessions else "not_available_for_gene_or_organism",
        "query": query,
        "accessions": sorted(accessions),
        "note": None,
    }


def transcript_prefix_rank(accession: str) -> int:
    """Prefer curated coding/noncoding RefSeq before predicted records."""
    for rank, prefix in enumerate(("NM_", "NR_", "XM_", "XR_")):
        if accession.startswith(prefix):
            return rank
    return 99


def transcript_is_partial(record: dict[str, Any]) -> bool:
    return "partial" in str(record.get("title", "")).casefold()


def automatic_transcript_primary_rank(record: dict[str, Any]) -> tuple[Any, ...]:
    """Biological/record-quality rank before the deterministic accession tie-break."""
    return (
        not bool(record.get("is_refseq_select", False)),
        transcript_prefix_rank(str(record.get("accession", ""))),
        transcript_is_partial(record),
        -int(record.get("length", 0)),
    )


def automatic_transcript_rank(record: dict[str, Any]) -> tuple[Any, ...]:
    return automatic_transcript_primary_rank(record) + (str(record.get("accession", "")),)


def _transcript_summary(record: dict[str, Any]) -> dict[str, Any]:
    return {
        "accession": record.get("accession"),
        "length_nt": int(record.get("length", 0)),
        "title": record.get("title", ""),
        "refseq_class": str(record.get("accession", ""))[:3],
        "partial": transcript_is_partial(record),
        "refseq_select": bool(record.get("is_refseq_select", False)),
    }


def hcrprobeforge_cache_root() -> Path:
    """Return the root used by transcript-sequence and annotation caches.

    ``HCRPROBEFORGE_CACHE_DIR`` is authoritative when set. Otherwise, direct
    runs use the selected project-level ``<project>/cache/<species>/<workflow>``
    tree.
    Batch runners set the same location explicitly so all front ends share one
    cache layout.
    """
    override = os.getenv("HCRPROBEFORGE_CACHE_DIR")
    if override:
        return Path(override).expanduser()
    species = safe_name(os.getenv("HCRPROBEFORGE_SPECIES", "default")) or "default"
    workflow = os.getenv("HCRPROBEFORGE_WORKFLOW", "design")
    # The short names remain useful internally, but cache folders use the same
    # readable names as run folders.  Import lazily to avoid the core/reference
    # module import cycle.
    try:
        from .references import species_component, workflow_directory

        species = species_component(species)
        workflow = workflow_directory(workflow)
    except ImportError:  # pragma: no cover - only relevant during bootstrap
        workflow = safe_name(workflow) or "workflow"
    project_root = os.getenv("HCRPROBEFORGE_PROJECT_ROOT")
    if project_root:
        return Path(project_root).expanduser().resolve() / "cache" / species / workflow
    return Path.cwd() / "hcr_results" / "cache" / species / workflow


def transcript_record_cache_paths(accession: str) -> tuple[Path, Path]:
    """Return metadata and FASTA cache paths for one exact accession-version."""
    root = hcrprobeforge_cache_root() / "transcripts" / safe_name(accession)
    return root / "record.json", root / "sequence.fa"


def _parse_fasta_text(
    text: str,
    *,
    expected_accession: str | None = None,
) -> tuple[str, str]:
    """Validate FASTA text and return (header_without_>, normalized_sequence)."""
    lines = [line.strip() for line in text.splitlines() if line.strip()]
    if not lines or not lines[0].startswith(">"):
        raise ValueError("response does not begin with a FASTA header")
    header = lines[0][1:].strip()
    if expected_accession:
        accession_pattern = rf"(?<![A-Za-z0-9_.]){re.escape(expected_accession)}(?![A-Za-z0-9_.])"
        if re.search(accession_pattern, header) is None:
            raise ValueError(
                f"FASTA header does not contain requested accession {expected_accession}: {header!r}"
            )
    if any(line.startswith(">") for line in lines[1:]):
        raise ValueError("FASTA response unexpectedly contains multiple records")
    sequence = "".join(lines[1:])
    sequence = sequence.upper().replace("U", "T")
    if not sequence or re.search(r"[^ACGTN]", sequence):
        raise ValueError("FASTA contains an empty or unsupported nucleotide sequence")
    return header, sequence


def _title_from_fasta_header(header: str, accession: str) -> str:
    title = str(header or "").strip()
    if title.startswith(accession):
        title = title[len(accession):].lstrip(" |")
    return title or str(header or accession)


def _organism_text_matches(text: str, expected: str) -> bool:
    """Return whether an NCBI title/organism field names *expected* organism."""
    expected_words = re.findall(r"[a-z0-9]+", str(expected).casefold())
    observed_words = re.findall(r"[a-z0-9]+", str(text).casefold())
    if not expected_words or not observed_words:
        return False
    # NCBI FASTA titles normally start with the binomial name, while some
    # records put it in square brackets.  Requiring the complete ordered word
    # sequence avoids accepting a related species with one shared word.
    width = len(expected_words)
    return any(observed_words[i : i + width] == expected_words for i in range(len(observed_words) - width + 1))


def _contains_known_alternate_organism(text: str, expected: str) -> bool:
    """Return whether *text* explicitly names another built-in species.

    NCBI titles are not required to contain a taxon name, so the absence of
    the selected name is not itself an error.  If a title does name one of
    HCRProbeForge's other public presets, however, accepting it would allow a
    stale cross-species cache entry to pass just because a separate metadata
    field still contains the selected organism.
    """
    from . import references

    expected_key = str(expected).casefold().strip()
    for preset in references.supported_species():
        if preset.scientific_name.casefold() == expected_key:
            continue
        if _organism_text_matches(text, preset.scientific_name):
            return True
    return False


def validate_transcript_organism(record: dict[str, Any], args: argparse.Namespace) -> None:
    """Reject an exact transcript whose NCBI organism disagrees with the run.

    Explicit accessions are otherwise valid identifiers on their own, so a
    manifest can accidentally combine an accession from one species with a
    different selected genome.  NCBI's organism/title metadata is the source
    of truth and is checked before a FASTA is written, cached, or designed.
    """
    expected = str(getattr(args, "organism", "") or "").strip()
    if not expected or record.get("source") == "user_fasta":
        return
    accession = str(record.get("accession") or "unknown")
    organism_text = str(record.get("organism") or "").strip()
    title_text = str(record.get("title") or "").strip()
    if organism_text and not _organism_text_matches(organism_text, expected):
        title = title_text or "the NCBI record"
        raise RuntimeError(
            f"Transcript accession {accession} does not belong to the selected organism "
            f"{expected!r}. NCBI reported: {title!r}. Choose an accession from the selected species."
        )
    if title_text and _contains_known_alternate_organism(title_text, expected):
        title = str(record.get("title") or "the NCBI record")
        raise RuntimeError(
            f"Transcript accession {accession} does not belong to the selected organism "
            f"{expected!r}. NCBI reported: {title!r}. Choose an accession from the selected species."
        )
    # Store the validated organism in newly written cache metadata.  Existing
    # caches without this field remain supported because title matching above
    # is also used for them.
    record["organism"] = expected


def _update_resolution_after_sequence(record: dict[str, Any]) -> None:
    resolution = record.get("transcript_resolution")
    if not isinstance(resolution, dict):
        return
    resolution["selected_accession"] = record.get("accession")
    candidates = resolution.get("transcript_candidates")
    if isinstance(candidates, list):
        summary = _transcript_summary(record)
        replaced = False
        for index, candidate in enumerate(candidates):
            if isinstance(candidate, dict) and candidate.get("accession") == record.get("accession"):
                candidates[index] = summary
                replaced = True
        if not replaced and len(candidates) <= 1:
            resolution["transcript_candidates"] = [summary]


def _set_transcript_cache_metadata(
    record: dict[str, Any],
    status: str,
    record_json: Path | None,
    fasta_path: Path | None,
) -> None:
    record["transcript_cache_status"] = status
    resolution = record.get("transcript_resolution")
    if isinstance(resolution, dict):
        resolution["transcript_cache"] = {
            "status": status,
            "record_json": str(record_json) if record_json else None,
            "fasta": str(fasta_path) if fasta_path else None,
        }


def load_cached_transcript_record(accession: str) -> dict[str, Any] | None:
    """Load and integrity-check a cached exact transcript record and FASTA."""
    if not re.fullmatch(r"(?:NM|NR|XM|XR)_\d+\.\d+", str(accession)):
        return None
    record_json, fasta_path = transcript_record_cache_paths(accession)
    if not record_json.exists() or not fasta_path.exists():
        return None
    try:
        metadata = json.loads(record_json.read_text())
        if not isinstance(metadata, dict) or metadata.get("accession") != accession:
            raise ValueError("cached accession metadata does not match the requested accession")
        fasta_text = fasta_path.read_text()
        header, sequence = _parse_fasta_text(fasta_text, expected_accession=accession)
        cached_length = int(metadata.get("length_nt", -1))
        if cached_length != len(sequence):
            raise ValueError(
                f"cached length mismatch: metadata={cached_length}, FASTA={len(sequence)}"
            )
        expected_sha = str(metadata.get("fasta_sha256") or "")
        actual_sha = hashlib.sha256(fasta_text.encode()).hexdigest()
        if expected_sha and expected_sha != actual_sha:
            raise ValueError("cached FASTA checksum does not match record metadata")
        return {
            "uid": metadata.get("uid"),
            "accession": accession,
            "title": str(metadata.get("title") or _title_from_fasta_header(header, accession)),
            "length": len(sequence),
            "source": "ncbi_transcript_cache",
            "_cached_fasta_path": str(fasta_path),
            "_transcript_cache_record_path": str(record_json),
            "transcript_cache_status": "hit",
            "organism": metadata.get("organism"),
        }
    except (OSError, ValueError, TypeError, json.JSONDecodeError) as exc:
        print(
            f"WARNING: ignoring invalid transcript cache for {accession}: {exc}",
            file=sys.stderr,
        )
        return None


def write_transcript_cache(
    record: dict[str, Any],
    fasta_text: str,
) -> tuple[Path, Path] | None:
    """Atomically cache one exact NCBI transcript FASTA plus compact metadata."""
    accession = str(record.get("accession") or "")
    if not re.fullmatch(r"(?:NM|NR|XM|XR)_\d+\.\d+", accession):
        return None
    try:
        header, sequence = _parse_fasta_text(fasta_text, expected_accession=accession)
        record_json, fasta_path = transcript_record_cache_paths(accession)
        record_json.parent.mkdir(parents=True, exist_ok=True)
        metadata = {
            "schema_version": 1,
            "accession": accession,
            "uid": record.get("uid"),
            "title": str(record.get("title") or _title_from_fasta_header(header, accession)),
            "organism": record.get("organism"),
            "length_nt": len(sequence),
            "fasta_sha256": hashlib.sha256(fasta_text.encode()).hexdigest(),
            "cached_at_utc": utc_now_iso(),
            "source": "NCBI EFetch",
        }
        fasta_tmp = fasta_path.with_name(fasta_path.name + f".tmp.{os.getpid()}")
        json_tmp = record_json.with_name(record_json.name + f".tmp.{os.getpid()}")
        fasta_tmp.write_text(fasta_text)
        json_tmp.write_text(
            json.dumps(metadata, indent=2, sort_keys=True, ensure_ascii=False) + "\n",
            encoding="utf-8",
        )
        os.replace(fasta_tmp, fasta_path)
        os.replace(json_tmp, record_json)
        return record_json, fasta_path
    except (OSError, ValueError, TypeError) as exc:
        print(
            f"WARNING: transcript FASTA for {accession} could not be cached: {exc}",
            file=sys.stderr,
        )
        return None


def fetch_ncbi_fasta_text(
    session: requests.Session,
    record: dict[str, Any],
    args: argparse.Namespace,
) -> tuple[str, str, str]:
    """Fetch and semantically validate an exact transcript FASTA with retries."""
    accession = str(record.get("accession") or "")
    identifier = accession or str(record.get("uid") or "")
    params = base_params(args) | {
        "db": "nuccore",
        "id": identifier,
        "rettype": "fasta",
        "retmode": "text",
    }
    last_error: Exception | None = None
    for attempt in range(4):
        try:
            response = ncbi_get(session, "efetch.fcgi", params)
            text = response.text.strip() + "\n"
            header, sequence = _parse_fasta_text(
                text,
                expected_accession=accession if accession else None,
            )
            return text, header, sequence
        except (RetryableNCBIError, ValueError) as exc:
            last_error = exc
            if attempt == 3:
                break
            delay = 2**attempt
            print(
                f"WARNING: NCBI FASTA response for {accession or identifier} was incomplete "
                f"or transient ({exc}); retrying in {delay}s ({attempt + 1}/4).",
                file=sys.stderr,
            )
            time.sleep(delay)
    raise RetryableNCBIError(
        f"NCBI did not return a valid FASTA for {accession or identifier} after retries: {last_error}"
    )


def resolve_transcripts(
    session: requests.Session, args: argparse.Namespace
) -> list[dict[str, Any]]:
    """Resolve explicit accessions or automatically choose a representative transcript.

    Exact versioned accessions use the shared transcript cache when available.
    On a cache miss, the exact accession is carried directly to EFetch rather
    than performing a redundant ESearch/ESummary round trip.  Versionless
    accessions still require NCBI resolution to determine the current version.
    """
    if args.accession:
        requested = str(args.accession)
        if "." in requested:
            cached = load_cached_transcript_record(requested)
            if cached is not None:
                record = cached
                selection_reason = (
                    "User supplied an exact versioned --accession; the verified local transcript "
                    "cache was reused and automatic transcript selection was bypassed."
                )
                record["transcript_resolution"] = {
                    "mode": "explicit_accession",
                    "query_gene_symbol": args.gene,
                    "organism": args.organism,
                    "requested_accession": requested,
                    "selected_accession": requested,
                    "transcript_selection_policy": "explicit_accession",
                    "selection_reason": selection_reason,
                    "transcript_candidates": [_transcript_summary(record)],
                }
                record_json = Path(str(record.get("_transcript_cache_record_path")))
                fasta_path = Path(str(record.get("_cached_fasta_path")))
                _set_transcript_cache_metadata(record, "hit", record_json, fasta_path)
                return [record]

            # The accession version is already authoritative.  Avoid a redundant
            # ESearch/ESummary lookup and validate it when the exact FASTA is
            # fetched below.  This is especially important for B2-B5 catalogue
            # runs, which normally hit the cache populated by the first channel.
            record = {
                "uid": None,
                "accession": requested,
                "title": "",
                "length": 0,
                "source": "ncbi_exact_accession",
            }
            selection_reason = (
                "User supplied an exact versioned --accession; automatic transcript selection "
                "was bypassed and the exact accession will be validated by NCBI EFetch."
            )
            record["transcript_resolution"] = {
                "mode": "explicit_accession",
                "query_gene_symbol": args.gene,
                "organism": args.organism,
                "requested_accession": requested,
                "selected_accession": requested,
                "transcript_selection_policy": "explicit_accession",
                "selection_reason": selection_reason,
                "transcript_candidates": [_transcript_summary(record)],
            }
            record_json, fasta_path = transcript_record_cache_paths(requested)
            _set_transcript_cache_metadata(record, "miss", record_json, fasta_path)
            return [record]

        # Versionless accessions still need a lookup to discover the current
        # NCBI accession-version before sequence caching can be authoritative.
        ids = search_ids(session, "nuccore", f'{requested}[accn]', args, retmax=5)
        records = nucleotide_summaries(session, ids, args)
        matches = [
            record
            for record in records
            if str(record["accession"]).split(".", 1)[0] == requested
        ]
        matches.sort(
            key=lambda record: -int(str(record["accession"]).rsplit(".", 1)[1])
            if "." in str(record["accession"]) and str(record["accession"]).rsplit(".", 1)[1].isdigit()
            else 0
        )
        if not matches:
            raise RuntimeError(f"Could not resolve accession {args.accession!r} at NCBI.")
        record = matches[0]
        selection_reason = (
            f"User supplied versionless accession {requested}; NCBI resolved the current version "
            f"{record['accession']}."
        )
        record["transcript_resolution"] = {
            "mode": "explicit_accession",
            "query_gene_symbol": args.gene,
            "organism": args.organism,
            "requested_accession": requested,
            "selected_accession": record["accession"],
            "transcript_selection_policy": "explicit_accession",
            "selection_reason": selection_reason,
            "transcript_candidates": [_transcript_summary(record)],
        }
        record_json, fasta_path = transcript_record_cache_paths(record["accession"])
        _set_transcript_cache_metadata(record, "miss", record_json, fasta_path)
        return [record]

    gene_id, gene_resolution = resolve_gene_record(session, args.gene, args.organism, args)
    records = nucleotide_summaries(session, linked_refseq_rna_ids(session, gene_id, args), args)
    refseq_rnas = [
        record
        for record in records
        if str(record.get("accession", "")).startswith(("NM_", "XM_", "NR_", "XR_"))
    ]
    if not refseq_rnas:
        raise RuntimeError(
            f"No linked RefSeq RNA was found for {args.gene} in {args.organism}. "
            "For species with limited annotation, retry with the gene's NCBI/RefSeq accession or provide a FASTA file."
        )

    select_lookup = {
        "status": "not_checked",
        "query": None,
        "accessions": [],
        "note": "RefSeq Select lookup is only needed for automatic single-transcript selection with multiple candidates.",
    }
    select_accessions: set[str] = set()
    if len(refseq_rnas) > 1 and args.transcript_policy == "auto" and not args.all_transcripts:
        select_accessions, select_lookup = refseq_select_accessions(
            session,
            gene_resolution.get("selected_gene_symbol") or args.gene,
            args.organism,
            args,
        )
    for record in refseq_rnas:
        record["is_refseq_select"] = str(record.get("accession", "")) in select_accessions

    automatic_ranked = sorted(refseq_rnas, key=automatic_transcript_rank)
    transcript_ambiguity_note = None
    if args.all_transcripts:
        selected_records = automatic_ranked
        selection_reason = "--all-transcripts requested every linked RefSeq RNA, in automatic preference order."
    elif args.transcript_policy == "require-accession" and len(refseq_rnas) > 1:
        accessions = ", ".join(record["accession"] for record in automatic_ranked)
        raise RuntimeError(
            f"Multiple RefSeq transcripts are linked to {args.gene}: {accessions}. "
            "Supply --accession or use --transcript-policy auto."
        )
    elif args.transcript_policy == "longest":
        selected_records = [
            sorted(
                refseq_rnas,
                key=lambda record: (-int(record.get("length", 0)), transcript_prefix_rank(record["accession"]), record["accession"]),
            )[0]
        ]
        selection_reason = "Longest linked RefSeq RNA; accession class used only as a tie-breaker."
    elif args.transcript_policy == "interactive" and len(refseq_rnas) > 1:
        choice = _prompt_number(
            f"Multiple RefSeq transcripts are linked to {args.gene} (GeneID {gene_id}):",
            [
                f"{record['accession']} | {record['length']} nt | "
                f"{'partial' if transcript_is_partial(record) else 'complete/unspecified'} | {record['title']}"
                for record in automatic_ranked
            ],
        )
        if choice is None:
            selected_records = [automatic_ranked[0]]
            selection_reason = (
                "Interactive selection was requested, but stdin was non-interactive; "
                "the automatic deterministic policy was used."
            )
        else:
            selected_records = [automatic_ranked[choice]]
            selection_reason = "Selected interactively in the current invocation."
    else:
        selected = automatic_ranked[0]
        tied = [
            record
            for record in automatic_ranked
            if automatic_transcript_primary_rank(record) == automatic_transcript_primary_rank(selected)
        ]
        if len(tied) > 1:
            choice = _prompt_number(
                f"Multiple RefSeq transcripts remain tied after automatic ranking for {args.gene} (GeneID {gene_id}):",
                [
                    f"{record['accession']} | {record['length']} nt | {record['title']}"
                    for record in tied
                ],
            )
            if choice is not None:
                selected = tied[choice]
                transcript_ambiguity_note = "A genuinely tied automatic transcript choice was resolved interactively."
            else:
                transcript_ambiguity_note = (
                    "Several transcripts remained tied after RefSeq Select/class/completeness/length ranking in "
                    f"non-interactive mode; {selected['accession']} was selected by deterministic accession order."
                )
                print(f"WARNING: {transcript_ambiguity_note}", file=sys.stderr)
        selected_records = [selected]
        class_label = {
            "NM_": "curated protein-coding RefSeq",
            "NR_": "curated noncoding RefSeq",
            "XM_": "predicted protein-coding RefSeq",
            "XR_": "predicted noncoding RefSeq",
        }.get(selected["accession"][:3], "RefSeq RNA")
        completeness = "non-partial record" if not transcript_is_partial(selected) else "partial record"
        select_reason = "NCBI RefSeq Select; " if selected.get("is_refseq_select") else ""
        selection_reason = (
            f"Automatic deterministic policy: {select_reason}preferred {class_label}, then {completeness}, "
            "then the longest transcript within that class."
        )

    candidates = [_transcript_summary(record) for record in automatic_ranked]
    for record in selected_records:
        record["gene_id"] = gene_id
        record["resolved_gene_symbol"] = gene_resolution.get("selected_gene_symbol")
        record["transcript_resolution"] = {
            "mode": "gene_symbol_lookup",
            "query_gene_symbol": args.gene,
            "organism": args.organism,
            "gene_resolution": gene_resolution,
            "selected_gene_id": gene_id,
            "selected_accession": record["accession"],
            "transcript_selection_policy": "all" if args.all_transcripts else args.transcript_policy,
            "selection_reason": selection_reason,
            "transcript_ambiguity_note": transcript_ambiguity_note,
            "refseq_select_lookup": select_lookup,
            "transcript_candidate_count": len(candidates),
            "transcript_candidates": candidates,
        }
    return selected_records


def fetch_fasta(
    session: requests.Session,
    record: dict[str, Any],
    output: Path,
    args: argparse.Namespace,
) -> int:
    """Write the exact transcript FASTA, preferring the verified local cache."""
    accession = str(record.get("accession") or "")
    cached_path_text = str(record.get("_cached_fasta_path") or "")
    if cached_path_text:
        cached_path = Path(cached_path_text)
        try:
            fasta_text = cached_path.read_text()
            header, sequence = _parse_fasta_text(
                fasta_text,
                expected_accession=accession if accession else None,
            )
            expected_length = int(record.get("length", 0) or 0)
            if expected_length and expected_length != len(sequence):
                raise ValueError(
                    f"cached FASTA length {len(sequence)} does not match record length {expected_length}"
                )
            record["length"] = len(sequence)
            if not record.get("title"):
                record["title"] = _title_from_fasta_header(header, accession)
            validate_transcript_organism(record, args)
            output.write_text(fasta_text)
            record_json = Path(str(record.get("_transcript_cache_record_path") or ""))
            _set_transcript_cache_metadata(record, "hit", record_json, cached_path)
            _update_resolution_after_sequence(record)
            return len(sequence)
        except (OSError, ValueError) as exc:
            print(
                f"WARNING: cached transcript FASTA for {accession} could not be used; "
                f"falling back to NCBI: {exc}",
                file=sys.stderr,
            )
            record.pop("_cached_fasta_path", None)

    fasta_text, header, sequence = fetch_ncbi_fasta_text(session, record, args)
    previous_length = int(record.get("length", 0) or 0)
    if previous_length and previous_length != len(sequence):
        raise RuntimeError(
            f"NCBI FASTA length mismatch for {accession}: summary={previous_length}, "
            f"FASTA={len(sequence)}."
        )
    record["length"] = len(sequence)
    if not record.get("title"):
        record["title"] = _title_from_fasta_header(header, accession)
    validate_transcript_organism(record, args)
    output.write_text(fasta_text)

    cached = write_transcript_cache(record, fasta_text)
    if cached is not None:
        record_json, fasta_path = cached
        record["_cached_fasta_path"] = str(fasta_path)
        record["_transcript_cache_record_path"] = str(record_json)
        _set_transcript_cache_metadata(record, "miss_fetched_cached", record_json, fasta_path)
    else:
        record_json, fasta_path = transcript_record_cache_paths(accession)
        _set_transcript_cache_metadata(record, "miss_fetched_cache_write_failed", record_json, fasta_path)
    _update_resolution_after_sequence(record)
    return len(sequence)


# -----------------------------------------------------------------------------
# Transcript annotation (reporting / plotting only)
# -----------------------------------------------------------------------------


def _xml_local_name(tag: str) -> str:
    """Return an XML local tag name without a namespace prefix."""
    return str(tag).rsplit("}", 1)[-1]


def _xml_direct_text(node: ET.Element, names: set[str]) -> str:
    """Return text from the first direct child whose local name is in *names*."""
    for child in list(node):
        if _xml_local_name(child.tag) in names:
            return (child.text or child.attrib.get("value") or "").strip()
    return ""


def _truthy_xml_text(text: str) -> bool:
    return str(text or "").strip().casefold() in {"true", "1", "yes"}


def transcript_annotation_cache_path(accession: str) -> Path:
    """Return the shared per-accession annotation cache path.

    The cache is intentionally outside individual B1-B5 output directories so
    the all-channel helper can fetch an exact RefSeq annotation once and reuse
    it for the remaining channels. Helper/manifest workflows keep it under the
    catalogue root; HCRPROBEFORGE_CACHE_DIR remains an explicit override.
    """
    root = hcrprobeforge_cache_root()
    return root / "transcript_annotations" / f"{safe_name(accession)}.json"


def _annotation_unavailable(
    accession: str,
    sequence_length: int,
    reason: str,
    *,
    source: str = "NCBI RefSeq/GenBank feature table (EFetch)",
) -> dict[str, Any]:
    return {
        "status": "unavailable",
        "accession": accession,
        "transcript_length_nt": int(sequence_length),
        "coordinate_system": "1-based inclusive transcript coordinates",
        "source": source,
        "reason": reason,
        "regions": [],
    }


def parse_genbank_transcript_annotation(
    xml_text: str,
    *,
    accession: str,
    sequence_length: int,
    record_partial: bool = False,
) -> dict[str, Any]:
    """Parse an exact RefSeq transcript CDS from NCBI GenBank XML.

    The feature table is used only for annotation/reporting.  It does not alter
    candidate generation, QC, or selection.  Coordinates are normalized to the
    same 1-based inclusive transcript convention used by HCRProbeForge maps.
    """
    try:
        root = ET.fromstring(xml_text)
    except ET.ParseError as exc:
        raise RuntimeError(f"Could not parse NCBI GenBank XML for {accession}: {exc}") from exc

    sequence_nodes = [
        element
        for element in root.iter()
        if _xml_local_name(element.tag) in {"GBSeq", "INSDSeq"}
    ]
    if _xml_local_name(root.tag) in {"GBSeq", "INSDSeq"} and root not in sequence_nodes:
        sequence_nodes.insert(0, root)

    matched: ET.Element | None = None
    seen_accessions: list[str] = []
    for seq_node in sequence_nodes:
        acc = _xml_direct_text(
            seq_node,
            {"GBSeq_accession-version", "INSDSeq_accession-version"},
        )
        if acc:
            seen_accessions.append(acc)
        if acc == accession:
            matched = seq_node
            break

    if matched is None:
        raise RuntimeError(
            f"NCBI annotation accession mismatch for {accession}; XML contained "
            f"{', '.join(seen_accessions) if seen_accessions else 'no accession-version field'}."
        )

    xml_length_text = _xml_direct_text(matched, {"GBSeq_length", "INSDSeq_length"})
    try:
        xml_length = int(xml_length_text)
    except (TypeError, ValueError):
        raise RuntimeError(f"NCBI annotation for {accession} did not contain a valid transcript length.")

    if xml_length != int(sequence_length):
        raise RuntimeError(
            f"NCBI annotation length mismatch for {accession}: annotation={xml_length}, "
            f"downloaded FASTA={sequence_length}."
        )

    cds_features: list[dict[str, Any]] = []
    for feature in matched.iter():
        if _xml_local_name(feature.tag) not in {"GBFeature", "INSDFeature"}:
            continue
        key = _xml_direct_text(feature, {"GBFeature_key", "INSDFeature_key"})
        if key != "CDS":
            continue

        location = _xml_direct_text(feature, {"GBFeature_location", "INSDFeature_location"})
        partial5 = _truthy_xml_text(
            _xml_direct_text(feature, {"GBFeature_partial5", "INSDFeature_partial5"})
        ) or "<" in location
        partial3 = _truthy_xml_text(
            _xml_direct_text(feature, {"GBFeature_partial3", "INSDFeature_partial3"})
        ) or ">" in location

        intervals: list[tuple[int, int]] = []
        for interval in feature.iter():
            if _xml_local_name(interval.tag) not in {"GBInterval", "INSDInterval"}:
                continue
            start_text = _xml_direct_text(interval, {"GBInterval_from", "INSDInterval_from"})
            end_text = _xml_direct_text(interval, {"GBInterval_to", "INSDInterval_to"})
            point_text = _xml_direct_text(interval, {"GBInterval_point", "INSDInterval_point"})
            try:
                if start_text and end_text:
                    a, b = int(start_text), int(end_text)
                    intervals.append((min(a, b), max(a, b)))
                elif point_text:
                    point = int(point_text)
                    intervals.append((point, point))
            except ValueError:
                continue

        if not intervals and location:
            for a_text, b_text in re.findall(r"<?(\d+)\.\.>?(\d+)", location):
                a, b = int(a_text), int(b_text)
                intervals.append((min(a, b), max(a, b)))
            if not intervals:
                points = [int(value) for value in re.findall(r"\d+", location)]
                if len(points) == 1:
                    intervals.append((points[0], points[0]))

        intervals = sorted(set(intervals))
        if not intervals:
            continue
        cds_start = min(start for start, _ in intervals)
        cds_end = max(end for _, end in intervals)
        if cds_start < 1 or cds_end > xml_length or cds_start > cds_end:
            continue
        cds_features.append(
            {
                "start": cds_start,
                "end": cds_end,
                "intervals": [{"start": a, "end": b} for a, b in intervals],
                "location": location,
                "partial_5prime": bool(partial5),
                "partial_3prime": bool(partial3),
                "total_interval_nt": sum(b - a + 1 for a, b in intervals),
            }
        )

    retrieved_at = utc_now_iso()
    source = "NCBI RefSeq/GenBank feature table (EFetch)"
    base: dict[str, Any] = {
        "accession": accession,
        "transcript_length_nt": xml_length,
        "coordinate_system": "1-based inclusive transcript coordinates",
        "source": source,
        "retrieved_at_utc": retrieved_at,
        "cds_feature_count": len(cds_features),
    }

    if not cds_features:
        if accession.startswith(("NR_", "XR_")):
            return base | {
                "status": "noncoding_no_cds",
                "transcript_type": "noncoding",
                "cds": None,
                "five_prime_utr": None,
                "three_prime_utr": None,
                "regions": [
                    {"kind": "noncoding", "start": 1, "end": xml_length, "label": "non-coding transcript"}
                ],
                "note": "No CDS feature was present; UTR labels are not assigned to a non-coding transcript.",
            }
        return base | {
            "status": "cds_unavailable",
            "transcript_type": "unknown_or_untranslated",
            "cds": None,
            "five_prime_utr": None,
            "three_prime_utr": None,
            "regions": [],
            "note": "No usable CDS feature was present in the exact transcript record; the map remains unsegmented.",
        }

    # RefSeq mRNAs normally contain one CDS.  If several are present, use the
    # longest interval span deterministically and retain the count in metadata.
    cds_features.sort(
        key=lambda item: (
            -int(item["total_interval_nt"]),
            int(item["start"]),
            int(item["end"]),
        )
    )
    cds = dict(cds_features[0])
    partial5 = bool(cds["partial_5prime"])
    partial3 = bool(cds["partial_3prime"])
    if record_partial and not (partial5 or partial3):
        # The title marks the transcript partial but the feature table does not
        # specify which end.  Do not mislabel terminal sequence as a true UTR.
        partial5 = True
        partial3 = True
        cds["partial_5prime"] = True
        cds["partial_3prime"] = True

    five_utr = None
    three_utr = None
    regions: list[dict[str, Any]] = []
    if cds["start"] > 1:
        if partial5:
            regions.append(
                {
                    "kind": "partial_or_unannotated",
                    "start": 1,
                    "end": cds["start"] - 1,
                    "label": "5′ terminal sequence (partial/unannotated)",
                }
            )
        else:
            five_utr = {"start": 1, "end": cds["start"] - 1}
            regions.append({"kind": "5UTR", **five_utr, "label": "5′ UTR"})

    cds_intervals = [(int(item["start"]), int(item["end"])) for item in cds.get("intervals", [])]
    if not cds_intervals:
        cds_intervals = [(int(cds["start"]), int(cds["end"]))]
    for interval_index, (interval_start, interval_end) in enumerate(cds_intervals):
        if interval_index > 0:
            previous_end = cds_intervals[interval_index - 1][1]
            if interval_start > previous_end + 1:
                regions.append(
                    {
                        "kind": "partial_or_unannotated",
                        "start": previous_end + 1,
                        "end": interval_start - 1,
                        "label": "internal non-CDS / unannotated sequence",
                    }
                )
        regions.append({"kind": "CDS", "start": interval_start, "end": interval_end, "label": "CDS"})

    if cds["end"] < xml_length:
        if partial3:
            regions.append(
                {
                    "kind": "partial_or_unannotated",
                    "start": cds["end"] + 1,
                    "end": xml_length,
                    "label": "3′ terminal sequence (partial/unannotated)",
                }
            )
        else:
            three_utr = {"start": cds["end"] + 1, "end": xml_length}
            regions.append({"kind": "3UTR", **three_utr, "label": "3′ UTR"})

    note_parts: list[str] = []
    if len(cds_features) > 1:
        note_parts.append(
            f"{len(cds_features)} CDS features were present; the longest CDS feature was used for visualization."
        )
    if partial5 or partial3:
        note_parts.append(
            "The CDS/transcript is marked partial; uncertain terminal sequence is not labelled as a UTR."
        )

    return base | {
        "status": "resolved_partial" if (partial5 or partial3) else "resolved",
        "transcript_type": "protein_coding",
        "cds": {
            "start": int(cds["start"]),
            "end": int(cds["end"]),
            "intervals": cds["intervals"],
            "location": cds["location"],
            "partial_5prime": partial5,
            "partial_3prime": partial3,
        },
        "five_prime_utr": five_utr,
        "three_prime_utr": three_utr,
        "regions": regions,
        "note": " ".join(note_parts) if note_parts else None,
    }


def resolve_transcript_annotation(
    session: requests.Session | None,
    record: dict[str, Any],
    args: argparse.Namespace,
    sequence_length: int,
) -> dict[str, Any]:
    """Return exact-transcript 5′ UTR/CDS/3′ UTR annotation when available.

    Annotation failure is deliberately non-fatal: probe design must continue if
    NCBI feature metadata are temporarily unavailable.  Exact accession and
    transcript-length checks prevent a stale or mismatched annotation from being
    drawn on the map.
    """
    accession = str(record.get("accession") or "")
    if record.get("source") == "ncbi_premrna" and isinstance(record.get("premrna_target"), dict):
        target = record["premrna_target"]
        return {
            "status": "resolved_premrna",
            "mode": target.get("mode") or "pre-mRNA-intronic",
            "selection_mode": target.get("selection_mode") or "intronic",
            "accession": accession,
            "transcript_length_nt": int(sequence_length),
            "regions": target.get("regions", []),
            "source": target.get("source"),
            "note": target.get("note"),
        }
    if record.get("source") == "user_fasta" or session is None:
        return _annotation_unavailable(
            accession,
            sequence_length,
            "Custom FASTA input has no verified RefSeq feature annotation.",
            source="user FASTA",
        )
    if not accession:
        return _annotation_unavailable(accession, sequence_length, "No exact transcript accession was available.")

    cache_path = transcript_annotation_cache_path(accession)
    try:
        if cache_path.exists():
            cached = json.loads(cache_path.read_text())
            if (
                cached.get("accession") == accession
                and int(cached.get("transcript_length_nt", -1)) == int(sequence_length)
            ):
                cached["cache_status"] = "hit"
                cached["cache_path"] = str(cache_path)
                return cached
    except (OSError, ValueError, TypeError, json.JSONDecodeError):
        pass

    # Use the exact versioned accession rather than a numeric UID. This keeps
    # the requested RefSeq version authoritative even when NCBI has a newer
    # sequence revision for the same accession stem.
    params = base_params(args) | {
        "db": "nuccore",
        "id": accession,
        "rettype": "gb",
        "retmode": "xml",
    }
    try:
        response = ncbi_get(session, "efetch.fcgi", params)
        annotation = parse_genbank_transcript_annotation(
            response.text,
            accession=accession,
            sequence_length=sequence_length,
            record_partial=transcript_is_partial(record),
        )
        annotation["cache_status"] = "miss_fetched"
        annotation["cache_path"] = str(cache_path)
        try:
            cache_path.parent.mkdir(parents=True, exist_ok=True)
            write_json(cache_path, annotation)
        except OSError as exc:
            annotation["cache_write_warning"] = str(exc)
        return annotation
    except (RuntimeError, OSError, requests.RequestException, ET.ParseError) as exc:
        print(
            f"WARNING: transcript annotation for {accession} could not be resolved; "
            f"probe design will continue without UTR/CDS map segmentation: {exc}",
            file=sys.stderr,
        )
        result = _annotation_unavailable(accession, sequence_length, str(exc))
        result["cache_status"] = "unavailable"
        result["cache_path"] = str(cache_path)
        return result


def _interval_overlap_nt(start_a: int, end_a: int, start_b: int, end_b: int) -> int:
    return max(0, min(int(end_a), int(end_b)) - max(int(start_a), int(start_b)) + 1)


def annotate_probe_with_transcript_region(
    row: dict[str, Any],
    annotation: dict[str, Any] | None,
) -> dict[str, Any]:
    """Add transcript-region overlap fields to a probe/candidate dictionary."""
    start = int(row["start"])
    end = int(row.get("end", start + int(row["length"]) - 1))
    row["end"] = end

    overlaps = {
        "5UTR": 0,
        "CDS": 0,
        "3UTR": 0,
        "noncoding": 0,
        "partial_or_unannotated": 0,
    }
    for region in (annotation or {}).get("regions", []) or []:
        kind = str(region.get("kind", ""))
        if kind not in overlaps:
            continue
        overlaps[kind] += _interval_overlap_nt(start, end, int(region["start"]), int(region["end"]))

    row["five_prime_utr_overlap_nt"] = overlaps["5UTR"]
    row["cds_overlap_nt"] = overlaps["CDS"]
    row["three_prime_utr_overlap_nt"] = overlaps["3UTR"]
    row["noncoding_overlap_nt"] = overlaps["noncoding"]
    row["unannotated_overlap_nt"] = overlaps["partial_or_unannotated"]

    premrna_regions = [
        region for region in (annotation or {}).get("regions", []) or []
        if str(region.get("kind", "")) in {"exon", "intron"}
    ]
    row.setdefault("premrna_intron_number", "")
    row.setdefault("premrna_genomic_start", "")
    row.setdefault("premrna_genomic_end", "")

    positive = [kind for kind in ("5UTR", "CDS", "3UTR") if overlaps[kind] > 0]
    if overlaps["noncoding"] > 0 and not positive:
        region_label = "noncoding"
    elif len(positive) == 1 and overlaps["partial_or_unannotated"] == 0:
        region_label = positive[0]
    elif len(positive) >= 2:
        region_label = "-".join(positive) + "_boundary"
    elif positive and overlaps["partial_or_unannotated"] > 0:
        region_label = positive[0] + "-partial_boundary"
    elif overlaps["partial_or_unannotated"] > 0:
        region_label = "partial_or_unannotated"
    else:
        region_label = "unannotated"
    row["transcript_region"] = region_label
    if premrna_regions:
        fully_contained = [
            region for region in premrna_regions
            if start >= int(region["start"]) and end <= int(region["end"])
        ]
        intron_hits = [
            region for region in premrna_regions
            if str(region.get("kind")) == "intron"
            and _interval_overlap_nt(start, end, int(region["start"]), int(region["end"])) > 0
        ]
        if len(fully_contained) == 1:
            region = fully_contained[0]
            if str(region.get("kind")) == "intron":
                row["premrna_intron_number"] = int(region["intron_number"])
                row["premrna_genomic_start"] = int(region.get("genomic_start", 0) or 0)
                row["premrna_genomic_end"] = int(region.get("genomic_end", 0) or 0)
                row["transcript_region"] = f"intron_{int(region['intron_number'])}"
            else:
                row["transcript_region"] = f"exon_{int(region.get('exon_number', 0) or 0)}"
        elif intron_hits:
            row["premrna_intron_number"] = "boundary"
            row["transcript_region"] = "intron_boundary"
        else:
            row["premrna_intron_number"] = "separator"
            row["transcript_region"] = "intron_separator"
    return row


def annotate_probe_rows(
    rows: list[dict[str, Any]],
    annotation: dict[str, Any] | None,
) -> list[dict[str, Any]]:
    for row in rows:
        annotate_probe_with_transcript_region(row, annotation)
    return rows


def selected_probe_region_counts(rows: list[dict[str, Any]]) -> dict[str, int]:
    counts: dict[str, int] = {}
    for row in rows:
        key = str(row.get("transcript_region") or "unannotated")
        counts[key] = counts.get(key, 0) + 1
    return dict(sorted(counts.items()))


def load_sibling_transcript_annotation(path: Path) -> dict[str, Any] | None:
    """Best-effort annotation discovery for --plot-only regeneration."""
    try:
        candidates = list(path.expanduser().resolve().parent.glob("*_transcript_resolution.json"))
        if not candidates:
            return None
        metadata_path = max(candidates, key=lambda item: item.stat().st_mtime_ns)
        data = json.loads(metadata_path.read_text())
        annotation = data.get("transcript_annotation")
        return annotation if isinstance(annotation, dict) else None
    except (OSError, ValueError, TypeError, json.JSONDecodeError):
        return None


def wrap_fasta_sequence(sequence: str, width: int = 80) -> str:
    """Return a FASTA-wrapped sequence string."""
    return "\n".join(sequence[i : i + width] for i in range(0, len(sequence), width))


def parse_fasta_records(path: Path) -> list[dict[str, str]]:
    """Parse a FASTA file and return records with header and cleaned sequence.

    HCRProbeForge currently treats --fasta as a single-target input mode. Multi-record
    FASTA files are rejected so users do not accidentally design only the first record.
    """
    fasta_path = path.expanduser().resolve()
    if not fasta_path.exists():
        raise RuntimeError(f"FASTA file does not exist: {fasta_path}")
    if not fasta_path.is_file():
        raise RuntimeError(f"FASTA path is not a file: {fasta_path}")

    records: list[dict[str, str]] = []
    header: str | None = None
    seq_parts: list[str] = []

    for line_number, raw_line in enumerate(fasta_path.read_text().splitlines(), start=1):
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                records.append({"header": header, "sequence": "".join(seq_parts)})
            header = line[1:].strip() or fasta_path.stem
            seq_parts = []
        else:
            if header is None:
                raise RuntimeError(
                    f"Invalid FASTA file {fasta_path}: sequence appears before a header at line {line_number}."
                )
            seq_parts.append(re.sub(r"\s+", "", line))

    if header is not None:
        records.append({"header": header, "sequence": "".join(seq_parts)})

    if not records:
        raise RuntimeError(f"No FASTA records found in {fasta_path}.")

    for index, record in enumerate(records, start=1):
        sequence = record["sequence"].upper()
        if not sequence:
            raise RuntimeError(f"FASTA record {index} in {fasta_path} has an empty sequence.")
        if re.search(r"[^ACGTUN]", sequence):
            bad_chars = "".join(sorted(set(re.findall(r"[^ACGTUN]", sequence))))
            raise RuntimeError(
                f"FASTA record {index} in {fasta_path} contains unsupported bases: {bad_chars!r}. "
                "Use A/C/G/T/U/N only. U is converted to T."
            )
        record["sequence"] = sequence.replace("U", "T")

    return records


def resolve_user_fasta_record(args: argparse.Namespace) -> list[dict[str, Any]]:
    """Resolve a user-supplied single-record FASTA into the internal target record format."""
    fasta_path = args.fasta.expanduser().resolve()
    records = parse_fasta_records(fasta_path)
    if len(records) != 1:
        raise RuntimeError(
            f"{fasta_path} contains {len(records)} FASTA records. "
            "This version of HCRProbeForge accepts one record per --fasta run. "
            "Split the file and run each target separately. Multi-record FASTA batch mode can be added later."
        )

    fasta_record = records[0]
    target_name = safe_name(fasta_path.stem) or "custom_target"
    sequence = fasta_record["sequence"]
    return [
        {
            "uid": None,
            "accession": target_name,
            "title": f"user-supplied FASTA: {fasta_path.name}; record: {fasta_record['header']}",
            "length": len(sequence),
            "source": "user_fasta",
            "input_fasta": str(fasta_path),
            "fasta_header": fasta_record["header"],
            "_target_label": target_name,
            "_sequence": sequence,
        }
    ]


def target_label_for_record(args: argparse.Namespace, record: dict[str, Any]) -> str:
    """Return the filesystem-safe technical target label used for files and designProbes --targetName."""
    if record.get("source") == "user_fasta":
        return safe_name(record.get("_target_label") or record.get("accession") or "custom_target")
    if record.get("source") == "ncbi_premrna":
        target = record.get("premrna_target") or {}
        if str(target.get("selection_mode") or "").casefold() == "whole":
            suffix = "premrna_whole"
        elif str(target.get("distribution") or "") == "selected_introns":
            selected = target.get("selected_introns") or []
            numbers = ",".join(str(int(number)) for number in selected)
            suffix = f"premrna_introns{numbers}" if numbers else "premrna_introns"
        else:
            suffix = "premrna_introns"
        # Keep the comma-separated intron list readable in the folder name;
        # the gene and accession portions have already been normalized.
        prefix = safe_name(f"{args.gene or 'target'}_{record['accession']}")
        return f"{prefix}_{suffix}"
    return safe_name(f"{args.gene or 'target'}_{record['accession']}")


def record_target_type(record: dict[str, Any]) -> str:
    """Return the stable report label for the actual biological target."""
    if record.get("source") != "ncbi_premrna":
        return "mature"
    target = record.get("premrna_target") or {}
    return "pre-mRNA-whole" if str(target.get("selection_mode") or "").casefold() == "whole" else "pre-mRNA"


def write_user_fasta(record: dict[str, Any], output: Path) -> int:
    """Write a normalised single-record FASTA for a user-supplied target sequence."""
    sequence = str(record.get("_sequence", "")).upper().replace("U", "T")
    if not sequence or re.search(r"[^ACGTN]", sequence):
        raise RuntimeError(f"Invalid user-supplied FASTA sequence for {record.get('accession', 'custom target')}.")
    header = safe_name(str(record.get("accession") or "custom_target")) or "custom_target"
    original_header = str(record.get("fasta_header") or "")
    title = f"{header} user_supplied_fasta original_header={original_header}".strip()
    output.write_text(f">{title}\n{wrap_fasta_sequence(sequence)}\n")
    return len(sequence)


def prepare_fasta_for_record(
    session: requests.Session | None,
    record: dict[str, Any],
    output: Path,
    args: argparse.Namespace,
) -> int:
    """Write the FASTA used by designProbes and return its sequence length."""
    if record.get("source") == "user_fasta":
        return write_user_fasta(record, output)
    if record.get("source") == "ncbi_premrna":
        sequence = str(record.get("_sequence") or "").upper().replace("U", "T")
        if not sequence or re.search(r"[^ACGTN]", sequence):
            raise RuntimeError("The prepared Pre-mRNA target contains an invalid base.")
        header = safe_name(str(record.get("accession") or "target")) or "target"
        output.write_text(f">{header}_Pre-mRNA\n{wrap_fasta_sequence(sequence)}\n", encoding="utf-8")
        return len(sequence)
    if session is None:
        raise RuntimeError("Internal error: NCBI transcript retrieval requires an active requests session.")
    return fetch_fasta(session, record, output, args)


def build_designprobes_command(
    executable: str,
    fasta: Path,
    target: str,
    probes_tsv: Path,
    idt_tsv: Path,
    args: argparse.Namespace,
) -> list[str]:
    """Build the exact designProbes command used by the pipeline."""
    command = [
        executable,
        str(fasta),
        "--species",
        args.species,
        "--channel",
        args.channel,
        "--targetName",
        designprobes_target_name(target),
        "--tileSize",
        str(args.tile_size),
        "--minGC",
        str(args.min_gc),
        "--maxGC",
        str(args.max_gc),
        "--minGibbs",
        str(args.min_gibbs),
        "--maxGibbs",
        str(args.max_gibbs),
        "--targetGibbs",
        str(args.target_gibbs),
        "--maxRunMismatches",
        str(args.max_run_mismatches),
        "--maxProbes",
        str(args.max_probes),
        "--num-hits-allowed",
        str(args.num_hits_allowed),
        "--output",
        str(probes_tsv),
        "--idt",
        str(idt_tsv),
    ]
    if args.index:
        command += ["--index", args.index]
    if args.no_genomemask:
        command.append("--no-genomemask")
    if args.dtm_filter:
        command += ["--dTmFilter", "--dTmMax", str(args.dtm_max)]
    return command


def make_run_metadata(
    *,
    fasta: Path,
    output_dir: Path,
    target: str,
    record: dict[str, Any],
    sequence_length: int,
    probes_tsv: Path,
    idt_tsv: Path,
    log_file: Path,
    run_parameters_json: Path,
    command: list[str],
    args: argparse.Namespace,
) -> dict[str, Any]:
    """Create a reproducibility record for a designProbes run."""
    return {
        "status": "started",
        "started_at_utc": utc_now_iso(),
        "completed_at_utc": None,
        "duration_seconds": None,
        "return_code": None,
        "pipeline": {
            "script": str(Path(__file__).resolve()),
            "script_sha256": sha256_file(Path(__file__).resolve()),
            "argv": sanitized_argv(list(getattr(args, "_pipeline_argv", sys.argv[:]))),
            "python_executable": sys.executable,
            "python_version": sys.version,
            "platform": platform.platform(),
            "package_versions": {
                "hcrprobeforge": __version__,
                "hcrprobedesign": package_version("hcrprobedesign"),
                "primer3-py": package_version("primer3-py"),
                "pandas": package_version("pandas"),
                "openpyxl": package_version("openpyxl"),
                "requests": package_version("requests"),
                "bowtie2": executable_version("bowtie2"),
                "bowtie2-build": executable_version("bowtie2-build"),
            },
        },
        "input": {
            "target_name": target,
            "gene_symbol": args.gene,
            "display_label": args.gene or record.get("accession"),
            "species": args.species,
            "organism": args.organism,
            "reference": getattr(args, "reference", None),
            "fasta_path": str(fasta),
            "fasta_sha256": sha256_file(fasta) if fasta.exists() else None,
            "record": {
                "source": record.get("source", "ncbi"),
                "target_type": record_target_type(record),
                "uid": record.get("uid"),
                "accession": record.get("accession"),
                "title": record.get("title", ""),
                "input_fasta": record.get("input_fasta"),
                "fasta_header": record.get("fasta_header"),
                "ncbi_reported_length_nt": record.get("length") if record.get("source", "ncbi") == "ncbi" else None,
                "target_sequence_length_nt": sequence_length,
                "downloaded_fasta_length_nt": sequence_length if record.get("source", "ncbi") == "ncbi" else None,
                "transcript_resolution": record.get("transcript_resolution"),
                "transcript_annotation": record.get("transcript_annotation"),
                "premrna_target": record.get("premrna_target"),
                "transcript_resolution_json": record.get("transcript_resolution_json"),
            },
        },
        "arguments": sanitized_args(args),
        "curation_context": {
            "phase": getattr(args, "_curation_phase", None),
            "plan": getattr(args, "_curation_plan", None),
        },
        "designprobes": {
            "executable": command[0],
            "command": command,
            "hcrprobedesign_version": package_version("hcrprobedesign"),
        },
        "outputs": {
            "output_dir": str(output_dir),
            "probes_tsv": str(probes_tsv),
            "idt_tsv": str(idt_tsv),
            "log_file": str(log_file),
            "run_parameters_json": str(run_parameters_json),
        },
    }


def write_log_header(log_file: Path, run_metadata: dict[str, Any]) -> None:
    """Write reproducibility metadata at the top of the designProbes log."""
    log_file.write_text(
        "# HCRProbeForge run parameters\n"
        + json.dumps(json_ready(run_metadata), indent=2, sort_keys=True)
        + "\n\n# designProbes output\n"
    )


def run_design_probes(
    fasta: Path,
    output_dir: Path,
    target: str,
    args: argparse.Namespace,
    record: dict[str, Any],
    sequence_length: int,
) -> tuple[Path, Path, Path, Path, dict[str, Any]]:
    executable = shutil.which("designProbes")
    if not executable:
        raise RuntimeError("designProbes was not found. Install it with: pip install hcrprobedesign")

    probes_tsv = output_dir / f"{target}_probes.tsv"
    idt_tsv = output_dir / f"{target}_IDT.tsv"
    log_file = output_dir / f"{target}.log"
    run_parameters_json = output_dir / f"{target}_run_parameters.json"

    command = build_designprobes_command(executable, fasta, target, probes_tsv, idt_tsv, args)
    run_metadata = make_run_metadata(
        fasta=fasta,
        output_dir=output_dir,
        target=target,
        record=record,
        sequence_length=sequence_length,
        probes_tsv=probes_tsv,
        idt_tsv=idt_tsv,
        log_file=log_file,
        run_parameters_json=run_parameters_json,
        command=command,
        args=args,
    )

    # Write a started-state JSON immediately, then run designProbes.
    # Do not pre-write or prepend metadata to the main .log file: HCRProbeDesign/designProbes
    # may create that log internally, and users generally expect it to remain the native
    # designProbes log. Reproducibility metadata is stored in *_run_parameters.json and
    # copied into the summary JSON instead.
    write_json(run_parameters_json, run_metadata)
    start_time = time.time()

    # Preserve the original package-style log: do not prepend metadata and do not
    # append helper sections for captured stdout/stderr. The run parameters are
    # stored only in *_run_parameters.json and the summary JSON. Redirect stderr
    # to <target>.log, matching the earlier behavior, while leaving
    # designProbes stdout out of the log so lines such as printed Namespace(...)
    # do not pollute the package log.
    raise_if_cancelled()
    phase = str(getattr(args, "_curation_phase", "initial") or "initial")
    notify_progress(
        "designprobes",
        f"Running designProbes and Bowtie2 specificity screening ({phase} pass)…",
    )
    with log_file.open("w") as log_handle:
        process = subprocess.Popen(
            command,
            cwd=output_dir,
            stdout=subprocess.DEVNULL,
            stderr=log_handle,
            text=True,
            start_new_session=(os.name != "nt"),
            creationflags=(subprocess.CREATE_NEW_PROCESS_GROUP if os.name == "nt" else 0),
        )
        try:
            while True:
                try:
                    return_code = process.wait(timeout=0.25)
                    break
                except subprocess.TimeoutExpired:
                    try:
                        raise_if_cancelled()
                    except RunCancelled:
                        _stop_design_process(process)
                        try:
                            process.wait(timeout=3)
                        except subprocess.TimeoutExpired:
                            _stop_design_process(process, force=True)
                            process.wait(timeout=3)
                        raise
        except RunCancelled:
            run_metadata["completed_at_utc"] = utc_now_iso()
            run_metadata["duration_seconds"] = round(time.time() - start_time, 3)
            run_metadata["return_code"] = None
            run_metadata["status"] = "cancelled"
            write_json(run_parameters_json, run_metadata)
            raise

    completed_at = utc_now_iso()
    run_metadata["completed_at_utc"] = completed_at
    run_metadata["duration_seconds"] = round(time.time() - start_time, 3)
    run_metadata["return_code"] = return_code
    run_metadata["status"] = "completed" if return_code == 0 else "failed"

    write_json(run_parameters_json, run_metadata)

    if return_code != 0:
        raise subprocess.CalledProcessError(return_code, command)
    notify_progress(
        "designprobes",
        f"Completed designProbes ({phase} pass); reading candidate probes…",
    )
    return probes_tsv, idt_tsv, log_file, run_parameters_json, run_metadata


# -----------------------------------------------------------------------------
# Final HCR oligo-structure QC
# -----------------------------------------------------------------------------

# Dependencies are imported lazily so --help and --plot-only keep working even in
# environments where the QC dependencies are not installed. A normal design run
# uses this QC by default after designProbes creates the IDT oligo table.
_pd = None
_primer3 = None


def require_qc_dependencies():
    """Load pandas and primer3-py only when final-oligo QC is requested."""
    global _pd, _primer3
    if _pd is None:
        try:
            import pandas as pandas_module
        except ImportError as exc:  # pragma: no cover - environment-specific
            raise RuntimeError(
                "Final oligo QC requires pandas and openpyxl. Install them with: "
                "python -m pip install pandas openpyxl, or rerun with --no-oligo-qc."
            ) from exc
        _pd = pandas_module
    if _primer3 is None:
        try:
            import primer3 as primer3_module
        except ImportError as exc:  # pragma: no cover - environment-specific
            raise RuntimeError(
                "Final oligo QC requires primer3-py. Install it with: "
                "python -m pip install primer3-py, or rerun with --no-oligo-qc."
            ) from exc
        _primer3 = primer3_module
    return _pd, _primer3


# Thresholds are tuned for final HCR detection oligos. Hairpin Tm alone is not
# treated as a hard failure; self-dimers and same-pair P1/P2 heterodimers are
# weighted more heavily than unrelated cross-pool heterodimers.
HAIRPIN_REVIEW_TM = 50.0
HAIRPIN_STRONG_TM = 60.0
HAIRPIN_REVIEW_DG = -5.0
HAIRPIN_STRONG_DG = -8.0
HAIRPIN_STRONG_TM_DG = -7.0

SELF_DIMER_REVIEW_DG = -9.0
SELF_DIMER_STRONG_DG = -12.0
SELF_DIMER_PRACTICAL_REJECT_DG = -15.0

SAME_PAIR_REVIEW_DG = -10.0
SAME_PAIR_STRONG_DG = -12.0
SAME_PAIR_STRONG_TM = 37.0
SAME_PAIR_STRONG_TM_DG = -11.0

CROSS_DIMER_REVIEW_DG = -10.5
CROSS_DIMER_STRONG_DG = -13.0
CROSS_DIMER_PRACTICAL_REJECT_DG = -15.0

FLAG_RANK = {"PASS": 0, "REVIEW": 1, "STRONG_REVIEW": 2}


def clean_oligo_seq(seq: Any) -> str:
    """Return uppercase A/C/G/T sequence; strip whitespace and non-DNA chars."""
    pd, _ = require_qc_dependencies()
    if pd.isna(seq):
        return ""
    seq = str(seq).upper().replace(" ", "").replace("\n", "").replace("\r", "")
    seq = re.sub(r"[^ACGT]", "", seq)
    return seq


def normalize_colname(name: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", str(name).strip().lower()).strip("_")


def read_oligo_table(path: Path):
    """Read CSV/TSV/TXT/XLSX and return a name/sequence oligo table.

    Accepted inputs:
      * IDT-style tables with name and sequence columns.
      * Headerless two-column files: name, sequence.
      * Probe-pair tables with P1 and P2 columns; these are expanded into
        one oligo row per P1/P2 sequence.
    """
    pd, _ = require_qc_dependencies()
    suffix = path.suffix.lower()
    if suffix in {".xlsx", ".xls"}:
        df = pd.read_excel(path)
    else:
        # Auto-detect comma vs tab. This handles normal IDT CSV and TSV files.
        df = pd.read_csv(path, sep=None, engine="python")

    original_cols = list(df.columns)
    normalized_to_original = {normalize_colname(c): c for c in original_cols}

    name_candidates = ["name", "oligo_name", "oligo", "id", "probe_name"]
    sequence_candidates = ["sequence", "seq", "oligo_sequence", "bases"]

    name_col = next((normalized_to_original[c] for c in name_candidates if c in normalized_to_original), None)
    seq_col = next((normalized_to_original[c] for c in sequence_candidates if c in normalized_to_original), None)

    # Support selected-pairs tables that contain one P1 and one P2 sequence per row.
    if seq_col is None:
        p1_col = normalized_to_original.get("p1")
        p2_col = normalized_to_original.get("p2")
        if p1_col and p2_col:
            long_rows = []
            for idx, row in df.iterrows():
                if name_col:
                    base_name = str(row[name_col]).strip()
                elif "probe" in normalized_to_original:
                    base_name = str(row[normalized_to_original["probe"]]).strip()
                elif "start" in normalized_to_original:
                    start_value = str(row[normalized_to_original["start"]]).strip()
                    base_name = f"probe_{start_value}"
                else:
                    base_name = f"probe_{idx + 1}"
                long_rows.append({"name": f"{base_name}_P1", "sequence": row[p1_col]})
                long_rows.append({"name": f"{base_name}_P2", "sequence": row[p2_col]})
            df = pd.DataFrame(long_rows)
            name_col = "name"
            seq_col = "sequence"

    # Support headerless IDT-style two-column files if header detection failed.
    if name_col is None or seq_col is None:
        if suffix not in {".xlsx", ".xls"}:
            raw = pd.read_csv(path, sep=None, engine="python", header=None)
            if raw.shape[1] >= 2:
                candidate = raw.iloc[:, :2].copy()
                candidate.columns = ["name", "sequence"]
                cleaned = candidate["sequence"].map(clean_oligo_seq)
                # Treat as headerless if most non-empty sequence entries look like DNA.
                dna_like = cleaned.str.len().gt(0).mean()
                if dna_like >= 0.8:
                    df = candidate
                    name_col = "name"
                    seq_col = "sequence"

    if name_col is None or seq_col is None:
        raise RuntimeError(
            f"Missing required oligo QC columns in {path}. Found columns: {original_cols}. "
            "Expected at least name and sequence, or a probe-pair table with P1 and P2 columns."
        )

    if name_col != "name":
        df = df.rename(columns={name_col: "name"})
    if seq_col != "sequence":
        df = df.rename(columns={seq_col: "sequence"})

    df = df.copy()
    df["name"] = df["name"].astype(str).str.strip()
    df["sequence"] = df["sequence"].map(clean_oligo_seq)

    empty_seq = df["sequence"].eq("")
    if empty_seq.any():
        bad_names = ", ".join(df.loc[empty_seq, "name"].head(5).tolist())
        raise RuntimeError(f"Empty/invalid sequence after cleaning for: {bad_names}")

    duplicate_names = df["name"].duplicated(keep=False)
    if duplicate_names.any():
        dupes = ", ".join(sorted(df.loc[duplicate_names, "name"].unique())[:10])
        raise RuntimeError(f"Duplicate oligo names found. Names must be unique. Examples: {dupes}")

    return df[["name", "sequence"]]


def primer3_metric(result: Any, attr: str) -> Any:
    if not getattr(result, "structure_found", False):
        return None
    return getattr(result, attr, None)


def dg_kcal(result: Any) -> Any:
    """primer3-py reports deltaG in cal/mol; convert to kcal/mol."""
    if not getattr(result, "structure_found", False):
        return None
    return result.dg / 1000.0


def call_flag_hairpin(hp: Any) -> Tuple[str, str, bool]:
    """Return flag, reason, practical_reject boolean for a hairpin result."""
    if not getattr(hp, "structure_found", False):
        return "PASS", "no hairpin found", False

    tm = float(hp.tm)
    dg = float(dg_kcal(hp))

    if dg <= HAIRPIN_STRONG_DG:
        return "STRONG_REVIEW", f"hairpin deltaG <= {HAIRPIN_STRONG_DG:g} kcal/mol", False
    if tm >= HAIRPIN_STRONG_TM and dg <= HAIRPIN_STRONG_TM_DG:
        return (
            "STRONG_REVIEW",
            f"hairpin Tm >= {HAIRPIN_STRONG_TM:g} degC and deltaG <= {HAIRPIN_STRONG_TM_DG:g} kcal/mol",
            False,
        )
    if tm >= HAIRPIN_REVIEW_TM:
        return "REVIEW", f"hairpin Tm >= {HAIRPIN_REVIEW_TM:g} degC", False
    if dg <= HAIRPIN_REVIEW_DG:
        return "REVIEW", f"hairpin deltaG <= {HAIRPIN_REVIEW_DG:g} kcal/mol", False
    return "PASS", "hairpin below review thresholds", False


def call_flag_self_dimer(result: Any) -> Tuple[str, str, bool]:
    """Return flag, reason, practical_reject boolean for a homodimer result."""
    if not getattr(result, "structure_found", False):
        return "PASS", "no self-dimer found", False

    dg = float(dg_kcal(result))

    if dg <= SELF_DIMER_PRACTICAL_REJECT_DG:
        return (
            "STRONG_REVIEW",
            f"self-dimer deltaG <= {SELF_DIMER_PRACTICAL_REJECT_DG:g} kcal/mol; practical reject if alternatives exist",
            True,
        )
    if dg <= SELF_DIMER_STRONG_DG:
        return "STRONG_REVIEW", f"self-dimer deltaG <= {SELF_DIMER_STRONG_DG:g} kcal/mol", False
    if dg <= SELF_DIMER_REVIEW_DG:
        return "REVIEW", f"self-dimer deltaG <= {SELF_DIMER_REVIEW_DG:g} kcal/mol", False
    return "PASS", "self-dimer below review thresholds", False


def parse_probe_name(name: str) -> Tuple[str, str]:
    """Return pair_key and arm label from common HCR oligo names."""
    text = str(name).strip()
    match = re.match(r"^(?P<base>.*?)(?:[_\-\s]+)(?P<arm>P1|P2|ODD|EVEN)$", text, flags=re.IGNORECASE)
    if not match:
        return text, ""

    base = match.group("base")
    arm_raw = match.group("arm").upper()
    if arm_raw == "ODD":
        arm = "P1"
    elif arm_raw == "EVEN":
        arm = "P2"
    else:
        arm = arm_raw
    return base, arm


def classify_heterodimer(name_a: str, name_b: str) -> dict[str, str]:
    key_a, arm_a = parse_probe_name(name_a)
    key_b, arm_b = parse_probe_name(name_b)

    if key_a == key_b and {arm_a, arm_b} == {"P1", "P2"}:
        dimer_class = "same_pair_P1_P2"
    elif key_a == key_b:
        dimer_class = "same_probe_region"
    else:
        dimer_class = "cross_pair"

    return {
        "pair_key_1": key_a,
        "arm_1": arm_a,
        "pair_key_2": key_b,
        "arm_2": arm_b,
        "dimer_class": dimer_class,
    }


def call_flag_heterodimer(result: Any, dimer_class: str) -> Tuple[str, str, bool]:
    """Return flag, reason, practical_reject boolean for a heterodimer result."""
    if not getattr(result, "structure_found", False):
        return "PASS", "no heterodimer found", False

    dg = float(dg_kcal(result))
    tm = float(result.tm)

    if dimer_class in {"same_pair_P1_P2", "same_probe_region"}:
        if dg <= SELF_DIMER_PRACTICAL_REJECT_DG:
            return (
                "STRONG_REVIEW",
                f"same-pair heterodimer deltaG <= {SELF_DIMER_PRACTICAL_REJECT_DG:g} kcal/mol; practical reject if alternatives exist",
                True,
            )
        if dg <= SAME_PAIR_STRONG_DG:
            return "STRONG_REVIEW", f"same-pair heterodimer deltaG <= {SAME_PAIR_STRONG_DG:g} kcal/mol", False
        if tm >= SAME_PAIR_STRONG_TM and dg <= SAME_PAIR_STRONG_TM_DG:
            return (
                "STRONG_REVIEW",
                f"same-pair heterodimer Tm >= {SAME_PAIR_STRONG_TM:g} degC and deltaG <= {SAME_PAIR_STRONG_TM_DG:g} kcal/mol",
                False,
            )
        if dg <= SAME_PAIR_REVIEW_DG:
            return "REVIEW", f"same-pair heterodimer deltaG <= {SAME_PAIR_REVIEW_DG:g} kcal/mol", False
        return "PASS", "same-pair heterodimer below review thresholds", False

    if dg <= CROSS_DIMER_PRACTICAL_REJECT_DG:
        return (
            "STRONG_REVIEW",
            f"cross-pool heterodimer deltaG <= {CROSS_DIMER_PRACTICAL_REJECT_DG:g} kcal/mol; practical reject if recurrent",
            True,
        )
    if dg <= CROSS_DIMER_STRONG_DG:
        return "STRONG_REVIEW", f"cross-pool heterodimer deltaG <= {CROSS_DIMER_STRONG_DG:g} kcal/mol", False
    if dg <= CROSS_DIMER_REVIEW_DG:
        return "REVIEW", f"cross-pool heterodimer deltaG <= {CROSS_DIMER_REVIEW_DG:g} kcal/mol", False
    return "PASS", "cross-pool heterodimer below review thresholds", False


def flag_rank(flag: str) -> int:
    return FLAG_RANK.get(flag, -1)


def make_single_qc(df):
    pd, primer3 = require_qc_dependencies()
    rows = []
    for _, row in df.iterrows():
        name = row["name"]
        seq = row["sequence"]

        pair_key, arm = parse_probe_name(name)
        hp = primer3.calc_hairpin(seq)
        self_dim = primer3.calc_homodimer(seq)

        hp_flag, hp_reason, hp_reject = call_flag_hairpin(hp)
        self_flag, self_reason, self_reject = call_flag_self_dimer(self_dim)

        rows.append(
            {
                "name": name,
                "pair_key": pair_key,
                "arm": arm,
                "length": len(seq),
                "sequence": seq,
                "hairpin_found": hp.structure_found,
                "hairpin_tm": primer3_metric(hp, "tm"),
                "hairpin_dg_kcal": dg_kcal(hp),
                "hairpin_flag": hp_flag,
                "hairpin_reason": hp_reason,
                "hairpin_practical_reject": hp_reject,
                "self_dimer_found": self_dim.structure_found,
                "self_dimer_tm": primer3_metric(self_dim, "tm"),
                "self_dimer_dg_kcal": dg_kcal(self_dim),
                "self_dimer_flag": self_flag,
                "self_dimer_reason": self_reason,
                "self_dimer_practical_reject": self_reject,
                "single_max_flag_rank": max(flag_rank(hp_flag), flag_rank(self_flag)),
            }
        )
    return pd.DataFrame(rows)


def make_pair_qc(df):
    pd, primer3 = require_qc_dependencies()
    rows = []
    for (_, a), (_, b) in itertools.combinations(df.iterrows(), 2):
        name_a = a["name"]
        name_b = b["name"]

        result = primer3.calc_heterodimer(a["sequence"], b["sequence"])
        class_info = classify_heterodimer(name_a, name_b)
        dim_flag, dim_reason, dim_reject = call_flag_heterodimer(result, class_info["dimer_class"])

        rows.append(
            {
                "name_1": name_a,
                "name_2": name_b,
                **class_info,
                "heterodimer_found": result.structure_found,
                "heterodimer_tm": primer3_metric(result, "tm"),
                "heterodimer_dg_kcal": dg_kcal(result),
                "heterodimer_flag": dim_flag,
                "heterodimer_reason": dim_reason,
                "heterodimer_practical_reject": dim_reject,
                "heterodimer_flag_rank": flag_rank(dim_flag),
            }
        )
    return pd.DataFrame(rows)


def make_qc_thresholds_table():
    pd, _ = require_qc_dependencies()
    rows = [
        ["hairpin", "REVIEW", f"Tm >= {HAIRPIN_REVIEW_TM:g} degC OR deltaG <= {HAIRPIN_REVIEW_DG:g} kcal/mol"],
        ["hairpin", "STRONG_REVIEW", f"deltaG <= {HAIRPIN_STRONG_DG:g} kcal/mol OR Tm >= {HAIRPIN_STRONG_TM:g} degC with deltaG <= {HAIRPIN_STRONG_TM_DG:g} kcal/mol"],
        ["self_dimer", "REVIEW", f"deltaG <= {SELF_DIMER_REVIEW_DG:g} kcal/mol"],
        ["self_dimer", "STRONG_REVIEW", f"deltaG <= {SELF_DIMER_STRONG_DG:g} kcal/mol"],
        ["self_dimer", "practical_reject_note", f"deltaG <= {SELF_DIMER_PRACTICAL_REJECT_DG:g} kcal/mol"],
        ["same_pair_P1_P2", "REVIEW", f"deltaG <= {SAME_PAIR_REVIEW_DG:g} kcal/mol"],
        ["same_pair_P1_P2", "STRONG_REVIEW", f"deltaG <= {SAME_PAIR_STRONG_DG:g} kcal/mol OR Tm >= {SAME_PAIR_STRONG_TM:g} degC with deltaG <= {SAME_PAIR_STRONG_TM_DG:g} kcal/mol"],
        ["cross_pair", "REVIEW", f"deltaG <= {CROSS_DIMER_REVIEW_DG:g} kcal/mol"],
        ["cross_pair", "STRONG_REVIEW", f"deltaG <= {CROSS_DIMER_STRONG_DG:g} kcal/mol"],
        ["cross_pair", "practical_reject_note", f"deltaG <= {CROSS_DIMER_PRACTICAL_REJECT_DG:g} kcal/mol"],
    ]
    return pd.DataFrame(rows, columns=["qc_type", "flag", "rule"])


def make_qc_summary(single_qc, pair_qc, input_path: Path):
    pd, _ = require_qc_dependencies()
    rows = [
        ["input_file", str(input_path)],
        ["total_oligos", int(len(single_qc))],
        ["total_pairwise_heterodimers", int(len(pair_qc))],
        ["hairpin_PASS", int((single_qc["hairpin_flag"] == "PASS").sum())],
        ["hairpin_REVIEW", int((single_qc["hairpin_flag"] == "REVIEW").sum())],
        ["hairpin_STRONG_REVIEW", int((single_qc["hairpin_flag"] == "STRONG_REVIEW").sum())],
        ["self_dimer_PASS", int((single_qc["self_dimer_flag"] == "PASS").sum())],
        ["self_dimer_REVIEW", int((single_qc["self_dimer_flag"] == "REVIEW").sum())],
        ["self_dimer_STRONG_REVIEW", int((single_qc["self_dimer_flag"] == "STRONG_REVIEW").sum())],
        ["heterodimer_PASS", int((pair_qc["heterodimer_flag"] == "PASS").sum())],
        ["heterodimer_REVIEW", int((pair_qc["heterodimer_flag"] == "REVIEW").sum())],
        ["heterodimer_STRONG_REVIEW", int((pair_qc["heterodimer_flag"] == "STRONG_REVIEW").sum())],
        ["same_pair_P1_P2_heterodimers_flagged", int(((pair_qc["dimer_class"] == "same_pair_P1_P2") & (pair_qc["heterodimer_flag"] != "PASS")).sum())],
        ["cross_pair_heterodimers_flagged", int(((pair_qc["dimer_class"] == "cross_pair") & (pair_qc["heterodimer_flag"] != "PASS")).sum())],
        ["single_oligo_practical_reject_notes", int((single_qc["hairpin_practical_reject"] | single_qc["self_dimer_practical_reject"]).sum())],
        ["heterodimer_practical_reject_notes", int(pair_qc["heterodimer_practical_reject"].sum())],
    ]
    return pd.DataFrame(rows, columns=["metric", "value"])


def write_qc_workbook(output_path: Path, single_qc, pair_qc, summary) -> None:
    flagged_single = single_qc[
        (single_qc["hairpin_flag"] != "PASS") | (single_qc["self_dimer_flag"] != "PASS")
    ].copy()
    flagged_single = flagged_single.sort_values(
        by=["single_max_flag_rank", "self_dimer_practical_reject", "hairpin_practical_reject", "self_dimer_dg_kcal", "hairpin_dg_kcal"],
        ascending=[False, False, False, True, True],
        na_position="last",
    )

    flagged_pairs = pair_qc[pair_qc["heterodimer_flag"] != "PASS"].copy()
    flagged_pairs = flagged_pairs.sort_values(
        by=["heterodimer_flag_rank", "heterodimer_practical_reject", "dimer_class", "heterodimer_dg_kcal"],
        ascending=[False, False, True, True],
        na_position="last",
    )

    try:
        with _pd.ExcelWriter(output_path) as writer:
            summary.to_excel(writer, sheet_name="summary", index=False)
            make_qc_thresholds_table().to_excel(writer, sheet_name="thresholds", index=False)
            single_qc.drop(columns=["single_max_flag_rank"]).to_excel(writer, sheet_name="single_oligo_QC", index=False)
            pair_qc.drop(columns=["heterodimer_flag_rank"]).to_excel(writer, sheet_name="all_pair_heterodimers", index=False)
            flagged_single.drop(columns=["single_max_flag_rank"]).to_excel(writer, sheet_name="flagged_single_oligos", index=False)
            flagged_pairs.drop(columns=["heterodimer_flag_rank"]).to_excel(writer, sheet_name="flagged_heterodimers", index=False)
    except ImportError as exc:  # pragma: no cover - environment-specific
        raise RuntimeError(
            "Writing the QC workbook requires openpyxl. Install it with: "
            "python -m pip install openpyxl, or rerun with --no-oligo-qc."
        ) from exc


def qc_summary_to_dict(summary) -> dict[str, Any]:
    return {str(row["metric"]): json_ready(row["value"]) for _, row in summary.iterrows()}


def run_oligo_structure_qc(input_path: Path, output_path: Path) -> dict[str, Any]:
    """Run final-oligo hairpin/self-dimer/heterodimer QC and write xlsx."""
    if output_path.suffix.lower() not in {".xlsx", ".xlsm"}:
        raise RuntimeError("Oligo QC output must be an .xlsx file.")
    output_path.parent.mkdir(parents=True, exist_ok=True)

    df = read_oligo_table(input_path)
    single_qc = make_single_qc(df)
    pair_qc = make_pair_qc(df)
    summary = make_qc_summary(single_qc, pair_qc, input_path)
    write_qc_workbook(output_path, single_qc, pair_qc, summary)

    summary_dict = qc_summary_to_dict(summary)
    summary_dict["output_file"] = str(output_path)
    return summary_dict


def _looks_like_no_candidate_probe_output(text: str) -> bool:
    """Recognize small status-only files emitted when no probes are found.

    HCRProbeDesign releases do not all represent a zero-result tier
    identically.  Some create an empty file, some write only a header, and some
    write a short human-readable status line.  These forms are biologically
    empty results, not malformed design tables, during adaptive auto-curation.
    """
    normalized = " ".join(text.replace("\ufeff", "").lower().split())
    if not normalized:
        return True
    markers = (
        "no probes",
        "no probe",
        "no candidates",
        "no candidate",
        "0 probes",
        "zero probes",
        "no passing",
        "nothing to write",
    )
    return any(marker in normalized for marker in markers)


_PROBE_COORDINATE_SUFFIX = re.compile(r"^(?P<prefix>.*?)(?P<start>\d+)-(?P<end>\d+)$")


def normalize_probe_name_coordinates(name: Any, start: int, end: int, length: int) -> str:
    """Keep a coordinate suffix consistent with the package's inclusive fields.

    Some HCRProbeDesign releases write a name such as ``probe_101-152`` for a
    52-nt probe starting at 101.  The numeric fields use 1-based inclusive
    coordinates, so the corresponding suffix is ``101-152`` only when the
    upstream convention is inclusive and ``101-153`` when it is exclusive at
    the right edge.  Normalize only the unambiguous exclusive form; arbitrary
    user labels are preserved.
    """
    text = str(name or "").strip()
    match = _PROBE_COORDINATE_SUFFIX.fullmatch(text)
    if not match:
        return text
    try:
        name_start = int(match.group("start"))
        name_end = int(match.group("end"))
    except ValueError:
        return text
    if name_start != int(start) or name_end != int(end) + 1:
        return text
    return f"{match.group('prefix')}{int(start)}-{int(end)}"


def read_probe_rows(tsv: Path, *, allow_empty: bool = False) -> list[dict[str, Any]]:
    """Read an HCRProbeDesign probe table.

    In adaptive auto-curation, a successful designProbes run may legitimately
    produce a zero-byte, whitespace-only, header-only, or status-only table when
    a tier has no passing candidates.  ``allow_empty=True`` represents those
    outcomes as an empty list so later relaxed tiers can still run.  Non-empty
    tables with genuine data but missing required columns remain hard errors.
    """
    if not tsv.exists():
        raise RuntimeError(f"Expected HCRProbeDesign output was not created: {tsv}")

    raw_text = tsv.read_text(errors="replace")
    if not raw_text.strip():
        if allow_empty:
            return []
        raise RuntimeError(f"{tsv} is empty or lacks expected HCRProbeDesign columns.")

    # Remove blank lines and UTF-8 BOMs before parsing.  A few upstream versions
    # produce an otherwise valid header preceded by an empty line or BOM.
    cleaned_lines = [
        line.lstrip("\ufeff")
        for line in raw_text.splitlines()
        if line.strip()
    ]
    cleaned_text = "\n".join(cleaned_lines) + ("\n" if cleaned_lines else "")

    reader = csv.DictReader(cleaned_text.splitlines(), delimiter="\t")
    fieldnames = {str(name).strip() for name in (reader.fieldnames or []) if name is not None}
    rows = [
        row
        for row in reader
        if any(str(value or "").strip() for value in row.values())
    ]

    required = {"start", "length", "GC", "Tm", "dTm", "GibbsFE"}
    if not required.issubset(fieldnames):
        if allow_empty and _looks_like_no_candidate_probe_output(raw_text):
            return []
        raise RuntimeError(
            f"{tsv} is non-empty but lacks expected HCRProbeDesign columns: "
            + ", ".join(sorted(required - fieldnames))
        )
    if not rows:
        if allow_empty:
            return []
        raise RuntimeError(f"{tsv} contains headers but no probe rows.")

    for row in rows:
        row["start"] = int(row["start"])
        row["length"] = int(row["length"])
        row["end"] = row["start"] + row["length"] - 1
        if "name" in row:
            row["name"] = normalize_probe_name_coordinates(
                row["name"], row["start"], row["end"], row["length"]
            )
        for field in ("GC", "Tm", "dTm", "GibbsFE"):
            row[field] = float(row[field])
    return sorted(rows, key=lambda row: row["start"])


def _annotation_label_mode(
    *,
    kind: str,
    width_nt: float,
    label_width_nt: float,
    premrna_map: bool,
) -> str:
    """Choose the placement of one transcript-feature label.

    Labels are evaluated independently so a narrow feature cannot move every
    other annotation label into external callout lanes.  Intronic Pre-mRNA
    maps omit exon/intron labels that cannot fit in their own region; mature
    transcript feature labels may use an external callout only when their
    individual region is too narrow.
    """
    fits_inside = width_nt >= label_width_nt + 8.0
    if fits_inside:
        return "inside"
    if premrna_map and kind in {"exon", "intron"}:
        return "omit"
    return "outside"


_PROBE_LABEL_SIMILARITY_RATIO = 1.15


def _uniform_probe_label_mode(
    box_widths_px: list[float],
    label_widths_px: list[float],
) -> str | None:
    """Return one label mode for similarly sized probe boxes.

    A ``None`` result means the boxes differ materially in width and the
    caller may make an individual fit decision.  When widths are similar,
    every label uses the same inside/outside mode to avoid a visually mixed
    map.
    """
    if not box_widths_px or len(box_widths_px) != len(label_widths_px):
        return None
    if max(box_widths_px) / max(min(box_widths_px), 1e-9) > _PROBE_LABEL_SIMILARITY_RATIO:
        return None
    all_fit_inside = all(
        box_width >= label_width + 5.0
        for box_width, label_width in zip(box_widths_px, label_widths_px)
    )
    return "inside" if all_fit_inside else "outside"


def _grouped_probe_label_modes(
    box_widths_px: list[float],
    label_widths_px: list[float],
    *,
    similarity_ratio: float = _PROBE_LABEL_SIMILARITY_RATIO,
) -> list[str | None]:
    """Choose one label mode for each group of similarly sized probe boxes.

    A map can contain several tile lengths. Applying one decision to the
    whole map makes a short tile alter the label placement of a much larger
    group, while deciding each box independently creates a mixed appearance
    among equal-length probes. Width groups are built from rendered pixel
    widths, so the policy follows the actual figure rather than a nucleotide
    length assumption. The default 15% width tolerance keeps adjacent
    curation tiers such as 50 and 52 nt together while separating materially
    shorter tiers such as 38 or 40 nt.
    """
    if len(box_widths_px) != len(label_widths_px):
        return [None] * len(box_widths_px)
    if not box_widths_px:
        return []
    if similarity_ratio <= 1.0:
        raise ValueError("similarity_ratio must be greater than 1")

    ordered = sorted(range(len(box_widths_px)), key=lambda index: box_widths_px[index])
    groups: list[list[int]] = []
    current: list[int] = []
    current_min = 0.0
    for index in ordered:
        width = max(float(box_widths_px[index]), 1e-9)
        if not current:
            current = [index]
            current_min = width
        elif width / current_min <= similarity_ratio:
            current.append(index)
        else:
            groups.append(current)
            current = [index]
            current_min = width
    if current:
        groups.append(current)

    modes: list[str | None] = [None] * len(box_widths_px)
    for group in groups:
        mode = _uniform_probe_label_mode(
            [box_widths_px[index] for index in group],
            [label_widths_px[index] for index in group],
        )
        for index in group:
            modes[index] = mode
    return modes


def plot_probe_map(
    rows: list[dict[str, Any]],
    transcript_length: int,
    title: str,
    output: Path,
    *,
    theme: str = "pastel",
    color_by: str = "gc",
    dpi: int = 300,
    show_labels: bool = True,
    transcript_annotation: dict[str, Any] | None = None,
) -> None:
    """Draw a publication-style map for mature or full Pre-mRNA targets."""
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.cm import ScalarMappable
        from matplotlib.colors import LinearSegmentedColormap, Normalize
        from matplotlib.patches import Patch, Rectangle
    except ImportError as exc:
        raise RuntimeError("Plotting requires matplotlib: pip install matplotlib") from exc

    if transcript_length < 1:
        raise ValueError("Transcript length must be positive.")
    if not rows:
        raise ValueError("No probe rows were supplied for plotting.")

    themes = {
        "pastel": {
            "background": "#FCFCFD",
            "text": "#334E5C",
            "muted": "#6F818A",
            "track": "#DDE8ED",
            "grid": "#D9E2E7",
            "probe_edge": "#8FA7B3",
            "palette": ["#A8DADC", "#B8E0D2", "#CDB4DB", "#F6D6AD", "#F4B6C2"],
        },
        "minimal": {
            "background": "#FFFFFF",
            "text": "#2F3E46",
            "muted": "#66747C",
            "track": "#D9E0E4",
            "grid": "#E5EAED",
            "probe_edge": "#8FA7B3",
            "palette": ["#B7C9D6", "#BFD8C1", "#D2C4DE", "#E6CFB3", "#DDBFC6"],
        },
    }
    style = themes[theme]

    probe_y = 0.90
    transcript_y = 0.31
    figure_width = max(10.0, min(16.0, transcript_length / 105.0))
    annotation_regions = (transcript_annotation or {}).get("regions", []) or []
    premrna_map = any(str(region.get("kind")) in {"exon", "intron"} for region in annotation_regions)
    figure_height = 4.8 if annotation_regions else 3.9

    fig, ax = plt.subplots(figsize=(figure_width, figure_height), facecolor=style["background"])
    ax.set_facecolor(style["background"])
    ax.set_xlim(0, transcript_length + 1)
    ax.set_ylim(0.0, 1.48 if show_labels else 1.24)

    values: list[float] = []
    scalar_mappable = None
    if color_by == "gc":
        values = [float(row["GC"]) for row in rows]
        legend_label = "Probe GC (%)"
    elif color_by == "dtm":
        values = [float(row["dTm"]) for row in rows]
        legend_label = "Probe-pair dTm (°C)"
    else:
        legend_label = ""

    if values:
        low, high = min(values), max(values)
        if low == high:
            low -= 0.5
            high += 0.5
        cmap = LinearSegmentedColormap.from_list(
            "hcr_pastel", ["#F6D6AD", "#B8E0D2", "#A8DADC", "#CDB4DB"]
        )
        norm = Normalize(vmin=low, vmax=high)
        scalar_mappable = ScalarMappable(norm=norm, cmap=cmap)
        scalar_mappable.set_array([])

    # Render the target as one uninterrupted coordinate backbone. Mature
    # targets use UTR/CDS regions; Pre-mRNA targets use exon/intron regions.
    # Coordinates are 1-based and inclusive; half-base boundaries preserve the
    # full border of a probe or annotation segment ending at transcript_length.
    transcript_left = 0.5
    transcript_right = transcript_length + 0.5
    region_colors = {
        "5UTR": "#E6E9EC",
        "CDS": "#C3CBD1",
        "3UTR": "#E6E1DB",
        "noncoding": "#D9DEDC",
        "partial_or_unannotated": "#ECEDEF",
        "exon": "#D6DCE2",
        "intron": "#CFE8D5",
        "separator": "#F2F4F5",
    }
    region_labels = {
        "5UTR": "5′ UTR",
        "CDS": "CDS",
        "3UTR": "3′ UTR",
        "noncoding": "non-coding transcript",
        "partial_or_unannotated": "partial / unannotated",
        "exon": "exon",
        "intron": "intron",
        "separator": "spacer",
    }
    # Use the rendered font width rather than a fixed nucleotide threshold.
    # This keeps labels adaptive across short and very long targets and avoids
    # placing a label inside a box that cannot actually contain it.
    label_width_cache: dict[tuple[str, float, str], float] = {}

    def label_width_nt(text: str, *, fontsize: float, fontweight: str = "normal") -> float:
        cache_key = (text, float(fontsize), fontweight)
        if cache_key in label_width_cache:
            return label_width_cache[cache_key]
        temporary = ax.text(
            0,
            0,
            text,
            fontsize=fontsize,
            fontweight=fontweight,
            alpha=0.0,
            clip_on=False,
        )
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        width_px = temporary.get_window_extent(renderer=renderer).width
        temporary.remove()
        x0 = ax.transData.transform((0, 0))[0]
        x1 = ax.transData.transform((1, 0))[0]
        pixels_per_nt = max(abs(x1 - x0), 1e-9)
        width_nt = width_px / pixels_per_nt
        label_width_cache[cache_key] = width_nt
        return width_nt

    if annotation_regions:
        track_height = 0.105
        legend_kinds: list[str] = []
        rendered_regions: list[dict[str, Any]] = []
        for region in annotation_regions:
            kind = str(region.get("kind", "partial_or_unannotated"))
            start = max(1, int(region.get("start", 1)))
            end = min(transcript_length, int(region.get("end", transcript_length)))
            if end < start:
                continue
            width = end - start + 1
            ax.add_patch(
                Rectangle(
                    (start - 0.5, transcript_y - track_height / 2),
                    width,
                    track_height,
                    linewidth=0.6,
                    edgecolor=style["probe_edge"],
                    facecolor=region_colors.get(kind, style["track"]),
                    zorder=2,
                )
            )
            if kind not in legend_kinds and kind != "separator":
                legend_kinds.append(kind)
            region_center = start + (width - 1) / 2
            if premrna_map and kind in {"exon", "intron"}:
                number_key = "exon_number" if kind == "exon" else "intron_number"
                number = region.get(number_key)
                # Use the informative numbered label.  If the complete label
                # cannot fit inside its own interval, the region is omitted
                # below rather than receiving an external callout.
                label = f"{kind} {int(number)}" if number not in {None, ""} else kind
            else:
                label = region_labels.get(kind, str(region.get("label") or kind))
            label_width = label_width_nt(
                label,
                fontsize=7.4,
                fontweight="bold" if kind == "CDS" else "normal",
            )
            # Annotation labels are evaluated per region.  One short exon or
            # intron must not force every other feature label outside the
            # transcript track.  Intronic Pre-mRNA maps intentionally omit a
            # label when the complete numbered name cannot fit in its own
            # region; external exon/intron callouts make dense maps unreadable.
            label_mode = _annotation_label_mode(
                kind=kind,
                width_nt=width,
                label_width_nt=label_width,
                premrna_map=premrna_map,
            )
            rendered_regions.append(
                {
                    "kind": kind,
                    "center": region_center,
                    "width": width,
                    "label": label,
                    "label_width": label_width,
                    "label_mode": label_mode,
                }
            )

        annotation_label_lanes = [0.45, 0.51, 0.57, 0.63, 0.69, 0.75]
        annotation_lane_right = [-float("inf")] * len(annotation_label_lanes)
        for item in rendered_regions:
            kind = str(item["kind"])
            label = str(item["label"])
            center = float(item["center"])
            label_width = float(item["label_width"])
            label_mode = str(item["label_mode"])
            if label_mode == "omit":
                continue
            if label_mode == "inside":
                ax.text(
                    center,
                    transcript_y,
                    label,
                    ha="center",
                    va="center",
                    fontsize=7.4,
                    fontweight="bold" if kind == "CDS" else "normal",
                    color=style["text"],
                    zorder=3,
                )
                continue

            lane = next(
                (
                    index
                    for index, right in enumerate(annotation_lane_right)
                    if center - label_width / 2.0 >= right + 2.0
                ),
                None,
            )
            if lane is None:
                # Keep the label visible even on extremely dense maps.  The
                # lane assignment remains deterministic and the leader line
                # still identifies the corresponding annotation box.
                lane = min(range(len(annotation_lane_right)), key=annotation_lane_right.__getitem__)
            label_y = annotation_label_lanes[lane]
            annotation_lane_right[lane] = center + label_width / 2.0
            ax.annotate(
                label,
                xy=(center, transcript_y + track_height / 2),
                xytext=(center, label_y),
                ha="center",
                va="center",
                fontsize=7.2,
                fontweight="bold" if kind == "CDS" else "normal",
                color=style["text"],
                arrowprops={
                    "arrowstyle": "-",
                    "linewidth": 0.55,
                    "color": style["probe_edge"],
                    "shrinkA": 1.5,
                    "shrinkB": 1.5,
                },
                clip_on=False,
                zorder=4,
            )
        region_counts = selected_probe_region_counts(annotate_probe_rows(rows, transcript_annotation))
        legend_handles = []
        for kind in legend_kinds:
            count = 0
            if kind == "5UTR":
                count = region_counts.get("5UTR", 0)
            elif kind == "CDS":
                count = region_counts.get("CDS", 0)
            elif kind == "3UTR":
                count = region_counts.get("3UTR", 0)
            elif kind == "noncoding":
                count = region_counts.get("noncoding", 0)
            elif kind == "intron":
                count = sum(
                    value for key, value in region_counts.items()
                    if str(key).startswith("intron_") and str(key) not in {"intron_boundary", "intron_separator"}
                )
            label = region_labels.get(kind, kind)
            if kind in {"5UTR", "CDS", "3UTR", "noncoding", "intron"}:
                unit = "probe" if count == 1 else "probes"
                label = f"{label} ({count} {unit})"
            legend_handles.append(Patch(facecolor=region_colors.get(kind, style["track"]), edgecolor=style["probe_edge"], label=label))
        if any("boundary" in str(row.get("transcript_region", "")) for row in rows):
            boundary_count = sum(1 for row in rows if "boundary" in str(row.get("transcript_region", "")))
            legend_handles.append(Patch(facecolor="none", edgecolor="none", label=f"boundary-spanning ({boundary_count} {'probe' if boundary_count == 1 else 'probes'})"))
        if legend_handles:
            ax.legend(
                handles=legend_handles,
                loc="upper right",
                bbox_to_anchor=(1.0, 1.01),
                frameon=False,
                fontsize=7.4,
                ncol=min(4, len(legend_handles)),
                handlelength=1.4,
                columnspacing=1.0,
            )
    else:
        ax.plot(
            [transcript_left, transcript_right],
            [transcript_y, transcript_y],
            color=style["track"],
            linewidth=11,
            solid_capstyle="round",
            zorder=2,
        )

    label_lanes = [1.04, 1.16, 1.28, 1.40]
    probe_number_fontsize = 7.2
    label_last_centers: list[tuple[float, float]] = [(-float("inf"), 0.0) for _ in label_lanes]
    fig.canvas.draw()
    probe_measurements: list[tuple[float, float, float]] = []
    for number, row in enumerate(rows, start=1):
        center = row["start"] + (row["length"] - 1) / 2
        box_left_px = ax.transData.transform((row["start"] - 0.5, probe_y))[0]
        box_right_px = ax.transData.transform((row["end"] + 0.5, probe_y))[0]
        box_width_px = abs(box_right_px - box_left_px)
        number_width_nt = label_width_nt(
            str(number), fontsize=probe_number_fontsize, fontweight="bold"
        )
        number_width_px = abs(
            ax.transData.transform((number_width_nt, probe_y))[0]
            - ax.transData.transform((0, probe_y))[0]
        )
        probe_measurements.append((center, box_width_px, number_width_px))
    box_widths = [measurement[1] for measurement in probe_measurements]
    number_widths = [measurement[2] for measurement in probe_measurements]
    probe_label_modes = _grouped_probe_label_modes(box_widths, number_widths)
    omitted_labels = 0
    for number, (row, measurement) in enumerate(zip(rows, probe_measurements), start=1):
        y = probe_y
        center, box_width_px, number_width_px = measurement
        if color_by == "gc":
            facecolor = scalar_mappable.to_rgba(float(row["GC"]))
        elif color_by == "dtm":
            facecolor = scalar_mappable.to_rgba(float(row["dTm"]))
        else:
            facecolor = style["palette"][(number - 1) % len(style["palette"])]

        patch = Rectangle(
            (row["start"] - 0.5, y - 0.13),
            row["length"],
            0.26,
            linewidth=0.8,
            edgecolor=style["probe_edge"],
            facecolor=facecolor,
            joinstyle="round",
            zorder=3,
        )
        ax.add_patch(patch)
        number_text = str(number)
        can_place_inside = box_width_px >= number_width_px + 5.0
        label_mode = probe_label_modes[number - 1] if number - 1 < len(probe_label_modes) else None
        place_inside = show_labels and (
            label_mode == "inside"
            or (label_mode is None and can_place_inside)
        )
        if place_inside:
            ax.text(
                center,
                y,
                number_text,
                ha="center",
                va="center",
                fontsize=probe_number_fontsize,
                fontweight="bold",
                color=style["text"],
                zorder=4,
            )
        elif show_labels:
            center_px = ax.transData.transform((center, y))[0]
            label_width_px = max(number_width_px, 8.0)
            chosen_lane = None
            for lane_index, (last_center, last_width) in enumerate(label_last_centers):
                if abs(center_px - last_center) >= (label_width_px + last_width) / 2.0 + 5.0:
                    chosen_lane = lane_index
                    break
            if chosen_lane is None:
                omitted_labels += 1
            else:
                label_last_centers[chosen_lane] = (center_px, label_width_px)
                label_y = label_lanes[chosen_lane]
                ax.plot([center, center], [y + 0.13, label_y - 0.035], color=style["probe_edge"], linewidth=0.45, zorder=3)
                ax.text(
                    center,
                    label_y,
                    str(number),
                    ha="center",
                    va="center",
                    fontsize=probe_number_fontsize,
                    fontweight="bold",
                    color=style["text"],
                    zorder=4,
                )

    first_start = min(row["start"] for row in rows)
    last_end = max(row["end"] for row in rows)
    selection_mode = str((transcript_annotation or {}).get("selection_mode") or "").casefold()
    target_description = (
        "full Pre-mRNA (intronic regions)"
        if premrna_map and selection_mode == "intronic"
        else "full Pre-mRNA"
        if premrna_map
        else "mature transcript"
    )
    subtitle = (
        f"{len(rows)} probe pairs   •   {transcript_length:,} nt {target_description}   •   "
        f"probe span {first_start:,}–{last_end:,} nt"
    )
    if omitted_labels:
        subtitle += f"   •   {omitted_labels} crowded label(s) omitted"

    # Keep a half-base of padding outside both transcript boundaries. The
    # padding prevents Matplotlib from clipping the stroke of terminal probes.
    ax.set_xlim(0, transcript_length + 1)
    ax.set_ylim(0.0, 1.48 if show_labels else 1.24)
    ax.set_yticks([])
    ax.set_xlabel("Position on full Pre-mRNA (nt)" if premrna_map else "Position on mature mRNA (nt)", fontsize=10, color=style["text"], labelpad=10)
    ax.set_title(title, loc="left", fontsize=15, fontweight="bold", color=style["text"], pad=27)
    ax.text(
        0.0,
        1.035,
        subtitle,
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=9.5,
        color=style["muted"],
    )
    ax.text(
        transcript_left,
        transcript_y - 0.16,
        "5′",
        ha="left",
        va="top",
        fontsize=9,
        color=style["muted"],
    )
    ax.text(
        transcript_right,
        transcript_y - 0.16,
        "3′",
        ha="right",
        va="top",
        fontsize=9,
        color=style["muted"],
    )

    ax.grid(False)
    ax.tick_params(
        axis="x",
        colors=style["muted"],
        labelsize=9,
        length=4,
        width=0.8,
        color=style["grid"],
    )
    # The plotting limits include non-coordinate padding at 0 and L+1; do not
    # expose those padding values as transcript-position tick labels.
    ax.set_xticks([tick for tick in ax.get_xticks() if 1 <= tick <= transcript_length])
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_color(style["grid"])

    if scalar_mappable is not None:
        colorbar = fig.colorbar(
            scalar_mappable,
            ax=ax,
            orientation="horizontal",
            fraction=0.055,
            pad=0.19,
            aspect=45,
        )
        colorbar.set_label(legend_label, fontsize=9, color=style["muted"], labelpad=5)
        colorbar.ax.tick_params(labelsize=8, colors=style["muted"], length=0)
        colorbar.outline.set_visible(False)

    fig.tight_layout()
    fig.savefig(output, dpi=dpi, bbox_inches="tight", facecolor=fig.get_facecolor())
    fig.savefig(output.with_suffix(".svg"), bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close(fig)


def write_summary(
    record: dict[str, Any],
    rows: list[dict[str, Any]],
    sequence_length: int,
    output: Path,
    run_metadata: dict[str, Any] | None = None,
) -> None:
    gaps = [rows[i]["start"] - rows[i - 1]["end"] - 1 for i in range(1, len(rows))]
    summary = {
        "gene_symbol": record.get("gene_symbol"),
        "accession": record["accession"],
        "title": record.get("title", ""),
        "target_type": record_target_type(record),
        "premrna_target": record.get("premrna_target"),
        "transcript_length_nt": sequence_length,
        "probe_pairs": len(rows),
        "design_result_status": "completed" if rows else "completed_no_candidates",
        "first_probe_start": rows[0]["start"] if rows else None,
        "last_probe_end": rows[-1]["end"] if rows else None,
        "largest_internal_gap_nt": max(gaps, default=0) if rows else None,
        "mean_GC_percent": round(sum(r["GC"] for r in rows) / len(rows), 2) if rows else None,
        "max_dTm_C": round(max(r["dTm"] for r in rows), 2) if rows else None,
        "transcript_resolution": record.get("transcript_resolution"),
        "transcript_annotation": record.get("transcript_annotation"),
        "selected_probe_region_counts": selected_probe_region_counts(rows),
    }
    if run_metadata is not None:
        summary["run_parameters"] = json_ready(run_metadata)
    write_json(output, summary)



# -----------------------------------------------------------------------------
# Smart auto-curation mode
# -----------------------------------------------------------------------------

AUTO_CURATE_PLANS: dict[str, list[dict[str, Any]]] = {
    "conservative": [
        {"tile_size": 52, "min_gc": 45.0, "max_gc": 55.0, "max_run_mismatches": 2, "min_gibbs": -70.0, "max_gibbs": -50.0, "target_gibbs": -60.0},
        {"tile_size": 52, "min_gc": 40.0, "max_gc": 60.0, "max_run_mismatches": 2, "min_gibbs": -70.0, "max_gibbs": -50.0, "target_gibbs": -60.0},
        {"tile_size": 50, "min_gc": 40.0, "max_gc": 60.0, "max_run_mismatches": 2, "min_gibbs": -70.0, "max_gibbs": -48.0, "target_gibbs": -58.0},
        {"tile_size": 48, "min_gc": 40.0, "max_gc": 60.0, "max_run_mismatches": 2, "min_gibbs": -70.0, "max_gibbs": -47.0, "target_gibbs": -56.0},
        {"tile_size": 46, "min_gc": 40.0, "max_gc": 65.0, "max_run_mismatches": 1, "min_gibbs": -70.0, "max_gibbs": -45.0, "target_gibbs": -55.0},
    ],
    "standard": [
        {"tile_size": 52, "min_gc": 45.0, "max_gc": 55.0, "max_run_mismatches": 2, "min_gibbs": -70.0, "max_gibbs": -50.0, "target_gibbs": -60.0},
        {"tile_size": 52, "min_gc": 40.0, "max_gc": 60.0, "max_run_mismatches": 2, "min_gibbs": -70.0, "max_gibbs": -50.0, "target_gibbs": -60.0},
        {"tile_size": 52, "min_gc": 35.0, "max_gc": 65.0, "max_run_mismatches": 2, "min_gibbs": -70.0, "max_gibbs": -50.0, "target_gibbs": -60.0},
        {"tile_size": 52, "min_gc": 30.0, "max_gc": 65.0, "max_run_mismatches": 1, "min_gibbs": -70.0, "max_gibbs": -50.0, "target_gibbs": -60.0},
        {"tile_size": 50, "min_gc": 40.0, "max_gc": 60.0, "max_run_mismatches": 2, "min_gibbs": -70.0, "max_gibbs": -48.0, "target_gibbs": -58.0},
        {"tile_size": 48, "min_gc": 40.0, "max_gc": 60.0, "max_run_mismatches": 2, "min_gibbs": -70.0, "max_gibbs": -47.0, "target_gibbs": -56.0},
        {"tile_size": 46, "min_gc": 40.0, "max_gc": 65.0, "max_run_mismatches": 1, "min_gibbs": -70.0, "max_gibbs": -45.0, "target_gibbs": -55.0},
        {"tile_size": 44, "min_gc": 40.0, "max_gc": 65.0, "max_run_mismatches": 1, "min_gibbs": -70.0, "max_gibbs": -43.0, "target_gibbs": -52.0},
        {"tile_size": 42, "min_gc": 40.0, "max_gc": 65.0, "max_run_mismatches": 1, "min_gibbs": -70.0, "max_gibbs": -40.0, "target_gibbs": -50.0},
    ],
    "deep": [
        {"tile_size": 52, "min_gc": 45.0, "max_gc": 55.0, "max_run_mismatches": 2, "min_gibbs": -70.0, "max_gibbs": -50.0, "target_gibbs": -60.0},
        {"tile_size": 52, "min_gc": 40.0, "max_gc": 60.0, "max_run_mismatches": 2, "min_gibbs": -70.0, "max_gibbs": -50.0, "target_gibbs": -60.0},
        {"tile_size": 52, "min_gc": 35.0, "max_gc": 65.0, "max_run_mismatches": 2, "min_gibbs": -70.0, "max_gibbs": -50.0, "target_gibbs": -60.0},
        {"tile_size": 52, "min_gc": 30.0, "max_gc": 65.0, "max_run_mismatches": 1, "min_gibbs": -70.0, "max_gibbs": -50.0, "target_gibbs": -60.0},
        {"tile_size": 50, "min_gc": 35.0, "max_gc": 65.0, "max_run_mismatches": 1, "min_gibbs": -70.0, "max_gibbs": -48.0, "target_gibbs": -58.0},
        {"tile_size": 48, "min_gc": 35.0, "max_gc": 65.0, "max_run_mismatches": 1, "min_gibbs": -70.0, "max_gibbs": -47.0, "target_gibbs": -56.0},
        {"tile_size": 46, "min_gc": 35.0, "max_gc": 65.0, "max_run_mismatches": 1, "min_gibbs": -70.0, "max_gibbs": -45.0, "target_gibbs": -55.0},
        {"tile_size": 44, "min_gc": 35.0, "max_gc": 65.0, "max_run_mismatches": 1, "min_gibbs": -70.0, "max_gibbs": -43.0, "target_gibbs": -52.0},
        {"tile_size": 42, "min_gc": 35.0, "max_gc": 65.0, "max_run_mismatches": 1, "min_gibbs": -70.0, "max_gibbs": -40.0, "target_gibbs": -50.0},
        {"tile_size": 40, "min_gc": 35.0, "max_gc": 65.0, "max_run_mismatches": 1, "min_gibbs": -70.0, "max_gibbs": -38.0, "target_gibbs": -46.0},
        {"tile_size": 38, "min_gc": 35.0, "max_gc": 65.0, "max_run_mismatches": 0, "min_gibbs": -70.0, "max_gibbs": -36.0, "target_gibbs": -44.0},
        {"tile_size": 36, "min_gc": 35.0, "max_gc": 65.0, "max_run_mismatches": 0, "min_gibbs": -70.0, "max_gibbs": -34.0, "target_gibbs": -42.0},
    ],
}



def designable_sequence_from_fasta(path: Path) -> str:
    """Return the concatenated sequence from a single-record FASTA for capacity estimates."""
    records = parse_fasta_records(path)
    if len(records) != 1:
        raise RuntimeError(f"Expected one FASTA record in {path}; found {len(records)}.")
    return records[0]["sequence"].upper().replace("U", "T")


def tile_capacity(sequence: str, tile_size: int) -> int:
    """Geometric non-overlapping capacity across contiguous unambiguous A/C/G/T segments."""
    if tile_size < 1:
        return 0
    return sum(len(segment) // tile_size for segment in re.findall(r"[ACGT]+", sequence.upper().replace("U", "T")))


def plan_min_tile_size(plan_name: str) -> int:
    return min(int(params["tile_size"]) for params in AUTO_CURATE_PLANS[plan_name])


def parameter_signature(params: dict[str, Any]) -> tuple[Any, ...]:
    """Stable identity for a design tier, independent of plan name or run number."""
    fields = (
        "tile_size",
        "min_gc",
        "max_gc",
        "max_run_mismatches",
        "min_gibbs",
        "max_gibbs",
        "target_gibbs",
    )
    return tuple(params.get(field) for field in fields)


def remaining_plan_tiers(
    plan_name: str, completed_signatures: set[tuple[Any, ...]]
) -> list[tuple[int, dict[str, Any]]]:
    """Return distinct, not-yet-run plan entries with their original plan index."""
    remaining: list[tuple[int, dict[str, Any]]] = []
    seen = set(completed_signatures)
    for plan_index, params in enumerate(AUTO_CURATE_PLANS[plan_name], start=1):
        signature = parameter_signature(params)
        if signature in seen:
            continue
        seen.add(signature)
        remaining.append((plan_index, params))
    return remaining


def reclassify_candidates(candidates: list[dict[str, Any]], stringency: str) -> None:
    """Reapply candidate acceptance policy without recalculating Primer3 metrics."""
    for candidate in candidates:
        reject_reason = pair_reject_reason(candidate, stringency)
        candidate["auto_qc_reject_reason"] = reject_reason
        candidate["auto_qc_status"] = "REJECT" if reject_reason else "PASS"
        candidate["auto_qc_stringency"] = stringency


def tile_size_counts(selected: list[dict[str, Any]]) -> dict[str, int]:
    counts: dict[str, int] = {}
    for candidate in selected:
        key = str(int(candidate.get("length", 0)))
        counts[key] = counts.get(key, 0) + 1
    return dict(sorted(counts.items(), key=lambda item: int(item[0]), reverse=True))


COVERAGE_MIN_CANONICAL_CAPACITY = 10
COVERAGE_BIN_COUNT = 10
COVERAGE_MIN_SPAN_FRACTION = 0.85
COVERAGE_MAX_GAP_MULTIPLIER = 2.5
COVERAGE_MIN_OCCUPIED_BINS = 8


def coverage_is_enabled(policy: str, canonical_capacity: int) -> bool:
    """Enable coverage balancing only for transcripts supporting >=10 canonical tiles.

    Short transcripts remain strictly count/QC driven, even if shorter deep-plan
    tiles could raise their apparent capacity above eight.
    """
    return policy == "balanced" and canonical_capacity >= COVERAGE_MIN_CANONICAL_CAPACITY


def selected_coverage_metrics(selected: list[dict[str, Any]], sequence_length: int) -> dict[str, Any]:
    ordered = sorted(selected, key=lambda cand: (int(cand["start"]), int(cand["end"])))
    count = len(ordered)
    if count == 0 or sequence_length < 1:
        return {
            "selected_probe_pairs": count,
            "probe_span_fraction": 0.0,
            "maximum_untargeted_gap_nt": sequence_length if sequence_length > 0 else 0,
            "ideal_spacing_nt": None,
            "maximum_gap_to_ideal_spacing_ratio": None,
            "occupied_bins": 0,
            "total_bins": COVERAGE_BIN_COUNT,
            "mean_ideal_position_deviation_nt": None,
        }

    first_start = int(ordered[0]["start"])
    last_end = int(ordered[-1]["end"])
    gaps = [max(0, first_start - 1)]
    gaps.extend(
        max(0, int(ordered[index]["start"]) - int(ordered[index - 1]["end"]) - 1)
        for index in range(1, count)
    )
    gaps.append(max(0, int(sequence_length) - last_end))
    ideal_spacing = float(sequence_length) / float(count)
    midpoints = [(float(cand["start"]) + float(cand["end"])) / 2.0 for cand in ordered]
    ideal_midpoints = [
        (index + 0.5) * float(sequence_length) / float(count)
        for index in range(count)
    ]
    occupied_bins = {
        min(COVERAGE_BIN_COUNT - 1, max(0, int((midpoint - 1.0) * COVERAGE_BIN_COUNT / sequence_length)))
        for midpoint in midpoints
    }
    maximum_gap = max(gaps, default=0)
    return {
        "selected_probe_pairs": count,
        "first_probe_start": first_start,
        "last_probe_end": last_end,
        "probe_span_nt": max(0, last_end - first_start + 1),
        "probe_span_fraction": round(max(0.0, min(1.0, (last_end - first_start + 1) / float(sequence_length))), 6),
        "maximum_untargeted_gap_nt": int(maximum_gap),
        "ideal_spacing_nt": round(ideal_spacing, 3),
        "maximum_gap_to_ideal_spacing_ratio": round(maximum_gap / ideal_spacing, 6) if ideal_spacing else None,
        "occupied_bins": len(occupied_bins),
        "total_bins": COVERAGE_BIN_COUNT,
        "mean_ideal_position_deviation_nt": round(
            sum(abs(actual - ideal) for actual, ideal in zip(midpoints, ideal_midpoints)) / count,
            3,
        ),
    }


def coverage_review_reasons(
    selected: list[dict[str, Any]],
    *,
    sequence_length: int,
    enabled: bool,
) -> list[str]:
    if not enabled or len(selected) < COVERAGE_MIN_CANONICAL_CAPACITY:
        return []
    metrics = selected_coverage_metrics(selected, sequence_length)
    reasons: list[str] = []
    span_fraction = float(metrics.get("probe_span_fraction") or 0.0)
    gap_ratio = float(metrics.get("maximum_gap_to_ideal_spacing_ratio") or 0.0)
    occupied_bins = int(metrics.get("occupied_bins") or 0)
    required_bins = min(COVERAGE_MIN_OCCUPIED_BINS, len(selected), COVERAGE_BIN_COUNT)
    if span_fraction < COVERAGE_MIN_SPAN_FRACTION:
        reasons.append(
            f"probe span fraction {span_fraction:.3f} is below {COVERAGE_MIN_SPAN_FRACTION:.2f}"
        )
    if gap_ratio > COVERAGE_MAX_GAP_MULTIPLIER:
        reasons.append(
            f"maximum untargeted gap is {gap_ratio:.2f}× ideal spacing, above {COVERAGE_MAX_GAP_MULTIPLIER:.1f}×"
        )
    if occupied_bins < required_bins:
        reasons.append(
            f"probes occupy {occupied_bins}/{COVERAGE_BIN_COUNT} transcript bins; at least {required_bins} are preferred"
        )
    return reasons


def selected_set_reaches_clean_goal(
    selected: list[dict[str, Any]],
    goal: int,
    *,
    sequence_length: int,
    coverage_enabled: bool,
) -> bool:
    return (
        goal > 0
        and len(selected) >= goal
        and selected_dimer_review_count(selected) == 0
        and not coverage_review_reasons(
            selected,
            sequence_length=sequence_length,
            enabled=coverage_enabled,
        )
    )


def clone_args_with(args: argparse.Namespace, **updates: Any) -> argparse.Namespace:
    """Shallow-copy parsed arguments and override selected fields."""
    data = dict(vars(args))
    data.update(updates)
    return argparse.Namespace(**data)


def tier_label(index: int, params: dict[str, Any], *, initial: bool = False) -> str:
    prefix = "run_00_initial" if initial else f"run_{index:02d}"
    return (
        f"{prefix}_t{int(params['tile_size'])}_"
        f"gc{int(params['min_gc'])}_{int(params['max_gc'])}_"
        f"mrm{int(params['max_run_mismatches'])}_"
        f"tg{str(params['target_gibbs']).replace('-', 'm').replace('.', 'p')}"
    )


def run_design_tier(
    *,
    fasta: Path,
    run_dir: Path,
    target: str,
    args: argparse.Namespace,
    record: dict[str, Any],
    sequence_length: int,
) -> dict[str, Any]:
    """Run one design tier and its final-oligo QC."""
    run_dir.mkdir(parents=True, exist_ok=True)
    probes_tsv, idt_tsv, log_file, run_parameters_json, run_metadata = run_design_probes(
        fasta,
        run_dir,
        target,
        args,
        record,
        sequence_length,
    )
    rows = read_probe_rows(probes_tsv, allow_empty=True)
    qc_input = idt_tsv
    eligible_idt_csv: Path | None = None
    if record.get("source") == "ncbi_premrna":
        rows = filter_premrna_candidates(rows, record)
        run_metadata.setdefault("outputs", {})["candidate_idt_tsv"] = str(idt_tsv)
        if rows:
            eligible_idt_csv = run_dir / f"{target}_eligible_IDT_order.csv"
            write_idt_csv(final_idt_rows(rows, target, args.channel), eligible_idt_csv)
            qc_input = eligible_idt_csv
            run_metadata.setdefault("outputs", {})["eligible_idt_order_csv"] = str(eligible_idt_csv)
        else:
            run_metadata.setdefault("outputs", {})["eligible_idt_order_csv"] = None
    qc_xlsx = run_dir / f"{target}_oligo_structure_QC.xlsx"

    if rows:
        notify_progress("qc", "Running final oligo structure QC on eligible probe candidates…")
        qc_report = run_oligo_structure_qc(qc_input, qc_xlsx)
        design_result_status = "completed"
        run_metadata.setdefault("outputs", {})["oligo_qc_xlsx"] = str(qc_xlsx)
        run_metadata.setdefault("outputs", {})["oligo_qc_input"] = str(qc_input)
    else:
        qc_report = {
            "status": "skipped",
            "reason": "design tier returned no candidate probe pairs",
            "output_file": None,
        }
        design_result_status = "completed_no_candidates"
        run_metadata.setdefault("outputs", {})["oligo_qc_xlsx"] = None
        print("    No candidate probe pairs were returned; oligo QC was skipped and curation will continue.")

    run_metadata["status"] = design_result_status
    run_metadata["design_result"] = {
        "status": design_result_status,
        "probe_pairs": len(rows),
        "oligo_qc_status": qc_report.get("status", "completed"),
        "oligo_qc_reason": qc_report.get("reason"),
    }
    run_metadata["final_oligo_qc"] = qc_report
    write_json(run_parameters_json, run_metadata)
    summary = run_dir / f"{target}_summary.json"
    write_summary(record, rows, sequence_length, summary, run_metadata)
    return {
        "status": design_result_status,
        "run_dir": run_dir,
        "probes_tsv": probes_tsv,
        "idt_tsv": idt_tsv,
        "eligible_idt_csv": eligible_idt_csv,
        "qc_input": qc_input if rows else None,
        "log_file": log_file,
        "run_parameters_json": run_parameters_json,
        "summary_json": summary,
        "qc_xlsx": qc_xlsx if rows else None,
        "qc_report": qc_report,
        "run_metadata": run_metadata,
        "rows": rows,
    }


def oligo_single_qc_for_pair(p1: str, p2: str, pair_key: str) -> tuple[Any, Any, dict[str, Any]]:
    """Return single-oligo QC, same-pair heterodimer QC, and summary flags for one pair."""
    pd, _ = require_qc_dependencies()
    oligos = pd.DataFrame(
        [
            {"name": f"{pair_key}_P1", "sequence": clean_oligo_seq(p1)},
            {"name": f"{pair_key}_P2", "sequence": clean_oligo_seq(p2)},
        ]
    )
    single_qc = make_single_qc(oligos)
    pair_qc = make_pair_qc(oligos)
    p1_row = single_qc.iloc[0]
    p2_row = single_qc.iloc[1]
    same_pair_rows = pair_qc[pair_qc["dimer_class"].isin(["same_pair_P1_P2", "same_probe_region"])]

    if same_pair_rows.empty:
        same_flag = "PASS"
        same_found = False
        same_dg = None
        same_tm = None
        same_reason = "no same-pair heterodimer row"
        same_practical = False
    else:
        same = same_pair_rows.iloc[0]
        same_flag = str(same["heterodimer_flag"])
        same_found = bool(same["heterodimer_found"])
        same_dg = same["heterodimer_dg_kcal"]
        same_tm = same["heterodimer_tm"]
        same_reason = str(same["heterodimer_reason"])
        same_practical = bool(same["heterodimer_practical_reject"])

    # Keep both the compact QC flags and the raw values that triggered those flags.
    # These fields are exported to candidate_master.tsv, backup_candidates.tsv,
    # rejected_candidates.tsv, and final_selected_pairs.tsv so manual review does
    # not require opening the per-tier QC workbooks just to see Tm/deltaG values.
    flags = {
        "p1_hairpin_flag": str(p1_row["hairpin_flag"]),
        "p1_hairpin_found": bool(p1_row["hairpin_found"]),
        "p1_hairpin_tm": p1_row["hairpin_tm"],
        "p1_hairpin_dg_kcal": p1_row["hairpin_dg_kcal"],
        "p1_hairpin_reason": str(p1_row["hairpin_reason"]),
        "p1_hairpin_practical_reject": bool(p1_row["hairpin_practical_reject"]),
        "p2_hairpin_flag": str(p2_row["hairpin_flag"]),
        "p2_hairpin_found": bool(p2_row["hairpin_found"]),
        "p2_hairpin_tm": p2_row["hairpin_tm"],
        "p2_hairpin_dg_kcal": p2_row["hairpin_dg_kcal"],
        "p2_hairpin_reason": str(p2_row["hairpin_reason"]),
        "p2_hairpin_practical_reject": bool(p2_row["hairpin_practical_reject"]),
        "p1_self_dimer_flag": str(p1_row["self_dimer_flag"]),
        "p1_self_dimer_found": bool(p1_row["self_dimer_found"]),
        "p1_self_dimer_tm": p1_row["self_dimer_tm"],
        "p1_self_dimer_dg_kcal": p1_row["self_dimer_dg_kcal"],
        "p1_self_dimer_reason": str(p1_row["self_dimer_reason"]),
        "p1_self_dimer_practical_reject": bool(p1_row["self_dimer_practical_reject"]),
        "p2_self_dimer_flag": str(p2_row["self_dimer_flag"]),
        "p2_self_dimer_found": bool(p2_row["self_dimer_found"]),
        "p2_self_dimer_tm": p2_row["self_dimer_tm"],
        "p2_self_dimer_dg_kcal": p2_row["self_dimer_dg_kcal"],
        "p2_self_dimer_reason": str(p2_row["self_dimer_reason"]),
        "p2_self_dimer_practical_reject": bool(p2_row["self_dimer_practical_reject"]),
        "same_pair_heterodimer_flag": same_flag,
        "same_pair_heterodimer_found": same_found,
        "same_pair_heterodimer_tm": same_tm,
        "same_pair_heterodimer_dg_kcal": same_dg,
        "same_pair_heterodimer_reason": same_reason,
        "same_pair_heterodimer_practical_reject": same_practical,
        # Backward-compatible combined practical-reject aliases used by reject logic.
        "p1_practical_reject": bool(p1_row["hairpin_practical_reject"] or p1_row["self_dimer_practical_reject"]),
        "p2_practical_reject": bool(p2_row["hairpin_practical_reject"] or p2_row["self_dimer_practical_reject"]),
        "same_pair_practical_reject": same_practical,
    }
    return single_qc, pair_qc, flags


def pair_reject_reason(flags: dict[str, Any], stringency: str) -> str:
    """Return empty string if candidate is acceptable under stringency; otherwise reason."""
    if flags.get("p1_practical_reject") or flags.get("p2_practical_reject") or flags.get("same_pair_practical_reject"):
        return "practical reject note in final-oligo QC"
    if flags.get("p1_self_dimer_flag") == "STRONG_REVIEW" or flags.get("p2_self_dimer_flag") == "STRONG_REVIEW":
        return "strong self-dimer"
    if flags.get("same_pair_heterodimer_flag") == "STRONG_REVIEW":
        return "strong same-pair P1/P2 heterodimer"
    if stringency in {"balanced", "strict"}:
        if flags.get("p1_hairpin_flag") == "STRONG_REVIEW" or flags.get("p2_hairpin_flag") == "STRONG_REVIEW":
            return "strong hairpin"
    if stringency == "strict":
        if flags.get("p1_self_dimer_flag") != "PASS" or flags.get("p2_self_dimer_flag") != "PASS":
            return "strict mode excludes self-dimer review"
        if flags.get("same_pair_heterodimer_flag") != "PASS":
            return "strict mode excludes same-pair heterodimer review"
    return ""


def optional_float(value: Any) -> float | None:
    """Return a finite float or None for missing/NaN-like QC values."""
    if value is None:
        return None
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    if math.isnan(number) or math.isinf(number):
        return None
    return number



def hairpin_penalty(flags: dict[str, Any], prefix: str) -> float:
    """Return a curation penalty for one oligo hairpin flag.

    Mild Tm-only hairpin REVIEW calls are intentionally penalised only lightly.
    In manual curation, these have been less concerning than self-dimers or
    same-pair P1/P2 dimers when hairpin deltaG is weak. Hairpin REVIEW caused
    by more negative deltaG receives a larger, but still non-dimer, penalty.
    """
    flag = str(flags.get(f"{prefix}_hairpin_flag", "PASS"))
    if flag == "PASS":
        return 0.0
    if flag == "STRONG_REVIEW":
        return 200.0

    dg = optional_float(flags.get(f"{prefix}_hairpin_dg_kcal"))
    reason = str(flags.get(f"{prefix}_hairpin_reason", "")).lower()
    if flag == "REVIEW":
        # Tm-only REVIEW with weak deltaG is common and was not treated as a
        # major concern in the manually curated sox9/atp1a1 sets.
        if "tm >=" in reason and (dg is None or dg > HAIRPIN_REVIEW_DG):
            return 2.0
        if dg is not None and dg <= HAIRPIN_REVIEW_DG:
            return 35.0 + max(0.0, HAIRPIN_REVIEW_DG - dg) * 15.0
        return 8.0
    return 0.0


def candidate_qc_penalty(flags: dict[str, Any]) -> float:
    """Penalty used for reporting candidate score; reject logic is separate.

    Selection now uses a lexicographic dimer-first objective under the hood.
    This scalar score is retained for review tables and backup ranking, and it
    follows the same philosophy: self-dimer and same-pair P1/P2 dimer REVIEW
    calls dominate mild Tm-only hairpin REVIEW calls.
    """
    penalty = 0.0
    penalty += hairpin_penalty(flags, "p1")
    penalty += hairpin_penalty(flags, "p2")

    for prefix in ("p1", "p2"):
        flag = str(flags.get(f"{prefix}_self_dimer_flag", "PASS"))
        dg = optional_float(flags.get(f"{prefix}_self_dimer_dg_kcal"))
        if flag == "REVIEW":
            penalty += 500.0 + dimer_review_severity(dg, SELF_DIMER_REVIEW_DG) * 120.0
        elif flag == "STRONG_REVIEW":
            penalty += 1800.0 + dimer_review_severity(dg, SELF_DIMER_STRONG_DG) * 160.0

    same_flag = str(flags.get("same_pair_heterodimer_flag", "PASS"))
    same_dg = optional_float(flags.get("same_pair_heterodimer_dg_kcal"))
    if same_flag == "REVIEW":
        penalty += 800.0 + dimer_review_severity(same_dg, SAME_PAIR_REVIEW_DG) * 160.0
    elif same_flag == "STRONG_REVIEW":
        penalty += 2600.0 + dimer_review_severity(same_dg, SAME_PAIR_STRONG_DG) * 200.0
    return penalty


def dimer_review_severity(dg: float | None, threshold: float) -> float:
    """How far a dimer deltaG is beyond its REVIEW threshold.

    primer3 reports more stable dimers as more negative deltaG. A value of
    -9.05 kcal/mol just beyond a -9 threshold is less concerning than -11.8.
    Missing deltaG values get a small non-zero severity if the caller already
    knows a flag was present.
    """
    if dg is None:
        return 0.05
    return max(0.0, float(threshold) - float(dg))


def hairpin_review_profile(cand: dict[str, Any]) -> tuple[int, float, int]:
    """Return hairpin strong count, deltaG severity, and mild Tm-only count."""
    strong = 0
    dg_severity = 0.0
    mild_tm_only = 0
    for prefix in ("p1", "p2"):
        flag = str(cand.get(f"{prefix}_hairpin_flag", "PASS"))
        dg = optional_float(cand.get(f"{prefix}_hairpin_dg_kcal"))
        reason = str(cand.get(f"{prefix}_hairpin_reason", "")).lower()
        if flag == "STRONG_REVIEW":
            strong += 1
            if dg is not None:
                dg_severity += max(0.0, HAIRPIN_STRONG_DG - dg)
        elif flag == "REVIEW":
            if dg is not None and dg <= HAIRPIN_REVIEW_DG:
                dg_severity += max(0.0, HAIRPIN_REVIEW_DG - dg)
            elif "tm >=" in reason:
                mild_tm_only += 1
            else:
                mild_tm_only += 1
    return strong, dg_severity, mild_tm_only


def dimer_review_profile(cand: dict[str, Any]) -> tuple[int, int, float]:
    """Return self-dimer count, same-pair dimer count, and severity."""
    self_count = 0
    same_count = 0
    severity = 0.0

    for prefix in ("p1", "p2"):
        flag = str(cand.get(f"{prefix}_self_dimer_flag", "PASS"))
        if flag != "PASS":
            self_count += 1
            threshold = SELF_DIMER_STRONG_DG if flag == "STRONG_REVIEW" else SELF_DIMER_REVIEW_DG
            severity += dimer_review_severity(optional_float(cand.get(f"{prefix}_self_dimer_dg_kcal")), threshold)

    same_flag = str(cand.get("same_pair_heterodimer_flag", "PASS"))
    if same_flag != "PASS":
        same_count += 1
        threshold = SAME_PAIR_STRONG_DG if same_flag == "STRONG_REVIEW" else SAME_PAIR_REVIEW_DG
        # Same-pair dimers are weighted slightly more because P1/P2 interaction
        # directly competes with target binding.
        severity += 1.5 * dimer_review_severity(optional_float(cand.get("same_pair_heterodimer_dg_kcal")), threshold)

    return self_count, same_count, severity


def candidate_technical_score(cand: dict[str, Any]) -> float:
    """Soft design-quality tie-breaker used after QC priorities."""
    gc = float(cand.get("GC", 50.0))
    dtm = float(cand.get("dTm", 0.0))
    gibbs = float(cand.get("GibbsFE", cand.get("source_target_gibbs", -60.0)))
    target_gibbs = float(cand.get("source_target_gibbs", -60.0))
    tile_size = int(cand.get("length", cand.get("source_tile_size", 52)))
    source_tier = int(cand.get("source_tier", 0))
    return (
        source_tier * 100.0
        + max(0, 52 - tile_size) * 3.0
        + abs(gc - 50.0) * 2.0
        + abs(gibbs - target_gibbs) * 1.5
        + max(0.0, dtm - 5.0) * 8.0
    )


def candidate_selection_vector(cand: dict[str, Any]) -> tuple[float, ...]:
    """Lexicographic objective for one candidate.

    Candidate selection maximizes probe count first. For sets with the same
    count, summed vectors are minimized in this order:
      1. any selected self-dimer/same-pair dimer REVIEW count,
      2. same-pair dimer REVIEW count,
      3. actual dimer deltaG severity,
      4. strong/severe hairpin burden,
      5. mild Tm-only hairpin count,
      6. conventional design-quality score.

    This prevents a mild hairpin Tm-only REVIEW, shorter tile, or later-tier
    candidate from outweighing an avoidable dimer issue.
    """
    self_count, same_count, dimer_severity = dimer_review_profile(cand)
    hairpin_strong, hairpin_dg_severity, mild_tm_only = hairpin_review_profile(cand)
    total_dimer_count = self_count + same_count
    return (
        float(total_dimer_count),
        float(same_count),
        round(float(dimer_severity), 6),
        float(hairpin_strong),
        round(float(hairpin_dg_severity), 6),
        float(mild_tm_only),
        round(candidate_technical_score(cand), 6),
        float(cand.get("auto_score", 0.0)),
    )


def add_objective_tuple(a: tuple[float, ...], b: tuple[float, ...]) -> tuple[float, ...]:
    return tuple(x + y for x, y in zip(a, b))


def zero_objective_tuple() -> tuple[float, ...]:
    return (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)


def coverage_position_penalty(
    cand: dict[str, Any],
    *,
    position_index: int,
    selected_count: int,
    sequence_length: int,
    coverage_enabled: bool,
) -> float:
    """Return normalized deviation from the ideal evenly spaced position.

    The penalty is zero when coverage balancing is disabled.  It is evaluated
    only after all explicit QC components in the selection vector, so coverage
    never trades a cleaner same-count set for one with more QC flags.
    """
    if not coverage_enabled or selected_count < 1 or sequence_length < 1:
        return 0.0
    midpoint = (float(cand["start"]) + float(cand["end"])) / 2.0
    ideal_midpoint = (position_index - 0.5) * float(sequence_length) / float(selected_count)
    return round(abs(midpoint - ideal_midpoint) / float(sequence_length), 9)


def positioned_candidate_selection_vector(
    cand: dict[str, Any],
    *,
    position_index: int,
    selected_count: int,
    sequence_length: int,
    coverage_enabled: bool,
) -> tuple[float, ...]:
    """Return the candidate vector with a set-position coverage component.

    Candidate count is still maximized before this vector is considered.  The
    first six components are the existing dimer/hairpin QC priorities.  The
    coverage component is inserted only after those QC priorities and before
    conventional technical tie-breakers.
    """
    base = candidate_selection_vector(cand)
    coverage = coverage_position_penalty(
        cand,
        position_index=position_index,
        selected_count=selected_count,
        sequence_length=sequence_length,
        coverage_enabled=coverage_enabled,
    )
    return base[:6] + (coverage,) + base[6:]


def zero_positioned_objective_tuple() -> tuple[float, ...]:
    return (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)


def selected_set_objective(
    selected: list[dict[str, Any]],
    *,
    sequence_length: int,
    coverage_enabled: bool,
) -> tuple[float, ...]:
    ordered = sorted(selected, key=lambda cand: (int(cand["start"]), int(cand["end"])))
    objective = zero_positioned_objective_tuple()
    selected_count = len(ordered)
    for position_index, cand in enumerate(ordered, start=1):
        objective = add_objective_tuple(
            objective,
            positioned_candidate_selection_vector(
                cand,
                position_index=position_index,
                selected_count=selected_count,
                sequence_length=sequence_length,
                coverage_enabled=coverage_enabled,
            ),
        )
    return objective


def selected_dimer_review_count(selected: list[dict[str, Any]]) -> int:
    """Count selected candidates with any self-dimer or same-pair dimer flag."""
    count = 0
    for cand in selected:
        self_count, same_count, _ = dimer_review_profile(cand)
        count += self_count + same_count
    return count


def selected_set_is_dimer_clean(selected: list[dict[str, Any]], target_count: int) -> bool:
    """Return True when the selected set reaches target_count and has no dimer REVIEW flags."""
    return len(selected) >= target_count and selected_dimer_review_count(selected) == 0


def make_candidate_from_row(
    row: dict[str, Any],
    *,
    source_run: str,
    source_tier: int,
    source_params: dict[str, Any],
    target_label: str,
    stringency: str,
) -> dict[str, Any]:
    """Convert one HCRProbeDesign row into a scored auto-curation candidate."""
    start = int(row["start"])
    length = int(row["length"])
    end = start + length - 1
    pair_key = safe_name(f"{target_label}_{start}-{end}_t{length}")
    _, _, flags = oligo_single_qc_for_pair(str(row["P1"]), str(row["P2"]), pair_key)
    reject_reason = pair_reject_reason(flags, stringency)

    gc = float(row.get("GC", 50.0))
    dtm = float(row.get("dTm", 0.0))
    gibbs = float(row.get("GibbsFE", source_params.get("target_gibbs", -60.0)))
    target_gibbs = float(source_params.get("target_gibbs", -60.0))
    tile_size = int(length)

    score = 0.0
    score += source_tier * 100.0
    score += max(0, 52 - tile_size) * 3.0
    score += abs(gc - 50.0) * 2.0
    score += abs(gibbs - target_gibbs) * 1.5
    score += max(0.0, dtm - 5.0) * 8.0
    score += candidate_qc_penalty(flags)

    candidate = dict(row)
    candidate.update(
        {
            "candidate_id": safe_name(f"{source_run}_{start}_{end}_{tile_size}"),
            "start": start,
            "end": end,
            "length": tile_size,
            "source_run": source_run,
            "source_tier": source_tier,
            "source_tile_size": int(source_params.get("tile_size", tile_size)),
            "source_min_gc": float(source_params.get("min_gc", 0.0)),
            "source_max_gc": float(source_params.get("max_gc", 100.0)),
            "source_target_gibbs": target_gibbs,
            "auto_qc_reject_reason": reject_reason,
            "auto_qc_status": "REJECT" if reject_reason else "PASS",
            "auto_qc_stringency": stringency,
            "auto_score": round(score, 4),
            **flags,
        }
    )
    return candidate


def deduplicate_candidates(candidates: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Keep the best exact candidate among duplicate intervals/sequences."""
    best: dict[tuple[Any, ...], dict[str, Any]] = {}
    for cand in candidates:
        key = (
            int(cand["start"]),
            int(cand["end"]),
            clean_oligo_seq(cand.get("P1", "")),
            clean_oligo_seq(cand.get("P2", "")),
        )
        if key not in best or candidate_selection_vector(cand) < candidate_selection_vector(best[key]):
            best[key] = cand
    return list(best.values())


def maximum_nonoverlapping_count(candidates: list[dict[str, Any]], target_count: int) -> int:
    """Return the maximum interval count, capped at target_count.

    For unweighted interval scheduling, choosing the earliest finishing
    compatible interval is count-optimal.
    """
    count = 0
    previous_end: int | None = None
    for cand in sorted(candidates, key=lambda item: (int(item["end"]), int(item["start"]))):
        start = int(cand["start"])
        end = int(cand["end"])
        if previous_end is None or start > previous_end:
            count += 1
            previous_end = end
            if count >= target_count:
                return int(target_count)
    return count


def select_best_nonoverlapping(
    candidates: list[dict[str, Any]],
    target_count: int,
    blocked_ids: set[str] | None = None,
    *,
    sequence_length: int,
    coverage_enabled: bool,
    blocked_intervals: Sequence[tuple[int, int]] | None = None,
) -> list[dict[str, Any]]:
    """Select up to target_count non-overlapping candidates.

    Dynamic programming maximizes probe count first. Among same-count solutions,
    it uses a lexicographic dimer-first objective rather than a single blended
    scalar. This better matches conservative manual curation: avoid selected
    self-dimer/same-pair dimer REVIEW flags if an alternative with the same
    probe count exists, even if that alternative is from a relaxed tier, shorter,
    or has only mild hairpin Tm-only REVIEW flags.
    """
    blocked_ids = blocked_ids or set()
    blocked_intervals = tuple(
        (int(start), int(end))
        for start, end in (blocked_intervals or ())
        if int(end) >= int(start)
    )
    usable = [
        cand
        for cand in candidates
        if (
            not cand.get("auto_qc_reject_reason")
            and cand.get("candidate_id") not in blocked_ids
            and not any(
                int(cand["start"]) <= blocked_end
                and int(cand["end"]) >= blocked_start
                for blocked_start, blocked_end in blocked_intervals
            )
        )
    ]
    usable.sort(key=lambda c: (int(c["end"]), int(c["start"]), candidate_selection_vector(c)))
    n = len(usable)
    kmax = maximum_nonoverlapping_count(usable, int(target_count))
    if n == 0 or kmax < 1:
        return []

    ends = [int(c["end"]) for c in usable]
    prev_counts = [bisect_right(ends, int(c["start"]) - 1) for c in usable]

    # The final count is known before QC/coverage optimization.  This preserves
    # the count-first rule and lets the j-th selected probe be compared with the
    # j-th ideal transcript position using an additive dynamic-programming cost.
    dp: list[list[tuple[tuple[float, ...], tuple[int, ...]] | None]] = [[None] * (kmax + 1) for _ in range(n + 1)]
    dp[0][0] = (zero_positioned_objective_tuple(), tuple())

    for i in range(1, n + 1):
        cand = usable[i - 1]
        pcount = prev_counts[i - 1]
        for k in range(0, kmax + 1):
            best_state = dp[i - 1][k]
            if k > 0 and dp[pcount][k - 1] is not None:
                prev_objective, prev_indices = dp[pcount][k - 1]
                cand_objective = positioned_candidate_selection_vector(
                    cand,
                    position_index=k,
                    selected_count=kmax,
                    sequence_length=sequence_length,
                    coverage_enabled=coverage_enabled,
                )
                take_state = (add_objective_tuple(prev_objective, cand_objective), prev_indices + (i - 1,))
                if best_state is None or take_state[0] < best_state[0]:
                    best_state = take_state
            dp[i][k] = best_state

    chosen_state = dp[n][kmax]
    if chosen_state is None:
        return []
    return [usable[i] for i in chosen_state[1]]


def _premrna_intron_regions(record: dict[str, Any]) -> list[dict[str, Any]]:
    annotation = record.get("transcript_annotation") or {}
    return [
        region for region in annotation.get("regions", []) or []
        if str(region.get("kind", "")) == "intron"
        and region.get("intron_number") is not None
        and bool(region.get("eligible", True))
    ]


def premrna_intron_capacity(
    record: dict[str, Any], sequence: str, tile_size: int
) -> int:
    """Return the non-overlapping geometric capacity of eligible introns.

    Pre-mRNA coordinates are 1-based and inclusive on the assembled genomic
    transcript.  Capacity is therefore calculated from each eligible intron
    sequence independently, using :func:`tile_capacity`; a tile cannot cross
    an intron boundary or an ambiguous base.  This deliberately counts
    non-overlapping tiles rather than all possible sliding-window starts.
    """
    if tile_size < 1:
        return 0
    normalized_sequence = str(sequence or "").upper().replace("U", "T")
    capacity = 0
    for region in _premrna_intron_regions(record):
        try:
            start = int(region["start"])
            end = int(region["end"])
        except (KeyError, TypeError, ValueError):
            continue
        if start < 1 or end < start:
            continue
        intron_sequence = normalized_sequence[start - 1 : end]
        if len(intron_sequence) != end - start + 1:
            continue
        capacity += tile_capacity(intron_sequence, tile_size)
    return capacity


def _premrna_is_whole(record: dict[str, Any]) -> bool:
    target = record.get("premrna_target") or record.get("transcript_annotation") or {}
    return (
        str(target.get("selection_mode") or "").casefold() == "whole"
        or str(target.get("mode") or "").casefold() == "pre-mrna-whole"
    )


def filter_premrna_candidates(
    candidates: list[dict[str, Any]],
    record: dict[str, Any],
) -> list[dict[str, Any]]:
    """Keep candidates valid for the selected genomic Pre-mRNA target."""
    if _premrna_is_whole(record):
        annotation = record.get("transcript_annotation") or {}
        usable: list[dict[str, Any]] = []
        for candidate in candidates:
            try:
                start = int(candidate["start"])
                end = int(candidate.get("end", start + int(candidate["length"]) - 1))
            except (KeyError, TypeError, ValueError):
                continue
            if start < 1 or end > int(record.get("length") or 0) or end < start:
                continue
            # Whole-target designs intentionally allow exon, intron, and
            # boundary-spanning pairs. Annotation is still applied so maps,
            # reports, and downstream review can distinguish those cases.
            annotate_probe_with_transcript_region(candidate, annotation)
            usable.append(candidate)
        return usable
    regions = _premrna_intron_regions(record)
    if not regions:
        return []
    usable: list[dict[str, Any]] = []
    for candidate in candidates:
        try:
            start = int(candidate["start"])
            end = int(candidate.get("end", start + int(candidate["length"]) - 1))
        except (KeyError, TypeError, ValueError):
            continue
        region = next(
            (
                item for item in regions
                if start >= int(item["start"]) and end <= int(item["end"])
            ),
            None,
        )
        if region is None:
            continue
        candidate["premrna_intron_number"] = int(region["intron_number"])
        candidate["premrna_genomic_start"] = int(region.get("genomic_start", 0) or 0)
        candidate["premrna_genomic_end"] = int(region.get("genomic_end", 0) or 0)
        candidate["transcript_region"] = f"intron_{int(region['intron_number'])}"
        usable.append(candidate)
    return usable


def select_premrna_candidates(
    candidates: list[dict[str, Any]],
    target_count: int,
    blocked_ids: set[str] | None,
    *,
    sequence_length: int,
    record: dict[str, Any],
) -> list[dict[str, Any]]:
    """Select a target-sized pool with an equal, capacity-aware intron quota.

    This is intentionally a separate selector.  Mature-transcript selection
    remains governed by the established single-backbone dynamic program.
    """
    blocked_ids = blocked_ids or set()
    annotation = record.get("transcript_annotation") or {}
    intron_regions = _premrna_intron_regions(record)
    if not intron_regions:
        return select_best_nonoverlapping(
            candidates,
            target_count,
            blocked_ids,
            sequence_length=sequence_length,
            coverage_enabled=False,
        )

    usable_by_intron: dict[int, list[dict[str, Any]]] = {int(region["intron_number"]): [] for region in intron_regions}
    for candidate in filter_premrna_candidates(candidates, record):
        if candidate.get("auto_qc_reject_reason") or candidate.get("candidate_id") in blocked_ids:
            continue
        for region in intron_regions:
            if candidate.get("premrna_intron_number") == int(region["intron_number"]):
                usable_by_intron[int(region["intron_number"])].append(candidate)
                break

    region_numbers = [int(region["intron_number"]) for region in intron_regions]
    if target_count < 1:
        return []
    base, remainder = divmod(int(target_count), len(region_numbers))
    quotas = {number: base + (index < remainder) for index, number in enumerate(region_numbers)}
    selected: list[dict[str, Any]] = []
    for number in region_numbers:
        quota = quotas[number]
        if quota < 1:
            continue
        region = next(region for region in intron_regions if int(region["intron_number"]) == number)
        selected.extend(
            select_best_nonoverlapping(
                usable_by_intron[number],
                quota,
                blocked_ids | {str(item.get("candidate_id")) for item in selected},
                sequence_length=max(1, int(region["end"]) - int(region["start"]) + 1),
                coverage_enabled=False,
                blocked_intervals=[
                    (int(item["start"]), int(item["end"])) for item in selected
                ],
            )
        )

    if len(selected) < target_count:
        selected_ids = {str(item.get("candidate_id")) for item in selected}
        eligible_pool = filter_premrna_candidates(candidates, record)
        selected.extend(
            select_best_nonoverlapping(
                eligible_pool,
                target_count - len(selected),
                blocked_ids | selected_ids,
                sequence_length=sequence_length,
                coverage_enabled=False,
                blocked_intervals=[
                    (int(item["start"]), int(item["end"])) for item in selected
                ],
            )
        )
    return sorted(selected[:target_count], key=lambda item: int(item["start"]))


def enforce_premrna_nonoverlap(
    candidates: list[dict[str, Any]], sequence_length: int
) -> list[dict[str, Any]]:
    """Apply interval packing to a filtered one-pass Pre-mRNA candidate table.

    HCRProbeDesign normally emits non-overlapping native selections, but a
    filtered genomic scope can expose overlapping rows from upstream output.
    Keep the one-pass candidate request unchanged while ensuring that the
    order-ready Pre-mRNA table never contains overlapping intervals. Smart
    curation uses :func:`select_premrna_candidates`, which applies the same
    invariant across intron quota and fallback pools.
    """
    if not candidates:
        return []
    return sorted(
        select_best_nonoverlapping(
            candidates,
            len(candidates),
            set(),
            sequence_length=sequence_length,
            coverage_enabled=False,
        ),
        key=lambda item: int(item["start"]),
    )


def final_idt_rows(selected: list[dict[str, Any]], target_label: str, channel: str) -> list[dict[str, str]]:
    rows = []
    for number, cand in enumerate(sorted(selected, key=lambda c: int(c["start"])), start=1):
        start = int(cand["start"])
        end = int(cand["end"])
        base = safe_name(f"{target_label}_{channel}_{number:02d}_{start}-{end}")
        rows.append({"name": f"{base}_P1", "sequence": clean_oligo_seq(cand["P1"])})
        rows.append({"name": f"{base}_P2", "sequence": clean_oligo_seq(cand["P2"])})
    return rows


def write_idt_csv(rows: list[dict[str, str]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["name", "sequence"])
        writer.writeheader()
        writer.writerows(rows)



AUTO_QC_VALUE_FIELDS = [
    "p1_hairpin_found",
    "p1_hairpin_tm",
    "p1_hairpin_dg_kcal",
    "p1_hairpin_reason",
    "p1_hairpin_practical_reject",
    "p2_hairpin_found",
    "p2_hairpin_tm",
    "p2_hairpin_dg_kcal",
    "p2_hairpin_reason",
    "p2_hairpin_practical_reject",
    "p1_self_dimer_found",
    "p1_self_dimer_tm",
    "p1_self_dimer_dg_kcal",
    "p1_self_dimer_reason",
    "p1_self_dimer_practical_reject",
    "p2_self_dimer_found",
    "p2_self_dimer_tm",
    "p2_self_dimer_dg_kcal",
    "p2_self_dimer_reason",
    "p2_self_dimer_practical_reject",
    "same_pair_heterodimer_found",
    "same_pair_heterodimer_tm",
    "same_pair_heterodimer_dg_kcal",
    "same_pair_heterodimer_reason",
    "same_pair_heterodimer_practical_reject",
]


def tsv_cell_value(value: Any) -> Any:
    """Return a review-friendly TSV cell value, with blank cells for missing QC values."""
    if value is None:
        return ""
    if isinstance(value, float):
        if math.isnan(value) or math.isinf(value):
            return ""
        return round(value, 4)
    # Handle numpy scalar floats without importing numpy as a dependency.
    if hasattr(value, "item"):
        try:
            return tsv_cell_value(value.item())
        except Exception:
            pass
    try:
        if value != value:  # NaN-like objects
            return ""
    except Exception:
        pass
    return value


def write_dict_rows_tsv(rows: list[dict[str, Any]], path: Path, fields: list[str], *, sort_key=None) -> None:
    """Write TSV rows with consistent QC-value formatting."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="	", extrasaction="ignore")
        writer.writeheader()
        ordered = sorted(rows, key=sort_key) if sort_key else rows
        for row in ordered:
            writer.writerow({field: tsv_cell_value(row.get(field, "")) for field in fields})

def write_selected_pairs_tsv(selected: list[dict[str, Any]], path: Path) -> None:
    base_fields = ["name", "probe", "start", "length", "P1", "P2", "channel", "GC", "Tm", "dTm", "GibbsFE"]
    extra_fields = [
        "end",
        "source_run",
        "source_tier",
        "source_tile_size",
        "source_min_gc",
        "source_max_gc",
        "source_target_gibbs",
        "auto_score",
        "auto_qc_status",
        "auto_qc_reject_reason",
        "p1_hairpin_flag",
        "p2_hairpin_flag",
        "p1_self_dimer_flag",
        "p2_self_dimer_flag",
        "same_pair_heterodimer_flag",
        "transcript_region",
        "five_prime_utr_overlap_nt",
        "cds_overlap_nt",
        "three_prime_utr_overlap_nt",
        "noncoding_overlap_nt",
        "unannotated_overlap_nt",
        "premrna_intron_number",
        "premrna_genomic_start",
        "premrna_genomic_end",
    ]
    fields = base_fields + extra_fields + AUTO_QC_VALUE_FIELDS
    write_dict_rows_tsv(selected, path, fields, sort_key=lambda c: int(c["start"]))


def write_candidates_tsv(candidates: list[dict[str, Any]], path: Path) -> None:
    fields = [
        "candidate_id", "name", "start", "end", "length", "channel", "GC", "Tm", "dTm", "GibbsFE",
        "source_run", "source_tier", "source_tile_size", "source_min_gc", "source_max_gc", "source_target_gibbs",
        "auto_score", "auto_qc_status", "auto_qc_reject_reason",
        "p1_hairpin_flag", "p2_hairpin_flag", "p1_self_dimer_flag", "p2_self_dimer_flag", "same_pair_heterodimer_flag",
        "transcript_region", "five_prime_utr_overlap_nt", "cds_overlap_nt", "three_prime_utr_overlap_nt",
        "noncoding_overlap_nt", "unannotated_overlap_nt", "premrna_intron_number", "premrna_genomic_start",
        "premrna_genomic_end",
        *AUTO_QC_VALUE_FIELDS,
        "P1", "P2", "probe",
    ]
    write_dict_rows_tsv(
        candidates,
        path,
        fields,
        sort_key=lambda c: (int(c.get("source_tier", 0)), int(c.get("start", 0))),
    )



def analyze_final_cross_dimers(idt_rows: list[dict[str, str]]) -> tuple[Any, list[str], list[dict[str, Any]]]:
    """Return final cross-pool heterodimer flags and candidate keys to consider blocking.

    Candidate-level QC already handles self-dimers and same-pair P1/P2 dimers.
    This final-pool pass looks only at cross-pair interactions among selected
    oligos. REVIEW cross-pool dimers are not necessarily failures, but if a
    same-count replacement exists the repair loop will prefer a cleaner pool.
    """
    pd, _ = require_qc_dependencies()
    if not idt_rows:
        return pd.DataFrame(), [], []
    df = pd.DataFrame(idt_rows)
    pair_qc = make_pair_qc(df)
    bad = pair_qc[
        (pair_qc["dimer_class"] == "cross_pair")
        & ((pair_qc["heterodimer_flag"] != "PASS") | (pair_qc["heterodimer_practical_reject"] == True))
    ].copy()
    blocked_pair_keys: list[str] = []
    details: list[dict[str, Any]] = []
    for _, row in bad.iterrows():
        key1 = str(row["pair_key_1"])
        key2 = str(row["pair_key_2"])
        details.append({key: json_ready(row[key]) for key in row.index})
        blocked_pair_keys.extend([key1, key2])
    return pair_qc, blocked_pair_keys, details


def final_cross_dimer_objective(pair_qc: Any) -> tuple[float, float, float, float]:
    """Return cross-pool dimer burden for a final selected set.

    Lower is better. Strong/practical-reject cross-pool flags dominate REVIEW
    flags, and actual deltaG severity is used as a tie-breaker.
    """
    if pair_qc is None or len(pair_qc) == 0:
        return (0.0, 0.0, 0.0, 0.0)
    cross = pair_qc[pair_qc["dimer_class"] == "cross_pair"].copy()
    if cross.empty:
        return (0.0, 0.0, 0.0, 0.0)
    practical = float(cross["heterodimer_practical_reject"].fillna(False).astype(bool).sum())
    strong = float((cross["heterodimer_flag"] == "STRONG_REVIEW").sum())
    review = float((cross["heterodimer_flag"] == "REVIEW").sum())
    severity = 0.0
    for _, row in cross[cross["heterodimer_flag"] != "PASS"].iterrows():
        dg = optional_float(row.get("heterodimer_dg_kcal"))
        threshold = CROSS_DIMER_STRONG_DG if row.get("heterodimer_flag") == "STRONG_REVIEW" else CROSS_DIMER_REVIEW_DG
        severity += dimer_review_severity(dg, threshold)
    return (practical, strong, review, round(severity, 6))


def selected_pair_key_to_candidate_id(selected: list[dict[str, Any]], target_label: str, channel: str) -> dict[str, str]:
    mapping = {}
    for number, cand in enumerate(sorted(selected, key=lambda c: int(c["start"])), start=1):
        start = int(cand["start"])
        end = int(cand["end"])
        base = safe_name(f"{target_label}_{channel}_{number:02d}_{start}-{end}")
        mapping[base] = str(cand["candidate_id"])
    return mapping


def write_auto_curate_report_md(path: Path, report: dict[str, Any]) -> None:
    capacity = report.get("capacity", {})
    adaptive = report.get("adaptive_curation", {})
    coverage = report.get("coverage", {})
    final_coverage = coverage.get("final_metrics", {})
    candidate_requests = report.get("candidate_requests", {})
    resolution = report.get("transcript_resolution") or {}
    select_lookup = resolution.get("refseq_select_lookup") or {}
    annotation = report.get("transcript_annotation") or resolution.get("transcript_annotation") or {}
    region_counts = report.get("selected_probe_region_counts") or {}
    cds = annotation.get("cds") or {}
    five_utr = annotation.get("five_prime_utr") or {}
    three_utr = annotation.get("three_prime_utr") or {}
    five_utr_text = f"{five_utr.get('start')}-{five_utr.get('end')}" if five_utr else "not annotated"
    cds_text = f"{cds.get('start')}-{cds.get('end')}" if cds else "not annotated"
    three_utr_text = f"{three_utr.get('start')}-{three_utr.get('end')}" if three_utr else "not annotated"
    lines = [
        f"# HCR auto-curation report: {report.get('display_target_name') or report.get('target_name', '')}",
        "",
        f"Status: **{report.get('status', '')}** ({report.get('status_basis', '')})",
        f"Requested target: **{report.get('target_probes', 0)}** probe pairs",
        f"Capacity-aware deep-plan goal: **{capacity.get('adaptive_goal', 0)}** probe pairs",
        f"Initial QC-passing candidates: **{report.get('initial_good_probe_pairs', 0)}**",
        f"Normal-phase selected probes: **{adaptive.get('normal_phase_probe_pairs', 0)}**",
        f"Final selected probes: **{report.get('final_probe_pairs', 0)}**",
        f"Final selected dimer review flags: **{report.get('final_dimer_review_flags', 0)}**",
        f"Selected tile sizes: **{report.get('selected_tile_size_counts', {})}**",
        "",
        "## Transcript selection",
        "",
        f"- Selected accession: `{report.get('accession', '')}`",
        f"- Policy: `{resolution.get('transcript_selection_policy', '')}`",
        f"- Reason: {resolution.get('selection_reason', '')}",
        f"- RefSeq Select lookup: `{select_lookup.get('status', 'not applicable')}`",
        f"- Transcript ambiguity note: {resolution.get('transcript_ambiguity_note') or 'none'}",
        "",
        "## Transcript annotation",
        "",
        f"- Annotation status: `{annotation.get('status', 'unavailable')}`",
        f"- Source: {annotation.get('source', 'unavailable')}",
        f"- Coordinate system: {annotation.get('coordinate_system', '1-based inclusive transcript coordinates')}",
        f"- 5′ UTR: {five_utr_text}",
        f"- CDS: {cds_text}",
        f"- 3′ UTR: {three_utr_text}",
        f"- Selected probe regions: **{region_counts}**",
        f"- Annotation note: {annotation.get('note') or annotation.get('reason') or 'none'}",
        "",
        "## Capacity and adaptive fallback",
        "",
        f"- Canonical 52-nt capacity: **{capacity.get('canonical_52nt_capacity', 0)}**",
        f"- Normal-plan minimum tile: **{capacity.get('normal_plan_min_tile_size', 0)} nt**",
        f"- Normal-plan geometric capacity: **{capacity.get('normal_plan_capacity', 0)}**",
        f"- Deep-plan minimum tile: **{capacity.get('deep_plan_min_tile_size', 0)} nt**",
        f"- Deep-plan geometric capacity: **{capacity.get('deep_plan_capacity', 0)}**",
        f"- Adaptive fallback triggered: **{adaptive.get('triggered', False)}**",
        f"- Trigger reason: {adaptive.get('trigger_reason') or 'not triggered'}",
        f"- Fallback stringency: `{adaptive.get('fallback_stringency') or 'not used'}`",
        f"- Permissive-only selected candidates: **{adaptive.get('permissive_only_selected', 0)}**",
        "",
        "## Candidate reservoir",
        "",
        f"- Initial design request: **{candidate_requests.get('initial_max_probes', 0)}** candidates",
        f"- Per-tier request: **{candidate_requests.get('per_tier_max_probes', 0)}** candidates",
        f"- Initial request rule: {candidate_requests.get('initial_request_rule', '')}",
        f"- Per-tier request rule: {candidate_requests.get('per_tier_request_rule', '')}",
        f"- User overrode initial request: **{candidate_requests.get('max_probes_user_override', False)}**",
        f"- User overrode per-tier request: **{candidate_requests.get('auto_curate_max_probes_user_override', False)}**",
        "",
        "## Transcript coverage",
        "",
        f"- Policy: `{coverage.get('policy', '')}`",
        f"- Enabled for selection: **{coverage.get('enabled_for_selection', False)}**",
        f"- Probe span fraction: **{final_coverage.get('probe_span_fraction', 0):.3f}**",
        f"- Maximum untargeted gap: **{final_coverage.get('maximum_untargeted_gap_nt', 0)} nt**",
        f"- Ideal spacing: **{final_coverage.get('ideal_spacing_nt')} nt**",
        f"- Maximum gap / ideal spacing: **{final_coverage.get('maximum_gap_to_ideal_spacing_ratio')}**",
        f"- Occupied transcript bins: **{final_coverage.get('occupied_bins', 0)}/{final_coverage.get('total_bins', 10)}**",
        f"- Remaining coverage review conditions: {coverage.get('final_review_reasons') or 'none'}",
        "",
        "## Runs",
    ]
    for item in report.get("runs", []):
        lines.append(
            f"- {item.get('run_name')}: phase={item.get('phase')}; plan={item.get('plan_name')}; "
            f"plan tier={item.get('plan_tier_index')}; status={item.get('status', 'completed')}; "
            f"stringency={item.get('qc_stringency')}; requested max={item.get('requested_max_probes', 0)}; "
            f"candidates={item.get('probe_pairs', 0)}; accepted={item.get('good_probe_pairs', 0)}; "
            f"oligo QC={item.get('oligo_qc_status', 'completed')}; dir `{item.get('run_dir')}`"
        )
    lines += ["", "## Final intervals", ""]
    for item in report.get("final_selected", []):
        lines.append(
            f"- {item.get('start')}-{item.get('end')} ({item.get('length')} nt) from "
            f"{item.get('source_run')} region={item.get('transcript_region', 'unannotated')} "
            f"score={item.get('auto_score')}"
        )
    if report.get("notes"):
        lines += ["", "## Notes", ""]
        for note in report["notes"]:
            lines.append(f"- {note}")
    path.write_text("\n".join(lines) + "\n")


def run_auto_curate_if_needed(
    *,
    fasta: Path,
    output_dir: Path,
    target: str,
    args: argparse.Namespace,
    record: dict[str, Any],
    sequence_length: int,
) -> None:
    """Run normal curation first, then automatically rescue low-yield targets.

    The normal phase uses the requested plan and QC stringency. If it does not
    reach the capacity-aware goal (or retains selected dimer reviews), existing
    candidates are reclassified under permissive QC and all remaining distinct
    tiers in the deep plan are considered, subject to --max-auto-runs.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    details_dir = output_dir / "details"
    details_dir.mkdir(parents=True, exist_ok=True)
    design_root = details_dir / "design_runs"
    design_root.mkdir(parents=True, exist_ok=True)
    target_count = int(args.target_probes)
    blocked_ids: set[str] = set()
    all_candidates: list[dict[str, Any]] = []
    run_records: list[dict[str, Any]] = []
    notes: list[str] = []
    completed_signatures: set[tuple[Any, ...]] = set()
    execution_index = 0
    normal_extra_runs = 0
    adaptive_extra_runs = 0

    transcript_sequence = designable_sequence_from_fasta(fasta)
    premrna_mode = record.get("source") == "ncbi_premrna"
    premrna_intronic = premrna_mode and not _premrna_is_whole(record)

    canonical_capacity = (
        premrna_intron_capacity(record, transcript_sequence, 52)
        if premrna_intronic
        else tile_capacity(transcript_sequence, 52)
    )
    normal_min_tile = plan_min_tile_size(args.auto_curate_plan)
    normal_capacity = (
        premrna_intron_capacity(record, transcript_sequence, normal_min_tile)
        if premrna_intronic
        else tile_capacity(transcript_sequence, normal_min_tile)
    )
    deep_min_tile = plan_min_tile_size("deep")
    deep_capacity = (
        premrna_intron_capacity(record, transcript_sequence, deep_min_tile)
        if premrna_intronic
        else tile_capacity(transcript_sequence, deep_min_tile)
    )
    normal_goal = min(target_count, normal_capacity)
    adaptive_goal = min(target_count, deep_capacity)
    canonical_goal = min(target_count, canonical_capacity)
    balance_coverage = coverage_is_enabled(args.coverage_policy, canonical_capacity)
    if premrna_intronic:
        # Pre-mRNA targets are intentionally balanced by intron quota. The
        # ordinary transcript-wide coverage heuristic would count exons even
        # though candidate selection is restricted to introns.
        balance_coverage = False
        notes.append("Pre-mRNA selection used equal, capacity-aware quotas across the selected introns; transcript-wide coverage balancing was not applied.")
    premrna_intron_count = sum(
        1 for region in (record.get("premrna_target") or {}).get("regions", [])
        if str(region.get("kind")) == "intron"
    ) if premrna_intronic else 0
    premrna_reservoir = (
        max(int(args.max_probes), target_count * max(2, min(4, premrna_intron_count)))
        if premrna_intronic else int(args.max_probes)
    )

    def select_candidates(pool: list[dict[str, Any]], count: int, blocked: set[str]) -> list[dict[str, Any]]:
        if premrna_intronic:
            return select_premrna_candidates(
                pool,
                count,
                blocked,
                sequence_length=sequence_length,
                record=record,
            )
        return select_best_nonoverlapping(
            pool,
            count,
            blocked,
            sequence_length=sequence_length,
            coverage_enabled=balance_coverage,
        )
    if args.coverage_policy == "balanced" and not balance_coverage and not premrna_intronic:
        notes.append(
            f"Coverage balancing was disabled because canonical 52-nt capacity is {canonical_capacity}, "
            f"below the short-transcript cutoff of {COVERAGE_MIN_CANONICAL_CAPACITY}. Count and QC remain primary."
        )

    initial_params = {
        "tile_size": args.tile_size,
        "min_gc": args.min_gc,
        "max_gc": args.max_gc,
        "max_run_mismatches": args.max_run_mismatches,
        "min_gibbs": args.min_gibbs,
        "max_gibbs": args.max_gibbs,
        "target_gibbs": args.target_gibbs,
    }
    completed_signatures.add(parameter_signature(initial_params))

    initial_args = clone_args_with(
        args,
        max_probes=premrna_reservoir,
        _curation_phase="initial",
        _curation_plan="initial",
    )
    initial_name = tier_label(0, initial_params, initial=True)
    notify_progress("curation", f"Starting auto-curation pass: {initial_name}…")
    print(f"Auto-curation: running initial design {initial_name}")
    initial_result = run_design_tier(
        fasta=fasta,
        run_dir=design_root / initial_name,
        target=target,
        args=initial_args,
        record=record,
        sequence_length=sequence_length,
    )
    initial_candidates = [
        make_candidate_from_row(
            row,
            source_run=initial_name,
            source_tier=0,
            source_params=initial_params,
            target_label=target,
            stringency=args.qc_stringency,
        )
        for row in initial_result["rows"]
    ]
    all_candidates.extend(initial_candidates)
    initial_good = [candidate for candidate in initial_candidates if not candidate.get("auto_qc_reject_reason")]
    if not initial_candidates:
        notes.append(
            f"Initial design tier {initial_name} completed with zero candidates; oligo QC was skipped and relaxed tiers remained eligible."
        )
    run_records.append(
        {
            "run_name": initial_name,
            "execution_index": 0,
            "phase": "initial",
            "plan_name": "initial",
            "plan_tier_index": 0,
            "run_dir": str(initial_result["run_dir"]),
            "status": initial_result.get("status", "completed"),
            "oligo_qc_status": initial_result.get("qc_report", {}).get("status", "completed"),
            "probe_pairs": len(initial_candidates),
            "good_probe_pairs": len(initial_good),
            "requested_max_probes": int(initial_args.max_probes),
            "qc_stringency": args.qc_stringency,
            "parameters": initial_params,
        }
    )

    selected = select_candidates(deduplicate_candidates(all_candidates), target_count, blocked_ids)
    initial_selected_snapshot = list(selected)
    initial_dimer_reviews = selected_dimer_review_count(selected)
    normal_needed = not selected_set_reaches_clean_goal(
        selected,
        normal_goal,
        sequence_length=sequence_length,
        coverage_enabled=balance_coverage,
    )

    if normal_goal == 0:
        notes.append(
            f"The transcript has no complete {normal_min_tile}-nt unambiguous segment; "
            "normal-plan capacity is zero. Deep fallback may still help if a shorter tier fits."
        )
    elif not normal_needed:
        print(
            f"Initial design reached the normal capacity-aware goal ({normal_goal}) with clean dimer QC; "
            "normal relaxed tiers were not needed."
        )
    else:
        print(
            f"Initial design selected {len(selected)} probe pair(s); normal capacity-aware goal is {normal_goal}. "
            f"Running distinct tiers from the {args.auto_curate_plan} plan."
        )
        for plan_index, params in remaining_plan_tiers(args.auto_curate_plan, completed_signatures):
            if normal_extra_runs >= int(args.max_auto_runs):
                notes.append(
                    f"Stopped normal auto-curation after {normal_extra_runs} distinct additional run(s) "
                    "because the per-phase --max-auto-runs limit was reached."
                )
                break
            if selected_set_reaches_clean_goal(
                selected,
                normal_goal,
                sequence_length=sequence_length,
                coverage_enabled=balance_coverage,
            ):
                break
            normal_extra_runs += 1
            execution_index += 1
            completed_signatures.add(parameter_signature(params))
            run_name = tier_label(execution_index, params)
            tier_args = clone_args_with(
                args,
                tile_size=int(params["tile_size"]),
                min_gc=float(params["min_gc"]),
                max_gc=float(params["max_gc"]),
                max_run_mismatches=int(params["max_run_mismatches"]),
                min_gibbs=float(params["min_gibbs"]),
                max_gibbs=float(params["max_gibbs"]),
                target_gibbs=float(params["target_gibbs"]),
                max_probes=max(premrna_reservoir, args.auto_curate_max_probes),
                _curation_phase="normal",
                _curation_plan=args.auto_curate_plan,
            )
            print(f"  Running normal tier {plan_index} as {run_name}")
            notify_progress("curation", f"Running auto-curation tier {run_name}…")
            try:
                result = run_design_tier(
                    fasta=fasta,
                    run_dir=design_root / run_name,
                    target=target,
                    args=tier_args,
                    record=record,
                    sequence_length=sequence_length,
                )
            except subprocess.CalledProcessError:
                notes.append(f"Normal tier {run_name} failed; see its log file.")
                continue
            tier_candidates = [
                make_candidate_from_row(
                    row,
                    source_run=run_name,
                    source_tier=plan_index,
                    source_params=params,
                    target_label=target,
                    stringency=args.qc_stringency,
                )
                for row in result["rows"]
            ]
            all_candidates.extend(tier_candidates)
            good = [candidate for candidate in tier_candidates if not candidate.get("auto_qc_reject_reason")]
            if not tier_candidates:
                notes.append(
                    f"Normal tier {run_name} completed with zero candidates; oligo QC was skipped and the next eligible tier was considered."
                )
            run_records.append(
                {
                    "run_name": run_name,
                    "execution_index": execution_index,
                    "phase_run_index": normal_extra_runs,
                    "phase": "normal",
                    "plan_name": args.auto_curate_plan,
                    "plan_tier_index": plan_index,
                    "run_dir": str(result["run_dir"]),
                    "status": result.get("status", "completed"),
                    "oligo_qc_status": result.get("qc_report", {}).get("status", "completed"),
                    "probe_pairs": len(tier_candidates),
                    "good_probe_pairs": len(good),
                    "requested_max_probes": int(tier_args.max_probes),
                    "qc_stringency": args.qc_stringency,
                    "parameters": params,
                }
            )
            selected = select_candidates(deduplicate_candidates(all_candidates), target_count, blocked_ids)
            print(
                f"    Pool supports {len(selected)} selected probe pair(s); "
                f"selected dimer reviews={selected_dimer_review_count(selected)}; "
                f"coverage reviews={len(coverage_review_reasons(selected, sequence_length=sequence_length, enabled=balance_coverage))}."
            )

    normal_candidates = deduplicate_candidates(all_candidates)
    normal_selected = select_candidates(normal_candidates, target_count, blocked_ids)
    normal_dimer_reviews = selected_dimer_review_count(normal_selected)
    normal_coverage_reasons = coverage_review_reasons(
        normal_selected,
        sequence_length=sequence_length,
        enabled=balance_coverage,
    )

    fallback_reasons: list[str] = []
    if len(normal_selected) < adaptive_goal:
        fallback_reasons.append(
            f"normal phase selected {len(normal_selected)}, below the deep-plan capacity-aware goal of {adaptive_goal}"
        )
    if normal_dimer_reviews > 0:
        fallback_reasons.append(f"normal selected set retained {normal_dimer_reviews} dimer review flag(s)")
    if normal_coverage_reasons:
        fallback_reasons.extend(f"coverage: {reason}" for reason in normal_coverage_reasons)
    adaptive_triggered = bool(fallback_reasons) and adaptive_goal > 0

    pre_adaptive_pairs_tsv = None
    pre_adaptive_idt_csv = None
    fallback_stringency = None
    if adaptive_triggered:
        trigger_reason = "; ".join(fallback_reasons)
        print(f"Adaptive deep fallback triggered: {trigger_reason}.")
        notes.append(f"Adaptive deep fallback triggered because {trigger_reason}.")
        pre_adaptive_pairs_tsv = details_dir / f"{target}_normal_phase_selected_pairs.tsv"
        pre_adaptive_idt_csv = details_dir / f"{target}_normal_phase_IDT_order.csv"
        annotate_probe_rows(normal_selected, record.get("transcript_annotation") or {})
        write_selected_pairs_tsv(normal_selected, pre_adaptive_pairs_tsv)
        write_idt_csv(final_idt_rows(normal_selected, target, args.channel), pre_adaptive_idt_csv)

        low_yield_or_dimer_rescue = len(normal_selected) < adaptive_goal or normal_dimer_reviews > 0
        fallback_stringency = "permissive" if low_yield_or_dimer_rescue else args.qc_stringency
        reclassify_candidates(all_candidates, fallback_stringency)
        selected = select_candidates(deduplicate_candidates(all_candidates), target_count, blocked_ids)
        if fallback_stringency == "permissive":
            notes.append(
                "Existing candidates were reclassified under permissive QC without recalculating or changing Primer3 metrics."
            )
        else:
            notes.append(
                "The adaptive fallback was triggered for coverage only, so the original QC stringency was retained to avoid introducing additional QC flags."
            )

        for plan_index, params in remaining_plan_tiers("deep", completed_signatures):
            if adaptive_extra_runs >= int(args.max_auto_runs):
                notes.append(
                    f"Adaptive deep fallback stopped after {adaptive_extra_runs} distinct fallback run(s) "
                    "because the per-phase --max-auto-runs limit was reached."
                )
                break
            if selected_set_reaches_clean_goal(
                selected,
                adaptive_goal,
                sequence_length=sequence_length,
                coverage_enabled=balance_coverage,
            ):
                break
            adaptive_extra_runs += 1
            execution_index += 1
            completed_signatures.add(parameter_signature(params))
            run_name = tier_label(execution_index, params)
            tier_args = clone_args_with(
                args,
                tile_size=int(params["tile_size"]),
                min_gc=float(params["min_gc"]),
                max_gc=float(params["max_gc"]),
                max_run_mismatches=int(params["max_run_mismatches"]),
                min_gibbs=float(params["min_gibbs"]),
                max_gibbs=float(params["max_gibbs"]),
                target_gibbs=float(params["target_gibbs"]),
                max_probes=max(premrna_reservoir, args.auto_curate_max_probes),
                qc_stringency=fallback_stringency,
                _curation_phase="adaptive_deep",
                _curation_plan="deep",
            )
            print(f"  Running adaptive deep tier {plan_index} as {run_name}")
            notify_progress("curation", f"Running adaptive auto-curation tier {run_name}…")
            try:
                result = run_design_tier(
                    fasta=fasta,
                    run_dir=design_root / run_name,
                    target=target,
                    args=tier_args,
                    record=record,
                    sequence_length=sequence_length,
                )
            except subprocess.CalledProcessError:
                notes.append(f"Adaptive deep tier {run_name} failed; see its log file.")
                continue
            tier_candidates = [
                make_candidate_from_row(
                    row,
                    source_run=run_name,
                    source_tier=plan_index,
                    source_params=params,
                    target_label=target,
                    stringency=fallback_stringency,
                )
                for row in result["rows"]
            ]
            all_candidates.extend(tier_candidates)
            good = [candidate for candidate in tier_candidates if not candidate.get("auto_qc_reject_reason")]
            if not tier_candidates:
                notes.append(
                    f"Adaptive deep tier {run_name} completed with zero candidates; oligo QC was skipped and the next eligible tier was considered."
                )
            run_records.append(
                {
                    "run_name": run_name,
                    "execution_index": execution_index,
                    "phase_run_index": adaptive_extra_runs,
                    "phase": "adaptive_deep",
                    "plan_name": "deep",
                    "plan_tier_index": plan_index,
                    "run_dir": str(result["run_dir"]),
                    "status": result.get("status", "completed"),
                    "oligo_qc_status": result.get("qc_report", {}).get("status", "completed"),
                    "probe_pairs": len(tier_candidates),
                    "good_probe_pairs": len(good),
                    "requested_max_probes": int(tier_args.max_probes),
                    "qc_stringency": fallback_stringency,
                    "parameters": params,
                }
            )
            selected = select_candidates(deduplicate_candidates(all_candidates), target_count, blocked_ids)
            print(
                f"    Adaptive pool supports {len(selected)} selected probe pair(s); "
                f"selected dimer reviews={selected_dimer_review_count(selected)}; "
                f"coverage reviews={len(coverage_review_reasons(selected, sequence_length=sequence_length, enabled=balance_coverage))}."
            )
    else:
        trigger_reason = None
        selected = normal_selected
        if adaptive_goal == 0:
            notes.append(
                f"Adaptive deep fallback was not run because no complete {deep_min_tile}-nt unambiguous segment exists."
            )
        else:
            notes.append(
                f"Adaptive deep fallback was not needed: normal phase reached {len(normal_selected)} selected probe pair(s) "
                f"against a capacity-aware goal of {adaptive_goal}, with {normal_dimer_reviews} dimer review flag(s)."
            )

    candidates = deduplicate_candidates(all_candidates)

    # Repair final set if cross-pool dimers appear, preserving probe count.
    for repair_round in range(1, 5):
        idt_rows = final_idt_rows(selected, target, args.channel)
        pair_qc, bad_keys, _ = analyze_final_cross_dimers(idt_rows)
        current_cross_objective = final_cross_dimer_objective(pair_qc)
        if not bad_keys:
            break
        key_to_id = selected_pair_key_to_candidate_id(selected, target, args.channel)
        candidate_by_id = {str(candidate["candidate_id"]): candidate for candidate in selected}
        mapped_ids = []
        for key in bad_keys:
            candidate_id = key_to_id.get(key)
            if candidate_id and candidate_id not in mapped_ids:
                mapped_ids.append(candidate_id)
        if not mapped_ids:
            notes.append("Final QC found cross-pool heterodimer flags but could not map them for automatic repair.")
            break

        worst_id = max(
            mapped_ids,
            key=lambda candidate_id: (
                candidate_selection_vector(candidate_by_id.get(candidate_id, {})),
                float(candidate_by_id.get(candidate_id, {}).get("auto_score", 0.0)),
            ),
        )
        trial_blocked = set(blocked_ids)
        trial_blocked.add(worst_id)
        new_selected = select_candidates(candidates, target_count, trial_blocked)
        if len(new_selected) < len(selected):
            notes.append("Stopping final-pool repair because replacement would reduce the selected probe count.")
            break

        new_pair_qc, _, _ = analyze_final_cross_dimers(final_idt_rows(new_selected, target, args.channel))
        new_cross_objective = final_cross_dimer_objective(new_pair_qc)
        if (
            new_cross_objective,
            selected_set_objective(
                new_selected,
                sequence_length=sequence_length,
                coverage_enabled=balance_coverage,
            ),
        ) < (
            current_cross_objective,
            selected_set_objective(
                selected,
                sequence_length=sequence_length,
                coverage_enabled=balance_coverage,
            ),
        ):
            blocked_ids = trial_blocked
            selected = new_selected
            notes.append(
                f"Repair round {repair_round}: blocked candidate {worst_id} and improved final cross-pool "
                f"heterodimer profile from {current_cross_objective} to {new_cross_objective}."
            )
        else:
            notes.append("Stopping final-pool repair because the attempted replacement did not improve the selected set.")
            break

    transcript_annotation = record.get("transcript_annotation") or {}
    annotate_probe_rows(candidates, transcript_annotation)
    annotate_probe_rows(selected, transcript_annotation)

    final_pairs_tsv = output_dir / f"{target}_final_selected_pairs.tsv"
    final_idt_csv = output_dir / f"{target}_final_IDT_order.csv"
    final_qc_xlsx = output_dir / f"{target}_final_oligo_structure_QC.xlsx"
    final_map_png = output_dir / f"{target}_final_probe_map.png"
    candidate_master_tsv = details_dir / f"{target}_candidate_master.tsv"
    backup_tsv = details_dir / f"{target}_backup_candidates.tsv"
    rejected_tsv = details_dir / f"{target}_rejected_candidates.tsv"

    write_selected_pairs_tsv(selected, final_pairs_tsv)
    idt_rows = final_idt_rows(selected, target, args.channel)
    write_idt_csv(idt_rows, final_idt_csv)

    if selected:
        final_qc_report = run_oligo_structure_qc(final_idt_csv, final_qc_xlsx)
        final_rows = read_probe_rows(final_pairs_tsv)
        plot_probe_map(
            final_rows,
            sequence_length,
            args.plot_title or f"{args.gene or record['accession']} final curated · {record['accession']}",
            final_map_png,
            theme=args.plot_theme,
            color_by=args.plot_color_by,
            dpi=args.plot_dpi,
            show_labels=not args.no_probe_labels,
            transcript_annotation=transcript_annotation,
        )
    else:
        final_qc_report = {
            "status": "skipped",
            "reason": "no QC-acceptable non-overlapping candidates were selected",
            "output_file": None,
        }
        notes.append("No QC-acceptable non-overlapping candidates were selected; final QC and map were skipped.")

    selected_ids = {str(candidate["candidate_id"]) for candidate in selected}
    candidate_ids_sorted = sorted(
        candidates,
        key=lambda candidate: (
            candidate.get("auto_qc_status") != "PASS",
            candidate_selection_vector(candidate),
            float(candidate.get("auto_score", 0.0)),
        ),
    )
    write_candidates_tsv(candidate_ids_sorted, candidate_master_tsv)
    write_candidates_tsv(
        [
            candidate
            for candidate in candidate_ids_sorted
            if candidate.get("auto_qc_status") == "PASS" and str(candidate.get("candidate_id")) not in selected_ids
        ],
        backup_tsv,
    )
    write_candidates_tsv(
        [candidate for candidate in candidate_ids_sorted if candidate.get("auto_qc_status") != "PASS"],
        rejected_tsv,
    )

    final_dimer_review_count = selected_dimer_review_count(selected)
    final_coverage_metrics = selected_coverage_metrics(selected, sequence_length)
    final_coverage_reasons = coverage_review_reasons(
        selected,
        sequence_length=sequence_length,
        enabled=balance_coverage,
    )
    if final_dimer_review_count:
        notes.append(
            f"Final selected set has {final_dimer_review_count} self-dimer/same-pair dimer review flag(s); "
            "no cleaner same-count set was found in the attempted tiers."
        )
    if final_coverage_reasons:
        notes.append(
            "Final selected set retains coverage review condition(s): "
            + "; ".join(final_coverage_reasons)
            + ". QC and probe count were not weakened to force coverage."
        )

    requested_target_reached = len(selected) >= target_count
    adaptive_goal_reached = adaptive_goal > 0 and len(selected) >= adaptive_goal
    effective_min_acceptable = min(
        int(args.min_acceptable_probes),
        max(1, canonical_goal if canonical_goal > 0 else adaptive_goal if adaptive_goal > 0 else 1),
    )
    if requested_target_reached:
        status = "success"
        status_basis = "requested target reached"
    elif adaptive_goal_reached:
        status = "success"
        status_basis = "capacity-aware deep-plan goal reached; requested target exceeds geometric capacity"
    elif not candidates:
        status = "no_designable_candidates"
        status_basis = "all attempted design tiers completed without returning candidate probe pairs"
    elif len(selected) >= effective_min_acceptable:
        status = "partial_success"
        status_basis = f"at least the effective minimum of {effective_min_acceptable} probe pairs was reached"
    else:
        status = "insufficient"
        status_basis = f"fewer than the effective minimum of {effective_min_acceptable} probe pairs were selected"

    if not requested_target_reached:
        notes.append(
            f"Final set contains {len(selected)} probe pairs; requested target was {target_count}; "
            f"capacity-aware deep-plan goal was {adaptive_goal}."
        )

    permissive_only_selected = 0
    if adaptive_triggered:
        permissive_only_selected = sum(
            1
            for candidate in selected
            if pair_reject_reason(candidate, args.qc_stringency)
            and not pair_reject_reason(candidate, "permissive")
        )

    report = {
        "status": status,
        "status_basis": status_basis,
        "target_name": target,
        "gene_symbol": record.get("gene_symbol"),
        "display_target_name": record.get("gene_symbol") or record.get("accession"),
        "species": args.species,
        "organism": args.organism,
        "reference": getattr(args, "reference", None),
        "target_type": record_target_type(record),
        "premrna_target": record.get("premrna_target"),
        "accession": record.get("accession"),
        "transcript_resolution": record.get("transcript_resolution"),
        "transcript_annotation": record.get("transcript_annotation"),
        "selected_probe_region_counts": selected_probe_region_counts(selected),
        "target_probes": target_count,
        "requested_target_reached": requested_target_reached,
        "min_acceptable_probes": args.min_acceptable_probes,
        "effective_min_acceptable_probes": effective_min_acceptable,
        "qc_stringency": args.qc_stringency,
        "auto_curate_plan": args.auto_curate_plan,
        "max_auto_runs": args.max_auto_runs,
        "max_auto_runs_semantics": "maximum distinct additional design runs per phase; the adaptive deep phase receives its own limit",
        "auto_curate_invoked": normal_needed or adaptive_triggered,
        "initial_good_probe_pairs": len(initial_good),
        "initial_selected_probe_pairs": len(initial_selected_snapshot),
        "initial_selected_dimer_review_flags": initial_dimer_reviews,
        "final_probe_pairs": len(selected),
        "final_dimer_review_flags": final_dimer_review_count,
        "selected_tile_size_counts": tile_size_counts(selected),
        "shortest_selected_tile_size": min((int(candidate["length"]) for candidate in selected), default=None),
        "candidate_pool_size": len(candidates),
        "candidate_requests": {
            "initial_request_rule": "4 x target-probes unless --max-probes is supplied",
            "per_tier_request_rule": "5 x target-probes unless --auto-curate-max-probes is supplied",
            "initial_max_probes": int(args.max_probes),
            "per_tier_max_probes": int(max(args.max_probes, args.auto_curate_max_probes)),
            "max_probes_user_override": bool(getattr(args, "_max_probes_user_override", False)),
            "auto_curate_max_probes_user_override": bool(
                getattr(args, "_auto_curate_max_probes_user_override", False)
            ),
            "note": "HCRProbeDesign may still export fewer candidates because it applies upstream filtering and non-overlap selection.",
        },
        "capacity": {
            "method": "sum floor(unambiguous ACGT segment length / tile size)",
            "canonical_52nt_capacity": canonical_capacity,
            "canonical_goal": canonical_goal,
            "normal_plan": args.auto_curate_plan,
            "normal_plan_min_tile_size": normal_min_tile,
            "normal_plan_capacity": normal_capacity,
            "normal_goal": normal_goal,
            "deep_plan_min_tile_size": deep_min_tile,
            "deep_plan_capacity": deep_capacity,
            "adaptive_goal": adaptive_goal,
            "note": "Geometric capacities are upper bounds before sequence, specificity, thermodynamic, and QC filters.",
        },
        "coverage": {
            "policy": args.coverage_policy,
            "enabled_for_selection": balance_coverage,
            "short_transcript_rule": (
                "not applicable to intronic Pre-mRNA; equal capacity-aware intron quotas are used"
                if premrna_intronic
                else f"disabled when canonical 52-nt capacity is below {COVERAGE_MIN_CANONICAL_CAPACITY}"
            ),
            "selection_priority": (
                "probe count first; all explicit dimer and hairpin QC components second; "
                "coverage third; conventional technical and tile-length tie-breakers last"
            ),
            "thresholds": {
                "minimum_probe_span_fraction": COVERAGE_MIN_SPAN_FRACTION,
                "maximum_gap_to_ideal_spacing_ratio": COVERAGE_MAX_GAP_MULTIPLIER,
                "minimum_occupied_bins": COVERAGE_MIN_OCCUPIED_BINS,
                "total_bins": COVERAGE_BIN_COUNT,
            },
            "initial_metrics": selected_coverage_metrics(initial_selected_snapshot, sequence_length),
            "normal_phase_metrics": selected_coverage_metrics(normal_selected, sequence_length),
            "normal_phase_review_reasons": normal_coverage_reasons,
            "final_metrics": final_coverage_metrics,
            "final_review_reasons": final_coverage_reasons,
        },
        "adaptive_curation": {
            "enabled_automatically_with_auto_curate": True,
            "triggered": adaptive_triggered,
            "trigger_reason": trigger_reason,
            "normal_phase_probe_pairs": len(normal_selected),
            "normal_phase_dimer_review_flags": normal_dimer_reviews,
            "fallback_plan": "deep" if adaptive_triggered else None,
            "fallback_stringency": fallback_stringency,
            "normal_distinct_runs_executed": normal_extra_runs,
            "adaptive_deep_distinct_runs_executed": adaptive_extra_runs,
            "total_distinct_additional_runs_executed": normal_extra_runs + adaptive_extra_runs,
            "permissive_only_selected": permissive_only_selected,
        },
        "runs": run_records,
        "final_oligo_qc": final_qc_report,
        "outputs": {
            "final_selected_pairs_tsv": str(final_pairs_tsv),
            "final_idt_order_csv": str(final_idt_csv),
            "final_oligo_qc_xlsx": str(final_qc_xlsx) if selected else None,
            "final_probe_map_png": str(final_map_png) if selected else None,
            "final_probe_map_svg": str(final_map_png.with_suffix(".svg")) if selected else None,
            "candidate_master_tsv": str(candidate_master_tsv),
            "backup_candidates_tsv": str(backup_tsv),
            "rejected_candidates_tsv": str(rejected_tsv),
            "normal_phase_selected_pairs_tsv": str(pre_adaptive_pairs_tsv) if pre_adaptive_pairs_tsv else None,
            "normal_phase_idt_order_csv": str(pre_adaptive_idt_csv) if pre_adaptive_idt_csv else None,
            "transcript_resolution_json": record.get("transcript_resolution_json"),
        },
        "final_selected": [
            {
                "start": int(candidate["start"]),
                "end": int(candidate["end"]),
                "length": int(candidate["length"]),
                "source_run": candidate.get("source_run"),
                "source_tier": candidate.get("source_tier"),
                "candidate_id": candidate.get("candidate_id"),
                "auto_qc_stringency": candidate.get("auto_qc_stringency"),
                "auto_score": candidate.get("auto_score"),
                "transcript_region": candidate.get("transcript_region"),
                "five_prime_utr_overlap_nt": candidate.get("five_prime_utr_overlap_nt", 0),
                "cds_overlap_nt": candidate.get("cds_overlap_nt", 0),
                "three_prime_utr_overlap_nt": candidate.get("three_prime_utr_overlap_nt", 0),
            }
            for candidate in sorted(selected, key=lambda candidate: int(candidate["start"]))
        ],
        "notes": notes,
    }
    report_json = details_dir / f"{target}_auto_curate_report.json"
    report_md = details_dir / f"{target}_auto_curate_report.md"
    write_json(report_json, report)
    write_auto_curate_report_md(report_md, report)

    print(f"  final probes: {final_pairs_tsv}")
    print(f"  final IDT:    {final_idt_csv}")
    if selected:
        print(f"  final QC:     {final_qc_xlsx}")
        print(f"  final map:    {final_map_png}")
    else:
        print("  final QC:     skipped; no selected probes")
        print("  final map:    skipped; no selected probes")
    print(f"  report:       {report_json}")
    notify_progress("output", "Writing final probe tables, maps, QC workbook, and metadata…")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=f"HCRProbeForge {__version__}: design, QC, auto-curate, plot, annotate, and export split-initiator HCR probe sets.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog=(
            "For exact NCBI RefSeq transcripts, HCRProbeForge caches the validated transcript FASTA/record "
            "and feature annotation by accession-version below the selected project cache/species/workflow "
            "directory. HCRPROBEFORGE_CACHE_DIR can explicitly override that location. "
            "Incomplete/transient NCBI JSON responses are schema-validated and retried. Mature-transcript "
            "annotation is reporting-only for mature targets; the opt-in generalized Pre-mRNA mode uses matching genomic "
            "annotation to assemble a full exon-plus-intron target before restricting selection to introns. Custom FASTA inputs remain unannotated."
        ),
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {SCRIPT_BUILD}",
        help="Show the installed HCRProbeForge version and exit",
    )
    parser.add_argument("gene", nargs="?", help="Gene symbol/display label, e.g. twist1; may be used alone or with --accession. Filesystem-unsafe characters are preserved in metadata/display text while output path components use the existing safe-name rule")

    retrieval = parser.add_argument_group("transcript retrieval")
    retrieval.add_argument(
        "--organism",
        default=None,
        help="NCBI organism name; automatically selected from --species for the built-in presets",
    )
    retrieval.add_argument("--accession", help="Use a specific RefSeq accession, alone or with a gene symbol")
    retrieval.add_argument("--fasta", type=Path, help="Use a user-supplied single-record FASTA instead of NCBI gene/accession lookup; target name is inferred from the FASTA filename")
    retrieval.add_argument(
        "--target-type",
        choices=["mature", "pre-mrna", "pre-mrna-whole"],
        default="mature",
        help="Target sequence type. Mature uses the existing transcript workflow. Pre-mRNA · intronic reconstructs the full genomic transcript and restricts probe selection to complete selected introns. Pre-mRNA · whole reconstructs the same full genomic transcript and permits probe selection across exons and introns.",
    )
    retrieval.add_argument(
        "--premrna-introns",
        help="Optional comma-separated 1-based intron numbers for --target-type pre-mrna. Omit this option to use every annotated intron with even capacity-aware distribution. This option does not apply to pre-mrna-whole.",
    )
    retrieval.add_argument("--all-transcripts", action="store_true", help="Design every linked RefSeq RNA transcript instead of choosing one representative transcript")
    retrieval.add_argument(
        "--transcript-policy",
        choices=["auto", "longest", "interactive", "require-accession"],
        default="auto",
        help=(
            "How gene-only input chooses among linked RefSeq transcripts: auto prefers RefSeq Select when "
            "available, then curated, non-partial, longer records; longest ignores accession class; interactive "
            "prompts in the same invocation; require-accession preserves strict manual selection"
        ),
    )
    retrieval.add_argument("--email", default=os.getenv("NCBI_EMAIL"), help="Email sent with NCBI requests")
    retrieval.add_argument("--api-key", default=os.getenv("NCBI_API_KEY"), help="Optional NCBI API key")

    design = parser.add_argument_group("probe design")
    design.add_argument(
        "--species",
        default="xtr",
        help="Built-in or user-created registered genome alias. The selected preset supplies the organism, assembly, and matching genomic index.",
    )
    design.add_argument("--index", help="Bowtie2 index prefix; overrides the registered species index")
    design.add_argument("--channel", default="B1", choices=["B1", "B2", "B3", "B4", "B5"], help="HCR initiator channel")
    design.add_argument("--min-gc", type=float, default=45.0, help="Minimum tile GC percentage")
    design.add_argument("--max-gc", type=float, default=55.0, help="Maximum tile GC percentage")
    design.add_argument("--min-gibbs", type=float, default=-70.0, help="Lower Gibbs free-energy bound")
    design.add_argument("--max-gibbs", type=float, default=-50.0, help="Upper Gibbs free-energy bound")
    design.add_argument("--target-gibbs", type=float, default=-60.0, help="Preferred Gibbs value for candidate ranking")
    design.add_argument("--max-run-mismatches", type=int, default=2, help="Tolerance used by the C/G-rich run filter")
    design.add_argument(
        "--max-probes",
        type=int,
        default=argparse.SUPPRESS,
        help=(
            "Maximum non-overlapping probe pairs requested from the initial designProbes run. "
            "Default: 4 x --target-probes in smart auto-curation mode; 30 in one-run mode"
        ),
    )
    design.add_argument("--tile-size", type=int, default=52, help="Target tile length in nucleotides")
    design.add_argument("--num-hits-allowed", type=int, default=1, help="Maximum genomic alignments per candidate")
    design.add_argument("--dtm-filter", action="store_true", help="Enable probe-arm dTm filtering")
    design.add_argument("--dtm-max", type=float, default=5.0, help="Maximum dTm when --dtm-filter is enabled")
    design.add_argument("--no-genomemask", action="store_true", help="Disable genomic specificity screening")

    plotting = parser.add_argument_group("plotting")
    plotting.add_argument("--plot-only", type=Path, help="Plot an existing HCRProbeDesign TSV and exit")
    plotting.add_argument("--transcript-length", type=int, help="Full mature-transcript length for --plot-only")
    plotting.add_argument("--plot-theme", choices=["pastel", "minimal"], default="pastel", help="Plot theme")
    plotting.add_argument("--plot-color-by", choices=["gc", "dtm", "order"], default="gc", help="Variable used to color probe blocks")
    plotting.add_argument("--plot-dpi", type=int, default=300, help="PNG resolution; SVG is always generated")
    plotting.add_argument("--no-probe-labels", action="store_true", help="Hide probe numbers")
    plotting.add_argument("--plot-title", help="Override the automatic plot title")

    qc = parser.add_argument_group("final oligo QC")
    qc.add_argument("--no-oligo-qc", action="store_true", help="Skip automatic primer3 final-oligo QC after probe design")
    qc.add_argument("--qc-only", type=Path, help="Run final-oligo QC on an existing IDT CSV/TSV/XLSX or selected-pairs TSV and exit")
    qc.add_argument("--qc-output", type=Path, help="Output .xlsx path for --qc-only; automatic runs use <target>_oligo_structure_QC.xlsx")

    auto = parser.add_argument_group("smart auto-curation")
    auto.add_argument(
        "--auto-curate-if-needed",
        action="store_true",
        help=(
            "Run capacity-aware smart curation: use the requested plan/stringency first, then automatically "
            "continue with remaining deep-plan tiers under permissive candidate acceptance when yield is below "
            "the transcript-specific geometric goal or selected dimer reviews remain. With the default balanced "
            "coverage policy, additional tiers also run when a sufficiently long transcript has poor probe span, "
            "large gaps, or too few occupied transcript bins"
        ),
    )
    auto.add_argument("--target-probes", type=int, default=20, help="Desired final number of probe pairs in smart auto-curation mode")
    auto.add_argument("--min-acceptable-probes", type=int, default=12, help="User minimum for partial success; automatically capped by the transcript's canonical capacity for short targets")
    auto.add_argument("--qc-stringency", choices=["strict", "balanced", "permissive"], default="balanced", help="QC stringency used to classify candidates during auto-curation")
    auto.add_argument("--auto-curate-plan", choices=sorted(AUTO_CURATE_PLANS), default="standard", help="Tiered design plan used after the initial run if more candidates are needed")
    auto.add_argument("--max-auto-runs", type=int, default=12, help="Maximum number of distinct additional design runs in each phase; the normal plan and an adaptive deep fallback each receive this limit, and duplicate initial-equivalent tiers do not count")
    auto.add_argument(
        "--auto-curate-max-probes",
        type=int,
        default=argparse.SUPPRESS,
        help=(
            "Maximum probes requested from designProbes for each additional tier. "
            "Default: 5 x --target-probes; effective tier request is never lower than the initial request"
        ),
    )
    auto.add_argument(
        "--coverage-policy",
        choices=["balanced", "qc-only"],
        default="balanced",
        help=(
            "Selection policy for same-count, same-QC candidate sets. balanced favors even transcript coverage "
            "and can trigger additional tiers for large gaps; qc-only preserves the QC/technical objective. "
            "Coverage is automatically disabled for transcripts with canonical 52-nt capacity below 10"
        ),
    )

    output = parser.add_argument_group("output")
    output.add_argument("--outdir", type=Path, default=Path("hcr_results"), help="Project output root; runs and cache are created below it")
    parser.add_argument("--_species-scoped-outdir", action="store_true", help=argparse.SUPPRESS)
    return parser


def resolve_dynamic_defaults(args: argparse.Namespace) -> None:
    """Resolve candidate-request defaults after target-probes is known."""
    args._max_probes_user_override = hasattr(args, "max_probes")
    args._auto_curate_max_probes_user_override = hasattr(args, "auto_curate_max_probes")
    if not hasattr(args, "max_probes"):
        args.max_probes = 4 * int(args.target_probes) if args.auto_curate_if_needed else 30
    if not hasattr(args, "auto_curate_max_probes"):
        args.auto_curate_max_probes = 5 * int(args.target_probes)


def validate_args(args: argparse.Namespace) -> None:
    if not 0 <= args.min_gc < args.max_gc <= 100:
        raise ValueError("Require 0 <= --min-gc < --max-gc <= 100.")
    if args.min_gibbs >= args.max_gibbs:
        raise ValueError("Require --min-gibbs < --max-gibbs.")
    if not args.min_gibbs <= args.target_gibbs <= args.max_gibbs:
        raise ValueError("--target-gibbs must lie between --min-gibbs and --max-gibbs.")
    for name in ("tile_size", "max_probes", "num_hits_allowed"):
        if getattr(args, name) < 1:
            raise ValueError(f"--{name.replace('_', '-')} must be at least 1.")
    if args.max_run_mismatches < 0:
        raise ValueError("--max-run-mismatches cannot be negative.")
    if args.dtm_max < 0:
        raise ValueError("--dtm-max cannot be negative.")
    if args.plot_dpi < 72:
        raise ValueError("--plot-dpi must be at least 72.")
    if args.transcript_length is not None and args.transcript_length < 1:
        raise ValueError("--transcript-length must be positive.")
    if args.plot_only and args.qc_only:
        raise ValueError("Use only one of --plot-only or --qc-only.")
    if args.qc_output is not None and args.qc_output.suffix.lower() not in {".xlsx", ".xlsm"}:
        raise ValueError("--qc-output must end in .xlsx or .xlsm.")
    if args.target_probes < 1:
        raise ValueError("--target-probes must be at least 1.")
    if args.min_acceptable_probes < 1:
        raise ValueError("--min-acceptable-probes must be at least 1.")
    if args.min_acceptable_probes > args.target_probes:
        raise ValueError("--min-acceptable-probes cannot be greater than --target-probes.")
    if args.max_auto_runs < 0:
        raise ValueError("--max-auto-runs cannot be negative.")
    if args.auto_curate_max_probes < 1:
        raise ValueError("--auto-curate-max-probes must be at least 1.")
    if args.auto_curate_if_needed and args.no_oligo_qc:
        raise ValueError("--auto-curate-if-needed requires oligo QC; remove --no-oligo-qc.")
    if not args.qc_only and not args.plot_only:
        if args.fasta and (args.gene or args.accession):
            raise ValueError("--fasta cannot be combined with a gene symbol or --accession.")
        if not (args.gene or args.accession or args.fasta):
            raise ValueError("Provide a gene symbol, --accession, or --fasta.")
    if args.fasta and args.all_transcripts:
        raise ValueError("--all-transcripts applies only to NCBI gene lookup, not --fasta.")
    if args.accession and args.all_transcripts:
        raise ValueError("--all-transcripts cannot be combined with an explicit --accession.")
    if args.fasta and args.transcript_policy != "auto":
        raise ValueError("--transcript-policy applies only to NCBI gene-symbol lookup, not --fasta.")
    if args.target_type in {"pre-mrna", "pre-mrna-whole"}:
        if args.fasta:
            raise ValueError("Pre-mRNA design requires an NCBI transcript and the matching registered genomic assembly; it cannot use --fasta.")
        if args.all_transcripts:
            raise ValueError("Pre-mRNA design currently requires one selected transcript; remove --all-transcripts or provide --accession.")
        if args.index:
            raise ValueError("Pre-mRNA design uses the registered species assembly; omit --index so annotation and specificity screening use the same reference.")
        if args.target_type == "pre-mrna":
            from .premrna import parse_intron_selection

            parse_intron_selection(args.premrna_introns)
        elif args.premrna_introns:
            raise ValueError("--premrna-introns applies only to --target-type pre-mrna, not pre-mrna-whole.")
    elif args.premrna_introns:
        raise ValueError("--premrna-introns requires --target-type pre-mrna.")


def _main_impl(argv: list[str] | None = None) -> int:
    """Run the scientific pipeline through the package or original CLI path."""
    raise_if_cancelled()
    notify_progress("starting", "Starting HCRProbeForge and resolving the requested target…")
    args = build_parser().parse_args(argv)
    args._pipeline_argv = sys.argv[:] if argv is None else ["hcrprobeforge", *argv]
    # Print the exact executable path and build so helper-script runs can be
    # diagnosed unambiguously when multiple copies exist on the same machine.
    print(f"HCRProbeForge build: {SCRIPT_BUILD}")
    print(f"HCRProbeForge script: {Path(__file__).resolve()}")
    resolve_species_arguments(args, list(argv) if argv is not None else sys.argv[1:])
    os.environ["HCRPROBEFORGE_SPECIES"] = str(args.species)
    workflow = "plot" if args.plot_only else "design"
    os.environ["HCRPROBEFORGE_WORKFLOW"] = workflow
    # Keep direct CLI transcript/annotation caches beside the selected project
    # rather than silently creating a second cache tree under the cwd.  An
    # explicit HCRPROBEFORGE_CACHE_DIR remains authoritative.
    if not os.getenv("HCRPROBEFORGE_CACHE_DIR"):
        from .references import species_cache_root

        os.environ["HCRPROBEFORGE_PROJECT_ROOT"] = str(
            species_cache_root(args.outdir, args.species, workflow).parents[2]
        )
    resolve_dynamic_defaults(args)
    validate_args(args)
    if not args.qc_only and not args.plot_only and not args.no_genomemask and args.index is None:
        from . import references

        if not references.registered_index_is_ready(args.species):
            raise RuntimeError(missing_genome_index_message(args.species))

    # QC-only mode should not create the default hcr_results/ directory.
    # If --qc-output is not provided, place the QC workbook beside the input file.
    if args.qc_only:
        qc_input = args.qc_only.expanduser().resolve()
        qc_output = (
            args.qc_output.expanduser().resolve()
            if args.qc_output
            else qc_input.with_name(f"{qc_input.stem}_oligo_structure_QC.xlsx")
        )
        qc_summary = run_oligo_structure_qc(qc_input, qc_output)
        print(qc_output)
        print(
            "QC summary: "
            f"hairpin strong={qc_summary.get('hairpin_STRONG_REVIEW', 0)}, "
            f"self-dimer strong={qc_summary.get('self_dimer_STRONG_REVIEW', 0)}, "
            f"heterodimer strong={qc_summary.get('heterodimer_STRONG_REVIEW', 0)}"
        )
        return 0

    args.outdir = args.outdir.expanduser().resolve()
    if not args._species_scoped_outdir:
        from .references import species_run_root

        workflow = "plot" if args.plot_only else "design"
        args.outdir = species_run_root(args.outdir, args.species, workflow)
    args.outdir.mkdir(parents=True, exist_ok=True)

    if args.plot_only:
        rows = read_probe_rows(args.plot_only)
        final_probe_end = max(row["end"] for row in rows)
        length = args.transcript_length or final_probe_end
        if length < final_probe_end:
            raise ValueError("--transcript-length is shorter than the final probe coordinate.")
        output = args.outdir / f"{args.plot_only.stem}_probe_map.png"
        plot_annotation = load_sibling_transcript_annotation(args.plot_only)
        if plot_annotation:
            annotate_probe_rows(rows, plot_annotation)
        plot_probe_map(
            rows,
            length,
            args.plot_title or args.plot_only.stem,
            output,
            theme=args.plot_theme,
            color_by=args.plot_color_by,
            dpi=args.plot_dpi,
            show_labels=not args.no_probe_labels,
            transcript_annotation=plot_annotation,
        )
        print(output)
        return 0

    if args.fasta:
        session = None
        records = resolve_user_fasta_record(args)
    else:
        if requests is None:
            raise RuntimeError(
                "NCBI-backed design requires requests. Install the package dependencies with "
                "python -m pip install hcrprobeforge."
            )
        session = requests.Session()
        session.headers.update({"User-Agent": "hcrprobeforge"})
        notify_progress("transcript", "Resolving the selected NCBI transcript…")
        records = resolve_transcripts(session, args)

    for record in records:
        if args.target_type in {"pre-mrna", "pre-mrna-whole"}:
            from . import premrna

            notify_progress(
                "premrna",
                "Preparing the full genomic Pre-mRNA target and intron map…",
            )
            premrna.prepare_target(session, record, args)
        # Preserve the biological/display label separately from filesystem-safe target IDs.
        record["gene_symbol"] = args.gene
        label = target_label_for_record(args, record)
        workdir = args.outdir / label
        workdir.mkdir(parents=True, exist_ok=True)
        details_dir = workdir / "details"
        details_dir.mkdir(parents=True, exist_ok=True)
        fasta = details_dir / f"{label}.fa"
        sequence_length = prepare_fasta_for_record(session, record, fasta, args)

        transcript_annotation = resolve_transcript_annotation(session, record, args, sequence_length)
        record["transcript_annotation"] = transcript_annotation

        if record.get("source") == "user_fasta" and not record.get("transcript_resolution"):
            record["transcript_resolution"] = {
                "mode": "user_fasta",
                "selected_accession": record.get("accession"),
                "transcript_selection_policy": "user_fasta",
                "selection_reason": "The user supplied a single-record FASTA; NCBI resolution was bypassed.",
                "input_fasta": record.get("input_fasta"),
            }
        if record.get("transcript_resolution") is not None:
            record["transcript_resolution"]["reference"] = args.reference
            record["transcript_resolution"]["species"] = args.species
            record["transcript_resolution"]["scientific_name"] = args.organism
            record["transcript_resolution"]["target_type"] = record_target_type(record)
            if record.get("source") == "ncbi_premrna":
                record["transcript_resolution"]["premrna_target"] = record.get("premrna_target")
            record["transcript_resolution"]["transcript_annotation"] = transcript_annotation
        resolution_json = details_dir / f"{label}_transcript_resolution.json"
        write_json(resolution_json, record.get("transcript_resolution") or {"transcript_annotation": transcript_annotation})
        record["transcript_resolution_json"] = str(resolution_json)

        if record.get("source") == "user_fasta":
            source_label = "user FASTA"
        elif record.get("transcript_cache_status") == "hit":
            source_label = "local NCBI transcript cache"
        else:
            source_label = "NCBI"
        print(f"Selected {record['accession']} ({sequence_length} nt, {source_label}): {record['title']}")
        if record.get("transcript_resolution"):
            print(f"  selection: {record['transcript_resolution'].get('selection_reason', '')}")
            print(f"  resolution metadata: {resolution_json}")
        annotation_status = transcript_annotation.get("status", "unavailable")
        if annotation_status == "resolved_premrna":
            intron_count = sum(
                1 for region in (transcript_annotation.get("regions") or [])
                if str(region.get("kind")) == "intron"
            )
            if str(transcript_annotation.get("selection_mode") or "") == "whole":
                print("  target annotation: whole Pre-mRNA; exons and introns are eligible")
            else:
                print(
                    "  target annotation: intronic Pre-mRNA; "
                    f"{intron_count} selected intron region(s) with even quota selection"
                )
        elif annotation_status in {"resolved", "resolved_partial"}:
            cds = transcript_annotation.get("cds") or {}
            five = transcript_annotation.get("five_prime_utr")
            three = transcript_annotation.get("three_prime_utr")
            five_text = f"{five.get('start')}-{five.get('end')}" if five else "not annotated"
            cds_text = f"{cds.get('start')}-{cds.get('end')}" if cds else "not annotated"
            three_text = f"{three.get('start')}-{three.get('end')}" if three else "not annotated"
            print(
                "  transcript annotation: "
                f"{annotation_status}; 5′ UTR={five_text}, CDS={cds_text}, 3′ UTR={three_text}"
            )
        else:
            print(f"  transcript annotation: {annotation_status}")
        if args.auto_curate_if_needed:
            notify_progress(
                "design",
                "Running designProbes, Bowtie2 specificity screening, and probe curation…",
            )
            run_auto_curate_if_needed(
                fasta=fasta,
                output_dir=workdir,
                target=label,
                args=args,
                record=record,
                sequence_length=sequence_length,
            )
            continue

        notify_progress("design", "Running designProbes and Bowtie2 specificity screening…")
        probes_tsv, idt_tsv, log_file, run_parameters_json, run_metadata = run_design_probes(fasta, details_dir, label, args, record, sequence_length)
        rows = read_probe_rows(probes_tsv)
        candidate_idt_tsv = idt_tsv
        eligible_idt_csv: Path | None = None
        if record.get("source") == "ncbi_premrna":
            rows = filter_premrna_candidates(rows, record)
            rows = enforce_premrna_nonoverlap(rows, sequence_length)
            # Keep the native upstream table for diagnostics, but never expose
            # it as the order table for an intronic target: it can contain
            # exonic and boundary-spanning candidates that are removed above.
            candidate_idt_tsv = details_dir / f"{label}_candidate_IDT.tsv"
            os.replace(idt_tsv, candidate_idt_tsv)
            if rows:
                eligible_idt_csv = workdir / f"{label}_eligible_IDT_order.csv"
                write_idt_csv(final_idt_rows(rows, label, args.channel), eligible_idt_csv)
            outputs = run_metadata.setdefault("outputs", {})
            outputs.pop("idt_tsv", None)
            outputs["candidate_idt_tsv"] = str(candidate_idt_tsv)
            outputs["eligible_idt_order_csv"] = str(eligible_idt_csv) if eligible_idt_csv else None
        annotate_probe_rows(rows, transcript_annotation)
        map_png = workdir / f"{label}_probe_map.png"
        if rows:
            plot_probe_map(
                rows,
                sequence_length,
                args.plot_title
                or f"{args.gene or record['accession']}  ·  {record['accession']}",
                map_png,
                theme=args.plot_theme,
                color_by=args.plot_color_by,
                dpi=args.plot_dpi,
                show_labels=not args.no_probe_labels,
                transcript_annotation=transcript_annotation,
            )
        else:
            map_png = None
            print("  No eligible probe pairs remained after Pre-mRNA intron filtering; probe map skipped.")
        qc_xlsx = workdir / f"{label}_oligo_structure_QC.xlsx"
        if rows and not args.no_oligo_qc:
            qc_input = eligible_idt_csv or candidate_idt_tsv
            qc_report = run_oligo_structure_qc(qc_input, qc_xlsx)
            run_metadata.setdefault("outputs", {})["oligo_qc_xlsx"] = str(qc_xlsx)
            run_metadata.setdefault("outputs", {})["oligo_qc_input"] = str(qc_input)
            run_metadata["final_oligo_qc"] = qc_report
            write_json(run_parameters_json, run_metadata)
        elif not rows:
            qc_report = {"status": "skipped", "reason": "no eligible probe pairs after target filtering"}
            run_metadata["final_oligo_qc"] = qc_report
            write_json(run_parameters_json, run_metadata)
        else:
            qc_report = None
            run_metadata["final_oligo_qc"] = {"status": "skipped", "reason": "--no-oligo-qc"}
            write_json(run_parameters_json, run_metadata)

        summary = details_dir / f"{label}_summary.json"
        write_summary(record, rows, sequence_length, summary, run_metadata)

        print(f"  FASTA:   {fasta}")
        print(f"  probes:  {probes_tsv}")
        if record.get("source") == "ncbi_premrna":
            print(f"  candidate IDT: {candidate_idt_tsv}")
            print(f"  eligible IDT:  {eligible_idt_csv}")
        else:
            print(f"  IDT:     {idt_tsv}")
        print(f"  log:     {log_file}")
        print(f"  map:     {map_png}")
        if qc_report is not None:
            print(f"  QC:      {qc_xlsx}")
        print(f"  summary: {summary}")
        print(f"  params:  {run_parameters_json}")

    return 0


def _snapshot_filesystem(root: Path) -> set[Path]:
    """Record a tree before a run so newly created failure artifacts can roll back."""
    root = Path(root).expanduser().resolve()
    if not root.is_dir():
        return set()
    return {root, *root.rglob("*")}


def _remove_new_filesystem_entries(root: Path, before: set[Path]) -> None:
    """Remove only entries created below *root* after a run started.

    The rollback is deliberately conservative: existing user results are not
    removed, while partial FASTA/JSON/plot/QC files and empty directories from
    a failed new run are cleaned up.  Failure logs and manifest status records
    are owned by the batch runner and remain available for diagnosis.
    """
    root = Path(root).expanduser().resolve()
    if not root.exists():
        return
    entries = sorted(
        (path for path in root.rglob("*") if path not in before),
        key=lambda path: len(path.parts),
        reverse=True,
    )
    for path in entries:
        try:
            if path.is_symlink() or path.is_file():
                path.unlink()
            elif path.is_dir() and not any(path.iterdir()):
                path.rmdir()
        except OSError:
            # Leave an entry alone if another process is using it or the user
            # changed it while the run was being rolled back.
            continue
    try:
        if root not in before and root.is_dir() and not any(root.iterdir()):
            root.rmdir()
    except OSError:
        pass


def _backup_existing_filesystem_contents(roots: list[Path]) -> tuple[Path, list[tuple[Path, Path]]]:
    """Back up regular files below output roots for failure restoration.

    A filesystem-entry snapshot removes new artifacts but cannot restore an
    existing result file that a rerun overwrote. Only the visible output root
    is backed up by :func:`main`; reference assets and transcript caches have
    their own transactional writers and are not copied here.
    """
    backup_root = Path(tempfile.mkdtemp(prefix="hcrprobeforge-output-rollback-"))
    backups: list[tuple[Path, Path]] = []
    try:
        for root_index, root in enumerate(roots):
            root = Path(root).expanduser().resolve()
            if not root.is_dir():
                continue
            for path in root.rglob("*"):
                if not path.is_file() or path.is_symlink():
                    continue
                relative = path.relative_to(root)
                backup = backup_root / str(root_index) / relative
                backup.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(path, backup)
                backups.append((path, backup))
    except Exception:
        shutil.rmtree(backup_root, ignore_errors=True)
        raise
    return backup_root, backups


def _restore_existing_filesystem_contents(
    backup_root: Path,
    backups: list[tuple[Path, Path]],
) -> None:
    """Restore files saved by :func:`_backup_existing_filesystem_contents`."""
    try:
        for original, backup in backups:
            original.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(backup, original)
    finally:
        shutil.rmtree(backup_root, ignore_errors=True)


def _discard_existing_filesystem_backup(backup_root: Path) -> None:
    shutil.rmtree(backup_root, ignore_errors=True)


def _cleanup_roots_for_argv(argv: list[str]) -> list[Path]:
    """Resolve output/cache trees without performing network or design work."""
    try:
        args = build_parser().parse_args(argv)
    except SystemExit:
        return []

    if args.qc_only:
        input_path = args.qc_only.expanduser().resolve()
        output_path = (
            args.qc_output.expanduser().resolve()
            if args.qc_output
            else input_path.with_name(f"{input_path.stem}_oligo_structure_QC.xlsx")
        )
        return [output_path.parent]

    try:
        resolve_species_arguments(args, argv)
    except (RuntimeError, OSError, ValueError):
        # Species/input validation failed before the pipeline can create an
        # output tree.  There is nothing safe to roll back in this case.
        return []
    workflow = "plot" if args.plot_only else "design"
    if args._species_scoped_outdir:
        output_root = args.outdir.expanduser().resolve()
    else:
        from .references import species_run_root

        output_root = species_run_root(args.outdir, args.species, workflow)
    roots = [output_root]
    cache_override = os.getenv("HCRPROBEFORGE_CACHE_DIR")
    if cache_override:
        roots.append(Path(cache_override).expanduser().resolve())
    else:
        from .references import species_cache_root

        roots.append(species_cache_root(args.outdir, args.species, workflow))
    return roots


def main(argv: list[str] | None = None) -> int:
    """Run one pipeline invocation with conservative failure rollback."""
    values = list(argv) if argv is not None else sys.argv[1:]
    roots = _cleanup_roots_for_argv(values)
    snapshots = [(root, _snapshot_filesystem(root)) for root in roots]
    output_backup_root, output_backups = _backup_existing_filesystem_contents(roots[:1])
    try:
        # Preserve the original argv metadata when the public API was called
        # without an explicit argument list.
        result = _main_impl(None if argv is None else values)
        _discard_existing_filesystem_backup(output_backup_root)
        return result
    except Exception:
        for root, before in snapshots:
            _remove_new_filesystem_entries(root, before)
        _restore_existing_filesystem_contents(output_backup_root, output_backups)
        raise


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except subprocess.CalledProcessError as exc:
        command = [str(item) for item in exc.cmd] if isinstance(exc.cmd, (list, tuple)) else [str(exc.cmd)]
        print(design_process_failure_message(int(exc.returncode or 1), command), file=sys.stderr)
        raise SystemExit(exc.returncode)
    except RetryableNCBIError as exc:
        print(f"ERROR: temporary NCBI failure after retries: {exc}", file=sys.stderr)
        raise SystemExit(75)
    except (RuntimeError, OSError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
