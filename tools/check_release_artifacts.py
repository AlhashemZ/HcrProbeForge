#!/usr/bin/env python3
"""Check wheel and sdist contents for accidental local/generated artifacts.

This checker is intentionally dependency-free so it can be run before the
development extras are installed.  It checks archive member paths and scans
text metadata for common absolute local paths.  It is a release hygiene check,
not a substitute for ``twine check`` or a clean-environment install test.
"""

from __future__ import annotations

import argparse
import re
import sys
import tarfile
import zipfile
from pathlib import PurePosixPath
from pathlib import Path


FORBIDDEN_PARTS = {
    ".DS_Store",
    "__MACOSX",
    "build",
    "dist",
    "release-dist",
    "hcr_results",
    "index_metadata",
}
FORBIDDEN_SUFFIXES = (
    ".sqlite",
    ".sqlite-shm",
    ".sqlite-wal",
    ".bt2",
    ".bt2l",
)
ABSOLUTE_PATH_PATTERNS = (
    re.compile(r"/(?:Users|home|workspace|private/var)/(?!<|\$|Name/)[^\s'\"<>]+"),
    re.compile(r"[A-Za-z]:[\\/](?:Users|home|workspace)[\\/](?!<|\$|Name[\\/])[^\s'\"<>]+"),
)
TEXT_SUFFIXES = {
    ".c",
    ".cfg",
    ".cff",
    ".csv",
    ".ini",
    ".json",
    ".md",
    ".py",
    ".rst",
    ".toml",
    ".tsv",
    ".txt",
    ".yaml",
    ".yml",
}


def path_issues(name: str) -> list[str]:
    """Return release-hygiene issues for one archive member path."""
    issues: list[str] = []
    normalized = name.replace("\\", "/")
    path = PurePosixPath(normalized)
    parts = path.parts
    if normalized.startswith("/") or (parts and parts[0].endswith(":")):
        issues.append("absolute member path")
    if ".." in parts:
        issues.append("parent traversal member path")
    if any(part == "._" or part.startswith("._") for part in parts):
        issues.append("macOS resource-fork metadata")
    if any(part in FORBIDDEN_PARTS for part in parts):
        issues.append("generated or machine-local directory/file")
    if path.name in {".DS_Store", "Thumbs.db"}:
        issues.append("machine metadata file")
    if path.suffix.lower() in FORBIDDEN_SUFFIXES:
        issues.append("reference/index/database artifact")
    return issues


def _text_issues(name: str, data: bytes) -> list[str]:
    if PurePosixPath(name).suffix.lower() not in TEXT_SUFFIXES:
        return []
    try:
        text = data.decode("utf-8")
    except UnicodeDecodeError:
        return []
    if any(pattern.search(text) for pattern in ABSOLUTE_PATH_PATTERNS):
        return ["absolute local path in text file"]
    return []


def check_zip(path: Path) -> list[str]:
    issues: list[str] = []
    with zipfile.ZipFile(path) as archive:
        for member in archive.infolist():
            name = member.filename
            issues.extend(f"{name}: {issue}" for issue in path_issues(name))
            if not member.is_dir():
                issues.extend(f"{name}: {issue}" for issue in _text_issues(name, archive.read(member)))
    return issues


def check_tar(path: Path) -> list[str]:
    issues: list[str] = []
    with tarfile.open(path, "r:*") as archive:
        for member in archive.getmembers():
            name = member.name
            issues.extend(f"{name}: {issue}" for issue in path_issues(name))
            if member.isfile():
                handle = archive.extractfile(member)
                if handle is not None:
                    issues.extend(f"{name}: {issue}" for issue in _text_issues(name, handle.read()))
    return issues


def check_artifact(path: Path) -> list[str]:
    if not path.is_file():
        return [f"{path}: not a regular file"]
    if path.suffix == ".whl" or path.name.endswith(".zip"):
        return check_zip(path)
    if path.name.endswith((".tar.gz", ".tar.bz2", ".tar.xz", ".tgz")):
        return check_tar(path)
    return [f"{path}: unsupported artifact type"]


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("artifacts", nargs="+", type=Path)
    args = parser.parse_args(argv)

    issues: list[str] = []
    for artifact in args.artifacts:
        artifact_issues = check_artifact(artifact)
        if artifact_issues:
            issues.extend(f"{artifact}: {issue}" for issue in artifact_issues)
        else:
            print(f"OK {artifact}")
    if issues:
        for issue in issues:
            print(f"ERROR {issue}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
