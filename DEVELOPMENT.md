# Maintainer and release guide

This document contains repository and release procedures that are intentionally
kept out of the end-user README.

## Versioning

Update the version consistently in:

- `pyproject.toml`;
- `VERSION`;
- `src/hcrprobeforge/core.py` (`__version__` and CLI build string);
- `CITATION.cff`;
- `CHANGELOG.md`; and
- `VERSION_HISTORY.md`; and
- `VALIDATION_REPORT.md`.

Use a new version whenever the wheel/sdist contents or public behavior change.
The package metadata and `hcrprobeforge --version` must report the same value.

## Validation

From a clean checkout, run:

```bash
python -m compileall -q src tests
PYTHONPATH=src python -m unittest discover -s tests -p 'test_*.py'
python -m ruff check src tests
```

On Windows PowerShell, set the checkout's source directory for the test
process before running the equivalent command:

```powershell
$env:PYTHONPATH = "src"
python -m unittest discover -s tests -p 'test_*.py'
```

Build into a fresh directory so stale artifacts cannot be mistaken for the
current release:

```bash
rm -rf release-dist
mkdir release-dist
python -m pip wheel --no-deps --no-build-isolation --wheel-dir release-dist .
python -m build --sdist --outdir release-dist .
python -m twine check release-dist/*
python tools/check_release_artifacts.py release-dist/*
```

On Windows PowerShell, use the equivalent clean-directory commands:

```powershell
Remove-Item -Recurse -Force release-dist -ErrorAction SilentlyContinue
New-Item -ItemType Directory release-dist | Out-Null
python -m pip wheel --no-deps --no-build-isolation --wheel-dir release-dist .
python -m build --sdist --outdir release-dist .
python -m twine check (Get-ChildItem release-dist)
python tools/check_release_artifacts.py (Get-ChildItem release-dist | ForEach-Object { $_.FullName })
```

If `python -m build` or `twine` is not installed, install the release tools in
the isolated development environment before publishing. Inspect the wheel and
source archive contents. They must not contain generated caches, downloaded
genomes, SQLite databases, local result files, `.DS_Store`, `__MACOSX`, or
absolute local paths. The artifact checker performs these portable structural
checks; review its output rather than treating a successful build alone as
proof that the release contents are clean.

Install the wheel into a clean environment and verify:

```bash
python -m pip install --force-reinstall release-dist/hcrprobeforge-*.whl
hcrprobeforge --version
hcrprobeforge --help
```

The HCRProbeDesign Python dependency and its command providers are installed by
the wheel's dependency set. Bowtie2 remains an external deployment prerequisite
and is not bundled; network-backed NCBI checks and reference indexes are also
not bundled.

## Release bundle

A GitHub release bundle should use a release-specific top-level directory. The
source tree is the GitHub repository content; the PyPI upload files are kept in
the separate `pypi/` directory so the wheel can be attached to the GitHub
release under the same release name:

```text
release-bundle-<release>/
├── github/
│   ├── hcrprobeforge-<release>/       # clean repository source tree
│   └── hcrprobeforge-<release>-source.zip
├── pypi/
│   ├── hcrprobeforge-<release>-py3-none-any.whl
│   └── hcrprobeforge-<release>.tar.gz
├── RELEASE_NOTES.md
└── SHA256SUMS
```

The README's release-installation instructions point users to the wheel in the
GitHub release assets or to `pypi/` in a downloaded bundle. Generate checksums
from the final files only. On Linux use:

```bash
sha256sum release-dist/* > SHA256SUMS
```

On macOS use:

```bash
shasum -a 256 release-dist/* > SHA256SUMS
```

On Windows PowerShell use:

```powershell
Get-FileHash release-dist/* -Algorithm SHA256
```

Do not include virtual environments, build directories, test output, user
reference data, or transient logs. Before publishing, compare the wheel and
sdist version, inspect `PKG-INFO`, and verify that `README.md` renders without
developer-only instructions.

## Publication checklist

- [ ] Tests and static checks pass.
- [ ] Wheel and sdist build from a clean tree.
- [ ] `twine check` passes.
- [ ] Wheel/sdist versions match all source metadata.
- [ ] README commands match the installed entry points.
- [ ] Pre-mRNA output and accession-version behavior are covered by tests.
- [ ] `CITATION.cff`, changelog, and version history are updated.
- [ ] `VALIDATION_REPORT.md` records the current release and checks.
- [ ] License and third-party notices are included in the sdist.
- [ ] `tools/check_release_artifacts.py` passes for the final wheel and sdist.
- [ ] Checksum file is generated from final artifacts.
- [ ] A clean-environment smoke test has been completed.
