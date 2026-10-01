# HCRProbeForge validation report

## Release

This report describes the reproducible checks for HCRProbeForge 1.3.7. It
contains no user paths, reference genomes, transcript sequences, result files,
or machine-specific state.

## Changes covered

The release keeps the established mature-transcript design pipeline and adds
regression coverage for the following correctness and release-safety fixes:

- annotation SQLite databases are closed cleanly after atomic creation and can
  be validated/reopened immediately;
- Pre-mRNA one-pass output filters exonic, boundary-spanning, and
  non-selected-intron candidates, then packs the eligible intervals without
  overlap before final oligo QC and order-table export;
- ready local references work with arbitrary registered aliases during
  Pre-mRNA preparation;
- explicitly versioned accessions cannot resolve to a different accession
  version;
- external compressed genome/annotation inputs are not deleted by cache
  normalization;
- reference files replaced during a failed index build are restored from a
  rollback backup;
- existing files in a visible output tree are restored if a failed rerun
  overwrote them;
- Pre-mRNA quota and fallback selection cannot reselect an interval overlapping
  a pair already selected from another quota pass;
- intronic Pre-mRNA geometric capacities count non-overlapping tiles separately
  within each eligible intron and contiguous unambiguous sequence segment;
- unambiguous upstream exclusive coordinate suffixes are normalized to the
  package's inclusive `start`/`end` convention;
- unresolved manifest genes/accessions remain failed rows and invalid preflight
  rows return nonzero instead of being silently treated as species skips; and
- similarly sized map boxes use one label placement mode, with dynamic lanes
  for labels that do not fit inside; probe boxes within roughly 15% of one
  another in rendered width are grouped together, including adjacent adaptive
  tile tiers, while materially different groups may use different modes.
  Mature feature labels are evaluated independently, and narrow Pre-mRNA
  exon/intron labels are omitted rather than emitted as external callouts;
- the webapp explains that normal runs require a selected-species genome index,
  identifies **Build a genome index** as the preparation step, and warns that
  disabling genome masking removes off-target screening; and
- the release bundle has a dependency-free archive hygiene checker for wheel
  and sdist contents.

## Automated checks

Run these commands from the source checkout:

```bash
python -m compileall -q src tests
PYTHONPATH=src python -m unittest discover -s tests -p 'test_*.py'
python -m ruff check src tests
python -m pip wheel --no-deps --no-build-isolation --wheel-dir dist .
python -m build --sdist --outdir dist .
python -m twine check dist/*
python tools/check_release_artifacts.py dist/*
```

The unit tests cover GFF3 transcript/exon linking, forward/reverse intron
parsing and selection, transcript accession matching, SQLite schema validation
and atomic cleanup, compressed-reference handling, failed-build restoration,
reference download validation, custom-species lifecycle behavior, and webapp
security/form rendering. They also cover intronic interval non-overlap,
coordinate-label normalization, grouped probe-label placement across adjacent
short-tile tiers, and manifest failure classification.

The release build must contain one version across project metadata, the CLI,
and the webapp; declare `Requires-Python: >=3.10`; include the MIT license,
public documentation, package assets, and tests; and exclude generated
caches, reference genomes, SQLite files, result folders, absolute local paths,
`.DS_Store`, and `__MACOSX` metadata.

## Deployment checks

Before publishing, install the wheel into a clean environment and run:

```bash
hcrprobeforge --version
hcrprobeforge --help
hcrprobeforge index list
```

For a full deployment smoke test, run one small mature-transcript design, one
Pre-mRNA design against a verified matching genome/GFF3/index, one all-channel
design, one plot, one QC run, and one local or NCBI-backed index operation.
These checks require the HCRProbeDesign dependency/executables installed by
pip, Bowtie2 installed separately, Primer3 support, and any required NCBI
access; those external programs and reference indexes are not bundled in this
source distribution.

The biological correctness of a Pre-mRNA result depends on pairing the genome
FASTA, GFF3 annotation, assembly accession, transcript model, and registered
Bowtie2 index from the same assembly. HCRProbeForge records that provenance
and validates the annotation database structure before reuse.
