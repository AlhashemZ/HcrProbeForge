# Contributing to HCRProbeForge

Thank you for helping improve HCRProbeForge. Changes that affect probe
selection, transcript orientation, annotation mapping, reference management,
or output semantics should include a focused regression test and an update to
the user documentation when the public behavior changes.

## Development setup

Create an environment and install the checkout in editable mode:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e '.[dev]'
```

The runtime dependency `hcrprobedesign` also supplies the external
`designProbes` and `buildGenomeIndex` commands used by full design/reference
workflows. Tests that exercise those workflows may need the corresponding
HCRProbeDesign installation and Bowtie2 on `PATH`.

## Before opening a pull request

Run the canonical validation sequence in
[`DEVELOPMENT.md`](DEVELOPMENT.md). It covers compilation, the unit-test suite,
Ruff, wheel and sdist builds, metadata validation, and release-artifact
inspection. If a tool is unavailable, report that explicitly rather than
silently omitting the check.

## Scientific changes

Keep mature-transcript behavior and file contracts stable unless a change is
deliberately planned and documented. For Pre-mRNA changes, test at least:

- forward- and reverse-strand transcript models;
- exon/intron coordinate mapping and transcript orientation;
- selected-intron filtering and boundary exclusion;
- whole-mode eligibility when exons and boundaries are intended to be valid;
- exact accession-version matching; and
- annotation-database creation, validation, reuse, and cancellation behavior.

Tests should use small synthetic FASTA/GFF3 fixtures where possible. Do not
commit downloaded genomes, annotation databases, generated result folders,
credentials, or machine-specific paths.

## Documentation changes

The README is for general users. Keep release engineering, private paths,
internal module names, and maintainer procedures in this file or
`DEVELOPMENT.md`. Commands in public documentation must correspond to an
installed entry point or an explicitly documented Python invocation.

When output behavior changes, document which files are candidates, which files
are order-ready, whether non-overlap selection has happened, and which table is
used for final QC.

## Issue reports

Include the HCRProbeForge version, Python version, operating system, command,
species/assembly metadata, and a redacted traceback. For Pre-mRNA reports,
include the transcript accession-version and the genome/GFF3 assembly identity.
