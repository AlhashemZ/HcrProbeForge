# Third-party notices

HCRProbeForge is distributed under the MIT License. The package does not
relicense its dependencies: each dependency keeps the copyright and license
terms supplied by its own authors. This file is a release summary, not a
replacement for the license files shipped by or installed with those projects.

## Python packages

The runtime dependency set declared in `pyproject.toml` includes:

| Package | Role in HCRProbeForge |
|---|---|
| `hcrprobedesign` | HCR probe candidate generation and reference-index helpers; also provides the `designProbes` and `buildGenomeIndex` commands. |
| `primer3-py` | Primer3-backed oligo hairpin, self-dimer, and heterodimer checks. |
| `requests` | NCBI and other HTTP requests. |
| `pandas` | Tabular input/output and manifest processing. |
| `openpyxl` | Excel manifest support. |
| `matplotlib` | Probe-map rendering. |
| `biopython` | Sequence and annotation utilities. |
| `pysam` | Compressed genomic/reference-file support. |
| `beautifulsoup4` | HTML/report parsing support. |
| `pyyaml` | Configuration and metadata parsing support. |

The exact license and copyright notices for the installed versions are
available in each distribution's `METADATA` file and license files. Check
those records when redistributing a locked environment because dependency
license metadata and versions can change after this release. In particular,
review the installed `hcrprobedesign`, `primer3-py`, and Primer3 notices rather
than relying on a time-dependent statement in this document.

## External executable

Bowtie2 is not bundled. It is an external executable used for genomic
specificity screening and index construction; obtain and distribute it under
the license terms supplied by the Bowtie2 project; the Bowtie2 distribution
used for this release is GPLv3-licensed. Verify the exact installed release
before redistribution.

## Biological data and HCR reagents

Reference sequences and annotations retrieved from NCBI/RefSeq remain subject
to NCBI data-use, attribution, and source-database terms. Keep the assembly
accession, annotation release, transcript accession-version, retrieval date,
and any local modifications with published designs.

HCR v3 initiator halves and amplifier/channel sequences are biological reagent
sequences supplied for use with the relevant HCR v3 system. Their use may be
subject to Molecular Instruments' documentation, reagent, trademark, or other
intellectual-property terms. HCRProbeForge does not grant rights to those
sequences or reagents; users are responsible for confirming the terms that
apply to their intended use.
