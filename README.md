# HCRProbeForge

<p align="center">
  <img src="https://raw.githubusercontent.com/AlhashemZ/hcrprobeforge/main/src/hcrprobeforge/assets/hcrprobeforge-logo.svg" alt="HCRProbeForge" width="520">
</p>

HCRProbeForge is a cross-platform package for designing, checking, curating,
plotting, and exporting HCR v3 split-initiator probe sets. It provides a
command-line interface and a local browser interface. The CLI, helper commands,
and webapp use the same scientific pipeline and file naming conventions, but
their default selection policies differ: direct `hcrprobeforge` designs are
one-pass by default, while the all-channel helper, manifest helper, and webapp
enable smart auto-curation by default.

The package accepts a gene symbol, a RefSeq transcript accession, one FASTA
record, or a gene/accession manifest. It retrieves or reuses exact NCBI
transcripts, generates channel-specific candidates with HCRProbeDesign and
Bowtie2, runs Primer3-based oligo QC, performs the configured curation and
selection, draws probe maps, and writes candidate and selected order tables
with metadata. Smart-curated outputs and one-pass Pre-mRNA outputs include a
selected non-overlapping order table. A one-pass mature-transcript run
deliberately preserves the native upstream candidate table and does not apply
final non-overlap selection.

The webapp defaults to a loopback-only local server. Sequences, probe tables,
and results stay on the computer running HCRProbeForge. Browser state-changing
requests use a per-server capability token and same-origin checks; non-loopback
binding requires an explicit access token. NCBI is contacted only when a
workflow needs a transcript, annotation, or genome assembly.

## Contents

- [What the package does](#what-the-package-does)
- [Installation](#installation)
  - [Choose a working directory](#choose-a-working-directory)
  - [macOS](#macos)
  - [Windows](#windows)
  - [Linux](#linux)
  - [Start the webapp later](#start-the-webapp-later)
  - [Installing from a release wheel](#installing-from-a-release-wheel)
  - [Install from the GitHub source with environment.yml](#install-from-the-github-source-with-environmentyml)
  - [Installation without conda](#installation-without-conda)
  - [Installing from a source checkout](#installing-from-a-source-checkout)
  - [Updating without reinstalling dependencies](#updating-without-reinstalling-dependencies)
- [Verify the installation](#verify-the-installation)
- [Species and genome indexes](#species-and-genome-indexes)
  - [Built-in species presets](#built-in-species-presets)
  - [Add a custom species preset](#add-a-custom-species-preset)
  - [Build a preset index](#build-a-preset-index)
  - [Use another assembly](#use-another-assembly)
  - [Build from a local genome FASTA](#build-from-a-local-genome-fasta)
  - [Where indexes and metadata are stored](#where-indexes-and-metadata-are-stored)
- [Launch the webapp](#launch-the-webapp)
  - [The five workflows](#the-five-workflows)
  - [Input inspection and progress](#input-inspection-and-progress)
  - [Webapp workflow instructions](#webapp-workflow-instructions)
  - [Webapp security and networking](#webapp-security-and-networking)
- [Manifest input and resume contract](#manifest-input-and-resume-contract)
- [CLI quick start](#cli-quick-start)
  - [Design one target](#design-one-target)
  - [Design a Pre-mRNA target](#design-a-pre-mrna-target)
  - [Design all linked transcripts](#design-all-linked-transcripts)
  - [Design selected channels](#design-selected-channels)
  - [Design a manifest](#design-a-manifest)
  - [Plot an existing table](#plot-an-existing-table)
  - [Run QC on an existing table](#run-qc-on-an-existing-table)
- [Complete CLI option reference](#complete-cli-option-reference)
- [Option interactions and reproducible command patterns](#option-interactions-and-reproducible-command-patterns)
- [How transcript selection works](#how-transcript-selection-works)
- [How the design pipeline works](#how-the-design-pipeline-works)
- [HCR channels and final-order conventions](#hcr-channels-and-final-order-conventions)
- [Auto-curation technical details](#auto-curation-technical-details)
- [Candidate ranking and auto_score](#candidate-ranking-and-auto_score)
- [Final oligo QC thresholds](#final-oligo-qc-thresholds)
- [Annotation and coordinate conventions](#annotation-and-coordinate-conventions)
- [Output layout](#output-layout)
- [Understanding the output files](#understanding-the-output-files)
- [Caching and reproducibility](#caching-and-reproducibility)
- [Troubleshooting](#troubleshooting)
- [Architecture and data flow](#architecture-and-data-flow)
- [Failure safety and reproducibility](#failure-safety-and-reproducibility)
- [Preserving provenance and getting help](#preserving-provenance-and-getting-help)
- [License and citation](#license-and-citation)
- [Additional project documents](#additional-project-documents)

## What the package does

For a normal NCBI-backed design, HCRProbeForge:

1. Validates the input and species/reference pairing.
2. Resolves a gene symbol to a NCBI Gene record and linked RefSeq RNA
   records, or validates the exact accession supplied by the user.
3. Retrieves the exact transcript FASTA and feature annotation, or validates
   a user-supplied single-record FASTA.
4. Reuses a valid versioned local transcript cache.
5. Calls HCRProbeDesign with the requested channel, tile, thermodynamic,
   genomic-specificity, and candidate-count settings.
6. Runs final-oligo Primer3 QC on the actual initiator-bearing P1/P2 oligos,
   unless --no-oligo-qc is requested.
7. If smart curation is enabled, applies the non-overlap selection,
   adaptive-tier, coverage, and cross-pair review rules. A normal one-pass
   mature-transcript run preserves the native candidate table. A normal
   one-pass Pre-mRNA run filters candidates to the requested genomic scope and
   enforces interval non-overlap before writing its order-ready table.
8. Writes the appropriate candidate or selected-pair tables, IDT order files,
   maps, QC workbooks, reports, logs, and run metadata.

The scientific pipeline is the same from the CLI and webapp. The webapp
collects options, starts the same functions, shows progress, and presents
the current run's files. Defaults are entry-point-specific: direct
`hcrprobeforge` is one-pass unless smart curation is requested, whereas
`hcrprobeforge-all-channels`, `hcrprobeforge-manifest`, and the webapp select
smart curation by default. To obtain comparable output behavior, specify
`--auto-curate-if-needed` or the one-pass setting explicitly.

## Installation

Use a dedicated conda environment for scientific executables, then install
the Python package with pip. The commands use python -m pip so pip belongs to
the active interpreter.

### Choose a working directory

Choose any directory where you want HCRProbeForge to save project results.
The paths below are examples only; replace them with a directory that suits
your project. The directory is a project root: HCRProbeForge creates its
`runs/` and `cache/` folders inside it. A descriptive name makes it clear that
the directory will contain all outputs from the project. For example:

```bash
mkdir -p "$HOME/hcrprobeforge_projects/hcr_probe_design_project"
cd "$HOME/hcrprobeforge_projects/hcr_probe_design_project"
```

The webapp uses the results directory selected in its form. CLI commands use
the path supplied to `--outdir`; if it is relative, it is resolved from the
current working directory. If the webapp or CLI receives
`$HOME/hcrprobeforge_projects/hcr_probe_design_project`, all run, cache, plot,
QC, and manifest workflow outputs are nested below that directory. The webapp
also maintains a small `index_metadata/` registry copy under the parent of the
selected project root (`<project parent>/index_metadata/`), rather than inside
the project root. Environment variables can relocate reference assets or a
direct-CLI cache elsewhere. A folder named `hcr_results` is still supported as
the default, but it is not automatically added to a user-selected path.

The first genome-masked design also requires a ready registered index. Build
the selected species index before the first design, or deliberately use
`--no-genomemask` when genomic specificity screening is not wanted. For the
default preset, the usual first setup is:

```bash
hcrprobeforge index fetch --species xtr --threads 4
```

Index downloads and builds are assembly-dependent. As a planning estimate,
allow at least a few GB of free disk for a small reference and tens of GB for
large vertebrate assemblies, including temporary files; small indexes may
finish in minutes, while large downloads/builds can take tens of minutes or
hours and may need roughly 8 GB or more of available RAM. Check the assembly
provider's size and the machine's available resources before starting. The
package does not bundle genome indexes.

### macOS

Install Miniforge once. The following commands select the correct Apple
Silicon or Intel installer automatically:

```bash
curl -L "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-MacOSX-$(uname -m).sh" -o Miniforge3.sh
bash Miniforge3.sh -b -p "$HOME/miniforge3"
macos_shell="${SHELL##*/}"; macos_shell="${macos_shell:-bash}"
"$HOME/miniforge3/bin/conda" init "$macos_shell"
source "$HOME/miniforge3/bin/activate"
```

Open a new Terminal window, change to your chosen working directory, and
create the environment:

```bash
cd "$HOME/hcrprobeforge_projects/hcr_probe_design_project"
conda create -n hcrprobe -c conda-forge -c bioconda python=3.11 pip bowtie2 pysam -y
conda activate hcrprobe
python -m pip install hcrprobeforge
hcrprobeforge-web
```

Start the webapp in later sessions with:

```bash
conda activate hcrprobe
cd "$HOME/hcrprobeforge_projects/hcr_probe_design_project"
hcrprobeforge-web
```

macOS keeps its native browser and folder-opening behavior.

### Windows

HCRProbeForge can run from native Windows Python when Python, Bowtie2, and the
HCRProbeDesign command-line tools are installed and available on `PATH`. WSL2
with Ubuntu is recommended because the Conda/bioconda installation below
provides the external scientific executables most consistently. The WSL2
instructions are not a separate HCRProbeForge mode; they run the same package
inside Linux on Windows.

For the recommended WSL2 setup, install Ubuntu from the Microsoft Store, or
run this in an elevated PowerShell window:

```powershell
wsl --install -d Ubuntu
```

Restart if prompted, open Ubuntu, create the Linux user, and run the following
inside the Ubuntu terminal:

```bash
curl -L "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-$(uname -m).sh" -o Miniforge3.sh
bash Miniforge3.sh -b -p "$HOME/miniforge3"
"$HOME/miniforge3/bin/conda" init bash
source "$HOME/miniforge3/bin/activate"
conda create -n hcrprobe -c conda-forge -c bioconda python=3.11 pip bowtie2 pysam -y
conda activate hcrprobe
python -m pip install hcrprobeforge
mkdir -p "$HOME/hcrprobeforge_projects/hcr_probe_design_project"
cd "$HOME/hcrprobeforge_projects/hcr_probe_design_project"
hcrprobeforge-web
```

HCRProbeForge detects WSL and opens the webapp in the Windows browser through
`explorer.exe`. Folder buttons use `wslpath` plus Windows Explorer. The URL
is always printed if manual opening is needed.

Ubuntu stores Linux files inside the WSL distribution. `explorer.exe .` opens
the current folder in Windows Explorer; the same files can be reached at
`\\wsl.localhost\Ubuntu\home\<your-linux-user>` (the older
`\\wsl$\Ubuntu\home\<your-linux-user>` form also works). Windows drives are mounted under
`/mnt`, so `C:\Users\Name\Documents` is normally
`/mnt/c/Users/Name/Documents`. Large reference files are usually faster under
the Linux filesystem in `~/` than under `/mnt/c`.

Start the webapp in later Ubuntu sessions with:

```bash
conda activate hcrprobe
cd "$HOME/hcrprobeforge_projects/hcr_probe_design_project"
hcrprobeforge-web
```

For a native PowerShell installation, install Python 3.10 or newer, create and
activate a virtual environment, and install HCRProbeForge with
`python -m pip install hcrprobeforge`. Pip installs the Python-side
`hcrprobedesign` dependency, including `designProbes` and
`buildGenomeIndex`. Install the external Bowtie2 executables separately and
put them on `PATH`:

```powershell
py -3 -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip hcrprobeforge
where.exe bowtie2
where.exe bowtie2-build
where.exe designProbes
where.exe buildGenomeIndex
```

If the installed HCRProbeDesign or Bowtie2 distribution does not provide a
working native Windows executable, use the WSL2 route above.

### Linux

```bash
curl -L "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-$(uname -m).sh" -o Miniforge3.sh
bash Miniforge3.sh -b -p "$HOME/miniforge3"
"$HOME/miniforge3/bin/conda" init bash
source "$HOME/miniforge3/bin/activate"
mkdir -p "$HOME/hcrprobeforge_projects/hcr_probe_design_project"
cd "$HOME/hcrprobeforge_projects/hcr_probe_design_project"
conda create -n hcrprobe -c conda-forge -c bioconda python=3.11 pip bowtie2 pysam -y
conda activate hcrprobe
python -m pip install hcrprobeforge
hcrprobeforge-web
```

Desktop Linux uses `xdg-open` or `gio` for folders; headless systems can use the
printed URL. To upgrade later, use the cross-platform instructions in
[Updating without reinstalling dependencies](#updating-without-reinstalling-dependencies).

### Start the webapp later

After installation, start HCRProbeForge in any later session with:

```bash
conda activate hcrprobe
cd "$HOME/hcrprobeforge_projects/hcr_probe_design_project"
hcrprobeforge-web
```

The command opens the local webapp in the default browser when the operating
system supports that behavior and always prints the URL in the terminal.
`hcrprobeforge-web` is the HCRProbeForge launcher. Installing
`hcrprobeforge` also installs the Python package `hcrprobedesign`, which
provides the scientific `designProbes` and `buildGenomeIndex` tools. Bowtie2
(`bowtie2` and `bowtie2-build`) remains an external executable and is supplied
by the Conda environment in the examples above.

### Install from the GitHub source with `environment.yml`

This option creates a Conda environment and installs the local checkout
through pip. It is not a separate Conda package:

```bash
git clone https://github.com/AlhashemZ/hcrprobeforge.git
cd hcrprobeforge
conda env create -f environment.yml
conda activate hcrprobe
hcrprobeforge-web
```

After updating the checkout:

```bash
cd hcrprobeforge
git pull
conda activate hcrprobe
python -m pip install --upgrade --editable .
```

### Installation without conda

Conda is recommended because it supplies Bowtie2 consistently, but a Python
virtual environment also works if `bowtie2` and `bowtie2-build` are already on
`PATH`. Pip installs HCRProbeDesign as an HCRProbeForge dependency, so its
`designProbes` and `buildGenomeIndex` commands are installed with the package:

```bash
python3 -m venv .venv
source .venv/bin/activate                    # macOS/Linux
python -m pip install --upgrade pip
python -m pip install .
```

On Windows PowerShell, activate with
`.\.venv\Scripts\Activate.ps1`. If the shell blocks activation, use a
per-user execution-policy change or launch the environment's Python directly.
This route does not install Bowtie2; install it separately and verify it with
`where.exe` (Windows) or `command -v` (macOS/Linux). If pip reports that a
dependency must be built locally, follow that dependency's build instructions
for the selected operating system.

### Installing from a release wheel

Download a wheel from [PyPI](https://pypi.org/project/hcrprobeforge/) or the
GitHub [Releases](https://github.com/AlhashemZ/hcrprobeforge/releases) page.
If you downloaded the prepared publication bundle, unzip it and use the wheel
from its `pypi/` directory:

```bash
cd /path/to/release-bundle-<release>
python -m pip install --upgrade ./pypi/hcrprobeforge-<release>-py3-none-any.whl
```

GitHub Releases also publish the wheel as a release asset. If you downloaded
that asset directly, change to the directory containing the downloaded wheel:

```bash
cd /path/to/downloaded-wheel
python -m pip install --upgrade ./hcrprobeforge-<release>-py3-none-any.whl
```

If the active conda environment already contains all required Python and
external dependencies:

```bash
python -m pip install --upgrade --no-deps ./hcrprobeforge-<release>-py3-none-any.whl
```

Use `--no-deps` only after checking the environment. It updates HCRProbeForge
without asking pip to resolve dependencies.

### Installing from a source checkout

From an unpacked source directory:

```bash
conda activate hcrprobe
python -m pip install --upgrade .
```

This installs the checkout as a normal package. If you are modifying the
source tree, use `--editable` instead; restart the CLI process or webapp after
changing code.

### Updating without reinstalling dependencies

Pip's default dependency strategy is already `only-if-needed`; the explicit
flag below documents that choice but is not required. To update a normal
installation:

```bash
python -m pip install --upgrade --upgrade-strategy only-if-needed hcrprobeforge
```

For a local wheel:

```bash
python -m pip install --upgrade --upgrade-strategy only-if-needed \
  ./hcrprobeforge-<release>-py3-none-any.whl
```

Pip replaces HCRProbeForge and keeps already-satisfied dependency versions; it
changes a dependency only when the installed version does not satisfy the
package requirement. Repository build and validation procedures are documented
in [`DEVELOPMENT.md`](DEVELOPMENT.md).

## Verify the installation

```bash
hcrprobeforge --version
hcrprobeforge --help
hcrprobeforge index --help
hcrprobeforge-web --help
hcrprobeforge-all-channels --help
hcrprobeforge-manifest --help
```

Then check `designProbes`, `buildGenomeIndex`, `bowtie2`, and `bowtie2-build`
with `command -v` on macOS/Linux or `where.exe` on Windows.

Primer3 is installed as the Python package `primer3-py`, not as a separate
executable. Check it with:

```bash
python -c "import primer3; print('primer3-py is available')"
```

## Species and genome indexes

Genomic masking is enabled by default. Use an index built from the same
organism and genome reference selected for transcript retrieval. If no
matching index is registered, build one first or provide an explicit Bowtie2
prefix with `--index`. `--no-genomemask` disables specificity screening and
### Built-in species presets

The package includes six built-in organism presets. Xenopus tropicalis uses the
HCRProbeDesign index alias `xtr`, while user-facing output uses
`Xenopus_tropicalis`. Preset aliases are short technical names accepted by the
CLI and index registry; project results, cache folders, and index-metadata
folders use the readable display name (sanitized only for filesystem safety).
`xtr` is the only built-in alias for Xenopus tropicalis; build or rebuild that
reference with `xtr`.

| CLI/webapp alias | Organism sent to NCBI | Assembly | NCBI accession | Results component |
|---|---|---|---|---|
| xtr | Xenopus tropicalis | UCB_Xtro_10.0 | GCF_000004195.4 | Xenopus_tropicalis |
| xla | Xenopus laevis | Xenopus_laevis_v10.1 | GCF_017654675.1 | Xenopus_laevis |
| zebrafish | Danio rerio | GRCz12ab | GCF_052040795.1 | Zebrafish |
| mouse | Mus musculus | GRCm39 | GCF_000001635.27 | Mouse |
| chicken | Gallus gallus | bGalGal1.mat.broiler.GRCg7b | GCF_016699485.2 | Chicken |
| human | Homo sapiens | GRCh38.p14 | GCF_000001405.40 | Human |

These are package defaults, not a permanent restriction. Record the assembly
and accession from reference.json with published probes. If a newer assembly
has consistent annotation and transcript support, add it with an explicit
assembly accession and a new index alias.

### Add a custom species preset

Organisms that are not built into the package can be added without changing
the scientific design pipeline. A custom preset has three required values:

- a stable alias used by the CLI and index registry;
- a readable display name shown in the webapp; and
- the exact NCBI organism name used for transcript and accession retrieval.

The assembly name is optional at preset-creation time. The NCBI assembly
accession is required for an NCBI download unless a local genome FASTA is
supplied. Provide one of those two genome sources before building the genome
index. The transcript side remains NCBI-backed through the organism name; the
local FASTA is used only to build the genome index.

From the CLI:

```bash
hcrprobeforge species add \
  --alias octopus_v1 \
  --display-name "Octopus vulgaris" \
  --organism "Octopus vulgaris" \
  --assembly-accession GCF_<ASSEMBLY_ID>.<VERSION> \
  --assembly-name Octopus_v1

# Or save a local genome FASTA with the preset:
hcrprobeforge species add \
  --alias octopus_local \
  --display-name "Octopus vulgaris" \
  --organism "Octopus vulgaris" \
  --assembly-name Octopus_local \
  --genome-fasta /path/to/octopus.fna

hcrprobeforge species list
```

Replace `GCF_<ASSEMBLY_ID>.<VERSION>` with a real versioned NCBI assembly
accession before running the command.

Then build the index using the existing index workflow:

```bash
hcrprobeforge index fetch --species octopus_v1 --threads 4
hcrprobeforge index build --species octopus_local --fasta /path/to/octopus.fna --threads 4
```

In the webapp, choose **Add a new species…** at the bottom of the Species
menu. The dialog marks required and optional fields explicitly and provides
tooltips for the alias, organism, assembly, and FASTA choices. An uploaded
genome FASTA is copied into the HCRProbeForge reference-data directory; the
original file is not modified. After the species is added, it remains in the
Species menu and can be selected in **Build a genome index**. The saved FASTA
is used automatically when the local FASTA field is left blank. The NCBI
organism field alone is not a genome source: it identifies the organism for
transcript retrieval, but does not select a reproducible assembly. The dialog
shows a red required marker beside the accession and explains the local-FASTA
exception. Supply a versioned NCBI assembly accession or a genome FASTA before
building the index.
The webapp exposes an adjacent **×** control for every custom preset. Clicking
it opens a confirmation dialog rather than a browser alert. After the
confirmation checkbox is selected, deletion removes the preset from the
species registry, its saved custom FASTA, HCRProbeForge reference metadata,
and the exact unshared Bowtie2 index files registered for that custom preset.
It also removes obsolete HCRProbeDesign registrations and any now-empty alias
folder.
The dialog previews the paths and file count before confirmation. An index
shared with another installed reference, or a prefix outside the managed
HCRProbeDesign data directory, is retained for safety. Built-in presets are
protected and their indexes are never changed by custom-species deletion.

The Species menu renders one option per custom preset using its display name.
If a successful index build has a separate assembly alias in its metadata, that
alias is not added as a second user-facing species option.

Custom presets are stored in:

```text
~/.hcrprobeforge/references/species_presets.json
```

The location follows `HCRPROBEFORGE_REFERENCE_DIR` when that environment
variable is set. A custom preset does not bypass species validation: the
selected organism is still passed to NCBI, and exact accessions must belong to
that organism. If a transcript is unavailable, the normal validation error is
returned before scientific output or transcript cache data are retained.

### Build a preset index

The index command resolves and validates the NCBI Assembly record, downloads
the genomic FASTA and matching GFF3 annotation, passes the uncompressed FASTA
to HCRProbeDesign's buildGenomeIndex, and registers the resulting Bowtie2
prefix. The annotation SQLite lookup used by Pre-mRNA is optional because many
users only need mature-transcript designs:

```bash
hcrprobeforge index fetch --species xtr --threads 4
hcrprobeforge index fetch --species xla --threads 4
hcrprobeforge index fetch --species zebrafish --threads 4
hcrprobeforge index fetch --species mouse --threads 4
hcrprobeforge index fetch --species chicken --threads 4
hcrprobeforge index fetch --species human --threads 4
hcrprobeforge index list

# Also prepare the annotation database now for future Pre-mRNA/intron designs
hcrprobeforge index fetch --species xtr --with-annotation-db --threads 4
```

The examples use four build threads, which is also the default. Choose a value
appropriate for the available CPU cores and memory.

A ready reference is reused. Use `--force` to rebuild it. The equivalent
webapp workflow is the Build index tab; leave assembly fields blank to use
the preset. When you click Run, the webapp asks whether you will need
Pre-mRNA/intron design. Selecting the option builds the annotation database
now and may take several minutes; leaving it clear builds only the genome
index and defers SQLite construction until the first Pre-mRNA run. Its result
page reports the full index prefix and metadata JSON.

Assembly downloads use the FTP directory returned by NCBI Assembly. The FASTA
filename is derived from that directory and normalized to NCBI's underscore
convention, so assembly names such as `Release 6 plus ISO1 MT` do not produce
invalid `%20` filename requests. The downloader validates the requested
accession in the ESummary response, retries transient HTTP failures, writes to
a temporary file, verifies that the compressed content can be decompressed,
and only then starts `buildGenomeIndex`. If the usual filename is not present,
it checks the Assembly directory listing for another `*_genomic.fna.gz` file.
If no usable file can be found, the error names the accession and directs you
to verify the Assembly record or use the local-FASTA index workflow. Successful
builds retain the uncompressed `genome.fna` and `annotation.gff` when NCBI
provides them. `annotation.sqlite` is retained only when
`--with-annotation-db` was selected; otherwise it is deferred until the first
Pre-mRNA run. Temporary compressed downloads and their `.gz` copies are not
retained. A genome index can still be registered when NCBI temporarily has no
annotation, but Pre-mRNA will report that the matching annotation must be
supplied or downloaded before design.

### Use another assembly

First verify that the desired NCBI assembly has appropriate annotation and
transcript support. Then provide its accession, optionally its expected name,
and preferably a new alias:

```bash
hcrprobeforge index fetch \
  --species zebrafish \
  --assembly-accession GCF_<ALTERNATE_ASSEMBLY_ID>.<VERSION> \
  --assembly-name <EXPECTED_ASSEMBLY_NAME> \
  --index-alias zebrafish_alternate \
  --threads 4
```

Replace the accession and expected name with a real compatible assembly. The
example intentionally uses placeholders rather than implying that an older or
built-in assembly is newer.

The same pattern works for either Xenopus preset, mouse, chicken, zebrafish,
or human. A new alias keeps an existing reference available. Reuse an alias
with `--force` only when intentionally replacing that registration.
`hcrprobeforge index list` shows registered
references and their status.

The webapp exposes these fields under Build index:

- NCBI assembly accession selects the exact record;
- Expected assembly name prevents an accidental mismatch;
- Index alias is the HCRProbeDesign registration name;
- Build threads controls index-building parallelism;
- Rebuild existing index corresponds to `--force`.

### Build from a local genome FASTA

```bash
hcrprobeforge index build \
  --species local_xtr \
  --fasta /path/to/genome.fna \
  --annotation /path/to/matching/annotation.gff3 \
  --assembly MyAssembly \
  --threads 4
```

Use a matching organism during design:

```bash
hcrprobeforge sox9 --species local_xtr \
  --organism "Xenopus tropicalis" --channel B3
```

If a matching GFF3 was supplied, the same registered local reference can be
used for Pre-mRNA design. The accession is still retrieved and validated from
NCBI using the organism name; the local genome and annotation provide the
genomic target and index:

```bash
hcrprobeforge sox9 --species local_xtr \
  --organism "Xenopus tropicalis" --target-type pre-mrna --channel B3
```

Or use a direct Bowtie2 prefix:

```bash
hcrprobeforge sox9 --species xtr \
  --index /path/to/index/prefix --channel B3
```

The prefix is the path before the Bowtie2 .1.bt2 or .1.bt2l suffix.

For a local build, `--annotation` is optional when you only need mature
transcript designs. Add `--with-annotation-db` only when the supplied GFF3 is
from the same assembly and you want to prepare Pre-mRNA lookup during the
index build:

```bash
hcrprobeforge index build \
  --species local_xtr \
  --fasta /path/to/genome.fna \
  --annotation /path/to/matching/annotation.gff3 \
  --with-annotation-db \
  --assembly MyAssembly \
  --threads 4
```

If you leave the option off, the first Pre-mRNA run builds the database from
the matching annotation. Existing ready references and compatible local
reference assets are reused before a new download is attempted.

The index and species commands use separate option groups. The assembly flags
are intentionally command-specific: `index fetch` uses
`--assembly-accession` to identify an NCBI record and `--assembly-name` to
verify its returned display name, while `index build` uses `--assembly` as the
display name recorded for a local FASTA build. They are not interchangeable.

| Command | Purpose |
|---|---|
| `hcrprobeforge index list` | List managed reference metadata and readiness. |
| `hcrprobeforge species list` | List built-in and user-created species presets. |
| `hcrprobeforge species add` | Create a persistent custom species preset. |
| `hcrprobeforge index fetch` | Download and build the selected NCBI assembly. |
| `hcrprobeforge index build` | Build an index from a local genome FASTA. |

| Command | Option | Meaning |
|---|---|---|
| `index fetch` | `--assembly-accession ACCESSION` | Exact versioned NCBI `GCF_` or `GCA_` assembly accession; defaults to the preset accession. |
| `index fetch` | `--assembly-name NAME` | Expected NCBI assembly name; the build stops if NCBI returns a different name. |
| `index fetch` | `--index-alias ALIAS` | Stable HCRProbeDesign alias for the registered index. Use a new alias for an alternate assembly. |
| `index build` | `--assembly NAME` | Display name recorded for a local FASTA build. |
| `species add` | `--genome-fasta FILE` | Local genome FASTA saved with a custom species preset. |
| `index fetch` / `index build` | `--threads N` | Number of index-building threads; default 4. |
| `index fetch` / `index build` | `--force` | Rebuild and replace an existing registration. |
| `index fetch` / `index build` | `--with-annotation-db` | Build the optional SQLite genomic-annotation lookup during index construction. Without it, construction is deferred until the first Pre-mRNA run. |
| `index fetch` | `--email EMAIL` / `--api-key KEY` | NCBI request credentials; `NCBI_EMAIL` and `NCBI_API_KEY` are used by default. |

Run `hcrprobeforge index fetch --help` or `hcrprobeforge index build --help` to
see the installed command's complete syntax.

### Where indexes and metadata are stored

The actual registered Bowtie2 files belong to HCRProbeDesign. Its default
data directory is:

```text
~/.hcrprobedesign/
```

Set `HCRPROBEDESIGN_DATA_DIR` to relocate it. The exact registered prefix is
printed by the index command and webapp result page.

HCRProbeForge reference records are stored by default in:

```text
~/.hcrprobeforge/references/<species>/<assembly>/reference.json
```

Set `HCRPROBEFORGE_REFERENCE_DIR` to relocate this registry. The webapp also
writes a metadata copy under the parent of the selected project root. If the
project root is `$HOME/hcrprobeforge_projects/hcr_probe_design_project`, the
phrase
`<project parent>` means `$HOME/hcrprobeforge_projects`:

```text
<project parent>/index_metadata/<display-name component>/<assembly>/reference.json
```

HCRProbeForge reference assets are kept separately from Bowtie2 indexes:

```text
~/.hcrprobeforge/references/<species>/<assembly>/
├── reference.json
├── genome.fna
├── annotation.gff          # when available or supplied
└── annotation.sqlite       # optional; built at index time or first Pre-mRNA run
```

Reusable project caches are kept below the selected project root rather than
inside the hidden reference directory:

```text
<project>/cache/<species>/<workflow>/
├── transcripts/<accession>/             # mature transcript cache
└── premrna/
    ├── models/<assembly>/               # parsed transcript/exon models
    └── targets/<assembly>/<accession>/  # assembled genomic targets by mode
        ├── intronic_all.fa/.json
        └── whole_all.fa/.json
```

The annotation database is a persistent reference asset because it is shared
by many targets; model and sequence caches are project workflow data. A
complete database is identified by its schema marker, final tables, indexes,
and SQLite structural validation. The normal readiness check is deliberately
fast; a full integrity audit is available for diagnostic validation. An
interrupted or incomplete `annotation.sqlite` is rebuilt atomically before a
Pre-mRNA design proceeds.

The reference-asset directory follows `HCRPROBEFORGE_REFERENCE_DIR`.
Compressed `.gz` files are used only as download/input sources and are not
retained beside these uncompressed assets.

The **Build a genome index** tab also provides **Open index folder**, which
opens HCRProbeDesign's persistent data directory in the operating system's
file browser. The actual Bowtie2 files are normally below:

```text
~/.hcrprobedesign/indices/<index alias>/
```

This location follows `HCRPROBEDESIGN_DATA_DIR` when that environment variable
is set. The webapp rechecks the Bowtie2 files while it is open, so deleting an
index folder changes the species badge from ready to not registered. A design
or manifest started afterward stops before transcript lookup and explains how
to rebuild the missing index. Closing or reloading a browser tab does not stop
the local server, so a design or index build cannot be interrupted by a tab
lifecycle event. Stop the server with Ctrl-C in its terminal or use an
explicit authenticated shutdown request after all jobs finish. Species
presets cannot be added or removed while a design, QC, manifest, plot, or
index job is queued or running.

The species component is derived from the display name, while the assembly
component identifies the built reference. Short index aliases remain in the
HCRProbeDesign registry and Bowtie2 prefix; they are not used as project
species folder names. `index_metadata` is a metadata directory, not the genome
index directory.

When a custom species is removed from the Species menu, the confirmation
dialog lists the saved FASTA, HCRProbeForge metadata, and exact Bowtie2 files
that will be deleted, plus obsolete HCRProbeDesign entries and empty alias
folders. Only unshared prefixes below the HCRProbeDesign data root are
eligible. Built-in preset indexes, shared custom indexes, and out-of-root
prefixes are retained. Config files are updated atomically with a backup.
If a previous installation left a custom alias registered without a complete
index, the next custom build repairs that incomplete registration
automatically.

## Launch the webapp

Activate the environment containing HCRProbeForge, HCRProbeDesign, Bowtie2, and
the Python dependencies:

```bash
hcrprobeforge-web
```

Open http://127.0.0.1:8766/. If that port is already in use, HCRProbeForge
automatically tries the next available ports and prints the actual URL. You
can also choose a starting port explicitly:

```bash
hcrprobeforge-web --port 8780
hcrprobeforge-web --port 8780 --no-browser
```

The browser opens automatically unless `--no-browser` is supplied. The default
host is local-only:

```bash
hcrprobeforge-web --host 127.0.0.1 --port 8780
```

For an intentional network bind, provide a long access token and use the
printed access URL. The access token is required whenever `--host` is not a
loopback address:

```bash
export SERVER_TOKEN='<LONG_RANDOM_TOKEN>'
hcrprobeforge-web --host 0.0.0.0 --auth-token "$SERVER_TOKEN" --no-browser
```

`<LONG_RANDOM_TOKEN>` is a placeholder; replace it with a newly generated,
private value. The default loopback server does not require a user-supplied
token for its browser URL. For a loopback server, stop it with Ctrl-C in the
terminal; the authenticated `/shutdown` request below is for a server started
with `--auth-token`.

### The five workflows

The setup page uses tabs:

1. Design one target: a gene, optional accession, or one FASTA record; one or all
   channels.
2. Design a gene list: CSV, TSV, XLSX, XLSM, or text input.
3. Plot a probe table: an existing HCRProbeDesign TSV/text table to PNG and SVG.
4. QC an oligo table: an existing IDT/selected-pairs table to an Excel workbook.
5. Build a genome index: retrieve/register an NCBI assembly or build from a local FASTA.

Species and organism are near the top so transcript retrieval and genomic
specificity use the same reference. Advanced settings contain the same
transcript, design, curation, QC, and plotting options as the CLI.

When a preset already has a registered index on the computer, the selected
species field shows an **Index ready on this computer** badge. If
the badge says **Index not registered yet**, build the preset index first or
select an explicit Bowtie2 prefix in Advanced settings/CLI.

Before designing, normal runs need a genome index for the selected species. If
the badge says the index is not registered, use **Build a genome index** first.
Turning off genome masking skips this requirement but is not recommended,
because probes will not be screened for off-target matches.

The setup page reads the finished reference metadata and refreshes both badges
when it opens, when the browser window regains focus, and when a page is
restored from browser history. Setup rendering and the live status endpoint
use the committed metadata/path state rather than opening the SQLite file, so
returning after an index build is immediate even for large annotation files.
The design worker still performs full SQLite validation before querying
annotations. Dynamic setup, progress, result, and status responses are marked
non-cacheable, so returning with **Start another run** does not require a
manual browser refresh after an index build. If a browser has an existing
HCRProbeForge page open, close that tab and start the installed
`hcrprobeforge-web` process again after an upgrade; already-running Python
processes do not load updated package code.

The How to start guide can be reopened from the setup page. It explains
species selection, index preparation, designing, plotting, QC, and results.

### Input inspection and progress

Inspect input checks the selected workflow before a long operation starts. It
does not contact NCBI, build an index, run `designProbes`, or write scientific
output. The checks are intentionally the same shape checks used when the real
workflow is started:

| Workflow | What Inspect input checks |
|---|---|
| Design one target | Species/organism pairing, gene/accession/FASTA exclusivity, one FASTA record, supported bases, accession syntax, transcript-policy conflicts, and design/QC option values. |
| Design a gene list | Manifest file type, readable rows, recognized gene/accession/channel columns, structural row validity, duplicate jobs, and directory-name collisions. An optional `channel` column must contain B1-B5. |
| Plot a probe table | Probe-table file type, required HCRProbeDesign columns, numeric coordinates/metrics, non-empty rows, and any supplied transcript-length bound. |
| QC an oligo table | File type, readable table format, required name/sequence or P1/P2 columns, non-empty rows, and sequence validity through the same QC reader used by the run. |
| Build a genome index | Built-in or custom species preset, build-thread value, assembly accession syntax, local FASTA availability, and the selected preset’s assembly details. NCBI confirms that an accession/name pair exists during the build itself. |

Inspection cannot prove that a gene symbol or accession exists for an
organism without making the NCBI request. Manifest rows are therefore checked
again during execution. A true species mismatch is recorded as
`skipped_species_mismatch`; a mistyped gene or unresolved accession is recorded
as `failed`. Neither starts a design or retains new design/cache files, and
valid rows continue to run. Structural errors are shown before the run so they
can be corrected first.

Long jobs show the active phase and progress. Cancel run stops the current
download, index build, or external design process when possible, prevents the
next channel or manifest job from starting, and returns to setup.

For an index build, the progress phases identify source retrieval, annotation
database preparation when selected, Bowtie2 index construction, registration,
and final metadata validation. The optional annotation database is committed
atomically before the reference is marked ready. Therefore a **Pre-mRNA
database ready** badge means the registered reference metadata points to the
finished database; the first Pre-mRNA design still performs a read-only schema
validation before use.

### Webapp workflow instructions

For Design one target:

1. Select the species and check the index-status badge.
2. If no index is registered, use Build index first.
3. Choose one channel or All channels (B1-B5).
4. Enter a gene, an accession, or upload one FASTA record.
5. Optionally give the workflow a readable folder name.
6. Inspect input, adjust Advanced settings only if needed, and run.

For Design a gene list:

1. Select the species and channels.
2. Upload CSV, TSV, XLSX, XLSM, or text.
3. Use gene_symbol and accession headers, or one gene per line; accession may
   be blank.
4. Optionally add a `channel` column with B1-B5 values for multiplex rows;
   leave it blank when the global channel selection should apply.
5. Inspect the manifest and fix any reported structural error or directory
   collision. Species mismatches are reported as skips; unresolved genes and
   accessions are reported as failures.

For Plot a probe table, choose a pastel or minimal theme and select whether
probe blocks are coloured by GC, dTm, or order. A sibling
`*_transcript_resolution.json` may provide annotation metadata, but a
standalone table still needs `--transcript-length` when its full sequence
length cannot be inferred from the table. For QC an oligo table, choose the
input file and optionally provide the workbook name. The QC workflow reads the
orderable sequences and writes its workbook in the QC run folder.

The webapp command itself has these options:

| Command | Meaning |
|---|---|
| `hcrprobeforge-web` | Start a local server on `127.0.0.1:8766` and open the browser. |
| `hcrprobeforge-web --host HOST --port PORT` | Choose a bind address and starting port. If the requested port is occupied, ports through `PORT + 9` are tried. Non-loopback hosts require `--auth-token`. |
| `hcrprobeforge-web --auth-token <LONG_RANDOM_TOKEN>` | Set the per-process capability token; required for non-loopback binding. The value is a placeholder and the access URL is printed in the terminal. |
| `hcrprobeforge-web --allowed-host HOST` | Add an allowed HTTP `Host` header name, such as a trusted machine name or IP; repeat the option for additional names. This is not a firewall or HTTPS setting. |
| `hcrprobeforge-web --no-browser` | Start the server without opening a browser window. |

The all-channel, manifest, index, and species helper commands are documented
in the corresponding sections of this README.

### Webapp security and networking

The default server listens only on `127.0.0.1`, so other computers cannot
connect to it. Every POST action also requires the per-process capability
token. The browser receives that token in a strict cookie, and the webapp
checks the request Host and, when supplied, the Origin. These protections stop
a separate webpage from silently starting a run, deleting a custom species,
opening local folders, or shutting down the server.

When `--host` is a non-loopback address, `--auth-token` is mandatory; remote
access is not enabled merely by changing the bind address. The
terminal prints an access URL containing the token once; open that URL only in
the intended browser and do not share it. The built-in server uses HTTP, not
HTTPS, so do not expose it to an untrusted network. For remote use, prefer an
SSH tunnel or an HTTPS reverse proxy in front of the loopback server. Use
`--allowed-host` to restrict which HTTP `Host` names may be accepted.

Closing, reloading, or navigating away from a tab never sends a shutdown
request. The server remains available for other tabs and active worker threads.
An explicit shutdown request is accepted only while the job queue is idle;
while a job is queued, running, or cancelling it returns a busy response and
leaves the server running. For the default loopback server, pressing `Ctrl-C`
in the terminal is the simplest shutdown method. The server still creates an
internal per-process token for browser POST actions, but it is not printed
because the browser's strict cookie carries it. If you need an HTTP shutdown
request, start the server with an explicit token so the same token can be sent
from the shell:

```bash
export SERVER_TOKEN='<LONG_RANDOM_TOKEN>'
hcrprobeforge-web --auth-token "$SERVER_TOKEN" --no-browser &

curl -X POST \
  -H "X-HCRProbeForge-Token: $SERVER_TOKEN" \
  http://127.0.0.1:8766/shutdown
```

Replace the placeholder with the private token used to start the authenticated
server. Cancellation
is cooperative: the worker stops at its cancellation checks and may need time
to finish an external subprocess. It is not a guarantee that every external
process stops immediately.

The port fallback changes only where the unchanged webapp listens. It does not
alter transcript retrieval, genome selection, Bowtie2 indexing, probe
generation, QC, curation, selection, output paths chosen by the user, or any
other scientific behavior.

The webapp is a front end, not a second implementation. Before a submission it
disables fields that belong to another tab, validates the values that can be
checked locally, and then calls the same core or batch functions as the CLI.
The Inspect input button does not download an assembly or run `designProbes`.
It checks the shape of the supplied input so simple mistakes can be corrected
before a long job begins.

## Manifest input and resume contract

A manifest is a list of independent target requests. It can be a CSV, TSV,
XLSX, XLSM, or plain UTF-8 text file. The parser ignores completely blank lines
and lines beginning with `#`. A row containing other fields but no gene label
is reported as invalid. A header is optional. When a header is present, these
names are recognized:

| Purpose | Accepted header names |
|---|---|
| Gene label | `gene_symbol`, `gene`, `gene_name`, or `gene name` |
| Transcript | `accession`, `transcript`, `transcript_accession`, or `refseq_accession` |
| HCR channel | `channel`, `channels`, `hcr_channel`, or `hcr_channels` |

The default command-level species is `xtr`; pass `--species` for another
preset. The default channel selection is `B1,B2,B3,B4,B5`, and a blank channel
field expands to that selection. The all-channel helper uses the same five
channel default. Plain-text fields may be separated by whitespace, commas, or
literal tab characters; CSV and TSV files should use their normal delimiters.
The manifest helper enables smart auto-curation by default, while a direct
`hcrprobeforge` invocation is one-pass unless `--auto-curate-if-needed` is
given. The webapp's Smart auto-curation control is selected by default; clear
it for a one-pass design.

Without a header, the first column is the gene label and the optional second
column is interpreted as either an exact accession or a channel (`B1`-`B5`).
Headerless three-column rows use `gene_symbol`, `accession`, `channel` order.
For example, all of these forms are valid:

```text
# gene only; the selected --channels value applies
sox9

# gene plus an exact transcript
sox9 NM_001016853.2

# gene plus a per-row HCR channel
foxj1 B1

# gene, exact transcript, and per-row HCR channel
pkd2 NM_001088456.1 B2
```

A gene-only row is resolved using the selected transcript policy; a row with
an accession is validated as that exact transcript. An explicit accession is
never expanded by `--all-transcripts`. When a `channel` column is present, its
single `B1`-`B5` value overrides the command's `--channels` setting for that
row. Rows with a blank channel continue to use `--channels`, so existing
manifests behave exactly as before.

Before design starts, preflight assigns each row a safe path component while
retaining the original biological label in metadata and status tables. It
reports rows with missing required values, malformed accessions, invalid
channel values, exact duplicate jobs, an accession used under multiple labels,
and path collisions. Completely blank lines are ignored rather than reported.
Ordinary invalid rows are excluded from execution but do not stop other valid
rows from running. Output-directory collisions are different: because they
could make two jobs overwrite one another, a collision stops the entire
manifest before any design starts. If no structurally valid jobs remain, the
manifest exits with an error without launching a design.
The per-row channel is part of the duplicate key: the same gene/accession can
legitimately occur once for B1 and once for B2. For example, `a/b` and `a_b`
both normalize to `a_b`; that collision stops the manifest before any design
is launched. A blank-channel row that expands to B1-B5 cannot overlap an
explicit row for the same target/channel, because that would make two rows
write the same result tree. Exact duplicate gene/accession/channel rows are
represented once in the executable job list. The first occurrence is run;
later occurrences are reported as duplicate bookkeeping rows and are not
executed again. Duplicates do not receive a separate `.skipped` marker because
they are not independent jobs; their status and reason remain in
`manifest_preflight.tsv`.

Each job is tracked independently. A successful job receives a `.done` marker;
an ordinary failed job receives a `.failed` marker and a log containing the
captured output. A true species-mismatch row receives a `.skipped` marker and a
status-table record, but no design or transcript-cache output. Explicit
per-row channel jobs have channel-specific marker and log names, so the same
gene/accession can be resumed independently for different channels. A
completed marker stores a signature of the HCRProbeForge version, species,
channels, smart-mode default, and forwarded design/QC/curation options. A later
invocation reuses a completed job only when that signature and the stored
species metadata both match. Changing an option causes that job to run again.
Completed markers without a parameter signature are rebuilt once so they can
acquire the current format. A temporary NCBI failure uses exit code 75. The
all-channel and manifest helpers retry that row or channel up to
`HCRPROBEFORGE_CHANNEL_ATTEMPTS` times (default 3); after the final attempt,
the manifest records a failed row and the helper returns a nonzero result.
Parsing, design, QC, and species-validation failures are not relabelled as
network failures.

The selected species is applied to both gene-symbol lookup and exact
transcript validation. An accession whose NCBI record belongs to another
organism is recorded as `skipped_species_mismatch`, produces a `.skipped`
marker, and does not make the manifest fail. A gene typo, unresolved accession,
or invalid transcript response is recorded as `failed` with a `.failed` marker
and makes the manifest return a nonzero status. In either case the row is
omitted from design and its new transcript/cache/output files are removed;
other manifest rows continue. The manifest command reports these outcomes
separately. Structural preflight-invalid rows also make the command return
nonzero; duplicate bookkeeping rows and true species-mismatch skips do not.
The
`manifest_preflight.tsv`, `run_status.tsv`, and outcome tables show the row and
the reason, while manifest logs retain the NCBI diagnostic for review.

Manifest status and result bookkeeping intentionally remain after a failure so
the user can identify and resume the failed row. Direct single-target, plot,
QC, and index failures are transactional and do not retain partial scientific
outputs or a misleading empty output directory.

## CLI quick start

The design examples below assume that the matching species index has already
been built, as described in [Species and genome indexes](#species-and-genome-indexes).
The default design performs genomic masking and stops if no matching index is
registered. Use `--no-genomemask` only when you deliberately want to skip that
specificity check. If `--species` is omitted, the CLI uses the built-in `xtr`
species preset; specify it explicitly for any other organism.

### Design one target

Gene-only input performs automatic transcript selection:

```bash
hcrprobeforge sox9 --species xtr --channel B3
```

This direct CLI command uses one-pass mode. Its native `<target>_IDT.tsv`
output is an upstream candidate table, not the final selected non-overlapping
order set. For a selected order-ready result, enable smart curation explicitly:

```bash
hcrprobeforge sox9 --species xtr --channel B3 \
  --auto-curate-if-needed --target-probes 20 \
  --outdir "$HOME/hcrprobeforge_projects/hcr_probe_design_project"
```

For a reproducible exact transcript:

```bash
hcrprobeforge sox9 --accession NM_001016853.2 \
  --species xtr --channel B3 --auto-curate-if-needed \
  --target-probes 20 \
  --outdir "$HOME/hcrprobeforge_projects/hcr_probe_design_project"
```

For a custom target sequence:

```bash
hcrprobeforge --fasta target.fa --species xtr \
  --channel B3 --outdir "$HOME/hcrprobeforge_projects/hcr_probe_design_project"
```

FASTA input must contain one record and cannot be combined with gene or
accession.

### Design a Pre-mRNA target

Pre-mRNA is an explicit target mode for any preset with a matching registered
genome and genomic annotation. It resolves one exact RefSeq transcript, then
uses that assembly's annotation and FASTA to reconstruct the full genomic
transcript in transcript 5′→3′ order: every exon and every intron is present
in the target sequence. On a minus-strand model, each genomic segment is
reverse-complemented and ordered in transcript orientation.

Use `--target-type pre-mrna` when probes must be intron-specific. It permits
only complete pairs inside selected introns, so exonic and intron-boundary
candidates are excluded before QC and curation. Use `pre-mrna-whole` when the
design should cover the entire unspliced transcript: exons, introns, and
boundary-spanning candidates are eligible. Both modes use the same genomic
target and annotation; the mode changes only candidate eligibility and
selection policy.

Leave `--premrna-introns` blank to use every annotated intron. When smart
curation is enabled, the Pre-mRNA selector distributes the requested final
pairs as evenly as sequence capacity allows across those introns. With smart
curation disabled, the one-pass design filters candidates to the selected
introns and enforces interval non-overlap before writing the order-ready table.
The one-pass path does not use `--target-probes` as a cap on the native
candidate request; use smart curation when you need a requested final count.
To target only particular introns, use a comma-separated 1-based list:

```bash
hcrprobeforge sox10 --species xtr \
  --target-type pre-mrna --premrna-introns 1,3,5 \
  --channel B3 --target-probes 20 --auto-curate-if-needed \
  --outdir "$HOME/hcrprobeforge_projects/hcr_probe_design_project"
```

For reproducibility, you may replace the gene-only input with the exact
versioned RefSeq accession resolved for the gene. Pre-mRNA design requires a
ready registered index and uses that registered assembly for both annotation
and genomic specificity screening; do not override it with `--index`. The
index workflow downloads and retains the uncompressed genome FASTA and GFF3
when an NCBI preset is built. If a ready index does not yet have the matching
reference assets, the first Pre-mRNA design downloads them once and caches the
indexed annotation and assembled target. The general smart-curation checkbox
remains optional in this mode: when enabled, curation distributes candidates
as evenly as capacity allows across selected introns; when disabled, the
ordinary one-pass design filters and packs eligible candidates after intron
filtering without targeting a final count. Packing means that candidates are
ordered by target coordinates and a candidate is retained only when it does
not overlap a pair already retained.

Intron numbering follows the transcript's 5′→3′ direction. Thus intron 1 is
the first intron encountered in transcript orientation, including for a
minus-strand genomic model. Pre-mRNA `start` and `end` coordinates are
1-based inclusive positions on the assembled full genomic-transcript target in
that same transcript orientation; mature-transcript coordinates refer to the
mature transcript sequence.

The output metadata records the assembly, selected transcript, intron numbers,
genomic coordinates, strand, target mode, cache paths, assembled target
length, and the distribution rule. For whole-target design:

```bash
hcrprobeforge sox10 --species xtr \
  --target-type pre-mrna-whole \
  --channel B3 --target-probes 20 --auto-curate-if-needed \
  --outdir "$HOME/hcrprobeforge_projects/hcr_probe_design_project"
```

The webapp exposes the two Pre-mRNA modes only for a selected species preset
and shows whether that species has a complete annotation database. If the
badge says the database will be built on the first Pre-mRNA run, the first run
may take time to scan the GFF3; later targets reuse the SQLite database and
the project cache. Mature-transcript mode remains the default and is
unchanged.

For a normal one-pass intronic run, the native HCRProbeDesign table is kept as
`<target>_candidate_IDT.tsv` under `details/` for diagnostics. It can contain
exonic or intron-boundary candidates from the upstream tool. The package
filters those candidates to complete selected introns, applies interval
non-overlap, writes
`<target>_eligible_IDT_order.csv`, and uses that eligible table for the normal
Pre-mRNA oligo QC. With `--auto-curate-if-needed`, every candidate reservoir is
filtered before curation and the final order file is
`<target>_final_IDT_order.csv`. In `pre-mrna-whole`, exonic, intronic, and
boundary-spanning candidates are eligible; the map and annotations identify
which regions each candidate overlaps.

### Design all linked transcripts

```bash
hcrprobeforge sox9 --species xtr \
  --all-transcripts --channel B3
```

Do not combine `--all-transcripts` with `--accession` or `--fasta`. The all-channel
helper forwards it to each requested channel:

```bash
hcrprobeforge-all-channels sox9 --species xtr --all-transcripts
```

### Design selected channels

```bash
hcrprobeforge-all-channels sox9 --species xtr
hcrprobeforge-all-channels sox9 --species xtr --channels B1,B3
hcrprobeforge-all-channels sox9 --accession NM_001016853.2 \
  --species xtr --channels B2,B4
```

The helper owns `--fasta`, `--outdir`, `--accession`, and `--channels`; pass
each of those to the helper itself, not as forwarded arguments. It also accepts
the design and curation options accepted by `hcrprobeforge`. Its default is all
five channels with smart auto-curation enabled.

### Design a manifest

```bash
hcrprobeforge-manifest genes.tsv "$HOME/hcrprobeforge_projects/hcr_probe_design_project" --species xtr
hcrprobeforge-manifest genes.csv "$HOME/hcrprobeforge_projects/hcr_probe_design_project" --species mouse --channels B1,B2,B3
hcrprobeforge-manifest genes.xlsx "$HOME/hcrprobeforge_projects/hcr_probe_design_project" --species human
hcrprobeforge-manifest multiplex.csv "$HOME/hcrprobeforge_projects/hcr_probe_design_project" --species xtr
```

Supported extensions are CSV, TSV, XLSX, XLSM, and plain text. Headered files
can use `gene_symbol`, `accession`, and optional `channel`. Headerless files may
contain one gene per line, gene plus accession, gene plus channel, or
gene/accession/channel in three-column order:

```text
sox9
sox9 NM_001016853.2
gnrhr2/nmi NM_001114076.1
foxj1 B1
pkd2 NM_001088456.1 B2
```

For a multiplex experiment, put the channel beside each target. This CSV
example avoids relying on visually invisible tab characters:

```csv
gene_symbol,accession,channel
foxj1,NM_<TRANSCRIPT_ID>.<VERSION>,B1
pkd2,NM_<TRANSCRIPT_ID>.<VERSION>,B2
arl13b,NM_<TRANSCRIPT_ID>.<VERSION>,B3
```

Replace the transcript placeholders with real RefSeq accessions. TSV files are
also supported; use actual tab separators rather than the two characters `\t`.

Run it once with the normal command. The optional per-row channel overrides
`--channels` for that row; rows without a channel use the global value. Exact
duplicate gene/accession/channel rows execute once: the first occurrence runs,
and later occurrences are recorded as duplicate bookkeeping rows without a
separate job marker.
Original biological labels remain in metadata and status tables; only path
components use the safe-name rule. Preflight detects invalid channel values
and path-normalization collisions before design starts. Species-
mismatched rows are recorded as skipped; mistyped genes and unresolved
accessions are recorded as failed. Neither starts a design, and valid rows
continue to run; the manifest exits nonzero when a row fails.

### Plot an existing table

```bash
hcrprobeforge --plot-only probes.tsv \
  --plot-theme minimal --plot-color-by dtm \
  --plot-dpi 300 --outdir "$HOME/hcrprobeforge_projects/hcr_probe_design_project"
```

Both PNG and SVG maps are written. The theme and color variable are separate:
minimal changes style and `--plot-color-by dtm` changes probe-block colors. Use
`--transcript-length` when the full length cannot be inferred.

### Run QC on an existing table

```bash
hcrprobeforge --qc-only final_IDT_order.csv \
  --qc-output "$HOME/hcrprobeforge_projects/hcr_probe_design_project/qc_review/final_oligo_structure_QC.xlsx"
```

Accepted inputs include IDT CSV/TSV/XLSX and selected-pairs tables. Without
`--qc-output`, CLI QC writes beside the input. Webapp QC writes in its QC run
directory.

## Complete CLI option reference

Run `hcrprobeforge --help` for the authoritative installed help.

### Input and transcript retrieval

| Option | Default | Meaning |
|---|---:|---|
| `gene` | none | Gene symbol or display label; unsafe characters stay in metadata. |
| `--organism` | preset | NCBI scientific name for a custom registered alias. |
| `--accession` | automatic | Exact RefSeq transcript accession. |
| `--fasta` | none | One FASTA record; cannot combine with gene/accession. |
| `--target-type` | mature | `mature` (default), `pre-mrna` (full target with intron-only eligibility), or `pre-mrna-whole` (full target with exon/intron eligibility); Pre-mRNA requires a matching registered genome/annotation. |
| `--premrna-introns` | all | Comma-separated 1-based intron numbers for `pre-mrna`; blank includes all introns in the full genomic target. Not used by `pre-mrna-whole`. |
| `--all-transcripts` | off | Design every linked RefSeq RNA transcript. |
| `--transcript-policy` | auto | auto, longest, interactive, or require-accession. |
| `--email` | `NCBI_EMAIL` | Email sent with NCBI requests. |
| `--api-key` | `NCBI_API_KEY` | Optional NCBI API key. |

### Reference and design

| Option | Default | Meaning |
|---|---:|---|
| `--species` | xtr | Registered species/index alias. |
| `--index` | registered alias | Explicit Bowtie2 prefix for mature-transcript designs. It cannot be used with `--target-type pre-mrna` or `pre-mrna-whole`, which must use the registered reference's matching genome and annotation. |
| `--channel` | B1 | B1, B2, B3, B4, or B5. |
| `--tile-size` | 52 | Target tile length in nucleotides. |
| `--min-gc` / `--max-gc` | 45 / 55 | Candidate GC bounds. |
| `--min-gibbs` / `--max-gibbs` | -70 / -50 | Gibbs free-energy bounds reported by HCRProbeDesign, in kcal/mol. |
| `--target-gibbs` | -60 | Preferred Gibbs free-energy value for candidate ranking, in kcal/mol. |
| `--max-run-mismatches` | 2 | Maximum tolerated mismatches in HCRProbeDesign's C/G-rich-run filter. The upstream run length is 7 nt and is not exposed as a HCRProbeForge option; this is not a Bowtie2 alignment-mismatch setting. |
| `--max-probes` | dynamic | Candidate-pair request sent to the initial `designProbes` run: 30 in one-pass mode, or the smart-mode initial request described in [How the design pipeline works](#how-the-design-pipeline-works). It is not the final selected count. |
| `--num-hits-allowed` | 1 | Maximum number of Bowtie2 genomic alignments accepted by HCRProbeDesign for each target tile. The intended genomic alignment consumes one reported hit, so `1` allows a unique on-target tile and rejects tiles with additional genomic matches. A tile that has zero accepted end-to-end genomic alignments, including a tile that spans an exon–exon junction with no contiguous genomic match, is not rejected by this maximum-hit setting alone; this option is not a positive on-target or junction-validity test. Screening uses the selected registered genome index; HCRProbeForge does not add a separate junction, paralog, or *X. laevis* homeolog screen. |
| `--dtm-filter` | off | Enable filtering by probe-arm dTm, the absolute difference between the two arm melting temperatures, in °C. |
| `--dtm-max` | 5 | Maximum dTm in °C when `--dtm-filter` is enabled. dTm is the absolute difference between the two arm Tm values. The auto-score separately adds a fixed penalty for dTm above 5 °C, even when a different `--dtm-max` is supplied. |
| `--no-genomemask` | off | Disable genomic specificity screening. |

### Plotting

| Option | Default | Meaning |
|---|---:|---|
| `--plot-only` | none | Plot an existing probe table and exit. |
| `--transcript-length` | inferred | Full sequence length for plotting. Supply it for a standalone table when it cannot be inferred from the table or sibling metadata. |
| `--plot-theme` | pastel | pastel or minimal. |
| `--plot-color-by` | gc | Color by GC percentage, probe-pair dTm in °C, or displayed table/order position (`order`). The order mode cycles through the map palette by row number; it is not a biological score. |
| `--plot-dpi` | 300 | PNG resolution; SVG is always generated. |
| `--no-probe-labels` | off | Hide probe numbers. |
| `--plot-title` | automatic | Replacement map title. |

### QC, curation, and output

| Option | Default | Meaning |
|---|---:|---|
| `--no-oligo-qc` | off | Skip final oligo QC. |
| `--qc-only` | none | QC an existing IDT/selected-pairs table. |
| `--qc-output` | automatic | XLSX/XLSM output for QC-only. |
| `--auto-curate-if-needed` | off | Run fallback tiers when needed. |
| `--target-probes` | 20 | Desired final pair count in smart mode. In one-pass mature and Pre-mRNA-intronic modes it does not cap the native candidate request or candidate table; one-pass Pre-mRNA still filters and packs eligible intronic candidates without overlap. Use `--max-probes` to control the one-pass candidate request. |
| `--min-acceptable-probes` | 12 | Minimum pair count used by smart-curation/reporting to label an undersized result as a useful `partial_success`. It does not change candidate generation, the `--target-probes` cap, or interval selection; the actual result status and counts are recorded in the report. |
| `--qc-stringency` | balanced | strict, balanced, or permissive. |
| `--auto-curate-plan` | standard | standard, conservative, or deep. |
| `--max-auto-runs` | 12 | Additional distinct runs per phase. |
| `--auto-curate-max-probes` | dynamic | Candidate request for added smart-curation tiers; the default is described in [How the design pipeline works](#how-the-design-pipeline-works). |
| `--coverage-policy` | balanced | balanced or qc-only. |
| `--outdir` | hcr_results | Project output root. |

The helper commands add `--channels` as a comma-separated subset of B1-B5; its
default is `B1,B2,B3,B4,B5`. The index command has separate options described
above.

## Option interactions and reproducible command patterns

The options are intentionally grouped by responsibility. A command can contain
options from several groups, but an option is acted on only by the workflow that
uses it. This makes it possible to keep one documented command line for a
project without accidentally changing a plot or QC run with a design-only
setting.

The following interactions are the ones most likely to matter in practice:

| Situation | Rule |
|---|---|
| `--species` without `--organism` | A built-in alias supplies the NCBI organism, assembly metadata, and registered HCRProbeDesign index alias. |
| Registered custom alias | Supply the matching `--organism`; the package will not guess the taxon from an arbitrary alias. |
| `--index` | Overrides the registered Bowtie2 prefix for a mature-transcript design. It does not change the NCBI organism used for transcript retrieval and is rejected for Pre-mRNA modes. |
| `--fasta` | Must be a single FASTA record and cannot be combined with a gene symbol or `--accession`. The selected species/index still controls genomic masking unless it is disabled. |
| `--all-transcripts` | Applies to gene-only input. It cannot be combined with `--accession`, because an explicit accession already identifies one transcript. |
| `--transcript-policy` | Applies only when a gene symbol needs automatic transcript selection. It has no effect on an explicit accession or custom FASTA. |
| `--target-type pre-mrna` or `pre-mrna-whole` | Requires one NCBI-selected transcript and a complete registered index with matching genomic assets; it cannot be combined with `--fasta`, `--all-transcripts`, or `--index`. `pre-mrna` restricts pairs to complete selected introns; `pre-mrna-whole` allows the full genomic target. Mature-pipeline QC and scientific parameters are unchanged. |
| `--premrna-introns` | Restricts Pre-mRNA to comma-separated 1-based intron numbers. Blank uses every annotated intron. The selector uses equal capacity-aware intron quotas. |
| `--auto-curate-if-needed` | Enables the capacity-aware curation search and requires final oligo QC; do not combine it with `--no-oligo-qc`. |
| `--target-probes` and `--max-probes` | `--target-probes` is the desired final pair count. `--max-probes` is the candidate request sent to the initial `designProbes` call. They are not interchangeable. |
| `--plot-only` and `--qc-only` | These are standalone modes. They read an existing table and do not run transcript retrieval or probe design. |
| `--plot-*` options | They affect plot-only output and the maps produced by a design run; they do not change candidate generation or selection. |
| `--qc-output` | It is used by QC-only mode. Automatic design QC remains beside the design output and uses the standard final-QC filename. |
| `--no-genomemask` | Deliberately skips genome specificity screening. It is useful for a custom sequence with no meaningful reference, but it removes an important design check. |
| `--outdir` | Is the project root. The package creates `runs/` and `cache/` below it and keeps species/workflow paths separate. |

For a reproducible design, pin the species, versioned transcript accession,
channel, and the design/curation options that differ from the defaults:

```bash
hcrprobeforge sox9 --accession NM_001016853.2 \
  --species xtr --channel B3 --auto-curate-if-needed \
  --target-probes 20 \
  --outdir "$HOME/hcrprobeforge_projects/hcr_probe_design_project"
```

For a one-pass design that must call `designProbes` only once, omit smart
curation and set the candidate request explicitly:

```bash
hcrprobeforge sox9 --accession NM_001016853.2 \
  --species xtr --channel B3 --max-probes 30 \
  --outdir "$HOME/hcrprobeforge_projects/hcr_probe_design_project"
```

For plotting and QC, the input table is the scientific input and the output
folder is presentation/bookkeeping output:

```bash
hcrprobeforge --plot-only selected_pairs.tsv --plot-color-by dtm \
  --plot-theme pastel --outdir "$HOME/hcrprobeforge_projects/hcr_probe_design_project"
hcrprobeforge --qc-only final_IDT_order.csv \
  --qc-output "$HOME/hcrprobeforge_projects/hcr_probe_design_project/qc_review.xlsx"
```

The installed command is authoritative if a future HCRProbeDesign release adds
an option. Use `hcrprobeforge --help`, `hcrprobeforge index fetch --help`, and
`hcrprobeforge-web --help` when copying commands between environments.

## How transcript selection works

With a gene and the default `--transcript-policy auto`, HCRProbeForge uses the
requested species and deterministic ranking:

1. exact active NCBI Gene match;
2. a record with linked RefSeq RNA;
3. RefSeq Select when available;
4. curated NM_/NR_ before predicted XM_/XR_;
5. complete before partial records;
6. longest remaining transcript;
7. accession as the final tie-break.

The selected accession, alternatives, GeneID, ranking fields, and reason are
saved in a file named `<target>_transcript_resolution.json` (where `<target>`
is replaced by the safe target name). `longest` ignores accession class,
interactive prompts in a terminal, and require-accession stops when a manual
accession is needed. The webapp supports the non-interactive policies.

Exact RefSeq sequence and feature data are validated against the accession
version. Annotation is reporting-only: it may add 5'-UTR, CDS, and 3'-UTR
tracks and region counts without changing candidate generation or selection.

Custom FASTA bypasses NCBI transcript lookup and annotation. It must contain
A/C/G/T/N; U is converted to T. Coordinates are 1-based within the FASTA
record.

## How the design pipeline works

1. HCRProbeDesign enumerates and filters candidate tiles and provides
   channel-specific HCR oligo sequences.
2. HCRProbeForge retains the candidate reservoir and calculates/report metrics
   used by QC and, when enabled, the selection pipeline.
3. Primer3 checks the final initiator-bearing P1/P2 oligos for hairpins,
   self-dimers, same-pair heterodimers, and pool-level heterodimers.
4. When smart auto-curation is enabled, it can request additional predefined
   tiers when the initial design is insufficient, dimer reviews remain, or
   balanced coverage is poor for a sufficiently long transcript.
5. Smart selection maximizes acceptable non-overlapping pair count first, then
   respects explicit QC priorities, coverage, and technical tie-breakers.
6. One-pass output, or smart-curated final output, and its metadata are written
   after the required QC processing succeeds. One-pass intronic Pre-mRNA output
   is filtered to complete selected introns and interval-packed so returned
   order-ready pairs do not overlap.

HCRProbeDesign reports `GibbsFE` as predicted target-binding free energy in
kcal/mol, using its documented 37 °C and 0.33 M salt-correction model. More
negative values indicate stronger predicted binding. HCRProbeDesign reports
`dTm` as the absolute melting-temperature difference between the two split
target-binding arms, in °C.

Smart mode normally requests four times the value of `--target-probes` initially
and five times that value from added tiers. These are candidate requests, not a
promise of final count.
Short transcripts may correctly yield fewer pairs because of geometric
capacity, sequence filtering, specificity, or QC.

The `_final_IDT_order.csv` file contains the actual selected initiator-bearing
oligos produced by smart curation. In a normal one-pass mature-transcript run,
the native `<target>_IDT.tsv` is the candidate order table and no non-overlap
selection has been applied. In a normal intronic Pre-mRNA run,
`<target>_eligible_IDT_order.csv` is the filtered, non-overlapping candidate
order table; the
native `<target>_candidate_IDT.tsv` is retained for diagnostics. Do not
reconstruct order sequences from target tiles.

## HCR channels and final-order conventions

Each selected probe pair is made from two orderable oligos. The 52-nt default
is the target-tile length supplied to HCRProbeDesign; it is not the length of
either final order oligo. HCRProbeDesign splits each target tile and adds the
channel-specific initiator and spacer structure before HCRProbeForge performs
final-oligo QC. In the standard 52-nt model, the target tile is split into two
25-nt target-binding arms with a 2-nt central gap (`25 + 2 + 25 = 52`); the
initiator and spacer sequences are then added to those arms. For the even tile
sizes used by the published fallback plans, each arm has
`(tile size − 2) / 2` nucleotides: 50 nt gives 24+2+24, 48 gives 23+2+23,
46 gives 22+2+22, 44 gives 21+2+21, 42 gives 20+2+20, 40 gives 19+2+19,
38 gives 18+2+18, and 36 gives 17+2+17. These shorter tiers are rescue
options, not equivalent substitutes for the standard 25-nt arms; shorter arms
can reduce target specificity and may have different Tm behavior. The final
QC and genomic screening rules still apply, but users should retain the tile
size, arm lengths, and HCRProbeDesign version with the design. The final order
file contains two rows per selected pair, normally named with `_P1` and `_P2`
suffixes.

The channel is part of the scientific design input. It must match the HCR
amplifier used in the experiment:

| Probe channel | P1/odd 20-nt channel sequence component | P2/even 20-nt channel sequence component | Matching amplifier |
|---|---|---|---|
| B1 | `GAGGAGGGCAGCAAACGGAA` | `TAGAAGAGTCTTCCTTTACG` | B1 H1/H2 |
| B2 | `CCTCGTAAATCCTCATCAAA` | `AAATCATCCAGTAAACCGCC` | B2 H1/H2 |
| B3 | `GTCCCTGCCTCTATATCTTT` | `TTCCACTCAACTTTAACCCG` | B3 H1/H2 |
| B4 | `CCTCAACCTACCTCCAACAA` | `ATTCTCACCATATTCGCTTC` | B4 H1/H2 |
| B5 | `CTCACTCCCAATCTCTATAA` | `AACTACCCTACAAATCCAAT` | B5 H1/H2 |

Each sequence in this table is a 20-nt channel-specific initiator-side
component that includes the 2-nt spacer used by the HCRProbeDesign channel
convention. The strings are shown 5′→3′ for their P1 or P2 attachment; they are
not complete order oligos, and they do not replace the P1/P2 sequences written
by HCRProbeDesign in the final IDT file.

The exact sequence and name columns in a final IDT CSV are the orderable
record. Candidate tables are useful for review, but a normal one-pass native
candidate table is not necessarily a non-overlapping final set. A design of 20
selected pairs normally produces 40 final oligos. Inspect the final QC
workbook and selected-pairs table before ordering, especially when a run used
permissive fallback tiers, a custom FASTA, or a Pre-mRNA target.

## Auto-curation technical details

Smart auto-curation is a capacity-aware search around the initial design
settings. It does not replace the candidate generator or silently change the
sequence calculations. It asks HCRProbeDesign for additional candidate
reservoirs, runs the same final-oligo QC on those candidates, and lets the
existing selector compare all acceptable candidates together.

### Initial and fallback phases

Phase 0 is always the exact initial request. No relaxation is applied to the
user's initial tile size, GC range, Gibbs range, C/G-run setting, or QC
stringency. Smart-mode candidate reservoir sizes use the initial and
per-tier defaults documented in the option reference unless the corresponding
candidate-count options are supplied. These are reservoir sizes, not promises
about the number of final pairs.

The normal phase follows the selected `--auto-curate-plan` and
`--qc-stringency`. A deep adaptive phase is started only when an applicable
condition indicates that more searching may help. The normal and deep phases
have separate `--max-auto-runs` limits. Duplicate parameter combinations are
skipped and do not consume a run slot.

The fallback decision considers three things:

1. The selected count is below the smaller of the requested target and the
   transcript-specific deep-plan geometric capacity.
2. The selected set still contains self-dimer or same-pair P1/P2 dimer review
   flags.
3. Balanced coverage is enabled and the selected set has an applicable span,
   gap, or occupied-bin review condition.

Low-yield or dimer rescue may use permissive candidate acceptance, but it still
rejects practical-reject structures and strong self- or same-pair dimers.
Coverage-only rescue keeps the original QC stringency; improved distribution
does not justify admitting a poorer-quality oligo. All runs, effective
parameters, triggers, and stopping reasons are recorded in the auto-curation
JSON and Markdown report.

### Geometric capacity and short-transcript protection

The geometric upper bound is calculated independently for every contiguous
unambiguous A/C/G/T segment:

```text
capacity(tile size) = sum floor(unambiguous segment length / tile size)
```

The report records capacity for the canonical 52-nt tile, the shortest tile in
the selected normal plan, and the 36-nt deep-plan minimum. These are only upper
bounds: GC, Gibbs, C/G-run filters, genomic specificity, candidate export,
Primer3 QC, and interval packing can reduce the achievable count.

Balanced coverage is automatically disabled when the canonical 52-nt capacity
is below ten pairs. Even when the policy is enabled, coverage review and
coverage-triggered fallback require a selected set of at least ten pairs.
Short transcripts are still allowed to use shorter tiles when that increases
the acceptable non-overlapping count, but they are not penalized for failing
long-transcript coverage criteria.

### Standard and deep plans

The standard plan uses these tiers, in order. The user's initial CLI/webapp
tile, GC, Gibbs, and C/G-run values are used exactly for the initial run. The
plan then considers each distinct standard tier, including the first row when
the user's initial values differ from it; an identical tier is not run twice.
Subsequent tiers are only considered when smart curation needs them.

| Tier | Tile (nt) | GC (%) | Max C/G-run mismatches | Gibbs range | Target Gibbs |
|---:|---:|---:|---:|---:|---:|
| 1 | 52 | 45–55 | 2 | −70 to −50 | −60 |
| 2 | 52 | 40–60 | 2 | −70 to −50 | −60 |
| 3 | 52 | 35–65 | 2 | −70 to −50 | −60 |
| 4 | 52 | 30–65 | 1 | −70 to −50 | −60 |
| 5 | 50 | 40–60 | 2 | −70 to −48 | −58 |
| 6 | 48 | 40–60 | 2 | −70 to −47 | −56 |
| 7 | 46 | 40–65 | 1 | −70 to −45 | −55 |
| 8 | 44 | 40–65 | 1 | −70 to −43 | −52 |
| 9 | 42 | 40–65 | 1 | −70 to −40 | −50 |

The conservative plan uses five tiers: 52 nt at 45–55% GC, 52 nt at 40–60%,
50 nt at 40–60% with a −70 to −48 Gibbs range, 48 nt at 40–60% with a −70 to
−47 range, and 46 nt at 40–65% with one allowed C/G-run mismatch. The deep
plan is explicit below; it broadens the 50–42 nt tiers to 35–65% GC and one
C/G-run mismatch, then adds shorter tiers. The exact plan used is written to
the report; the selector may retain a mix of tile sizes when that produces a
better valid set.

| Deep tier | Tile (nt) | GC (%) | Max C/G-run mismatches | Gibbs range | Target Gibbs |
|---:|---:|---:|---:|---:|---:|
| 1–4 | 52 | 45–55, 40–60, 35–65, 30–65 | 2, 2, 2, 1 | −70 to −50 | −60 |
| 5 | 50 | 35–65 | 1 | −70 to −48 | −58 |
| 6 | 48 | 35–65 | 1 | −70 to −47 | −56 |
| 7 | 46 | 35–65 | 1 | −70 to −45 | −55 |
| 8 | 44 | 35–65 | 1 | −70 to −43 | −52 |
| 9 | 42 | 35–65 | 1 | −70 to −40 | −50 |
| 10 | 40 | 35–65 | 1 | −70 to −38 | −46 |
| 11 | 38 | 35–65 | 0 | −70 to −36 | −44 |
| 12 | 36 | 35–65 | 0 | −70 to −34 | −42 |

### Balanced transcript coverage

Coverage is a tie-breaker after count and explicit QC. It never replaces a
cleaner set with a more heavily flagged set. For a selected set the pipeline
calculates:

- the fraction of the transcript spanned from the first probe start to the
  last probe end;
- the largest untargeted gap, including terminal gaps;
- ideal spacing (`transcript length / selected pair count`);
- the largest-gap-to-ideal-spacing ratio;
- the number of occupied bins out of ten transcript bins; and
- mean deviation from evenly spaced target positions.

For sets with at least ten probe pairs, balanced fallback can be triggered when
any of these conditions applies:

```text
probe span fraction < 0.85
maximum untargeted gap > 2.5 × ideal spacing
occupied transcript bins < 8 of 10
```

These are search triggers and report conditions, not hard biological failures.
If no equally clean alternative exists, the cleaner set is retained and the
remaining coverage issue is reported.

## Candidate ranking and auto_score

The final set is not chosen by taking the lowest `auto_score` values. The
selection objective is lexicographic and prioritizes practical quality:

1. maximize the number of non-overlapping selected pairs, up to
   `--target-probes`;
2. minimize selected self-dimer and same-pair P1/P2 dimer review burden;
3. minimize dimer deltaG severity;
4. minimize strong or energetically concerning hairpins;
5. minimize mild Tm-only hairpin reviews;
6. when balanced coverage is enabled, improve even transcript distribution
   among otherwise equivalent sets; and
7. use source tier, tile length, GC, Gibbs, dTm, and `auto_score` as technical
   tie-breakers.

`auto_score` is therefore a review aid and late tie-breaker. It is not a
thermodynamic measurement, Primer3 output, biological confidence score, or
pass/fail threshold. Lower is preferred only when the higher-priority
objectives are equivalent.

```text
auto_score =
    source_tier * 100
  + max(0, 52 - tile_size) * 3
  + abs(GC - 50) * 2
  + abs(GibbsFE - source_target_gibbs) * 1.5
  + max(0, dTm - 5) * 8
  + QC_penalty
```

The report's internal `source_tier` field is zero-based: it is `0` for the
initial run, which is user-facing Tier 1, `1` for the first relaxed tier (Tier
2), `2` for Tier 3, and so on. The source-tier penalty is therefore 0 for the
initial run, 100 for the first relaxed tier, 200 for the second, and so on. The
tile-size penalty is zero at 52 nt
and increases by three for every nucleotide below 52. This is a deliberately
simple technical tie-breaker that keeps the standard 52-nt design preferred
when higher-priority QC and coverage objectives are equivalent; it is not a
measured biological penalty or a claim that every nucleotide reduction has the
same experimental effect. GC is favored near 50%.
Gibbs is favored near the target for the source tier, and dTm above 5 °C adds
an increasing penalty. QC penalty components reflect the observed review
burden; they do not replace the explicit QC gate.

## Final oligo QC thresholds

Primer3 is run on the actual initiator-bearing P1 and P2 order sequences after
HCRProbeDesign constructs the channel-specific oligos. This is important:
the final order oligo can fold or dimerize differently from its transcript
target tile. The same QC logic is used in automatic design and the standalone
QC workflow.

HCRProbeForge calls `primer3-py` hairpin, homodimer, and heterodimer functions
without overriding their Primer3 thermodynamic defaults: 50 mM monovalent
salt, 1.5 mM divalent salt, 0.6 mM dNTP, 50 nM oligonucleotide concentration,
37 °C, and a maximum loop size of 30 bases. The QC report converts Primer3
deltaG values from cal/mol to kcal/mol and reports melting temperatures in °C.
Reproducibility records the installed `primer3-py` version; these are the
thermodynamic-analysis defaults, not a custom salt or oligo-concentration model
validated specifically for every HCR hybridization buffer. The package's QC
cutoffs are conservative review heuristics, not experimental guarantees; keep
the workbook and validate the final oligos under the intended assay conditions.

The following classifications are reported for each structure type:

| Structure | Review | Strong review | Practical reject |
|---|---|---|---|
| Hairpin | Tm ≥ 50 °C or ΔG ≤ −5 kcal/mol | ΔG ≤ −8 kcal/mol, or Tm ≥ 60 °C and ΔG ≤ −7 kcal/mol | not used as a standalone hard reject |
| Self-dimer | ΔG ≤ −9 kcal/mol | ΔG ≤ −12 kcal/mol | ΔG ≤ −15 kcal/mol |
| Same-pair P1/P2 heterodimer | ΔG ≤ −10 kcal/mol | ΔG ≤ −12 kcal/mol, or Tm ≥ 37 °C and ΔG ≤ −11 kcal/mol | not used as a standalone hard reject |
| Cross-pool heterodimer | ΔG ≤ −10.5 kcal/mol | ΔG ≤ −13 kcal/mol | ΔG ≤ −15 kcal/mol |

All stringencies exclude practical-reject structures and strong self- or
same-pair dimers. Balanced and strict modes also exclude strong hairpins.
Strict mode additionally excludes ordinary self- and same-pair review flags.
Permissive mode can retain strong-hairpin candidates during yield rescue, but
does not admit the hard dimer exclusions listed above. Cross-pool dimer reviews
are evaluated on the selected pool and can trigger replacement repair; they are
not silently treated as single-candidate hairpin or same-pair gates. A
cross-pool practical-reject threshold (ΔG ≤ −15 kcal/mol) is therefore a repair
trigger/review note rather than an unconditional candidate deletion; if no
count-preserving replacement is found, the remaining issue is retained in the
report. Mild Tm-only hairpin reviews are recorded and ranked below dimer
concerns.

The final QC workbook contains the structure calls and numerical values. A
QC-only run does not design probes; it reads the provided oligo table, checks
the orderable sequences, and writes a workbook to the requested QC folder.

## Annotation and coordinate conventions

Transcript coordinates are 1-based and inclusive. In mature-transcript mode,
probe `start` and `end` coordinates refer to the mature transcript sequence
used for design. In Pre-mRNA mode, they refer to the assembled full genomic
transcript in transcript orientation, as described above. The FASTA sequence
itself is normalized to uppercase DNA; U is converted to T and only A/C/G/T/N
are accepted.

NCBI feature annotation is reporting-only. When an exact RefSeq record is
available, the package can label 5′ UTR, CDS, 3′ UTR, partial/unannotated
terminal sequence, and unclassified non-CDS sequence on the map and in
selected-pair tables. It does not alter candidate generation, QC, non-overlap
selection, or order sequences. If annotation retrieval fails, design continues
with an unsegmented transcript map and the report records the reason.

Pre-mRNA is the separate opt-in exception to that reporting-only rule: it
reconstructs the full genomic transcript from the selected model, keeps
transcript-strand orientation (including reverse-complementing minus-strand
transcripts), and maps exons and introns onto the target. In `pre-mrna`
intronic mode, candidates must be completely contained in selected introns. In
`pre-mrna-whole` mode, candidates may be in exons, introns, or span a boundary.
The assembled target and transcript model are cached under the selected
project's `cache/<species>/<workflow>/premrna/` tree, while the
mature-transcript default remains reporting-only as described above.

Maps and metadata use the conventional `5′` and `3′` biological labels. JSON
files are UTF-8 encoded, so programs that parse them receive the same labels.
Map feature labels are placed inside their own annotation region whenever the
complete label fits. Mature-transcript labels may use an external callout when
their individual region is too narrow; Pre-mRNA exon/intron labels that do not
fit are omitted to keep dense maps readable. Probe-number labels use one
inside/outside mode for each group of boxes whose rendered widths are within
roughly 15% of one another, rather than mixing label placement within that
group. Clearly different tile-length groups may use different modes when
their available space differs.

Species matching is also part of input validation. For an NCBI-backed exact
accession, the title/organism metadata must agree with the selected species
before the sequence is written to the cache or passed to design. This prevents
a transcript from another organism being designed against the selected
genome. A FASTA supplied by the user is not assigned an NCBI organism; the
user is responsible for providing a sequence from the selected organism and
index.

## Output layout

The project root is the directory supplied to `--outdir` or entered in the
webapp's Results/project folder field. For example, if that value is
`$HOME/hcrprobeforge_projects/hcr_probe_design_project`, all run, cache, plot,
QC, and manifest outputs are nested below it. The webapp also keeps a small
`index_metadata/` reference-metadata copy under its parent directory, at
`<project parent>/index_metadata/`; that registry copy is described separately
below.

```text
hcr_probe_design_project/
├── cache/
│   └── <species>/
│       ├── individual_designs/                # transcript and target cache for single-target designs
│       ├── manifests/<manifest-name>/         # cache for one manifest run
│       ├── plot/                              # created only if plotting writes cache data
│       └── qc/                                # created only if QC writes cache data
└── runs/
    └── <species>/
        ├── individual_designs/
        │   ├── <target>/                     # one-target run
        │   └── <gene>/<gene>_B1/<target>/     # all-channel run
        ├── manifests/<manifest-name>/
        │   ├── results/<gene>/<gene>_B1/<target>/
        │   ├── status/
        │   │   ├── tables/*.tsv
        │   │   └── *.done, *.failed, or *.skipped
        │   └── logs/*.log
        ├── plot/<plot-name>/
        └── qc/<qc-name>/
```

`<species>` is the readable display-name component, such as
`Xenopus_tropicalis`; it is not necessarily the short HCRProbeDesign alias.
Important result files stay at the target root. Native design tables,
transcript-resolution data, run parameters, and other diagnostic files are
grouped under `details/`. Manifest status tables are grouped under
`status/tables/`, and manifest logs remain under `logs/`.

For a normal one-pass design, the target root contains the map, summary, and
QC workbook when applicable, while the native probe and IDT tables are under
`details/`. Smart-curated runs additionally place the selected-pairs table,
final IDT order file, final map, final QC workbook, and curation report at the
target root; its individual tier runs are under `details/design_runs/`. An
intronic Pre-mRNA one-pass run also places its filtered
`<target>_eligible_IDT_order.csv` at the target root; its unfiltered upstream
candidate order table stays under `details/`.

If a design, plot, QC, or index operation fails, newly created scientific
outputs and partial cache/index metadata are removed so an incomplete run is
not mistaken for a usable result. Files from an existing run that were
overwritten during the failed operation are restored. Manifest status tables,
failure markers, and logs are retained because they explain which row failed
and what to correct before rerunning.

## Understanding the output files

- `*_final_selected_pairs.tsv`: selected non-overlapping probe pairs from smart
  curation, with coordinates, metrics, QC, curation, and annotation fields.
- `*_final_IDT_order.csv`: order-oriented P1/P2 table with the selected sequences
  from smart curation.
- `*_final_probe_map.png` and `*_final_probe_map.svg`: transcript maps.
- `*_final_oligo_structure_QC.xlsx`: Primer3 structure-QC workbook.
- `*_auto_curate_report.json` and `.md`: capacity, tiers, fallback, coverage, QC,
  and selected-coordinate decisions.
- `*_candidate_master.tsv`: merged candidate reservoir.
- `*_backup_candidates.tsv` and `*_rejected_candidates.tsv`: alternatives and
  rejected candidates.
- `*_transcript_resolution.json`: accession selection and cache provenance.
- `*_probes.tsv` and `*_IDT.tsv`: native one-pass candidate tables under `details/`.
- `*_candidate_IDT.tsv`: native intronic Pre-mRNA candidates under `details/`;
  these are retained for diagnostics and are not the filtered order table.
- `*_eligible_IDT_order.csv`: intronic Pre-mRNA candidates that passed complete-
  intron eligibility, used as the one-pass QC input.
- `run_parameters.json`, `summary.json`, and logs: run provenance under `details/`.

Manifest status/tables include manifest_preflight.tsv, run_status.tsv,
channel_design_outcomes.tsv, and transcript_design_summary.tsv.

The webapp report shows final maps for design/plot workflows, selected-pair
files for design/manifest workflows, the full output directory, and a collapsible
console. It does not mix files from previous runs.

## Caching and reproducibility

Transcript FASTA, compact NCBI records, and feature annotation are cached by
exact accession-version. Direct CLI, helper, manifest, and webapp workflows
normally use one project cache tree under `<project>/cache/<species>`, with
the same readable workflow names used below `runs/`. A manifest cache is
grouped by manifest name. Plot and QC cache folders are not created unless
those workflows actually write cache data.

Override the cache location explicitly:

```bash
export HCRPROBEFORGE_CACHE_DIR=/path/to/shared/hcrprobeforge-cache
```

For direct CLI runs this environment variable replaces the project cache
location. The all-channel helper and webapp deliberately set their own
project-scoped cache path so channels and workflows cannot accidentally share
unrelated results; use the project directory or the workflow's explicit cache
setting for those front ends.

Windows PowerShell:

```powershell
$env:HCRPROBEFORGE_CACHE_DIR = "D:\data\hcrprobeforge-cache"
```

The supported environment variables are:

| Variable | Effect |
|---|---|
| `HCRPROBEFORGE_CACHE_DIR` | Override the direct CLI cache root. The all-channel helper and webapp normally set a project-scoped cache themselves. |
| `HCRPROBEFORGE_REFERENCE_DIR` | Relocate the HCRProbeForge reference registry, custom-species presets, and managed reference assets. |
| `HCRPROBEDESIGN_DATA_DIR` | Relocate HCRProbeDesign's registered Bowtie2 indexes. |
| `HCRPROBEFORGE_CHANNEL_ATTEMPTS` | Maximum retry attempts per helper/manifest channel for retryable NCBI failures; default `3`. |
| `NCBI_EMAIL` | Contact email used for NCBI requests when `--email` is omitted. |
| `NCBI_API_KEY` | Optional NCBI API key used when `--api-key` is omitted. |

If neither `--email` nor `NCBI_EMAIL` is supplied, HCRProbeForge omits the email
parameter; it does not invent or persist an address. Supplying an email is
recommended for responsible NCBI use, especially for larger jobs. An API key
is optional and only changes the NCBI request-rate allowance.

Gene-only input is intentionally time-dependent: it resolves against the
current NCBI Gene/RefSeq Select records. For reproducibility, pin an exact
versioned RefSeq accession and record the assembly accession and index alias.
Keep selected-pair and IDT files, `reference.json` or the index prefix,
transcript-resolution JSON, run parameters, curation report, exact input,
species, assembly, accession, and index alias with published designs. The run
parameters also record HCRProbeForge, HCRProbeDesign, Primer3, Bowtie2, and
relevant Python dependency versions when available.

## Troubleshooting

### designProbes or buildGenomeIndex was not found

Activate the environment containing HCRProbeDesign and check PATH with
command -v or where.exe. The Python package cannot execute a command invisible
to the environment that launched it.

### bowtie2 or bowtie2-build was not found

Install Bowtie2 with conda/mamba or add its directory to PATH, then recheck
from the same terminal.

### NCBI reports no Assembly record

Check that the accession is an NCBI Assembly accession and that its spelling
and version are correct. Provide the expected assembly name as an additional
safety check. An accession valid in another NCBI database may not be an
Assembly accession.

### NCBI requests fail or are slow

Set `NCBI_EMAIL` and optionally `NCBI_API_KEY`. Avoid many simultaneous runs.
Transient E-utilities and genome-download responses are retried; an
exact-accession cache or local FASTA avoids repeated downloads. Download
errors include the requested accession and the URL candidates that were
checked.

### An NCBI Assembly genome file is not found

The Assembly accession may be valid even when a guessed filename is not. The
package first uses the authoritative FTP directory from NCBI and derives the
filename from its normalized directory basename. It then tries a small set of
record-derived alternatives and checks the directory listing for a
`*_genomic.fna.gz` file. If all candidates fail, verify the accession in the
NCBI Assembly record and either retry later or download the genomic FASTA
yourself and use:

```bash
hcrprobeforge index build --species <alias> --fasta /path/to/genome.fna
```

For a local build, `--annotation` is optional for mature designs and should be
the GFF3 from the same assembly when Pre-mRNA will be used. The package copies
the uncompressed FASTA and annotation into its managed reference directory.
Use `--with-annotation-db` to build the SQLite lookup during index creation;
otherwise the first Pre-mRNA run builds it on demand. NCBI remains the
transcript source for NCBI-backed designs.

### The setup page shows an old index or Pre-mRNA status

The index and Pre-mRNA badges are derived from the registered reference
metadata, not from the transcript cache. A successful index build is recorded
only after its Bowtie2 files and any requested `annotation.sqlite` are
complete. The webapp refreshes the selected option automatically and all
dynamic responses use `no-store` cache headers. If a webapp process is still
running after an upgrade, stop it and launch the installed version again;
installing a new wheel does not replace code already loaded by a running
Python process.

If the badge remains pending after restarting, inspect the selected reference
directory and confirm that `reference.json` has `status: "ready"`,
`annotation_database_status: "ready"`, and an existing `annotation.sqlite`.
The design path validates the database schema before extracting a genomic
transcript and reports a specific rebuild or missing-annotation error when the
database is incomplete. Mature-transcript mode does not require this database.

### No registered index exists

For a design that uses genomic masking, HCRProbeForge reports the selected
species and tells you that its Bowtie2 index is missing. Build the matching
preset with `hcrprobeforge index fetch --species <alias>` or use the Build a
genome index tab, then rerun the design. You can also provide a valid explicit
`--index` prefix. A failed design is rolled back, so there is no failed-run
log to inspect; the error message is the diagnostic. Use `--no-genomemask`
only when you deliberately accept that genomic specificity screening is not
being performed.

### A gene symbol is not found for a custom species

Some non-standard organisms have incomplete or inconsistent NCBI Gene and
RefSeq links even when the sequence exists. If symbol lookup or linked RefSeq
RNA lookup fails, retry with the exact NCBI/RefSeq transcript or gene
accession. If NCBI does not provide a usable transcript record, use the FASTA
workflow and provide the transcript sequence directly. This fallback changes
only input resolution; it does not change probe generation, genomic masking,
oligo QC, or selection behavior.

### A manifest stops before design

Open `status/tables/manifest_preflight.tsv`. A manifest stops before any design
only when preflight finds an output-directory collision that could make two
jobs overwrite one another, or when no valid executable rows remain. Invalid
accessions and empty-gene rows are reported and excluded while other valid rows
continue. Fix the reported input or collision and rerun.

### All linked transcripts does not start

Do not combine `--all-transcripts` with `--accession` or `--fasta`. In the webapp,
choose Design one target, enter a gene, enable Design all linked transcripts in
Advanced settings, and leave accession blank.

### A plot still looks like pastel GC coloring

Pass both options:

```bash
hcrprobeforge --plot-only probes.tsv \
  --plot-theme minimal --plot-color-by dtm
```

The theme changes style; --plot-color-by dtm changes probe-block colors. The
webapp forwards both to the same plotting function.

### Fewer probe pairs are returned than requested

Read the curation report. Geometric capacity, sequence filtering, repeats,
specificity, or QC can limit the final set. Smart mode can try configured
fallback tiers but cannot create acceptable non-overlapping sequence that does
not exist.

### The webapp port is already in use

```bash
hcrprobeforge-web --port 8780
```

### An updated webapp page is not visible

Stop and relaunch the webapp after installing an update, then refresh the
browser. A running Python process continues using the package version that was
loaded when it started.

## Architecture and data flow

The package has one scientific pipeline and several ways to run it. The CLI,
all-channel helper, manifest helper, and webapp use the same transcript
resolution, candidate generation, QC, curation, and output rules. Their default
selection policies are different: direct `hcrprobeforge` is one-pass, while
the all-channel helper, manifest helper, and webapp use smart auto-curation by
default. The webapp collects and validates form values, starts the selected
workflow, shows its progress, and presents the files produced by that same
pipeline.

### Normal design data flow

The following order is important when interpreting a run:

1. The selected species is resolved to an organism name, assembly record, and
   HCRProbeDesign index alias. For example, the user-facing *Xenopus
   tropicalis* choice uses the registered `xtr` reference.
2. A gene symbol is resolved through NCBI Gene and linked RefSeq records, or an
   explicit accession is validated directly. A supplied FASTA bypasses NCBI
   transcript selection.
3. The selected exact transcript accession and organism metadata are checked
   before a transcript FASTA or annotation is accepted. A record naming a
   different built-in species is rejected before it can be cached or designed.
4. The exact versioned transcript sequence is retrieved from the project cache
   when it is valid, or downloaded from NCBI and then cached. Cache validation
   checks the accession-version, FASTA header, sequence alphabet, length, and
   stored SHA-256.
5. `designProbes` receives the same sequence, channel, species/index alias,
   tile, GC, Gibbs, C/G-run, dTm, specificity, and candidate-count arguments
   selected by the user. This is the candidate-generation boundary.
6. HCRProbeForge reads the native candidate output, attaches reporting fields,
   and runs Primer3 on the actual initiator-bearing P1/P2 order sequences.
7. Smart curation, when requested, combines eligible candidate reservoirs and
   selects a non-overlapping set using the documented QC and coverage rules.
   One-pass mode preserves the native candidates instead.
8. The relevant candidate or selected-pair tables, order table, map, QC
   workbook, reports, and provenance JSON are written after processing succeeds.

### Why the species choice affects both NCBI and Bowtie2

The species choice has two separate jobs that must agree. Its scientific name
is sent to NCBI so a gene symbol or accession resolves in the intended taxon;
its registered HCRProbeDesign alias selects the Bowtie2 genome index used for
genomic masking. The package records both values, plus the assembly and
assembly accession, in run metadata. This prevents a transcript from one
species from being silently screened against another species' genome.

The built-in preset is an organism/assembly pairing, not a promise that every
NCBI record is interchangeable. When another assembly is considered, check
that NCBI provides compatible RefSeq transcript and feature records first.
Fetch the assembly with a new index alias so the existing index remains
available:

```bash
hcrprobeforge index fetch --species xla \
  --assembly-accession GCF_<ALTERNATE_ASSEMBLY_ID>.<VERSION> \
  --assembly-name <EXPECTED_ASSEMBLY_NAME> \
  --index-alias xla_alternate
```

Replace the accession and expected name with a real compatible assembly; the
values above are placeholders.

The alias can then be supplied to a design as a registered species. The
transcript organism remains `Xenopus laevis`; the alternate assembly is
recorded in the installed-reference metadata. Never reuse an alias for a
different assembly unless replacement is intentional and `--force` is used.

### Webapp field behavior

The webapp is organized into five tabs. Each tab enables only the controls
that are meaningful for that operation, and controls for other workflows are
not submitted:

| Tab | Required input | Applicable optional controls |
|---|---|---|
| Design one target | species, channel, gene/accession/one-record FASTA | transcript policy, genome/index options, design thresholds, smart curation, plot and final-QC options |
| Design a gene list | species, channel, manifest file | transcript policy, design thresholds, smart curation, plot and final-QC options |
| Plot a probe table | existing probe TSV/text table | transcript length, plot theme, color, DPI, title, label visibility |
| QC an oligo table | existing IDT or selected-pairs table | QC workbook output path |
| Build a genome index | species preset and either a preset assembly accession or local genome FASTA | expected assembly name, alternate accession, alias, matching GFF3 for local builds, threads, rebuild |

Species and NCBI organism are shown before workflow-specific fields. A built-in
preset supplies its saved assembly information for index building. A custom
preset can save that information in the species dialog, but its NCBI organism
alone is not a genome source: provide either a versioned assembly accession or
a local genome FASTA before building. The expected assembly name, local FASTA,
and alias are optional controls when an alternative source or name is desired.
The **Index ready on this computer** status means that the registered Bowtie2
files are present on the computer; metadata alone is not treated as a usable
index.

The **Inspect input** button performs workflow-specific structural checks
without starting `designProbes`, downloading a genome, or writing a
scientific result. It is a convenience check, not a replacement for NCBI
validation or a completed index build. It checks the species/organism pairing,
the appropriate input file and required columns, FASTA record count and
sequence alphabet, manifest channel and path rules, plot coordinate bounds,
QC table readability, numeric ranges, allowed enumerated values, and
assembly-accession syntax. NCBI-backed gene and accession checks still happen
at execution time, where true species mismatches are skipped and unresolved
genes/accessions are failed individually.

### Manifest data flow and resume behavior

Manifest rows are normalized into a run record containing the original gene
label, an optional exact accession, an optional single HCR channel, and a safe
directory label. The original label remains in metadata and status
tables; only filesystem path components are normalized. Preflight rejects
path-normalization collisions before design. Exact duplicate
gene/accession/channel rows are represented by the first executable job; only
that first occurrence has a resume marker. The same target may be listed for
different channels.

For each valid row, the runner:

1. checks the stored done marker and validates that any existing output belongs
   to the selected species;
2. resolves a gene symbol within the selected organism, or validates an exact
   accession against the selected organism before writing transcript data;
3. chooses the row's channel when present, otherwise uses the global
   `--channels` selection;
4. runs the selected channel(s) through the same single-target pipeline;
5. records the effective accession, channel, output roots, status, and
   summaries; and
6. writes a done, failed, or skipped marker that allows the manifest history
   to be inspected and resumed.

An explicit accession is one exact transcript and is not expanded by
`--all-transcripts`. Gene-only rows may use the all-transcripts option. A
temporary NCBI failure uses exit code 75 and is eligible for the helper's
bounded retry. A true species mismatch is skipped after the organism check;
an unresolved gene or accession is a failed row. Design, QC, parsing, and
other failures remain visible as failed rows and do not abort later jobs. The
underlying design invocation is unchanged: each selected row/channel still receives the same sequence,
species, index, and design settings as a single-target run.

## Failure safety and reproducibility

The package uses best-effort rollback for a new direct design, plot, QC run, or
index build. Before work begins it records the existing filesystem entries and
backs up existing output files that may be replaced. If the operation raises
an error, newly created entries are removed and overwritten output files are
restored. This includes partial FASTA files, transcript annotations, candidate
tables, maps, QC workbooks, index metadata, and partial Bowtie2 files created
by a failed new build.

Rollback is scoped to the current operation and is not a replacement for
backups. Do not run two writers against the same project directory at the same
time; concurrent command-line runs are not globally locked or isolated.

Manifest status tables, failed markers, and manifest logs are intentionally
retained because they describe which row failed and allow a catalogue to be
resumed. They are bookkeeping, not scientific result files. A failed direct
design does not retain a per-target log, so the error message does not direct
the user to a nonexistent log. In particular, a missing species index is
reported by name with instructions to use the Build a genome index workflow.

The package does not silently fall back to a different species, assembly, or
index. If an exact NCBI accession's organism disagrees with the selected
organism, the run stops before caching or design. A user-supplied FASTA is
accepted as a custom sequence after sequence validation; its biological origin
must be checked by the user because it has no NCBI organism metadata.

## Preserving provenance and getting help

For a published probe set, keep the selected-pairs table, order file, final QC
workbook, `reference.json` or registered index prefix, transcript-resolution
JSON, run parameters, curation report, exact input, species, assembly,
transcript accession-version, and index alias together. This information makes
the design reproducible and allows a later reviewer to distinguish a changed
transcript, assembly, channel, or curation setting from a changed sequence.

When reporting a problem, include the HCRProbeForge version, Python version,
operating system, command, species and assembly metadata, and a redacted error
message. For Pre-mRNA problems, also include the transcript
accession-version, selected intron numbers, and the genome/GFF3 assembly
identity. Do not include NCBI credentials or private sequences in a public
issue.

## License and citation

HCRProbeForge is distributed under the MIT License; see `LICENSE`. The
`hcrprobedesign` and `primer3-py` Python packages are installed automatically
by pip as HCRProbeForge dependencies. Bowtie2 is an external executable and is
not bundled; install `bowtie2` and `bowtie2-build` separately, for example
through Conda or bioconda. Review the exact licenses of the versions installed
in your environment; the package's dependency notices identify the declared
dependencies and known external tools but do not replace their license files.
HCR v3 initiator sequences, amplifier reagents, and related protocol materials
may have separate provider terms or intellectual-property restrictions; users
must verify the terms that apply to their intended assay and supplier. NCBI
data and RefSeq records are subject to NCBI and database-provider usage terms;
retain the accession and assembly provenance with published designs. Cite
HCRProbeForge together with HCRProbeDesign, Bowtie2, Primer3, NCBI RefSeq, and
the exact genome assembly used for specificity screening. `CITATION.cff`
contains the software citation metadata for the current release.

## Additional project documents

The following repository documents provide complementary information:

- [`CONTRIBUTING.md`](CONTRIBUTING.md) describes contribution expectations and
  scientific regression-test requirements.
- [`DEVELOPMENT.md`](DEVELOPMENT.md) is the maintainer and release-validation
  guide, including clean builds and artifact checks.
- [`CHANGELOG.md`](CHANGELOG.md) lists user-visible changes by release.
- [`VERSION_HISTORY.md`](VERSION_HISTORY.md) provides a compact release history.
- [`VALIDATION_REPORT.md`](VALIDATION_REPORT.md) records the reproducible checks
  performed for the current release.
- [`THIRD_PARTY_NOTICES.md`](THIRD_PARTY_NOTICES.md) summarizes dependency,
  executable, data, and assay-material licensing considerations.
