# Changelog

## 1.3.7

- Fix probe-number label grouping when adjacent adaptive tile tiers are
  present, so similarly sized boxes such as 50 and 52 nt use one placement
  mode even when shorter 38- or 40-nt boxes are also on the map.
- Add a regression test for mixed short-tile Pre-mRNA maps and refresh the
  release validation metadata.

## 1.3.6

- Group probe-number labels by rendered box size so probes of similar size use
  one consistent inside or outside placement, while materially different
  tile-length groups can use different modes.
- Update the webapp's pre-designing notice to explain the normal genome-index
  requirement, the **Build a genome index** workflow, and the trade-off of
  disabling genome masking.
- Add regression tests for grouped probe-label placement and refresh release
  metadata and validation documentation.

## 1.3.5

- Correct Pre-mRNA intronic geometric-capacity reporting and adaptive goals to
  count non-overlapping tiles within each eligible intron and each contiguous
  unambiguous sequence segment, rather than counting every sliding-window
  start.
- Keep transcript-wide coverage balancing disabled for intronic Pre-mRNA for
  the documented reason—equal capacity-aware intron quotas—and no longer add a
  misleading short-transcript explanation to those reports.
- Clarify manifest stop conditions, one-pass `--target-probes` behavior,
  path-normalization collision wording, environment-variable formatting, and
  the webapp `index_metadata` location in the public README.

## 1.3.4

- Fix map annotation labels so each mature-transcript feature is placed inside
  its own region when it fits, instead of allowing one narrow feature to move
  every label outside the transcript line.
- Keep narrow Pre-mRNA exon and intron regions uncluttered by omitting labels
  that cannot fit inside their own region; labels that do fit remain inside the
  corresponding feature.
- Use one inside/outside probe-number mode for similarly sized probe boxes,
  while retaining per-box fallback for maps with materially different widths.
- Add plotting regression tests, a dependency-free release-artifact checker,
  and synchronized publication/development documentation.

## 1.3.3

- Fix a Pre-mRNA intronic auto-curation defect that could reselect a candidate
  overlapping an interval already chosen during intron-quota selection. All
  final intronic selections now enforce interval non-overlap across quota and
  fallback passes.
- Normalize unambiguous upstream exclusive coordinate suffixes so probe names
  agree with the package's 1-based inclusive `start`/`end` fields.
- Distinguish true species mismatches from unresolved genes/accessions in
  manifests; unresolved inputs now remain failed rows and return a nonzero
  status instead of being silently classified as skips. Invalid preflight rows
  also produce a nonzero result.
- Keep map label placement consistent for similarly sized annotation or probe
  boxes, with deterministic dynamic callout lanes when labels do not fit.
- Record HCRProbeForge, HCRProbeDesign, Bowtie2, and Bowtie2-build provenance
  alongside existing run metadata and refresh publication documentation.

## 1.3.2

- Make manifest resume markers parameter-aware so changed design, QC, curation,
  channel, species, or package-version settings cannot silently reuse stale
  completed jobs.
- Clarify installation prerequisites, project output roots, cache behavior,
  pre-mRNA coordinates, manifest parsing, QC terminology, index options, and
  webapp security in the public documentation.
- Replace placeholder and assembly examples with clearly labeled user-facing
  examples and correct the documented output layout.

## 1.3.1

- Fix persistent annotation-database locking so completed SQLite indexes can
  be validated and reused immediately.
- Make Pre-mRNA one-pass QC and order output use only candidates eligible for
  the selected intronic or whole-genomic target scope, while preserving the
  mature-transcript pipeline.
- Accept ready local reference metadata for arbitrary registered species
  aliases during Pre-mRNA preparation.
- Require exact accession-version matches when a version is explicitly
  requested; unversioned requests may still resolve to the current version.
- Preserve user-provided compressed genome/annotation inputs and restore
  replaced reference files if an index build fails after a partial commit;
  failed design reruns also restore overwritten files in the visible output
  tree.
- Replace the mixed internal/public README with a general-user guide and move
  maintainer procedures into contributor and development documentation.

## 1.3.0

- Make `xtr` the only built-in Xenopus tropicalis index alias and suppress
  stale duplicate installed-reference rows that describe a built-in assembly.
- Improve publication documentation, including the logo link, working-directory
  examples, Conda environment name, and regular webapp startup instructions.
- Keep the WSL browser/folder integration, protected built-in preset behavior,
  consistent map labels, and Pre-mRNA workflow fixes from the previous release.

## 1.2.32

- Add OS-aware WSL browser and Windows Explorer integration while preserving
  macOS's native browser and folder behavior.
- Make `xtr` the canonical Xenopus tropicalis alias and keep built-in preset
  discovery limited to the supported alias.
- Allow removal of built-in Bowtie2 index files without removing built-in
  presets or downloaded reference assets.
- Keep dynamic map labels but use consistent probe-number typography inside and
  above narrow boxes.
- Add `environment.yml`, remove Conda recipe publication instructions, and
  update the public installation/troubleshooting documentation.

## 1.2.31

- Removes the remaining full SQLite validation from setup-page rendering and the live `/index-status` badge endpoint. Both now use the atomically committed reference metadata and database path, eliminating the post-index return lag and stale badge state.
- Keeps full annotation-database structural validation in the design execution path, so the performance fix does not weaken scientific or database correctness checks.
- Adds regression tests for alternate installed references and live badge polling.

## 1.2.30

- Fix setup-page return latency after index builds by avoiding SQLite schema
  scans during the initial form render.
- Resolve Pre-mRNA readiness for custom and alternate assemblies through their
  preset/display metadata, refresh status on browser-history restore, and
  prevent stale setup/status pages from browser caches.
- Keep the Start another run link cache-busting and document the server-restart
  requirement when installing a new package version.

## 1.2.29

- Make the setup page return promptly after an index build by reusing one
  reference-metadata scan while rendering readiness badges.
- Refresh index and Pre-mRNA database status immediately when the setup page
  opens, and resolve alternate-species annotation status from the selected
  registered reference rather than a built-in preset.
- Start single-channel progress in the normal starting range instead of
  displaying the later 60% design phase before design has begun.
- Add MIT licensing and public-release metadata for GitHub, PyPI, and Conda
  packaging without changing scientific design behavior.

## 1.2.28

- Fix explicit multi-intron Pre-mRNA selections by sanitizing only the
  `designProbes` engine target name while retaining readable comma-separated
  intron scopes in project folders and reports.
- Restore adaptive exon/intron numbering on Pre-mRNA maps when the rendered
  interval has enough room, with collision-aware fallback labels.
- Use rendered text and box widths for probe-number placement so labels remain
  inside genuinely wide boxes and move above only when needed, without white
  label padding.
- Clarify the Design target guidance and update it dynamically for mature,
  intronic Pre-mRNA, and whole Pre-mRNA modes.
- Preserve mature-transcript scientific behavior and outputs.

## 1.2.27

- Keep probe numbers on sufficiently wide map boxes and use uncluttered
  callouts only for narrow boxes; remove unnecessary white label padding.
- Label Pre-mRNA map tracks as singular `exon` and `intron`, describe the
  target as “intronic regions”, and give explicit-intron output folders names
  such as `premrna_introns2,3`.
- Fix Pre-mRNA database readiness for registered alternate species aliases,
  including databases built for C. elegans, and make single-channel progress
  use the underlying phase progress directly.
- Write transcript-model caches as indented UTF-8 JSON, repair stale whole-mode
  eligibility metadata, and preserve the mature transcript length in
  Pre-mRNA metadata.
- Keep the mature-transcript pipeline and scientific candidate behavior
  unchanged; these are metadata, presentation, status, and navigation fixes.

## 1.2.26

- Validate annotation databases before reuse and automatically rebuild
  interrupted or incomplete SQLite builds instead of reporting them as ready.
- Keep Pre-mRNA model and assembled-target caches in the project-level
  `hcr_results/cache/<species>/<workflow>/premrna/` tree, and remove managed
  legacy cache artifacts after migration.
- Add a generalized whole-genomic Pre-mRNA target mode alongside mature and
  intronic modes. Whole mode allows candidates across exons, introns, and
  boundaries; intronic mode retains its existing complete-intron behavior.
- Add a live Pre-mRNA database badge, clearer target guidance, channel-aware
  monotonic all-channel progress, and collision-resistant Pre-mRNA maps.
- Preserve mature-transcript generation, QC, curation, and specificity
  behavior; add regression coverage for the new target and database paths.

## 1.2.25

- Make the genomic annotation SQLite database an explicit opt-in during index
  builds. The webapp now explains the trade-off before the build starts; mature
  transcript users can build only the Bowtie2 index, while the first later
  Pre-mRNA design builds the database when needed.
- Optimize annotation indexing with a single streamed GFF3 pass, in-memory
  transcript/alias resolution, batched SQLite inserts, deferred index creation,
  cancellation checks, atomic replacement, and byte-aware progress updates.
  The retained transcript/exon relationships and scientific target selection
  are unchanged.
- Reuse and migrate older `premrna/` genome and annotation assets into the
  canonical assembly directory, remove duplicate compressed/legacy copies, and
  avoid downloading assets that are already available.
- Keep mature-transcript behavior unchanged and add regression coverage for
  opt-in/deferred annotation builds, progress, sidecar cleanup, and legacy
  asset reuse.

## 1.2.24

- Generalize Pre-mRNA design beyond Xenopus tropicalis: the target now contains
  the complete genomic transcript (exons plus introns), while candidate
  selection remains restricted to complete selected introns.
- Make Pre-mRNA smart auto-curation optional, retain intron-balanced quotas
  when it is enabled, and add strand-aware exon/intron maps and metadata.
- Prepare genome FASTA and matching GFF3 annotation during NCBI index builds;
  retain only uncompressed reference assets, build a compact SQLite annotation
  lookup, and cache assembled Pre-mRNA targets for later runs. Existing ready
  indexes are still reusable and receive missing reference assets on their first
  Pre-mRNA design.
- Add local matching GFF3 support to local index builds, improve phase-aware
  webapp progress, and keep mature-transcript design behavior unchanged.

## 1.2.23

- Fix xtr10 Pre-mRNA readiness when a complete HCRProbeDesign index exists
  without an HCRProbeForge reference-metadata record. The built-in xtr10
  assembly metadata is now inferred only after the registered index files are
  verified, allowing existing indexes to be reused without rebuilding.
- Improve the missing-index diagnostic while keeping mature-transcript and
  all existing design behavior unchanged.

## 1.2.22

- Add an opt-in `Pre-mRNA` target mode for the `xtr10` Xenopus tropicalis
  preset. The selected RefSeq transcript is mapped to the same registered NCBI
  genomic assembly used for specificity screening, and its introns are
  extracted in transcript orientation.
- Design across all selected introns with equal, capacity-aware probe quotas,
  or restrict the run with comma-separated 1-based intron numbers. Introns are
  separated by `N` spacers so candidates cannot cross intron boundaries.
- Add compact conditional CLI/webapp controls, Pre-mRNA metadata, intron-aware
  probe tables, and intron-segment maps. Mature-transcript and FASTA workflows
  retain their existing design, QC, and selection behavior.

## 1.2.21

- Remove browser unload/page-hide shutdown requests so closing or reloading a
  tab cannot stop the server or interrupt active worker threads.
- Refuse explicit webapp shutdown and custom-species add/remove operations
  while a job is queued, running, or cancelling; the server remains available
  until active work is finished.
- Replace request-scoped staging directories with automatic temporary-folder
  cleanup and preserve uploaded inputs only below the output that records them.
- Drop the artificial Python upper bound, keeping the declared requirement at
  Python 3.10 or newer.
- Consolidate helper-command and webapp lifecycle guidance in README.md and
  remove the redundant standalone helper guide.

## 1.2.20

- Consolidate redundant global webapp CSS declarations while retaining the
  responsive media-query overrides and visual behavior.
- Correct the validation report so the current release scope is separated from
  the accumulated history of earlier scientific and UX changes.
- Remove the unsupported `scripts/` folder claim from the helper guide.
- Add PyPI-oriented project metadata, source-distribution inclusion rules, and
  documented build/check commands without changing scientific behavior.

## 1.2.19

- Add per-process webapp capability tokens, strict cookies, Host validation,
  and Origin checks for browser state-changing requests, protecting runs,
  species changes, folder actions, cancellation, and shutdown from unrelated
  webpages.
- Require an explicit access token for non-loopback webapp binding and document
  the trusted-network deployment contract without adding an authentication
  prompt to the normal local workflow.
- Retry the requested webapp port across a bounded range when another process
  already occupies it, then report the actual URL selected.
- Add HTTP-level security and occupied-port regression tests without changing
  scientific processing or pipeline parameters.

## 1.2.18

- Remove deleted custom aliases from HCRProbeDesign configuration and remove
  empty alias folders after exact Bowtie2 cleanup, so the same preset can be
  built again without `--force`.
- Repair stale custom registrations left by older releases before a build,
  including interrupted temporary index files, while protecting built-in and
  shared aliases.
- Extend the confirmation preview and regression coverage to include config
  entries, backups, alias folders, and stale-registration migration.

## 1.2.17

- Delete the exact unshared Bowtie2 index files associated with a removed
  custom species, along with its saved FASTA and HCRProbeForge reference
  metadata; built-in preset indexes remain protected.
- Add a pre-deletion preview to the webapp confirmation dialog showing the
  private files, index paths, and any shared or out-of-root files retained for
  safety.
- Add regression coverage for custom cleanup, built-in protection, and shared
  index preservation without changing scientific design behavior.

## 1.2.16

- Fix NCBI Assembly genome downloads for assembly names containing spaces or
  punctuation by deriving the FASTA filename from the authoritative FTP
  directory and applying NCBI's underscore-safe filename convention.
- Validate that every ESummary response matches the requested assembly
  accession instead of trusting the first Entrez search result.
- Retry transient genome-download responses, write downloads atomically, and
  reject empty or truncated files before index construction.
- Try a bounded set of valid NCBI filename candidates and, as a final
  fallback, discover a matching `*_genomic.fna.gz` link from the Assembly FTP
  directory listing.
- Improve custom-species and failed-download diagnostics, including the local
  FASTA fallback when an Assembly record cannot provide a usable genome file.
- Preserve probe generation, Primer3 QC, auto-curation, selection, annotation,
  plotting, and all existing scientific parameters.

## 1.2.15

- Keep one user-facing custom-species option per preset; internal assembly aliases no longer leak into the Species menu.
- Use the custom preset display name for project output, cache, and index-metadata folders while retaining aliases for internal index registration.
- Allow deletion of user-created presets whether or not their index is ready, with a confirmation dialog; built-in presets remain protected and existing Bowtie2 files are not silently deleted.
- Improve missing-annotation diagnostics for non-standard species by suggesting a gene accession or FASTA input after symbol/RefSeq lookup fails.
- Give the index-folder action a larger, evenly distributed inline layout and align advanced checkboxes to the same input-control row as the index rebuild control.

## 1.2.14

- Use the same vertical alignment model for advanced-settings checkboxes as the working index rebuild control.
- Add a red required marker to the custom-species NCBI assembly accession field and keep its local-FASTA exception explicit.
- Add an inline button for opening the HCRProbeDesign index data folder.
- Refresh index readiness from the running webapp so manually removed Bowtie2 files immediately change the species ribbon and are rejected before design.
- Close the local webapp server when its browser page is closed, while preserving normal internal navigation.
- Add accessible tooltip help to the custom-species removal control.

## 1.2.13

- Check for a missing genome index before NCBI gene lookup so users receive the correct build-index instruction.
- Clarify that custom presets need either an NCBI assembly accession or local genome FASTA before indexing.
- Replace the custom-preset removal button with a compact adjacent × control.
- Align advanced-settings checkboxes with neighboring form controls.

## 1.2.12

- Fixed custom-species index guidance when a preset has no assembly accession or local genome FASTA.
- Added safe removal of user-created presets that do not have a ready registered index.
- Fixed add-species tooltip clipping and aligned index/advanced-settings checkboxes.

## 1.2.11

- Add persistent user-created species presets through the CLI and webapp.
- Allow custom presets to use an NCBI assembly accession or an uploaded/local
  genome FASTA in the existing genome-index workflow.
- Keep custom species in the species menu after creation, with required and
  optional fields clearly marked and tooltip guidance for each reference input.
- Keep transcript lookup, species validation, caching, output organization, and
  the scientific probe-design pipeline unchanged.

## 1.2.10

- Add optional per-row HCR channel support to manifest workflows. Headered and
  headerless gene/accession lists can assign B1-B5 independently while rows
  without a channel continue to use the command-level channel selection.
- Skip species-incompatible manifest rows individually, record the reason in
  manifest status/preflight files, and continue valid rows without creating
  scientific outputs for the skipped targets.
- Extend input inspection to validate the input contract for design, manifest,
  plot, QC, and genome-index workflows before execution.
- Distribute workflow tabs evenly across the webapp tab card, retaining
  responsive layouts on narrower screens.
- Preserve the HCRProbeDesign, Primer3 QC, curation, selection, annotation,
  and plotting pipeline.

## 1.2.9

- Add Xenopus laevis as the second Xenopus preset, using NCBI's
  `Xenopus_laevis_v10.1` assembly (`GCF_017654675.1`). The preset is placed
  immediately below Xenopus tropicalis in the webapp and uses `xla` only as
  its internal/index alias.
- Keep the existing Xenopus tropicalis `xtr10` index compatibility while
  showing only `Xenopus tropicalis` in the species dropdown and writing
  readable species-specific output paths.
- Hide Smart auto-curation reliably outside design and manifest workflows;
  disabled hidden checkbox rows no longer leak through the shared checkbox
  display rule.
- Replace the generic missing-log `designProbes` error with a message that
  names the missing species index and directs users to the Build a genome
  index workflow. Failed direct runs continue to roll back their new files.
- Preserve the HCRProbeDesign command, candidate generation, Primer3 QC,
  curation, selection, annotation, and plotting behavior.

## 1.2.8

- Make workflow-specific visibility explicit in the webapp. Design and
  manifest settings, plot settings, and QC settings are now isolated to the
  workflows that use them; hidden controls are disabled before submission.
- Hide smart auto-curation outside design and manifest workflows, including
  genome-index builds.
- Tighten the logo canvas and header spacing so the supplied logo sits closer
  to the introductory copy without crowding it.
- Add a served SVG favicon using the circular logo mark to the setup,
  progress, result, and error pages.
- Preserve all workflow routing and scientific design behavior.

## 1.2.7

- Use the newly supplied publication logo artwork, including its Helvetica
  Bold typography, pastel palette, and updated probe-mark geometry.
- Enlarge and rebalance the webapp header so the logo sits comfortably beside
  the introduction at desktop and mobile widths.
- Center the input inspection action and place its feedback beneath the action
  when it is shown.
- Align checkbox controls, labels, and help markers on one readable baseline
  with consistent spacing and wrapping in every workflow and advanced-settings
  panel.
- Preserve all workflow routing and scientific design behavior.

## 1.2.6

- Use the supplied HCRProbeForge logo artwork in the webapp and package
  assets, with portable system fonts and the established pastel palette.
- Normalize form controls to consistent heights and grid alignment, including
  inline checkbox labels and help controls on desktop and mobile layouts.
- Rename the QC workflow tab to QC an oligo table for clarity.
- Remove empty target directories left behind when a batch job is rejected
  before scientific design, such as a cross-species accession mismatch.
- Roll back newly created scientific outputs, caches, and failed index metadata
  when a direct, webapp, batch, or index operation fails; existing results are
  preserved and manifest diagnostics remain available.
- Preserve the existing species guard and all scientific design behavior.

## 1.2.5

- Reject exact transcript accessions whose NCBI species metadata disagrees with
  the selected species before the sequence can be cached or designed. This
  prevents a mixed-species manifest from producing results under the wrong
  genome reference.
- Keep workflow directories consistent and readable: `individual_designs`,
  `manifests`, `plot`, and `qc` are used for new runs and caches.
- Stop creating empty plot and QC cache directories when those workflows do
  not write reusable cache data.
- Write transcript-annotation JSON as UTF-8 so labels such as `5′ UTR` and
  `3′ UTR` are readable in text editors; the annotation values and scientific
  behavior are unchanged.
- Scope manifest report discovery to the jobs and channels in the current run,
  including resumed jobs, and preserve mixed gene-only/exact-accession
  manifests.
- Improve workflow tabs, help text, layout, and technical documentation while
  keeping the design, QC, curation, selection, and plotting pipeline intact.

## 1.2.4

- Reorganized project output into `hcr_results/runs/<species>/...` and
  `hcr_results/cache/<species>/...`, with manifest status tables grouped under
  `status/tables/` and technical target files kept under `details/`.
- Fixed gene-list/all-transcript routing so batch jobs use the same core
  pipeline without relying on a removed output helper.
- Made direct CLI cache placement follow the selected `--outdir` while keeping
  `HCRPROBEFORGE_CACHE_DIR` as an explicit override.
- Added workflow tabs, workflow-specific report sections, explicit index-prefix
  and metadata paths, and clearer webapp validation.
- Rewrote the public README with complete installation, indexing, CLI, webapp,
  output, and troubleshooting guidance.
- Preserved the scientific candidate generation, Primer3 QC, adaptive tiers,
  selection, annotation, and plotting behavior.

## 1.2.3

- Added Human (`Homo sapiens`) as a fifth automatic species preset using the
  current NCBI RefSeq-compatible GRCh38.p14 assembly.
- Generalized the webapp launch guide and species help so they describe the
  package for all users without lab-specific wording.

## 1.2.2

- Fixed automatic NCBI Assembly lookup by resolving `GCF_`/`GCA_` accessions
  to Entrez Assembly UIDs before requesting ESummary records.
- Added bounded retry handling and clearer errors for transient NCBI failures.
- Kept `xtr10` as the HCRProbeDesign index alias while writing Xenopus output
  below `hcr_results/Xenopus_tropicalis/`.
- Updated the zebrafish preset to the current NCBI RefSeq reference assembly
  GRCz12ab (`GCF_052040795.1`); alternate assemblies remain available through
  the explicit assembly options.

## 1.2.1

- Added automatic NCBI genome-reference retrieval and Bowtie2 index building for
  Xenopus tropicalis (`xtr10`), Zebrafish, Mouse, and Chicken. Xenopus is the
  first preset and the CLI default so existing xtr10 installations continue to
  be the primary workflow.
- Restored xtr10 as a first-class preset using the existing UCB_Xtro_10.0 /
  GCF_000004195.4 reference metadata; it can be reused or rebuilt from the
  webapp and CLI like the other presets.
- Added alternate-assembly accession, expected-name, alias, and force-rebuild
  controls in the CLI and webapp.
- Added registered-reference discovery so alternate aliases can be selected for
  future transcript and probe-design runs.
- Scoped design output below hcr_results/<species> and recorded species,
  organism, assembly, accession, and index alias in metadata.
- Updated input inspection, webapp help, documentation, package metadata, and
  release validation records.
- Preserved the existing candidate generation, QC, curation, selection,
  annotation, plotting, and reporting pipeline.
