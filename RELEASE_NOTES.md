# HCRProbeForge 1.3.7

This release fixes probe-number placement on maps that contain adjacent
adaptive tile lengths. Probe boxes whose rendered widths are within roughly
15% of one another now share one inside-or-outside label mode, so 50- and
52-nt tiers are not split by unrelated 38- or 40-nt boxes. Materially
different box-size groups may still use different modes when their available
space differs.

The release also adds a regression test for the mixed-tier Pre-mRNA map case
and updates the README, changelog, version history, and validation report to
document the behavior. Mature-transcript design generation and selection are
unchanged.

## Validation

- 69 unit tests pass.
- The package compiles successfully.
- The 1.3.7 wheel and source distribution pass the release-artifact hygiene
  checker.
- A clean-environment wheel smoke test reports version 1.3.7 and exercises
  the CLI help and index-list commands.
