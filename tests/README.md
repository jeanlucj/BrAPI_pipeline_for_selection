# Tests

A [testthat](https://testthat.r-lib.org/) suite in three tiers, simplest first.

## Run

```bash
Rscript tests/run_tests.R                       # tiers 1 + 2 (offline, ~2 s)
RUN_LIVE_TESTS=true Rscript tests/run_tests.R   # also tier 3 (live network + downloads)
```

Or from R: `testthat::test_dir(here::here("tests","testthat"))`.

## Tiers

- **Tier 1 — `test-01-*.R` (offline, fast).** Pure helpers and logic on synthetic
  inputs: `.Gmatrix`/`std_grm`, the EM combiner, VCF parsing/thinning/merge/QC,
  pedigree-matrix loading, coverage-table parsing, trial selection (radius +
  `STUDY_NAMES` + multi-type, via a fake connection), phenotype/design extraction,
  Stage-1 BLUEs, the Stage-2 helpers, and `select_parents`.
- **Tier 2 — `test-02-models.R` (offline, heavier ~1 s).** Runs the real models on
  synthetic data: BGLR RKHS prediction + cross-validation, and a fully *mocked*
  `find_and_get_genotypes` exercising the multi-protocol EM-combine and
  single-protocol paths (synthetic VCFs, no network).
- **Tier 3 — `test-03-live.R` (network + downloads, opt-in).** Hits T3/Oat:
  connect, find NY trials, pull phenotypes, build a GRM from the Oat 3K protocol.
  Skipped unless `RUN_LIVE_TESTS` is set.

## How it works

`helper-setup.R` sources `code/` once, redirects `cache_path()`/`output_path()` to a
tempdir (so tests never touch `data/`/`output/`), and provides the synthetic-data
and mock-connection builders (`write_test_vcf`, `make_dosage`, `fake_conn`, …). Tests
pass `refresh = TRUE` so cached results are not reused.
