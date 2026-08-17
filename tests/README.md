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
  inputs: the `data/config/` list readers (`config_lines`/`config_traits`),
  `t3_login()`'s credential handling (stubbed connection, no network), the
  progress/timing helpers (silent and value-transparent when reporting is off),
  `.Gmatrix`/`std_grm`, the EM combiner, VCF parsing/thinning/merge/QC,
  pedigree-matrix loading + companion-file contracts, name→dbId resolution,
  coverage-table parsing, trial selection (radius + `TRAINING_TRIALS` + `TEST_TRIALS`
  roles + multi-type, via a fake connection), phenotype/design extraction with exact
  `trait_names` matching + `split_by_role`, Stage-1 BLUEs, the Stage-2 helpers, and
  `select_parents` (un-standardized index + two-block `breeders_output.csv`).
- **Tier 2 — `test-02-models.R` (offline, heavier).** Runs the real models on
  synthetic data: BGLR prediction + cross-validation, and a fully *mocked*
  `find_and_get_genotypes` exercising the multi-protocol EM-combine, single-protocol,
  prediction-target subsetting, pedigree-bridged test accessions, and training
  injection paths (synthetic VCFs, no network), and `test-02-partial-cache.R` proving
  step 4's intermediate caches are actually used (VCF parses are counted, not assumed)
  and change nothing about `G`. It ends with an **oracle** test built on
  `simulate_trials()`: a related population with known true breeding values, where both
  engines must recover those values on *held-out* lines and agree with each other. That
  is the regression that catches a kernel silently decoupled from the phenotypes — every
  shape stays correct while accuracy collapses to zero.
- **Tier 3 — `test-03-live.R` (network + downloads, opt-in).** Hits T3/Oat:
  connect, find NY trials, pull phenotypes, build a GRM from the Oat 3K protocol.
  Skipped unless `RUN_LIVE_TESTS` is set **and** `T3_USERNAME`/`T3_PASSWORD` are
  available (T3/Oat requires a login; see `.Renviron.example`).

## How it works

`helper-setup.R` sources `code/` once, redirects `cache_path()`/`output_path()` to a
tempdir (so tests never touch `data/`/`output/`), sets
`options(brapi.progress = FALSE)` so no test run is narrated, and provides the
synthetic-data builders — including `simulate_trials()`, which generates a **related**
population of inbred lines (founders → biparental crosses → doubled haploids), gives
every marker an effect, and derives phenotypes from the resulting genetic values, so
that the marker GRM is the correct kernel and Stage 2 has real signal to find (and a
known answer to be checked against) — and the synthetic-data
and mock-connection builders (`write_test_vcf`, `make_dosage`, `fake_conn`, …). Tests
pass `refresh = TRUE` so cached results are not reused.
