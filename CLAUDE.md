# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

A pipeline that pulls oat data from **T3/Oat** (`oat.triticeaetoolbox.org`) over
**BrAPI**, runs a two-stage genomic-prediction analysis, and ranks accessions as
parents for crossing for the **New York** environment. It is a
[workflowr](https://github.com/workflowr/workflowr) project.

Three living docs describe it; **keep them current when you change the pipeline**
(see `README.md` = usage, `DESIGN.md` = structure/data objects, `BACKGROUND.md` =
theory + design rationale, `code/README.md` = per-file function reference).

## Commands

```r
# Build the whole pipeline as the workflowr report (runs steps 1..7 in order):
workflowr::wflow_build("analysis/brapi_selection_pipeline.Rmd")

# Run/iterate on steps manually (sourcing order matters; see the Rmd's source chunk):
source(here::here("code", "config.R")); source(here::here("code", "01_connect.R")); ...
```

```bash
# Tests (testthat). Tiers 1+2 are offline; tier 3 hits the live server.
Rscript tests/run_tests.R
RUN_LIVE_TESTS=true Rscript tests/run_tests.R
```

There is no linter or formal build system. Beyond the test suite, verification is
done by running steps against the live server and sanity-checking outputs (counts,
matrix dims, EM convergence). The offline tests mock the BrAPI connection and use
synthetic data, so they do not exercise the real server — keep `tests/` in sync when
you change a step's signature or behavior (see `tests/README.md`). When changing a
step, re-run it with `refresh = TRUE` (or delete its `data/*.rds` cache) since every
step caches its result.

## Architecture (the big picture)

Linear, single-responsibility steps in `code/`, each consuming the previous step's
object and caching under `data/`:

```
config.R          parameters sourced by every step
01_connect        conn     T3/Oat BrAPI connection (anonymous; login fallback)
02_find_trials    trials   locations within radius OR explicit STUDY_NAMES -> studies
03_get_phenotypes pheno    /observationunits -> long phenotypes + field design + accessions
04_find_genotyping geno    coverage ranking -> per-protocol GRMs -> EM-combined GRM (+pedigree)
05_stage1_blues   blues    per-trial BLUEs (lme4), weight = 1/SE^2
06_stage2_*       gebv     BGLR genomic prediction (RKHS on the combined GRM) -> GEBVs
07_select         parents  selection index + relatedness flag -> output/selected_parents.csv
```

Shared: `grm_utils.R` (`.Gmatrix`, `std_grm`) and `em_covariance_combiner.R`
(`EMCovarianceCombiner`, copied verbatim from `../T3Predictathon2026`).

Things that require reading several files to grasp:
- **Two keying spaces.** Genotyping *coverage* is keyed by `germplasmDbId`;
  everything marker-/relationship-/BLUE-side is keyed by `germplasmName`. Step 04
  holds the dbid↔name map that bridges them — preserve this when editing 03/04/06.
- **`geno` shape drives Stage 2.** `find_and_get_genotypes()` returns a combined GRM
  `G` always, but a raw `markers` matrix **only** for the single-protocol case.
  `stage2_gblup()` uses RKHS on `G` by default and auto-downgrades marker-effect
  models (BRR/BayesB) to RKHS when `markers` is NULL.
- **`GENO_PROTOCOL_ID = NULL`** means "combine ALL covering protocols via the
  Wishart-EM combiner" (bridges = accessions on ≥2 platforms); an id means a single
  protocol. This single flag changes the whole genotyping path.

## Conventions (enforced; also in ~/.claude/CLAUDE.md)

- tidyverse + native pipe `|>`; `here::here()` for all paths (never `setwd()`/absolute);
  `here::i_am(...)` at the top of each script; non-tidyverse functions called as
  `pkg::fn()` rather than attaching the package.
- **Exception:** `httr` MUST be attached with `library(httr)` (done in
  `01_connect.R`) — BrAPI.R calls `timeout()` unqualified and needs it on the search
  path.

## Live-server gotchas (don't rediscover these)

- T3/Oat does **not** implement BrAPI genotyping endpoints; markers come from
  Breedbase VCF downloads. Prefer **archived** VCFs (`conn$vcf_archived`); on-the-fly
  `conn$vcf` stalls.
- `/studies` ignores a server-side `locationDbId` filter → studies are filtered
  client-side.
- GBS VCFs are large (hundreds of MB to multi-GB); large files are thinned to
  `TARGET_DENSITY` markers and undownloadable ones are skipped with a warning.
- QC drops high-missing **markers before accessions** (projects can mix genome
  annotations, leaving every accession ~50% missing otherwise).

## Cross-repo dependencies

- `../T3Predictathon2026/scripts/Analysis_Claude/` — source of the EM combiner and
  the VCF-thinning approach.
- `../BrAPI_pedigree_relmat/output/<id>/` — precomputed pedigree relationship
  matrices read via `PEDIGREE_DIR` (currently only small Oat groups exist; see
  `BACKGROUND.md`).
