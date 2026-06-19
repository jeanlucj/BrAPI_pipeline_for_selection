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
02_find_trials    trials   radius OR explicit TRAINING_TRIALS, + TEST_TRIALS, tagged by `role`
03_get_phenotypes pheno    /observationunits -> long phenotypes + field design + accessions; split_by_role()
04_find_genotyping geno    coverage(train+test) -> per-protocol GRMs -> EM-combined GRM (+pedigree), SUBSET to targets
05_stage1_blues   blues    per-trial BLUEs (lme4) on TRAINING trials, weight = 1/SE^2
06_stage2_*       gebv     BGLR genomic prediction (RKHS on the subset GRM) -> GEBVs for the targets
07_select         parents  selection index + relatedness flag -> output/selected_parents.csv
```

Shared: `grm_utils.R` (`.Gmatrix`, `std_grm`) and `em_covariance_combiner.R`
(`EMCovarianceCombiner`, copied verbatim from `../T3Predictathon2026`).

Things that require reading several files to grasp:
- **Two keying spaces.** Genotyping *coverage* is keyed by `germplasmDbId`;
  everything marker-/relationship-/BLUE-side is keyed by `germplasmName`. Step 04
  holds the dbid↔name map that bridges them — preserve this when editing 03/04/06.
  Pedigree group CSVs are keyed by `germplasmDbId`, so step 04 widens that map with
  the pedigree `germplasm_cache_<id>.rds` to name *non-phenotyped* pedigree
  accessions (needed for pedigree-bridge detection).
- **Prediction targets drive `G`'s size.** `find_and_get_genotypes(train_acc,
  test_acc, test_names)` builds the full EM-combined GRM (with bridges/pedigree) then
  **subsets it to training ∪ predictable-test** — bridges only inform the combine.
  Training accessions are force-kept: any with no genotype/pedigree are **injected**
  (diag = mean diag, off-diag 0). Test accessions absent from the GRM are silently
  dropped. Stage 1 trains on `split_by_role()`'s `train_pheno` only (test-trial
  phenotypes held out; a trial in both lists counts as training).
- **`geno` shape drives Stage 2.** `find_and_get_genotypes()` returns a combined GRM
  `G` always, but a raw `markers` matrix **only** for the single-protocol, no-pedigree,
  no-injection case. `stage2_gblup()` uses RKHS on `G` by default and auto-downgrades
  marker-effect models (BRR/BayesB) to RKHS when `markers` is NULL.
- **`GENO_PROTOCOL_ID = NULL`** means "combine ALL covering protocols via the
  Wishart-EM combiner" (bridges = accessions in ≥2 partials, counting each platform
  AND each pedigree group, so a single-platform accession that also sits in a
  pedigree group is kept); an id means a single protocol. This single flag changes
  the whole genotyping path.

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
  matrices read via `PEDIGREE_DIR`; step 04 also reads its companion
  `<id>_pedigree_groups.rds` (cheap per-group membership) and the sibling
  `../output/germplasm_cache_<id>.rds` (dbid→name for non-phenotyped accessions, for
  bridge detection). See `BACKGROUND.md`. **Step 04 depends on this output contract;
  if the sibling project changes any of it, update `04_find_genotyping.R`:**
  - group files named `<id>_group<N>.csv`, columns `germplasmDbId_i,
    germplasmDbId_j, relationship` (upper triangle + diagonal);
  - `<id>_pedigree_groups.rds` is a list whose element **named** `group<N>` (matching
    the CSV's `<N>`) holds that group's germplasmDbId vector;
  - `germplasm_cache_<id>.rds` is a **list of BrAPI germplasm records** (each with
    `$germplasmDbId`, `$germplasmName`), located in `dirname(PEDIGREE_DIR)`.
  Missing companions degrade gracefully (phenotyped-only pedigree, no pedigree
  bridges); a present-but-malformed companion fails loudly via
  `.pedigree_contract_error()` with a message pointing back here, rather than
  silently building a wrong/empty GRM. The offline tests pin the shape step 04
  *expects* (with synthetic fixtures) but don't observe the sibling repo's real
  output, so a genuine drift there is caught at runtime by those guards.
