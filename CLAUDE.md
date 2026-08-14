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
config.R          parameters sourced by every step; the long selection lists
                  (TRAINING_TRIALS, TEST_TRIALS, TEST_ACCESSIONS, TRAIT_NAMES +
                  weights/short names) come from data/config/*.txt
01_connect        conn     T3/Oat BrAPI connection (REQUIRES login; T3_USERNAME /
                           T3_PASSWORD from a gitignored .Renviron)
02_find_trials    trials   radius OR explicit TRAINING_TRIALS, + TEST_TRIALS, tagged by `role`
03_get_phenotypes pheno    /observationunits -> long phenotypes + field design + accessions; split_by_role()
04_find_genotyping geno    coverage(train+test) -> per-protocol GRMs -> EM-combined GRM (+pedigree), SUBSET to targets
05_stage1_blues   blues    per-trial BLUEs (lme4) on TRAINING trials, weight = 1/SE^2
06_stage2_*       gebv     BGLR genomic prediction (RKHS on the subset GRM) -> GEBVs for the targets
07_select         parents  selection index + relatedness flag -> output/selected_parents.csv
```

Shared: `grm_utils.R` (`.Gmatrix`, `std_grm`), `em_covariance_combiner.R`
(`EMCovarianceCombiner`, copied verbatim from `../T3Predictathon2026`), and
`progress.R` (console reporting).

Things that require reading several files to grasp:
- **All console reporting goes through `code/progress.R`** (`say`/`note`/`note_cache`
  for status, `step_start`/`step_done`/`print_timings` for per-step banners+timing,
  `pb_start`/`pb_tick`/`pb_done` and `pb_wrap` for bars), sourced via `config.R` and
  gated by `SHOW_PROGRESS` (default `interactive()`; `options(brapi.progress=)`
  overrides). So: bars are instrumented where the time actually goes (per study, per
  VCF, per imputed marker, per trait×study fit, per trait, per CV fold), cached steps
  announce themselves instead of returning silently, and everything is silent in
  batch/tests/`wflow_build()`. Two cli facts the wrappers exist to handle: a bar dies
  when the frame that created it returns (hence `.envir = parent.frame()` in
  `pb_start`), and cli shows only the innermost bar, so nested bars' labels must carry
  their outer context. `em_covariance_combiner.R` keeps its own `txtProgressBar` and is
  timed from the caller (`time_it`) rather than edited — it must stay verbatim.
- **Run inputs live outside the tracked source.** Credentials in `.Renviron`
  (gitignored; `.Renviron.example` is the tracked template) and the selection lists in
  `data/config/*.txt` — the one subfolder of `data/` that is *not* gitignored. Both are
  read by code, not by R startup: `connect_t3()` calls `readRenviron(here::here(".Renviron"))`
  itself (R only auto-reads `.Renviron` from its startup dir, so `wflow_build()` and
  `Rscript` from `analysis/` would otherwise see no credentials), and `config.R`'s
  `config_lines()` / `config_traits()` read the lists (one value per line; the trait
  file takes optional tab-separated weight and short-name columns, all-or-none per
  column; a missing file is an error, not an empty list).
- **Two keying spaces.** Genotyping *coverage* is keyed by `germplasmDbId`;
  everything marker-/relationship-/BLUE-side is keyed by `germplasmName`. Step 04
  holds the dbid↔name map that bridges them — preserve this when editing 03/04/06.
  Pedigree group CSVs are keyed by `germplasmDbId`, so step 04 widens that map with
  the pedigree `germplasm_cache_<id>.rds` to name *non-phenotyped* pedigree
  accessions (needed for pedigree-bridge detection).
- **VCF samples are matched by name, via a synonym map.** Coverage is resolved
  server-side by dbId, but VCF sample IDs are strings, and an accession can be
  genotyped under a *preliminary* line name later demoted to a SYNONYM — so a plain
  `germplasmName` match silently drops it. `code/synonyms.R` (`USE_SYNONYMS`) builds
  a cached alias→primary lookup (`T3_brapi_helpers::build_synonym_lookup` over
  `/search/germplasm`, sourced from the sibling repo when the package isn't
  installed) and `.vcf_to_dosage` canonicalizes every VCF sample name to the primary
  before matching/labeling. Keep this when editing the VCF→dosage path; the rest of
  the pipeline keys on primaries. NB: for the current T3/Oat targets this recovered
  ~0 extra matches — most disconnected accessions are genuinely absent from the
  archived VCFs (incl. the Oat 3K array, whose server "coverage" count far exceeds
  what the archived VCF actually contains), not hidden under a synonym.
- **Prediction targets drive `G`'s size.** `find_and_get_genotypes(train_acc,
  test_acc, test_names)` builds the full EM-combined GRM (with bridges/pedigree) then
  **subsets it to training ∪ predictable-test** — bridges only inform the combine.
  Training accessions are force-kept: any with no genotype/pedigree are **injected**
  (diag = mean diag, off-diag 0). Test accessions absent from the GRM are silently
  dropped. Stage 1 trains on `split_by_role()`'s `train_pheno` only (test-trial
  phenotypes held out; a trial in both lists counts as training).
- **Stage 2 is always kernel GBLUP.** `find_and_get_genotypes()` returns a combined
  GRM `G` always (plus a raw `markers` matrix only for the single-protocol,
  no-pedigree, no-injection case, used solely to build a kernel when `G` is absent).
  `stage2_gblup()` predicts on the relationship kernel via `MIXED_MODEL_ENGINE`
  (`"bglr"` or `"sommer"` REML GBLUP) — no marker-effect (BRR/BayesB-on-markers)
  path. sommer GEBVs cache to `gebv_sommer.rds`, BGLR to `gebv.rds`. **The "bglr"
  engine fits kernel GBLUP as ridge (BRR) on the relationship eigen-factor, NOT
  `model="RKHS"`: RKHS mishandles the 1/SE² weights and inflates GEBVs ~λ-fold —
  unweighted RKHS is fine, weighted is not (see `code/bglr_rkhs_vs_brr.R`). BGLR's
  `weights` are inverse-SDs (Var∝1/weights²), so pass `sqrt(1/SE²)`.**
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

- **Anonymous access is gone**; every run must log in. An unauthenticated call returns
  an *empty* result, not a 401, so a missing login looks like missing data. Two BrAPI
  `login()` traps `t3_login()` (`01_connect.R`) handles: with empty args it *prompts*
  (`readline`/`askpass`) and would hang `wflow_build()`; on a *wrong* password it only
  warns and returns normally, leaving `conn$auth_token` NULL — so check the token, not
  the return value.
- T3/Oat does **not** implement BrAPI genotyping endpoints; markers come from
  Breedbase VCF downloads. Prefer **archived** VCFs (`conn$vcf_archived`); on-the-fly
  `conn$vcf` stalls.
- `/studies` ignores a server-side `locationDbId` filter → studies are filtered
  client-side.
- GBS VCFs are large (hundreds of MB to multi-GB); large files are thinned to
  `TARGET_DENSITY` markers and undownloadable ones are skipped with a warning.
- Markers are keyed by **canonical CHROM_POS** (leading `chr` stripped), not the VCF
  `ID` column — T3 projects key inconsistently (SNP names vs `.`), which made identical
  markers look distinct and broke the multi-project merge. (`.vcf_to_dosage`.)
- QC drops high-missing **markers before accessions** (a protocol's projects can still
  mix genuinely different panels, leaving accessions ~50% missing otherwise).
- Each protocol's imputation + VanRaden GRM are estimated on a **panel** (its relevant
  accessions + non-relevant fillers up to `max(GRM_PANEL_MIN, n_relevant)`), then the
  GRM is subset back to the relevant accessions. Allele frequencies / imputation /
  diagonals estimated on only the few relevant accessions are unstable — diagonals
  become an artifact of who shares the panel (`.protocol_grm`; see
  `code/impute_diagnostic*.R`, BACKGROUND.md). Missing calls filled per `GRM_IMPUTE`
  (`.impute_glmnet` robust elastic-net default, or `.mean_impute`); `.impute_glmnet`
  retries each marker independently and fixes its CV folds, so one bad column can't
  abort it and the result is reproducible.

## Cross-repo dependencies

- `../T3_brapi_helpers/` — shared BrAPI helper package (also used by
  `T3Predictathon2026`). `code/synonyms.R` uses its `build_synonym_lookup` /
  `canonicalize_to_primary` / `get_synonyms_from_germplasm_names` for synonym
  canonicalization; it prefers the installed package but sources
  `../T3_brapi_helpers/R/brapi_germplasm.R` directly when the package isn't
  installed (an installed copy can lag the sibling repo). Degrades to exact-name
  matching if neither is present.
- `../T3Predictathon2026/scripts/Analysis_Claude/` — source of the EM combiner and
  the VCF-thinning approach (and `build_synonym_map.R`, the original of the synonym
  functions now generalized into `T3_brapi_helpers`).
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
