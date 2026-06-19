# BrAPI pipeline for selection

A pipeline that pulls oat data from **T3/Oat** (`oat.triticeaetoolbox.org`) over
**BrAPI**, runs a two-stage genomic-prediction analysis, and produces a ranked
list of accessions to use as **parents for crossing**, targeted at the **New York**
growing environment.

- **How to use it** → this file.
- **What it is / how it's structured** → [DESIGN.md](DESIGN.md).
- **Why it works the way it does (theory + decisions)** → [BACKGROUND.md](BACKGROUND.md).

This is a [workflowr](https://github.com/workflowr/workflowr) project: reusable
functions live in `code/`, the runnable report in `analysis/`, cached pulls in
`data/`, and results in `output/`.

## 1. Install dependencies

```r
remotes::install_github("TriticeaeToolbox/BrAPI.R")   # BrAPI wrapper
install.packages(c("tidyverse", "here", "geosphere",  # data + geography
                   "lme4", "BGLR", "vcfR",            # models + markers
                   "workflowr"))
```

R ≥ 4.5 is assumed. `httr` is attached automatically by `code/01_connect.R`
(BrAPI.R calls `timeout()` unqualified and needs it on the search path).

## 2. Configure the run — `code/config.R`

Everything is parameterised in one file. The settings you are most likely to change:

| Setting | Meaning |
|---------|---------|
| `DB_NAME` | BrAPI connection name (default `"T3/Oat"`). |
| `CENTER_LAT`, `CENTER_LON`, `RADIUS_KM` | Geographic search: trials within this radius of the center point (default Ithaca, NY / 500 km). |
| `TRAINING_TRIALS` | **If not NULL**, train on exactly these trials (by `studyName`) and ignore the radius search. NULL = use the radius search. Their accessions are always predicted. |
| `TEST_TRIALS` | Trials (by `studyName`) whose accessions you also want predicted; their phenotypes are **not** used for training. `NULL` = none. |
| `TEST_ACCESSIONS` | `germplasmName` vector to also predict (as many as possible). `NULL` = none. |
| `STUDY_TYPES` | Vector of allowable study types analyzed together (radius search). |
| `YEARS` | Seasons kept by the radius search. |
| `TRAIT_PATTERNS` | Case-insensitive substrings matched against trait names (e.g. `"yield"`); empty = all traits. |
| `GENO_PROTOCOL_ID` | `NULL` = use **all** covering genotyping protocols and EM-combine them; an id = a single protocol (also yields a raw marker matrix). |
| `TARGET_DENSITY` | Marker-thinning target for large VCFs (default 10000). |
| `PEDIGREE_DIR` | Folder of precomputed pedigree relationship matrices to stitch in; `NULL` disables. |
| `GRM_DF_MEAN`, `GRM_DF_STDEV` | EM degrees-of-freedom for marker GRMs: their effective-sample-size measure is re-centered on `GRM_DF_MEAN` with spread capped at `GRM_DF_STDEV`. |
| `PEDIGREE_DF` | Fixed EM degrees-of-freedom (weight) for the pedigree matrix; `GRM_DF_MEAN` vs this sets marker-vs-pedigree trust. |
| `BGLR_MODEL`, `BGLR_NITER`, `BGLR_BURNIN`, `SEED` | Genomic-prediction (BGLR) settings. |
| `MAX_MISSING`, `MIN_MAF` | Marker QC thresholds. |

### Choosing training trials
- **Geographic (default):** leave `TRAINING_TRIALS <- NULL`; trials are found within
  `RADIUS_KM` of the center point, of any type in `STUDY_TYPES`, grown in `YEARS`.
- **Explicit:** set `TRAINING_TRIALS <- c("Trial_A", "Trial_B", ...)`; exactly those
  trials train the model (type/year/radius ignored). Names not found on the server
  are reported with a warning.

### Choosing who to predict
By default the pipeline predicts only the `TRAINING_TRIALS` accessions. Widen the set
with `TEST_TRIALS` (predict accessions from those trials; their phenotypes are held out
of training) and/or `TEST_ACCESSIONS` (predict these `germplasmName`s). The prediction
set is the union; test accessions are predicted only if they end up in the relationship
matrix (genotype or pedigree), while training accessions are **always** predicted —
force-injected with a prior-only diagonal if they have neither genotypes nor pedigree.
A trial named in both `TRAINING_TRIALS` and `TEST_TRIALS` is treated as training.

## 3. Run it

The simplest path is to build the workflowr report, which runs every step in order:

```r
workflowr::wflow_build("analysis/brapi_selection_pipeline.Rmd")
```

Or run the steps yourself from R:

```r
library(tidyverse)
purrr::walk(c("config.R","01_connect.R","02_find_trials.R","03_get_phenotypes.R",
              "04_find_genotyping.R","05_stage1_blues.R",
              "06_stage2_genomic_prediction.R","07_select.R"),
            ~ source(here::here("code", .x)))

conn   <- connect_t3()                                  # 1. connect (public/anon)
trials <- find_ny_trials(conn)                          # 2. training + test trials
pheno  <- get_phenotypes(conn, trials$studyDbId)        # 3. phenotypes + design
sets   <- split_by_role(pheno, trials)                  #    train_pheno/train_acc/test_acc
geno   <- find_and_get_genotypes(conn, sets$train_acc,  # 4. prediction GRM
                                 sets$test_acc, TEST_ACCESSIONS)
blues  <- stage1_blues(sets$train_pheno)                # 5. per-trial BLUEs (training)
gebv   <- stage2_gblup(blues, geno)                     # 6. BGLR genomic prediction
parents<- select_parents(gebv, geno = geno)             # 7. ranked parents
```

Or just source the convenience driver, which runs steps 1–7 in order and leaves
all result objects in the global environment:

```r
source(here::here("code", "run_pipeline.R"))
```

Run it in an **interactive console** to watch the progress bars on the slow
phenotype/VCF downloads and the EM combine — `wflow_build()` renders in a
captured, non-interactive subprocess, so the bars don't display there. Steps are
cached, so set `PIPELINE_REFRESH <- TRUE` before sourcing to force fresh runs and
actually see the bars.

## 4. Outputs

- `output/selected_parents.csv` — every genotyped candidate ranked by the selection
  index, with per-trait GEBVs, a `selected` flag, and a `redundant_with` column that
  marks picks closely related to a better-ranked pick.
- The rendered report (`docs/brapi_selection_pipeline.html`) shows the trial map,
  trait summaries, genotyping coverage, Stage-1 diagnostics, cross-validated
  prediction accuracy, and the final parent list.

## 5. Caching & re-runs

Each data/compute step caches under `data/` (`ny_trials.rds`, `phenotypes.rds`,
`genotypes.rds`, `gebv.rds`) and downloaded VCFs under `data/vcf_cache/`. Re-runs
reuse the cache; pass `refresh = TRUE` to a step (or delete its cache file) to
recompute. **Delete the relevant cache after changing `config.R`**, otherwise a
step will return its stale cached result.

## 6. Notes & gotchas

- **Access** is anonymous/public by default. Set `login = TRUE` in `connect_t3()`
  (reads `T3_USER` / `T3_PASS` env vars) only if a call returns HTTP 401.
- **First run is slow**: pulling all NY-region trials' phenotypes takes ~30 min
  (network-bound), then it is cached.
- **GBS marker files are large** (hundreds of MB to multiple GB). Large files are
  thinned to ~`TARGET_DENSITY` markers; files that fail to download (e.g. multi-GB
  diversity panels) are skipped with a warning. For a fast run, set
  `GENO_PROTOCOL_ID` to a small protocol such as an Oat 3K SNP array.
- **Pedigree stitch** for T3/Oat currently engages rarely — see
  [BACKGROUND.md](BACKGROUND.md#pedigree-stitch).

## 7. Tests

A tiered [testthat](https://testthat.r-lib.org/) suite lives in `tests/`:

```bash
Rscript tests/run_tests.R                       # tiers 1 + 2 (offline, ~2 s)
RUN_LIVE_TESTS=true Rscript tests/run_tests.R   # also tier 3 (live network)
```

Tier 1 = pure offline unit tests; Tier 2 = real BGLR/EM on synthetic data with a
mocked connection; Tier 3 = live T3/Oat calls (opt-in). See `tests/README.md`.

---

See `code/README.md` for a per-file function reference.
