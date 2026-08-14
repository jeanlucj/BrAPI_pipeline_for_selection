# BrAPI pipeline for selection

A pipeline that pulls oat data from **T3/Oat** (`oat.triticeaetoolbox.org`) over
**BrAPI**, runs a two-stage genomic-prediction analysis, and produces a ranked
list of accessions to use as **parents for crossing**, targeted at the **New York**
growing environment.

- **How to use it** → this file.
- **What it is / how it's structured** → [DESIGN.md](DESIGN.md).
- **Why it works the way it does (theory + decisions)** → [BACKGROUND.md](BACKGROUND.md).
- **How to validate it, module by module** → [EVALUATION.md](EVALUATION.md), with
  [EVALUATION_CHECKLIST.md](EVALUATION_CHECKLIST.md) to tick off.

This is a [workflowr](https://github.com/workflowr/workflowr) project: reusable
functions live in `code/`, the runnable report in `analysis/`, cached pulls in
`data/`, and results in `output/`.

## 1. Install dependencies

```r
remotes::install_github("TriticeaeToolbox/BrAPI.R")   # BrAPI wrapper
install.packages(c("tidyverse", "here", "geosphere",  # data + geography
                   "lme4", "BGLR", "vcfR", "glmnet",  # models + markers
                   "cli", "workflowr"))               # progress reporting + site
```

R ≥ 4.5 is assumed. `httr` is attached automatically by `code/01_connect.R`
(BrAPI.R calls `timeout()` unqualified and needs it on the search path).

## 2. Set up your T3 credentials — `.Renviron`

T3/Oat **no longer allows anonymous access**, so the pipeline logs in on every run.
Copy the template and fill in your T3 account:

```bash
cp .Renviron.example .Renviron    # then edit: T3_USERNAME=..., T3_PASSWORD=...
```

`.Renviron` is gitignored — never commit it. `connect_t3()` reads it itself (via
`readRenviron(here::here(".Renviron"))`), so it works no matter which directory R
was started in, including inside `wflow_build()`. Missing credentials, or ones the
server rejects, stop the run at step 1 with an explicit message rather than
degrading to empty results.

## 3. Configure the run — `code/config.R` and `data/config/`

Everything is parameterised in `code/config.R`, except the four long selection lists,
which live as text files in `data/config/` (see below). The settings you are most
likely to change:

| Setting | Meaning |
|---------|---------|
| `DB_NAME` | BrAPI connection name (default `"T3/Oat"`). |
| `CENTER_LAT`, `CENTER_LON`, `RADIUS_KM` | Geographic search: trials within this radius of the center point (default Ithaca, NY / 500 km). |
| `TRAINING_TRIALS` | 📄 **If not NULL**, train on exactly these trials (by `studyName`) and ignore the radius search. NULL = use the radius search. Their accessions are always predicted. |
| `TEST_TRIALS` | 📄 Trials (by `studyName`) whose accessions you also want predicted; their phenotypes are **not** used for training. `NULL` = none. |
| `TEST_ACCESSIONS` | 📄 `germplasmName` vector to also predict (as many as possible). `NULL` = none. |
| `STUDY_TYPES` | Vector of allowable study types analyzed together (radius search). |
| `YEARS` | Seasons kept by the radius search. |
| `TRAIT_NAMES` | 📄 Exact `observationVariableName` strings to keep; empty = all traits. |
| `TRAIT_WEIGHTS` | 📄 Selection-index weight per trait (sign = direction, magnitude = importance/scale), named by `TRAIT_NAMES`. |
| `TRAIT_SHORT_NAMES` | 📄 Optional short labels (named by `TRAIT_NAMES`) used for the `breeders_output.csv` columns; `NULL` = full names. |
| `GENO_PROTOCOL_ID` | `NULL` = use **all** covering genotyping protocols and EM-combine them; an id = a single protocol (also yields a raw marker matrix). |
| `TARGET_DENSITY` | Marker-thinning target for large VCFs (default 10000). |
| `PEDIGREE_DIR` | Folder of precomputed pedigree relationship matrices to stitch in; `NULL` disables. |
| `GRM_DF_MEAN`, `GRM_DF_STDEV` | EM degrees-of-freedom for marker GRMs: their effective-sample-size measure is re-centered on `GRM_DF_MEAN` with spread capped at `GRM_DF_STDEV`. |
| `PEDIGREE_DF` | Fixed EM degrees-of-freedom (weight) for the pedigree matrix; `GRM_DF_MEAN` vs this sets marker-vs-pedigree trust. |
| `SE_FLOOR_FRAC` | Stage-1: floor on each BLUE's SE, as a fraction of the trial × trait response SD, so a (near-)perfect fit can't drive `1/SE²` to ∞ (default 0.01; 0 disables). Constant-response trials are dropped outright. |
| `MIXED_MODEL_ENGINE`, `BGLR_NITER`, `BGLR_BURNIN`, `SEED` | Stage-2 kernel-GBLUP settings: engine (`"bglr"`/`"sommer"`) and BGLR MCMC controls. |
| `MAX_MISSING`, `MIN_MAF` | Marker QC thresholds. |

### 📄 Selection lists — `data/config/*.txt`

The settings marked 📄 above are read from plain text files rather than typed into
`config.R`, because real lists run to hundreds of lines: retarget a run by swapping a
file, not by editing source. `data/config/` is the one subfolder of `data/` that is
**tracked in git**, so the lists that define a run are versioned with the analysis.

```r
TRAINING_TRIALS <- config_lines("Trial_Sel2026_Intersect20NoYld.txt")
TEST_ACCESSIONS <- config_lines("Acc_Sel2026.txt")
.traits <- config_traits("Trait_Sel2026.txt")
```

**Format:** one value per line. Blank lines and lines starting with `#` are ignored,
surrounding whitespace is trimmed, and duplicates are dropped with a message.
Everything else is used **verbatim** — a `germplasmName` like `APPLER|CIAV2680` keeps
its pipe. A file name is resolved inside `data/config/`; an absolute path also works.
A **missing file is an error**, not an empty list (a typo in `TRAINING_TRIALS` would
otherwise silently fall back to the radius search).

**The trait file** takes 1 to 3 **tab-separated** columns — name, index weight, short
label:

```
Grain yield - g/m2|CO_350:0000260	0.5	Yield g/m2
Freeze damage severity - 0-9 Rating|CO_350:0005001	-5	Freeze_damage 0-9
```

A column you use must be filled on *every* line (a half-filled weight column is an
error, since it would silently zero part of the selection index). Omit the weight
column for a plain all-traits prediction dump with no index; omit the short-name
column to use the full trait names.

To go back to an inline vector, just assign one — `TRAINING_TRIALS <- c("Trial_A",
"Trial_B")` — or `NULL` to disable the setting. `config.R` keeps both forms as
commented alternatives.

### Choosing training trials
- **Geographic (default):** leave `TRAINING_TRIALS <- NULL`; trials are found within
  `RADIUS_KM` of the center point, of any type in `STUDY_TYPES`, grown in `YEARS`.
- **Explicit (default):** `TRAINING_TRIALS <- config_lines("<file>.txt")` (or an
  inline `c("Trial_A", "Trial_B", ...)`); exactly those trials train the model
  (type/year/radius ignored). Names not found on the server are reported with a
  warning.

### Choosing who to predict
By default the pipeline predicts only the `TRAINING_TRIALS` accessions. Widen the set
with `TEST_TRIALS` (predict accessions from those trials; their phenotypes are held out
of training) and/or `TEST_ACCESSIONS` (predict these `germplasmName`s). The prediction
set is the union; test accessions are predicted only if they end up in the relationship
matrix (genotype or pedigree), while training accessions are **always** predicted —
force-injected with a prior-only diagonal if they have neither genotypes nor pedigree.
A trial named in both `TRAINING_TRIALS` and `TEST_TRIALS` is treated as training.

## 4. Run it

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

conn   <- connect_t3()                                  # 1. connect (logs in)
trials <- find_ny_trials(conn)                          # 2. training + test trials
pheno  <- get_phenotypes(conn, trials$studyDbId)        # 3. phenotypes + design
sets   <- split_by_role(pheno, trials)                  #    train_pheno/train_acc/test_acc
geno   <- find_and_get_genotypes(conn, sets$train_acc,  # 4. prediction GRM
                                 sets$test_acc, TEST_ACCESSIONS)
blues  <- stage1_blues(sets$train_pheno)                # 5. per-trial BLUEs (training)
gebv   <- stage2_gblup(blues, geno)                     # 6. BGLR genomic prediction
cv     <- map_dfr(unique(gebv$trait), ~ cv_accuracy(blues, geno, .x))
out    <- select_parents(gebv, cv = cv)                 # 7. breeders_output.csv
```

Or just source the convenience driver, which runs steps 1–7 in order and leaves
all result objects in the global environment:

```r
source(here::here("code", "run_pipeline.R"))
```

### Watching a run

Run the driver in an **interactive console**: each step announces itself, the slow
loops draw a progress bar with an ETA, and the run ends with a timing table.

```
── 4 Genotyping ────────────────────────────────────────────
ℹ Resolving genotyping coverage for 1 204 accessions ...
  3 covering protocol(s); using 3
⠙ protocol 66: VCF files 2/3 ■■■■■■■■■■■■■■        67% | ETA 1m
ℹ Building GRM for protocol 66 (1/3) ...
  panel: 812 relevant + 188 filler = 1000 accessions
  QC: 41 233 -> 17 004 markers, 1000 -> 998 accessions
⠹ protocol 66 (1/3): imputing markers 3412/17004 ■■■■■  20% | ETA 31m
ℹ EM-combining 4 partial covariance(s) over 1 118 accessions ...
✔ 4 Genotyping (12m 41.3s) -- protocols 66, 71, 83 | prediction GRM 812 x 812

── Run complete (24m 08s) ──────────────────────────────────
  1 Connect                   0.9s
  2 Trials                   12.4s   (cached)
  3 Phenotypes            6m 31.0s
  ...
```

Bars are drawn where the time actually goes: per study downloaded, per VCF file, per
marker imputed (the longest single loop in the pipeline), per trait×study BLUE, per
trait predicted, and per CV fold. A **cached step says so and returns immediately**,
so a fast re-run no longer looks like a hung one; set `PIPELINE_REFRESH <- TRUE`
before sourcing to force fresh runs.

Reporting is controlled by `SHOW_PROGRESS` in `config.R`, which defaults to
`interactive()` — `wflow_build()` renders in a captured, non-interactive subprocess
with no terminal to animate, and the test suite stays silent. Force it either way
with `options(brapi.progress = TRUE)` / `FALSE`.

## 5. Outputs

- `output/breeders_output.csv` — a spreadsheet-style file with two side-by-side blocks
  (separated by a blank column). **Block 1** (one row per predicted accession, sorted by
  selection index): `accession`, `In_Training` (1/0), one **un-standardized** GEBV
  column per trait, and the `Index` (= Σ `TRAIT_WEIGHTS` × GEBV). **Block 2** (one row
  per trait): short trait name, its index weight, and its cross-validated accuracy.
  With `TRAIT_WEIGHTS = NULL` there is no index: block 1 drops the `Index` column
  (rows sorted by accession) and block 2 drops the weight column — a plain
  all-traits prediction dump when no selection weights have been chosen.
- The rendered report (`docs/brapi_selection_pipeline.html`) shows the trial map,
  trait summaries, genotyping coverage, Stage-1 diagnostics, cross-validated
  prediction accuracy, and the final parent list.

## 6. Caching & re-runs

Each data/compute step caches under `data/` (`ny_trials.rds`, `phenotypes.rds`,
`genotypes.rds`, `gebv.rds` / `gebv_sommer.rds`) and downloaded VCFs under
`data/vcf_cache/`. Re-runs reuse the cache; pass `refresh = TRUE` to a step (or
delete its cache file) to recompute. **After changing `config.R` — or any
`data/config/*.txt` list — delete the caches of the steps it affects**; most steps
return their cached result verbatim.
*Exception:* Stage 2 (`stage2_gblup`) validates its GEBV cache against the request
and regenerates on its own when the relationship kernel's candidate set changes or a
new trait is requested, so switching engine / `G` / trait set yields fresh
predictions without re-downloading the upstream (expensive) trial and genotype data.

## 7. Notes & gotchas

- **Access requires a login** (section 2): `connect_t3()` authenticates with
  `T3_USERNAME` / `T3_PASSWORD` from `.Renviron`. `connect_t3(login = FALSE)` is for
  offline/mock use only — anonymous calls to T3/Oat come back empty, not with a 401.
  Note BrAPI's `conn$login()` neither errors on a wrong password (it only warns and
  leaves `conn$auth_token` empty) nor on empty arguments (it *prompts*, which would
  hang `wflow_build()`); `t3_login()` handles both.
- **First run is slow**: pulling all NY-region trials' phenotypes takes ~30 min
  (network-bound), then it is cached.
- **GBS marker files are large** (hundreds of MB to multiple GB). Large files are
  thinned to ~`TARGET_DENSITY` markers; files that fail to download (e.g. multi-GB
  diversity panels) are skipped with a warning. For a fast run, set
  `GENO_PROTOCOL_ID` to a small protocol such as an Oat 3K SNP array.
- **Pedigree stitch** for T3/Oat currently engages rarely — see
  [BACKGROUND.md](BACKGROUND.md#pedigree-stitch).

## 8. Tests

A tiered [testthat](https://testthat.r-lib.org/) suite lives in `tests/`:

```bash
Rscript tests/run_tests.R                       # tiers 1 + 2 (offline, ~2 s)
RUN_LIVE_TESTS=true Rscript tests/run_tests.R   # also tier 3 (live network)
```

Tier 1 = pure offline unit tests; Tier 2 = real BGLR/EM on synthetic data with a
mocked connection; Tier 3 = live T3/Oat calls (opt-in). See `tests/README.md`.

---

See `code/README.md` for a per-file function reference.
