# DESIGN

What this pipeline is and how it is put together. For usage see
[README.md](README.md); for the reasoning behind the methods see
[BACKGROUND.md](BACKGROUND.md).

## Purpose

Turn public T3/Oat data into a ranked list of parents for crossing, optimized for
the New York environment, with no manual data wrangling: discover the relevant
phenotyping trials, pull their phenotypes and the marker data for their accessions,
fit a two-stage genomic-prediction model, and rank candidates.

## Pipeline shape

A linear sequence of small, single-responsibility steps. Each consumes the previous
step's object and caches its own result under `data/`, so any step can be re-run in
isolation.

```
config.R ── parameters (sourced by every step)

01_connect          conn        connect to T3/Oat over BrAPI (anonymous)
02_find_trials      trials      training (radius or TRAINING_TRIALS) + TEST_TRIALS, tagged by role
03_get_phenotypes   pheno       observations + field design + accession list (split_by_role)
04_find_genotyping  geno        coverage -> per-protocol GRMs -> combined GRM, subset to targets
05_stage1_blues     blues       per-trial adjusted means (BLUEs) + weights (training trials only)
06_stage2_..._pred  gebv        BGLR genomic prediction -> GEBVs for the target accessions
07_select           parents     selection index + relatedness flag -> ranked CSV

analysis/brapi_selection_pipeline.Rmd   orchestrates 01..07 into a workflowr report
```

Shared helpers:
- `code/grm_utils.R` — `.Gmatrix()` (VanRaden GRM), `std_grm()` (mean-diagonal-1
  standardization), and `.effective_n()` (Galwey 2009 effective sample size, used to
  set each marker GRM's EM degrees of freedom), used by steps 4 and 6.
- `code/em_covariance_combiner.R` — the Wishart-EM combiner
  (`EMCovarianceCombiner`), copied verbatim from the sibling `T3Predictathon2026`
  project so this repo is self-contained.

## Module responsibilities

- **`02_find_trials.R`** — `find_ny_trials()`. Pulls `/locations` (GeoJSON
  coordinates → lon/lat, great-circle distance) and `/studies`. Training trials come
  from explicit `TRAINING_TRIALS` or the geographic search (`STUDY_TYPES` within
  `RADIUS_KM` over `YEARS`); `TEST_TRIALS` adds test trials. Returns one tibble with a
  `role` column (`"training"`/`"test"`; training wins on overlap). Filtering is
  client-side because `/studies` ignores a server-side location filter.
- **`03_get_phenotypes.R`** — `get_phenotypes()` (one call per study to
  `/observationunits?includeObservations=true` → long phenotypes, design, accession
  list) plus `split_by_role(pheno, trials)`, which partitions into `train_pheno`
  (training trials, feeds Stage 1) and `train_acc` / `test_acc` (distinct accessions
  per role).
- **`04_find_genotyping.R`** — `find_and_get_genotypes(conn, train_accessions,
  test_accessions, test_names)`. Ranks covering protocols/projects over **training +
  test** accessions (`TEST_ACCESSIONS` dbIds resolved via the germplasm cache);
  downloads archived VCFs (`conn$vcf_archived`, cached, skip-on-failure); thins large
  VCFs to `TARGET_DENSITY`; builds one standardized GRM per protocol; EM-combines them
  (with optional pedigree partials). The combined `G` is then **subset to the
  prediction targets** (training ∪ predictable-test) — bridges/extra accessions inform
  the combine but are dropped — and any training accession absent from `G` (no
  genotype, no pedigree) is **injected** (diagonal = mean diagonal, off-diagonals 0).
  Returns coverage tables, the protocols used, the subset `G`, and a raw marker matrix
  only for the single-protocol/no-pedigree/no-injection case. Bridges (used during the
  combine) are accessions in ≥2 partials, counting each protocol **and** each pedigree
  group; pedigree group membership/names come from the sibling project's
  `<id>_pedigree_groups.rds` and `germplasm_cache_<id>.rds`.
- **`05_stage1_blues.R`** — `stage1_blues()`. Per study × trait `lme4` mixed model
  (genotype fixed, design random) → BLUEs with weight = 1/SE²; falls back to `lm` /
  plot means when replication is absent.
- **`06_stage2_genomic_prediction.R`** — `stage2_gblup()`, `cv_accuracy()`. Removes
  environment main effects, weight-averages BLUEs per genotype, and fits BGLR (RKHS
  on `G` by default) to produce GEBVs for **every accession in `G`** (the prediction
  targets); target accessions without a training BLUE receive predicted-only GEBVs.
- **`07_select.R`** — `select_parents()`. Standardizes per-trait GEBVs into a
  selection index, ranks candidates, and flags picks closely related (via `G`) to a
  better-ranked pick.

## Key data objects

- **`trials`** — tibble: `studyDbId, studyName, studyType, locationDbId, year,
  trialName, locationName, longitude, latitude, distance_km, role` (training/test).
- **`pheno`** — list: `pheno` (long: study, germplasm, design cols, trait, value),
  `design`, `accessions` (germplasmDbId + germplasmName). `split_by_role(pheno,
  trials)` → `train_pheno`, `train_acc`, `test_acc`.
- **`geno`** — list: `protocols`, `projects` (coverage tables), `protocol_ids`,
  `G` (prediction GRM by germplasmName, subset to training ∪ predictable-test),
  `markers` (single-protocol/no-injection case only).
- **`blues`** — tibble: `trait, studyDbId, genotype, BLUE, SE, weight`.
- **`gebv`** — tibble: `trait, genotype, GEBV, phenotyped`.

Everything marker-/relationship-side is keyed by **germplasmName** (the VCF sample
labels and the BLUE genotype labels), while coverage ranking and the pedigree group
CSVs use **germplasmDbId**; `04` holds the dbid↔name mapping that bridges the two,
widened with the pedigree `germplasm_cache_<id>.rds` to cover non-phenotyped
accessions.

## Conventions

- **workflowr layout** + global R style: tidyverse, native pipe `|>`, `here::here()`
  for all paths, `here::i_am()` at the top of each script, non-tidyverse functions
  called as `pkg::fn()`.
- **Caching** is per-step under `data/`; results under `output/`. Steps take a
  `refresh` flag.
- **Self-contained**: external algorithms (EM combiner) are copied in, not sourced
  across repos; external *data* (pedigree matrices) is read from a configured path.
