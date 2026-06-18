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
02_find_trials      trials      locations within radius (or STUDY_NAMES) -> studies
03_get_phenotypes   pheno       observations + field design + accession list
04_find_genotyping  geno        coverage ranking -> per-protocol GRMs -> combined GRM
05_stage1_blues     blues       per-trial adjusted means (BLUEs) + weights
06_stage2_..._pred  gebv        BGLR genomic prediction -> GEBVs for all candidates
07_select           parents     selection index + relatedness flag -> ranked CSV

analysis/brapi_selection_pipeline.Rmd   orchestrates 01..07 into a workflowr report
```

Shared helpers:
- `code/grm_utils.R` — `.Gmatrix()` (VanRaden GRM) and `std_grm()` (mean-diagonal-1
  standardization), used by steps 4 and 6.
- `code/em_covariance_combiner.R` — the Wishart-EM combiner
  (`EMCovarianceCombiner`), copied verbatim from the sibling `T3Predictathon2026`
  project so this repo is self-contained.

## Module responsibilities

- **`02_find_trials.R`** — `find_ny_trials()`. Pulls `/locations` (GeoJSON
  coordinates → lon/lat, great-circle distance) and `/studies`. Two modes: explicit
  `STUDY_NAMES`, or geographic (`STUDY_TYPES` within `RADIUS_KM` over `YEARS`).
  Filtering is client-side because `/studies` ignores a server-side location filter.
- **`03_get_phenotypes.R`** — `get_phenotypes()`. One call per study to
  `/observationunits?includeObservations=true` yields observation values, germplasm,
  and field design (rep/block from `observationLevelRelationships`, row/col from
  position, entry type) together. Returns a long phenotype table, a design table, and
  the unique accession list (germplasmDbId + germplasmName).
- **`04_find_genotyping.R`** — `find_and_get_genotypes()`. Ranks covering protocols
  and projects (`conn$filter_geno_protocols` / `_projects`); downloads archived VCFs
  (`conn$vcf_archived`, cached, skip-on-failure); thins large VCFs to
  `TARGET_DENSITY`; builds one standardized GRM per protocol; EM-combines them (with
  optional pedigree partials) into one relationship matrix `G`. Returns coverage
  tables, the protocols used, `G`, and (single-protocol only) a raw marker matrix.
- **`05_stage1_blues.R`** — `stage1_blues()`. Per study × trait `lme4` mixed model
  (genotype fixed, design random) → BLUEs with weight = 1/SE²; falls back to `lm` /
  plot means when replication is absent.
- **`06_stage2_genomic_prediction.R`** — `stage2_gblup()`, `cv_accuracy()`. Removes
  environment main effects, weight-averages BLUEs per genotype, and fits BGLR (RKHS
  on `G` by default) to produce GEBVs for **every** genotyped candidate, including
  those never phenotyped.
- **`07_select.R`** — `select_parents()`. Standardizes per-trait GEBVs into a
  selection index, ranks candidates, and flags picks closely related (via `G`) to a
  better-ranked pick.

## Key data objects

- **`trials`** — tibble: `studyDbId, studyName, studyType, locationDbId, year,
  trialName, locationName, longitude, latitude, distance_km`.
- **`pheno`** — list: `pheno` (long: study, germplasm, design cols, trait, value),
  `design`, `accessions` (germplasmDbId + germplasmName).
- **`geno`** — list: `protocols`, `projects` (coverage tables), `protocol_ids`,
  `G` (combined GRM, by germplasmName), `markers` (single-protocol case only).
- **`blues`** — tibble: `trait, studyDbId, genotype, BLUE, SE, weight`.
- **`gebv`** — tibble: `trait, genotype, GEBV, phenotyped`.

Everything marker-/relationship-side is keyed by **germplasmName** (the VCF sample
labels and the BLUE genotype labels), while coverage ranking uses **germplasmDbId**;
`04` holds the dbid↔name mapping that bridges the two.

## Conventions

- **workflowr layout** + global R style: tidyverse, native pipe `|>`, `here::here()`
  for all paths, `here::i_am()` at the top of each script, non-tidyverse functions
  called as `pkg::fn()`.
- **Caching** is per-step under `data/`; results under `output/`. Steps take a
  `refresh` flag.
- **Self-contained**: external algorithms (EM combiner) are copied in, not sourced
  across repos; external *data* (pedigree matrices) is read from a configured path.
