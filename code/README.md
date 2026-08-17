# Code — BrAPI → T3/Oat selection pipeline

Reusable functions for a pipeline that pulls oat data from **T3/Oat** over
**BrAPI**, runs a two-stage genomic-prediction analysis, and ranks accessions to
use as parents for crossing in the **New York** environment. The orchestrating
report is `analysis/brapi_selection_pipeline.Rmd`.

## Steps

| File | Function | Does |
|------|----------|------|
| `progress.R` | `say()`, `note()`, `note_cache()`, `step_start()`/`step_done()`, `pb_start()`/`pb_tick()`/`pb_done()`, `pb_wrap()`, `time_it()`, `print_timings()` | Console status lines, cli progress bars and per-step timings. Sourced by `config.R`, so available in every step; gated by `SHOW_PROGRESS` (default `interactive()`), overridable with `options(brapi.progress = )` |
| `cache.R` | `cache_key()`, `cache_read()`, `cache_write()` | Request-keyed caching: a step reuses its `data/*.rds` only when it was built for the same request, else it reports what changed and rebuilds. The request is recorded in a sidecar `data/<name>.key.rds`, leaving the payload file's format untouched. Sourced by `config.R` |
| `config.R` | (parameters), `config_lines()`, `config_traits()` | Center point/radius, years, traits, protocol, BGLR & QC settings. The long selection lists (trials, test accessions, traits + optional tab-separated weights/short names) are read from `data/config/*.txt` |
| `01_connect.R` | `connect_t3()`, `t3_login()` | Connect to T3/Oat and log in (`T3_USERNAME`/`T3_PASSWORD` from `.Renviron`); fails fast on missing or rejected credentials |
| `02_find_trials.R` | `find_ny_trials()` | Training trials (radius or `TRAINING_TRIALS`) + `TEST_TRIALS`, tagged by `role` |
| `03_get_phenotypes.R` | `get_phenotypes()`, `split_by_role()` | Observations + field design per study; split into training phenotypes / training & test accessions |
| `04_find_genotyping.R` | `find_and_get_genotypes()` | Rank covering protocols/projects (train+test), download/thin VCFs, EM-combine per-protocol GRMs (+ pedigree), subset to prediction targets (inject ungenotyped training). Caches in four tiers — `data/vcf_cache/` (downloads), `data/partials/*.rds` (per-protocol GRM+dosage, and the pedigree partials), `data/combined_G.rds` (the EM combine), `data/genotypes.rds` (the final subset) — and prints a per-phase timing breakdown |
| `em_covariance_combiner.R` | `EMCovarianceCombiner()` | Wishart-EM combiner (copied from T3Predictathon2026) |
| `synonyms.R` | `build_alias_lookup()`, `canonicalize_to_primary()` | Cached alias→primary synonym map (wraps `T3_brapi_helpers`) so VCF samples carried under a synonym still match (`USE_SYNONYMS`) |
| `grm_utils.R` | `.Gmatrix()`, `std_grm()`, `.effective_n()`, `.impute_glmnet()`, `.mean_impute()` | VanRaden GRM, standardization, Galwey effective-sample-size (for EM df), and missing-marker imputation (robust per-marker elastic-net, or column mean), shared by steps 4 & 6 |
| `05_stage1_blues.R` | `stage1_blues()` | Per-trial BLUEs (lme4) with weights |
| `06_stage2_genomic_prediction.R` | `stage2_gblup()`, `cv_accuracy()` | Genomic prediction (RKHS on the combined GRM) → GEBVs; `engine = "bglr"` (Bayesian, default) or `"sommer"` (REML GBLUP, cross-check) |
| `07_select.R` | `select_parents()` | Un-standardized weighted selection index → two-block `breeders_output.csv` (per-accession GEBVs+index; per-trait weight+CV accuracy) |
| `evaluation.R` | `arm_evaluation()`/`disarm_evaluation()`, `eval_groups()`, `armed_functions()`, `peek()` | Evaluation-time tooling for `EVALUATION.md`: `debug()` a named group of functions, and print one-line health summaries of the objects flowing between steps. **Not sourced by the pipeline** — load it by hand, and disarm before any unattended run |
| `run_pipeline.R` | (script) | Interactive driver: sources + runs steps 1–7 in order so the download/EM progress bars are visible (set `PIPELINE_REFRESH <- TRUE` to force fresh runs) |

Each data-pulling/compute step caches its result under `data/` (e.g.
`ny_trials.rds`, `phenotypes.rds`, `genotypes.rds`, `gebv.rds`). Pass
`refresh = TRUE` (or delete the cache) to recompute.

## Dependencies

```r
remotes::install_github("TriticeaeToolbox/BrAPI.R")   # BrAPI wrapper
install.packages(c("tidyverse", "here", "geosphere",  # data + geo
                   "lme4", "BGLR", "vcfR"))            # models + markers
```

`httr` is attached by `01_connect.R` because BrAPI.R calls `timeout()` unqualified.

## Notes / gotchas (discovered against the live T3/Oat server)

- **`/studies` ignores a server-side `locationDbId` filter**, so studies are
  pulled in full and filtered client-side.
- Phenotypes + field design come from `/observationunits?includeObservations=true`.
  Pulling all ~58 NY-region trials takes a while on first run (it is cached).
- **Genotyping marker download**: T3's BrAPI genotyping endpoints are not
  implemented, so markers come from breedbase VCF downloads. On-the-fly
  generation (`conn$vcf`) is slow/unreliable, so the pipeline prefers
  pre-generated **archived** VCFs (`conn$vcf_archived`). GBS matrices —
  especially diversity panels — can be **multiple GB**; large files are thinned
  genome-wide to ~`TARGET_DENSITY` (10000) markers before reading, and any file
  that still fails to download is skipped with a warning.
- **Synonym canonicalization** (`USE_SYNONYMS`, `synonyms.R`): an accession may be
  genotyped under a *preliminary* line name later demoted to a synonym, so a VCF can
  carry it under a name ≠ our `germplasmName`. Step 4 builds a cached alias→primary
  lookup (`T3_brapi_helpers`) and relabels VCF samples to the primary before
  matching. (Empirically near-zero gain for the current T3/Oat targets — most
  disconnected accessions are simply absent from the archived VCFs — but it is the
  correct guard and is shared with the T3Predictathon project.)
- **Multi-platform combine**: `GENO_PROTOCOL_ID = NULL` builds one standardized
  GRM per covering protocol and **EM-combines** them into a single GRM
  (`em_covariance_combiner.R`), using accessions genotyped on ≥2 partials
  (platforms and/or pedigree groups) as bridges. Each partial's EM **degrees of
  freedom** weights how much it is trusted: marker GRMs use their effective number
  of independent samples (`grm_utils.R::.effective_n`, Galwey 2009) re-centered on
  `GRM_DF_MEAN` and capped at `GRM_DF_STDEV`; pedigree partials use the fixed
  `PEDIGREE_DF`. See `BACKGROUND.md` (degrees of freedom). Set `GENO_PROTOCOL_ID` to a single
  id (e.g. an **Oat 3K** array) for a faster single-platform run that also yields a
  raw marker matrix (so marker-effect BGLR models stay available). Downloaded VCFs
  are cached under `data/vcf_cache/`.
- **Per-protocol panel estimation**: each protocol's imputation and VanRaden GRM
  are estimated on a *panel* — all of our relevant accessions (targets + bridges)
  present in the protocol, **plus non-relevant accessions** from the same protocol
  up to `max(GRM_PANEL_MIN, n_relevant)` — then the GRM is subset back to the
  relevant accessions. Estimating allele frequencies / imputation / the GRM on only
  our handful of accessions makes the diagonals an artifact of who shares the panel
  (see `code/impute_diagnostic*.R`). Missing calls are filled by `GRM_IMPUTE`
  (`"glmnet"` = robust per-marker elastic-net, default; `"mean"` = column mean).
- **Prediction targets / subsetting**: protocol coverage is driven by the **training
  + test** accessions (`find_and_get_genotypes(train_acc, test_acc, test_names)`;
  `TEST_ACCESSIONS` dbIds resolved via the germplasm cache). After the combine the GRM
  is **subset to the prediction targets** (training ∪ predictable-test) — bridges/extra
  accessions only inform the combine — and any training accession with no
  genotype/pedigree is **injected** with a prior-only diagonal. `markers` is therefore
  kept only for the single-protocol, no-pedigree, no-injection case.
- **Pedigree stitch**: `PEDIGREE_DIR` points at the sibling `BrAPI_pedigree_relmat`
  project's `output/<id>` folder; group relationship matrices that overlap our
  accessions **or bridges** are added to the EM combine as extra partials (weighted
  by `PEDIGREE_DF`). A pedigree group counts as a partial for bridge detection, so
  an accession seen in only one protocol is kept when it also sits in a pedigree
  group. Set to `NULL` to disable; a missing folder is skipped. **Depends on the
  sibling project's output contract** (see `CLAUDE.md` → Cross-repo dependencies):
  `<id>_group<N>.csv` triplet files, the companion `<id>_pedigree_groups.rds`
  (per-group membership), and `germplasm_cache_<id>.rds` (dbid→name) one level above
  `PEDIGREE_DIR`. If those companions are absent it degrades to phenotyped-only
  pedigree treatment.
- VCF marker QC drops high-missing **markers before accessions**: T3 projects
  under one protocol can use different genome annotations (different marker
  names), so merging leaves every accession ~50% missing until the
  minority-annotation markers are removed.
