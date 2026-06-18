# Code — BrAPI → T3/Oat selection pipeline

Reusable functions for a pipeline that pulls oat data from **T3/Oat** over
**BrAPI**, runs a two-stage genomic-prediction analysis, and ranks accessions to
use as parents for crossing in the **New York** environment. The orchestrating
report is `analysis/brapi_selection_pipeline.Rmd`.

## Steps

| File | Function | Does |
|------|----------|------|
| `config.R` | (parameters) | Center point/radius, years, traits, protocol, BGLR & QC settings |
| `01_connect.R` | `connect_t3()` | Connect to T3/Oat (anonymous; optional login) |
| `02_find_trials.R` | `find_ny_trials()` | Locations within a radius of NY → phenotyping studies |
| `03_get_phenotypes.R` | `get_phenotypes()` | Observations + field design (rep/block/row/col) per study |
| `04_find_genotyping.R` | `find_and_get_genotypes()` | Rank covering protocols/projects, download/thin VCFs, EM-combine per-protocol GRMs (+ pedigree) → combined GRM |
| `em_covariance_combiner.R` | `EMCovarianceCombiner()` | Wishart-EM combiner (copied from T3Predictathon2026) |
| `grm_utils.R` | `.Gmatrix()`, `std_grm()` | VanRaden GRM + standardization, shared by steps 4 & 6 |
| `05_stage1_blues.R` | `stage1_blues()` | Per-trial BLUEs (lme4) with weights |
| `06_stage2_genomic_prediction.R` | `stage2_gblup()`, `cv_accuracy()` | BGLR genomic prediction (RKHS on the combined GRM) → GEBVs |
| `07_select.R` | `select_parents()` | Selection index, relatedness flag, ranked CSV |

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
- **Multi-platform combine**: `GENO_PROTOCOL_ID = NULL` builds one standardized
  GRM per covering protocol and **EM-combines** them into a single GRM
  (`em_covariance_combiner.R`), using accessions genotyped on ≥2 platforms as
  bridges. Set `GENO_PROTOCOL_ID` to a single id (e.g. an **Oat 3K** array) for a
  faster single-platform run that also yields a raw marker matrix (so
  marker-effect BGLR models stay available). Downloaded VCFs are cached under
  `data/vcf_cache/`.
- **Pedigree stitch**: `PEDIGREE_DIR` points at the sibling `BrAPI_pedigree_relmat`
  project's `output/<id>` folder; precomputed group relationship matrices that
  overlap our accessions are added to the EM combine as extra partials (weighted
  by `PEDIGREE_DF`). Set to `NULL` to disable; a missing folder is skipped.
  *Limitation for T3/Oat*: that project only precomputed matrices for its small
  pedigree groups (the ~27,800-accession main component exceeded its
  `max_relmat_size`), so the stitch rarely engages for Oat until it is rerun with
  a higher cap.
- VCF marker QC drops high-missing **markers before accessions**: T3 projects
  under one protocol can use different genome annotations (different marker
  names), so merging leaves every accession ~50% missing until the
  minority-annotation markers are removed.
