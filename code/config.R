# Central configuration for the BrAPI -> T3/Oat selection pipeline.
# Every step sources this file. Edit the values here to retarget the pipeline
# (different center point, region size, years, traits, or prediction settings).

library(tidyverse)
here::i_am("code/config.R")

# --- Database ----------------------------------------------------------------
DB_NAME <- "T3/Oat"

# --- Target environment: New York --------------------------------------------
# Center point (Ithaca, NY) and the radius that defines "NY + surrounding
# region". Trials at locations within RADIUS_KM of the center are kept.
CENTER_LAT <- 42.44
CENTER_LON <- -76.50
RADIUS_KM  <- 500

# --- Study selection ---------------------------------------------------------
# Explicit trials: a vector of studyName values. When NOT NULL the pipeline
# analyzes exactly these trials and ignores the geographic search (radius /
# years) below. NULL = use the location-radius search instead.
STUDY_NAMES <- NULL

# Allowable study types: a vector; trials of ANY of these types are analyzed
# together in the two-stage BLUE -> GBLUP process. (Used by the radius search;
# STUDY_NAMES, when set, selects by name regardless of type.)
STUDY_TYPES <- c("phenotyping_trial")   # T3/Oat studyType(s) for field trials
YEARS       <- 2015:2025                 # radius search: keep trials in these seasons

# --- Traits ------------------------------------------------------------------
# Case-insensitive substrings matched against the BrAPI observationVariableName
# (e.g. "Grain yield - g/m2|CO_350:0000260"). Empty vector = keep all traits.
TRAIT_PATTERNS <- c("yield")

# --- Genotyping --------------------------------------------------------------
# Which genotyping protocol(s) to use for markers.
#   NULL  = use ALL protocols that genotype our accessions and EM-combine their
#           genomic relationship matrices into one (Wishart-EM combiner). This is
#           the general, multi-platform default. Protocols whose archived VCF is
#           too large to download (e.g. multi-GB GBS diversity panels) are
#           skipped with a warning rather than aborting the run.
#   <id>  = use a single protocol only (faster; also yields a raw marker matrix
#           so marker-effect BGLR models stay available). Find ids in the
#           coverage table printed by step 4 (e.g. an Oat 3K SNP array).
GENO_PROTOCOL_ID <- NULL

# Marker thinning: when a downloaded VCF carries more than TARGET_DENSITY markers
# it is thinned genome-wide (keep every n-th marker, n = floor(n_markers /
# TARGET_DENSITY)) before reading, to bound memory on very large GBS files.
TARGET_DENSITY <- 10000

# --- Pedigree stitch ---------------------------------------------------------
# Folder of precomputed pedigree relationship matrices from the sibling
# BrAPI_pedigree_relmat project (sparse-triplet <id>_group<k>.csv files). When
# set, any group matrix overlapping our accessions is supplied to the EM combiner
# as an extra partial covariance to stitch otherwise-disjoint marker GRMs.
# Set to NULL to disable. A missing folder is skipped gracefully.
PEDIGREE_DIR <- here::here("..", "BrAPI_pedigree_relmat", "output", "T3_Oat")
PEDIGREE_DF  <- 30   # EM degrees-of-freedom weight given to pedigree partials

# --- Genomic prediction (Stage 2, BGLR) --------------------------------------
BGLR_MODEL  <- "RKHS"   # "RKHS" (genomic relationship kernel) or "BRR"/"BayesB"
BGLR_NITER  <- 12000
BGLR_BURNIN <- 2000
SEED        <- 1234

# --- Marker QC ---------------------------------------------------------------
MAX_MISSING <- 0.50   # drop markers / accessions with > this fraction missing
MIN_MAF     <- 0.01   # drop markers below this minor-allele frequency

# --- Paths (workflowr layout) ------------------------------------------------
cache_path  <- function(...) here::here("data", ...)    # cached BrAPI pulls
output_path <- function(...) here::here("output", ...)  # results / artifacts
