# Central configuration for the BrAPI -> T3/Oat selection pipeline.
# Every step sources this file. Edit the values here to retarget the pipeline
# (different center point, region size, years, traits, or prediction settings).

library(tidyverse)
here::i_am("code/config.R")

# --- Database ----------------------------------------------------------------
DB_NAME <- "T3/Oat"

# --- Refresh: force downloads from the database rather than use cache --------
PIPELINE_REFRESH <- FALSE

# --- Target environment: New York --------------------------------------------
# Center point (Ithaca, NY) and the radius that defines "NY + surrounding
# region". Trials at locations within RADIUS_KM of the center are kept.
CENTER_LAT <- 42.44
CENTER_LON <- -76.50
RADIUS_KM  <- 500

# --- Training-trial selection ------------------------------------------------
# Training trials: a vector of studyName values whose phenotypes TRAIN the
# genomic-prediction model. When NOT NULL the pipeline trains on exactly these
# trials and ignores the geographic search (radius / years) below. NULL = use
# the location-radius search instead. Accessions in the training trials are
# ALWAYS predicted (and force-included in the relationship matrix even if they
# have neither genotypes nor pedigree).
TRAINING_TRIALS <- NULL
TRAINING_TRIALS <- c("Cornell_WinterOatPeaIntercrop_2024_Ithaca",
                     "CU_2025_Ithaca_WOP_PLOT",
                     "CU_ARS_2026_WOP")

# Allowable study types: a vector; trials of ANY of these types are analyzed
# together in the two-stage BLUE -> GBLUP process. (Used by the radius search;
# TRAINING_TRIALS, when set, selects by name regardless of type.)
STUDY_TYPES <- c("phenotyping_trial", "Advanced Yield Trial",
                 "Preliminary Yield Trial", "Uniform Yield Trial",
                 "Variety Release Trial")   # T3/Oat studyType(s) for field trials
YEARS       <- 2015:2026                 # radius search: keep trials in these seasons

# --- Prediction targets (who to predict) -------------------------------------
# By default the pipeline predicts only the TRAINING_TRIALS accessions. Widen the
# prediction set with either or both of:
#   TEST_TRIALS     - a vector of studyName values; as many of their accessions as
#                     can be predicted (present in the relationship matrix).
#                     Their phenotypes are NOT used for training -- UNLESS a trial
#                     is also in TRAINING_TRIALS, in which case it trains.
#   TEST_ACCESSIONS - a vector of germplasmName values to predict (as many as can
#                     be). Their dbIds are resolved via the pedigree germplasm
#                     cache so their genotyping protocols are also downloaded.
# The prediction set is the union of TRAINING accessions, TEST_TRIALS accessions,
# and TEST_ACCESSIONS. NULL disables each.
TEST_TRIALS     <- NULL
TEST_ACCESSIONS <- NULL

# --- Traits ------------------------------------------------------------------
# Case-insensitive substrings matched against the BrAPI observationVariableName
# (e.g. "Grain yield - g/m2|CO_350:0000260"). Empty vector = keep all traits.
TRAIT_PATTERNS <- c("yield")
TRAIT_PATTERNS <- c()

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

# --- EM-combine degrees of freedom -------------------------------------------
# In the Wishart-EM combiner each partial covariance has a degrees-of-freedom
# (df = effective # of independent samples) that acts as its relative WEIGHT when
# the partials are merged. Marker GRMs get a df derived from their effective
# number of independent samples (Galwey 2009 measure on the GRM eigenvalues), but
# that measure only sets the GRMs' relative ordering: the values are then
# re-centered on GRM_DF_MEAN with a spread capped at GRM_DF_STDEV. PEDIGREE_DF is
# the fixed df for pedigree partials, so GRM_DF_MEAN vs PEDIGREE_DF controls how
# much more marker matrices are trusted than pedigree (here ~2x). Keep
# GRM_DF_MEAN comfortably larger than GRM_DF_STDEV so all dfs stay positive.
GRM_DF_MEAN  <- 60   # center of the marker-GRM dfs
GRM_DF_STDEV <- 15   # cap on the spread of the marker-GRM dfs
PEDIGREE_DF  <- 30   # fixed EM degrees-of-freedom weight for pedigree partials

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
