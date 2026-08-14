# Central configuration for the BrAPI -> T3/Oat selection pipeline.
# Every step sources this file. Edit the values here to retarget the pipeline
# (different center point, region size, years, traits, or prediction settings).

library(tidyverse)
here::i_am("code/config.R")

# Console progress/status/timing helpers (say, note, step_start, pb_start, ...) and
# the request-keyed cache (cache_key/cache_read/cache_write). Sourced here because
# every step sources config.R. cache.R uses say()/note(), so progress.R comes first.
source(here::here("code", "progress.R"))
source(here::here("code", "cache.R"))

# --- Progress reporting ------------------------------------------------------
# Whether the steps narrate themselves (status lines, progress bars, per-step
# timings). Only meaningful in a live console: Rscript, the test suite, and the
# non-interactive subprocess wflow_build() renders in have no terminal to animate,
# so the default is interactive(). Force either way with
# options(brapi.progress = TRUE/FALSE).
SHOW_PROGRESS <- interactive()

# --- Selection lists from data/config/ ---------------------------------------
# The long lists that define a run (which trials train, which accessions to
# predict, which traits to keep) live as plain text files in data/config/ rather
# than as literal vectors in this file: they run to hundreds of lines, and a run
# should be re-targeted by swapping a file, not by editing source. data/config/ is
# the one subfolder of data/ that is tracked in git.
#
# Format: one value per line. Blank lines and lines starting with # are ignored,
# and surrounding whitespace is trimmed. Values are used VERBATIM otherwise -- a
# germplasmName like "APPLER|CIAV2680" or a trait name like
# "Grain yield - g/m2|CO_350:0000260" keeps its pipe.

# Resolve a config file name to a path: bare names are looked up in data/config/,
# but an existing/absolute path is taken as given.
.config_path <- function(file) {
  if (file.exists(file)) return(file)
  here::here("data", "config", file)
}

# Read the non-blank, non-comment lines of a config file. read_lines() (not
# readLines()) because these files often lack a trailing newline.
.config_raw <- function(path) {
  if (!file.exists(path)) {
    stop("config file not found: ", path,
         "\n  (expected in data/config/; see data/README.md)", call. = FALSE)
  }
  readr::read_lines(path) |>
    stringr::str_trim() |>
    (\(x) x[nzchar(x) & !stringr::str_starts(x, "#")])()
}

#' Read a one-value-per-line config list (trials, accessions).
#'
#' @param file File name in data/config/, or a path.
#' @return Character vector of values, or NULL if the file holds none.
config_lines <- function(file) {
  path <- .config_path(file)
  x <- .config_raw(path)
  dup <- unique(x[duplicated(x)])
  if (length(dup)) {
    message("config_lines(", basename(path), "): dropping ", length(dup),
            " duplicate value(s), e.g. ", paste(head(dup, 3), collapse = ", "))
    x <- unique(x)
  }
  if (!length(x)) NULL else x
}

#' Read the trait config file: names, plus optional weights and short names.
#'
#' One trait per line, as 1 to 3 TAB-separated columns:
#'   <observationVariableName> [<index weight> [<short label>]]
#' A column that is present must be present on EVERY line -- a half-filled weight
#' column would silently zero out part of the selection index. A column absent
#' everywhere yields NULL, i.e. "no selection index" / "use the full names".
#'
#' @param file File name in data/config/, or a path.
#' @return list(names = chr, weights = named num or NULL, short = named chr or NULL).
config_traits <- function(file) {
  path <- .config_path(file)
  x <- .config_raw(path)
  if (!length(x)) return(list(names = c(), weights = NULL, short = NULL))

  cols <- stringr::str_split_fixed(x, "\t", 3) |>
    apply(2, stringr::str_trim)      # matrix: [line, column]
  dim(cols) <- c(length(x), 3)       # apply() drops dim when there is 1 line
  nms <- cols[, 1]

  # A column is either fully filled or fully empty; anything else is a mistake.
  .column <- function(j, what, parse = identity) {
    v <- cols[, j]
    if (!any(nzchar(v))) return(NULL)
    if (!all(nzchar(v))) {
      stop("config_traits(", basename(path), "): the ", what, " column (column ", j,
           ") is filled on some lines but not on: ",
           paste(nms[!nzchar(v)], collapse = ", "),
           "\n  Give every trait a ", what, ", or none of them.", call. = FALSE)
    }
    setNames(parse(v), nms)
  }
  .as_weight <- function(v) {
    w <- suppressWarnings(as.numeric(v))
    if (anyNA(w)) {
      stop("config_traits(", basename(path), "): non-numeric weight(s): ",
           paste(v[is.na(w)], collapse = ", "), call. = FALSE)
    }
    w
  }

  list(names   = nms,
       weights = .column(2, "weight", .as_weight),
       short   = .column(3, "short name"))
}

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
TRAINING_TRIALS <- config_lines("Trial_Sel2026_Intersect20NoYld.txt")
# TRAINING_TRIALS <- c("CU_2025_Ithaca_WOP_PLOT", "CU_ARS_2026_WOP")  # or inline
# TRAINING_TRIALS <- NULL                                # or the radius search

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
TEST_ACCESSIONS <- config_lines("Acc_Sel2026.txt")

# --- Traits ------------------------------------------------------------------
# All three come from one file, whose columns are TAB-separated (see
# config_traits() above):
#   col 1 TRAIT_NAMES       exact BrAPI observationVariableName strings to keep
#                           (matched exactly, not as substrings). Empty = all traits.
#   col 2 TRAIT_WEIGHTS     optional selection-index weight per trait (sign =
#                           direction; magnitude = importance/scale). Absent = no
#                           index, i.e. a plain all-traits prediction dump.
#   col 3 TRAIT_SHORT_NAMES optional short labels for the breeders_output.csv
#                           columns. Absent = use the full names.
.traits <- config_traits("Trait_Sel2026.txt")
TRAIT_NAMES       <- .traits$names
TRAIT_WEIGHTS     <- .traits$weights
TRAIT_SHORT_NAMES <- .traits$short
rm(.traits)

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

# Synonym canonicalization of VCF sample names. T3 accessions are sometimes
# genotyped under a *preliminary* line name that is later demoted to a SYNONYM when
# a final name is assigned, so an archived VCF can cover our accessions while its
# sample IDs are synonyms of -- not equal to -- our primary germplasmNames (this is
# exactly what loses ~all the Oat 3K array matches). When TRUE, step 4 builds an
# alias -> primary lookup for our accessions (T3_brapi_helpers::build_synonym_lookup
# over BrAPI /search/germplasm) and relabels VCF sample names to the primary name
# before matching. Cached to data/synonym_map.rds. FALSE = exact-name match only.
USE_SYNONYMS <- TRUE

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

# --- Stage-1 BLUEs -----------------------------------------------------------
# Floor on each BLUE's standard error, expressed as a fraction of that trial x
# trait's response SD, applied before forming the Stage-2 weight = 1/SE^2. A trial
# whose observations are (near-)constant within genotype -- e.g. duplicated records,
# or a degenerate design -- fits perfectly, driving the residual variance and hence
# every genotype SE toward 0 and the weight toward Inf, which lets one such trial
# dominate Stage 2. Flooring SE caps the weight without touching well-estimated
# BLUEs, whose SEs sit far above this fraction. A trial with NO response variation
# at all is dropped outright (no estimable genotypic signal), so this only backstops
# the partially-degenerate cases. Set to 0 to disable.
SE_FLOOR_FRAC <- 0.01

# --- Genomic prediction (Stage 2) --------------------------------------------
# Mixed-model engine for the kernel GBLUP. Prediction is ALWAYS kernel-based (a
# genomic/relationship matrix), never a marker-effect Bayesian-alphabet model, so
# the single-protocol and EM-combined cases are handled identically.
#   "bglr"   = BGLR Bayesian kernel GBLUP, fit as ridge regression (BRR) on the
#              relationship eigen-factor (MCMC; uses BGLR_NITER / BGLR_BURNIN).
#              NB: BGLR's model = "RKHS" is deliberately NOT used -- it mishandles
#              the 1/SE^2 weights and inflates GEBVs (see code/bglr_rkhs_vs_brr.R).
#   "sommer" = sommer REML GBLUP
MIXED_MODEL_ENGINE <- "bglr"
BGLR_NITER  <- 12000   # BGLR MCMC iterations (engine "bglr")
BGLR_BURNIN <- 2000
SEED        <- 1234

# --- Marker QC ---------------------------------------------------------------
MAX_MISSING <- 0.50   # drop markers / accessions with > this fraction missing
MIN_MAF     <- 0.01   # drop markers below this minor-allele frequency

# Imputation of the residual missing calls before a marker GRM is built.
#   "glmnet" = per-marker elastic-net regression on the most-correlated markers
#              (robust: a marker whose CV fails falls back to its mean); borrows
#              information across markers, so it needs a reasonably large panel to
#              be worthwhile (see GRM_PANEL_MIN).
#   "mean"   = fill each missing call with the marker mean (fast, no borrowing).
GRM_IMPUTE <- "glmnet"

# Per-protocol panel size for imputation + GRM estimation. VanRaden allele
# frequencies, imputation, and the MAF/missingness filters are all estimated from
# the accessions in the panel, so estimating them on only our handful of relevant
# accessions (targets + bridges) makes both the imputation and the GRM diagonals
# unstable (a self-relationship measured against a few panel-mates, not a stable
# reference). So we build each partial GRM on ALL its relevant accessions, then
# fill in with non-relevant accessions from the same protocol up to a panel of
# max(GRM_PANEL_MIN, n_relevant), and subset the GRM back to the relevant set
# afterwards. If n_relevant >= GRM_PANEL_MIN no filling happens.
GRM_PANEL_MIN <- 1000

# --- Paths (workflowr layout) ------------------------------------------------
cache_path  <- function(...) here::here("data", ...)    # cached BrAPI pulls
output_path <- function(...) here::here("output", ...)  # results / artifacts
