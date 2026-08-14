# Shared setup + synthetic-data builders for the test suite.
# testthat sources helper*.R into the test environment before running tests.

suppressMessages(library(tidyverse))

# Never narrate the test suite. SHOW_PROGRESS already defaults to interactive()
# (FALSE under Rscript), but the suite is also run from RStudio, so pin it.
options(brapi.progress = FALSE)

# Source the whole pipeline (offline: defines functions, does not connect).
for (f in c("config.R", "grm_utils.R", "em_covariance_combiner.R",
            "01_connect.R", "02_find_trials.R", "03_get_phenotypes.R",
            "04_find_genotyping.R", "05_stage1_blues.R",
            "06_stage2_genomic_prediction.R", "07_select.R")) {
  source(here::here("code", f))
}

# Redirect caches/outputs to a throwaway dir so tests never touch data/ or
# output/ in the repo. The pipeline functions are source()'d into the global
# environment (so is cache_path/output_path from config.R); we must override
# them THERE with `<<-`, not in this helper's environment, or the functions
# would keep resolving the real data/ paths.
.TEST_TMP <- file.path(tempdir(), "brapi_tests")
dir.create(.TEST_TMP, showWarnings = FALSE, recursive = TRUE)
cache_path  <<- function(...) { p <- file.path(.TEST_TMP, ...); dir.create(dirname(p), showWarnings = FALSE, recursive = TRUE); p }
output_path <<- function(...) { p <- file.path(.TEST_TMP, ...); dir.create(dirname(p), showWarnings = FALSE, recursive = TRUE); p }

# --- synthetic-data builders -------------------------------------------------

# Write a minimal vcfR-readable VCF from a dosage matrix (accessions x markers,
# values 0/1/2/NA). All IDs are "." so .vcf_to_dosage exercises its unique-ID fix.
write_test_vcf <- function(path, dosage) {
  gt_of <- c("0" = "0/0", "1" = "0/1", "2" = "1/1")
  hdr <- c("##fileformat=VCFv4.2",
           "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
           paste(c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT",
                   rownames(dosage)), collapse = "\t"))
  body <- vapply(seq_len(ncol(dosage)), function(j) {
    calls <- ifelse(is.na(dosage[, j]), "./.", gt_of[as.character(dosage[, j])])
    paste(c("1", j * 100L, ".", "A", "G", ".", "PASS", ".", "GT", calls), collapse = "\t")
  }, character(1))
  writeLines(c(hdr, body), path)
  path
}

# A reproducible accessions x markers dosage matrix.
make_dosage <- function(n_acc = 20, n_mar = 50, seed = 1, prefix = "L") {
  set.seed(seed)
  M <- matrix(sample(0:2, n_acc * n_mar, replace = TRUE), n_acc, n_mar)
  dimnames(M) <- list(paste0(prefix, seq_len(n_acc)),
                      paste0("m", seq_len(n_mar)))
  M
}

# Fake BrAPI connection whose $get(call, ...) returns canned locations/studies/
# observationunits. `locations`/`studies`/`units_by_study` are plain lists shaped
# like BrAPI $data elements.
fake_conn <- function(locations = list(), studies = list(), units_by_study = list()) {
  list(get = function(call, query = NULL, page = NULL, pageSize = NULL) {
    if (grepl("/locations", call)) return(list(data = locations))
    if (grepl("/studies",   call)) return(list(data = studies))
    if (grepl("/observationunits", call)) {
      sid <- as.character(query$studyDbId)
      return(list(data = units_by_study[[sid]] %||% list()))
    }
    stop("fake_conn: unexpected call ", call)
  })
}

loc_rec <- function(id, name, lon, lat, country = "USA") {
  list(locationDbId = id, locationName = name, countryCode = country,
       coordinates = list(geometry = list(coordinates = list(lon, lat, 100))))
}
study_rec <- function(id, name, type, loc, year) {
  list(studyDbId = id, studyName = name, studyType = type,
       locationDbId = loc, seasons = list(as.character(year)), trialName = name)
}
# One observation unit with embedded observations for the given traits/values.
unit_rec <- function(ou, gdbid, gname, rep, block, x, y, traits, values) {
  list(
    germplasmDbId = gdbid, germplasmName = gname, observationUnitDbId = ou,
    observationUnitPosition = list(
      positionCoordinateX = x, positionCoordinateY = y, entryType = "test",
      observationLevelRelationships = list(
        list(levelName = "rep",   levelCode = as.character(rep)),
        list(levelName = "block", levelCode = as.character(block)),
        list(levelName = "plot",  levelCode = ou))),
    observations = Map(function(t, v)
      list(observationVariableName = t, value = as.character(v)), traits, values))
}
