# Synonym canonicalization: relabel VCF sample names that are *synonyms* of our
# accessions back to the primary germplasmName the rest of the pipeline keys on.
#
# The generalized functions live in the sibling T3_brapi_helpers package
# (build_synonym_lookup / get_synonyms_from_germplasm_names / canonicalize_to_primary),
# shared with the T3Predictathon2026 project. The package may not be installed, and
# an installed copy can lag the sibling repo, so -- following
# T3Predictathon2026/scripts/Analysis_Claude/build_synonym_map.R -- we source the
# current helper file directly when the package is absent.

library(tidyverse)
here::i_am("code/synonyms.R")
source(here::here("code", "config.R"))

.t3_helpers_file <- here::here("..", "T3_brapi_helpers", "R", "brapi_germplasm.R")
if (requireNamespace("T3_brapi_helpers", quietly = TRUE)) {
  build_synonym_lookup           <- T3_brapi_helpers::build_synonym_lookup
  canonicalize_to_primary        <- T3_brapi_helpers::canonicalize_to_primary
  get_synonyms_from_germplasm_names <- T3_brapi_helpers::get_synonyms_from_germplasm_names
} else if (file.exists(.t3_helpers_file)) {
  source(.t3_helpers_file)                       # defines the three functions
}

# Fallback so .vcf_to_dosage still works (as a no-op) when neither the package nor
# the sibling repo is available: with an empty lookup this is the identity.
if (!exists("canonicalize_to_primary")) {
  canonicalize_to_primary <- function(sample_ids, alias_lookup) {
    if (length(alias_lookup) == 0) return(sample_ids)
    hit <- alias_lookup[sample_ids]
    ifelse(is.na(hit), sample_ids, unname(hit))
  }
}

# Whether synonym lookup is wired up (config on AND the helper available).
.synonyms_available <- function() {
  isTRUE(USE_SYNONYMS) && exists("build_synonym_lookup")
}

#' Cached alias -> primary lookup for our accessions.
#'
#' Builds (and caches to data/synonym_map.rds) a named character vector mapping
#' every known alias (each primary name to itself, plus each synonym to its
#' primary) for `names`, via BrAPI /search/germplasm. Returns an empty vector
#' (-> canonicalize_to_primary is a no-op) when USE_SYNONYMS is FALSE, the helper
#' is unavailable, or the BrAPI lookup fails -- so a synonym hiccup degrades to
#' exact-name matching rather than aborting the genotyping step.
build_alias_lookup <- function(conn, names, refresh = FALSE) {
  cache <- cache_path("synonym_map.rds")
  if (!refresh && file.exists(cache)) return(read_rds(cache))
  if (!.synonyms_available()) {
    if (isTRUE(USE_SYNONYMS))
      message("USE_SYNONYMS = TRUE but T3_brapi_helpers is unavailable ",
              "(install it or clone it next to this repo); using exact-name match.")
    return(character(0))
  }
  # A connection without a usable $search (e.g. the offline tests' fake_conn)
  # can't do a germplasm lookup; skip rather than spin through retries/backoff.
  if (!is.function(conn$search)) {
    message("BrAPI connection has no $search method; skipping synonym lookup.")
    return(character(0))
  }
  names <- unique(as.character(names))
  names <- names[!is.na(names) & nzchar(names)]
  lk <- tryCatch(
    build_synonym_lookup(names, brapi_connection = conn, verbose = TRUE),
    error = function(e) {
      warning("Synonym lookup failed (", conditionMessage(e),
              "); falling back to exact-name match.")
      character(0)
    })
  n_syn <- sum(names(lk) != unname(lk))
  message(sprintf("Synonym map: %d aliases over %d accessions (%d synonym aliases).",
                  length(lk), length(names), n_syn))
  write_rds(lk, cache)
  lk
}
