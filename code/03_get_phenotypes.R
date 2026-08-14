# Step 3: pull phenotypes + field design for a set of studies.
#
# /observationunits?includeObservations=true returns, in one call per study:
#   - the observation values (observationVariableName + value)
#   - the germplasm (germplasmDbId / germplasmName)
#   - the field design: rep & block from observationLevelRelationships, and
#     row / col from positionCoordinateY / positionCoordinateX, plus entryType.

library(tidyverse)
here::i_am("code/03_get_phenotypes.R")
source(here::here("code", "01_connect.R"))

# Pull the levelCode for a named design level (e.g. "rep", "block") out of an
# observationUnit's observationLevelRelationships list.
.level_code <- function(rels, name) {
  for (r in rels) {
    if (!is.null(r$levelName) && r$levelName == name) return(r$levelCode %||% NA_character_)
  }
  NA_character_
}

.study_cache_path <- function(sid) cache_path("pheno_cache", paste0(sid, ".rds"))

# One study's observation units, cached to data/pheno_cache/<studyDbId>.rds and
# stored UNFILTERED by trait: the trait filter is applied at read time in
# get_phenotypes(), so changing TRAIT_NAMES never costs a download. Mirrors the
# per-protocol VCF cache in step 4.
.study_observations <- function(conn, sid, refresh = FALSE) {
  path <- .study_cache_path(sid)
  if (!refresh && file.exists(path)) return(read_rds(path))

  resp <- conn$get("/observationunits",
                   query = list(studyDbId = sid, includeObservations = "true"),
                   page = "all", pageSize = 500)
  units <- resp$combined_data %||% resp$data
  rows <- map_dfr(units, function(u) {
    obs <- u$observations
    if (length(obs) == 0) return(tibble())
    pos  <- u$observationUnitPosition
    rels <- pos$observationLevelRelationships
    # one tibble per plot (observations vectorized) -- much faster than one
    # tibble per observation when a study has thousands of measurements.
    tibble(
      studyDbId     = sid,
      germplasmDbId = as.character(u$germplasmDbId %||% NA),
      germplasmName = u$germplasmName %||% NA_character_,
      obsUnitDbId   = as.character(u$observationUnitDbId %||% NA),
      rep       = .level_code(rels, "rep"),
      block     = .level_code(rels, "block"),
      row       = pos$positionCoordinateY %||% NA,
      col       = pos$positionCoordinateX %||% NA,
      entryType = pos$entryType %||% NA_character_,
      trait     = map_chr(obs, ~ .x$observationVariableName %||% NA_character_),
      value     = map_dbl(obs, ~ suppressWarnings(as.numeric(.x$value %||% NA)))
    )
  })
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  write_rds(rows, path)
  rows
}

#' Pull phenotypes and design for the given study ids.
#'
#' @param conn        BrAPI connection.
#' @param study_ids   character/numeric vector of studyDbIds.
#' @param trait_names exact observationVariableName strings to keep (empty = all
#'   traits).
#' @return list(pheno, design, accessions). Cached to data/phenotypes.rds, with the
#'   raw per-study pulls under data/pheno_cache/ (see .study_observations).
get_phenotypes <- function(conn, study_ids,
                           trait_names = TRAIT_NAMES, refresh = FALSE) {
  study_ids <- as.character(study_ids)
  key   <- cache_key(study_ids = study_ids, trait_names = trait_names)
  cache <- cache_path("phenotypes.rds")
  hit   <- cache_read(cache, key, refresh)
  if (!is.null(hit)) return(hit)

  # Assemble from the per-study caches, downloading only the studies we lack. So
  # adding a trial costs one study's download, not the whole set -- and changing
  # TRAIT_NAMES costs nothing at all, because the per-study files are unfiltered.
  cached <- map_lgl(study_ids, ~ file.exists(.study_cache_path(.x)) && !refresh)
  if (any(!cached)) {
    say("Downloading observations for ", sum(!cached), " of ", length(study_ids),
        " studies (network-bound; the first run takes a while) ...")
  }
  if (any(cached)) note(sum(cached), " study/studies already cached")
  rows <- map(study_ids, ~ .study_observations(conn, .x, refresh = refresh),
              .progress = pb_wrap("Phenotypes: studies")) |>
    list_rbind()

  if (length(trait_names) > 0 && nrow(rows) > 0) {
    kept <- filter(rows, trait %in% trait_names)   # exact observationVariableName match
    note(nrow(rows), " observations pulled; ", nrow(kept), " kept over the ",
         length(trait_names), " configured traits (",
         n_distinct(kept$trait), " of them present)")
    rows <- kept
  }

  pheno  <- rows |> filter(!is.na(value))
  design <- pheno |>
    distinct(studyDbId, obsUnitDbId, germplasmDbId, germplasmName,
             rep, block, row, col, entryType)
  accessions <- pheno |>
    distinct(germplasmDbId, germplasmName) |>
    filter(!is.na(germplasmDbId))

  out <- list(pheno = pheno, design = design, accessions = accessions)
  cache_write(cache, out, key)
  out
}

#' Split a phenotype pull into training vs test sets using trial roles.
#'
#' @param pheno  get_phenotypes() output (uses its `pheno` long table).
#' @param trials find_ny_trials() output with a `role` column ("training"/"test").
#' @return list(train_pheno, train_acc, test_acc):
#'   - train_pheno: long phenotypes from training trials only (feeds stage1_blues);
#'   - train_acc / test_acc: distinct (germplasmDbId, germplasmName) observed in the
#'     training / test trials respectively.
split_by_role <- function(pheno, trials) {
  train_ids <- trials$studyDbId[trials$role == "training"]
  test_ids  <- trials$studyDbId[trials$role == "test"]
  p <- pheno$pheno
  acc_in <- function(ids) {
    p |>
      filter(studyDbId %in% ids, !is.na(germplasmDbId)) |>
      distinct(germplasmDbId, germplasmName)
  }
  list(
    train_pheno = filter(p, studyDbId %in% train_ids),
    train_acc   = acc_in(train_ids),
    test_acc    = acc_in(test_ids)
  )
}
