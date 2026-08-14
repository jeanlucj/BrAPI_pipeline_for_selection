# Step 2: select the trials to analyze, tagged training vs test.
#
# Training trials (two modes):
#   * TRAINING_TRIALS set  -> train on exactly those trials (matched by studyName),
#     ignoring the geographic search.
#   * TRAINING_TRIALS NULL -> pull all locations (GeoJSON coordinates), keep those
#     within RADIUS_KM of the NY center point, and keep studies of any type in
#     STUDY_TYPES grown there in the requested YEARS.
# Test trials: studies matched by TEST_TRIALS studyName (their accessions are
# prediction targets; their phenotypes are not used for training). A trial named
# in both lists is treated as training. The returned tibble carries a `role`
# column ("training"/"test").
# Either way location metadata (name, lon/lat, distance_km) is joined on. (The
# /studies endpoint does not honor a server-side locationDbId filter, so
# filtering is done client-side.)

library(tidyverse)
here::i_am("code/02_find_trials.R")
source(here::here("code", "01_connect.R"))

#' Pull all locations as a tidy tibble with longitude / latitude extracted from
#' the GeoJSON `coordinates` (geometry$coordinates = c(lon, lat, alt)).
get_locations <- function(conn) {
  say("Fetching locations (all pages) ...")
  resp <- conn$get("/locations", page = "all", pageSize = 1000)
  dat  <- resp$combined_data %||% resp$data
  map_dfr(dat, function(l) {
    coords <- l$coordinates$geometry$coordinates
    tibble(
      locationDbId = as.character(l$locationDbId),
      locationName = l$locationName %||% NA_character_,
      countryCode  = l$countryCode %||% NA_character_,
      longitude    = if (length(coords) >= 1) as.numeric(coords[[1]]) else NA_real_,
      latitude     = if (length(coords) >= 2) as.numeric(coords[[2]]) else NA_real_
    )
  })
}

#' Pull all studies as a tidy tibble (one row per study).
get_studies <- function(conn) {
  say("Fetching studies (all pages; /studies ignores a server-side location filter) ...")
  resp <- conn$get("/studies", page = "all", pageSize = 1000)
  dat  <- resp$combined_data %||% resp$data
  map_dfr(dat, function(s) {
    tibble(
      studyDbId    = as.character(s$studyDbId),
      studyName    = s$studyName %||% NA_character_,
      studyType    = s$studyType %||% NA_character_,
      locationDbId = as.character(s$locationDbId %||% NA),
      year         = suppressWarnings(as.integer(first(unlist(s$seasons)))),
      trialName    = s$trialName %||% NA_character_
    )
  })
}

#' Select the trials to analyze (training + test).
#'
#' Training trials: if `training_trials` is not NULL, exactly those (matched by
#' studyName); otherwise trials of any type in `study_types` grown within
#' `radius_km` of the center point in `years`. Test trials: those matched by
#' `test_trials` that are not already training (training takes precedence).
#'
#' @return A tibble of studies (studyDbId, names, location, year, distance_km)
#'   with a `role` column ("training"/"test"). Cached to data/ny_trials.rds.
find_ny_trials <- function(conn,
                           center_lat = CENTER_LAT, center_lon = CENTER_LON,
                           radius_km = RADIUS_KM, years = YEARS,
                           study_types = STUDY_TYPES,
                           training_trials = TRAINING_TRIALS,
                           test_trials = TEST_TRIALS,
                           refresh = FALSE) {
  cache <- cache_path("ny_trials.rds")
  if (!refresh && file.exists(cache)) {
    note_cache(cache)
    return(read_rds(cache))
  }

  # Location metadata, with great-circle distance from the center point.
  locs <- get_locations(conn) |>
    filter(!is.na(latitude), !is.na(longitude)) |>
    mutate(distance_km = as.numeric(
      geosphere::distHaversine(cbind(longitude, latitude),
                               c(center_lon, center_lat))) / 1000)

  studies <- get_studies(conn)
  note(nrow(locs), " located sites, ", nrow(studies), " studies on the server")

  if (!is.null(training_trials)) {
    # Explicit training list: select by name, keep even if the location lacks
    # coordinates (left_join), regardless of type/year.
    train <- studies |>
      filter(studyName %in% training_trials) |>
      left_join(locs, by = "locationDbId") |>
      arrange(year, studyName)
    missing <- setdiff(training_trials, train$studyName)
    if (length(missing) > 0) {
      warning("TRAINING_TRIALS not found on the server: ",
              paste(missing, collapse = ", "))
    }
  } else {
    # Geographic search: trials of any allowable type within the radius / years.
    train <- studies |>
      filter(studyType %in% study_types, year %in% years) |>
      inner_join(filter(locs, distance_km <= radius_km), by = "locationDbId") |>
      arrange(distance_km, year)
  }
  train <- mutate(train, role = "training")

  # Test trials: by name, excluding any already selected as training.
  test <- tibble()
  if (!is.null(test_trials)) {
    test <- studies |>
      filter(studyName %in% test_trials, !studyDbId %in% train$studyDbId) |>
      left_join(locs, by = "locationDbId") |>
      arrange(year, studyName) |>
      mutate(role = "test")
    missing <- setdiff(test_trials, c(test$studyName, train$studyName))
    if (length(missing) > 0) {
      warning("TEST_TRIALS not found on the server: ",
              paste(missing, collapse = ", "))
    }
  }

  out <- bind_rows(train, test)
  write_rds(out, cache)
  out
}
