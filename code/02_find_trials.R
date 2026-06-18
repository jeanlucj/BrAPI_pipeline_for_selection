# Step 2: select the trials to analyze.
#
# Two modes:
#   * STUDY_NAMES set  -> analyze exactly those trials (matched by studyName),
#     ignoring the geographic search.
#   * STUDY_NAMES NULL -> pull all locations (GeoJSON coordinates), keep those
#     within RADIUS_KM of the NY center point, and keep studies of any type in
#     STUDY_TYPES grown there in the requested YEARS.
# Either way location metadata (name, lon/lat, distance_km) is joined on. (The
# /studies endpoint does not honor a server-side locationDbId filter, so
# filtering is done client-side.)

library(tidyverse)
here::i_am("code/02_find_trials.R")
source(here::here("code", "01_connect.R"))

#' Pull all locations as a tidy tibble with longitude / latitude extracted from
#' the GeoJSON `coordinates` (geometry$coordinates = c(lon, lat, alt)).
get_locations <- function(conn) {
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

#' Select the trials to analyze.
#'
#' If `study_names` is not NULL, returns exactly those trials (matched by
#' studyName), ignoring the geographic search. Otherwise returns trials of any
#' type in `study_types` grown within `radius_km` of the center point in `years`.
#'
#' @return A tibble of studies (studyDbId, names, location, year, distance_km).
#'   Cached to data/ny_trials.rds.
find_ny_trials <- function(conn,
                           center_lat = CENTER_LAT, center_lon = CENTER_LON,
                           radius_km = RADIUS_KM, years = YEARS,
                           study_types = STUDY_TYPES, study_names = STUDY_NAMES,
                           refresh = FALSE) {
  cache <- cache_path("ny_trials.rds")
  if (!refresh && file.exists(cache)) return(read_rds(cache))

  # Location metadata, with great-circle distance from the center point.
  locs <- get_locations(conn) |>
    filter(!is.na(latitude), !is.na(longitude)) |>
    mutate(distance_km = as.numeric(
      geosphere::distHaversine(cbind(longitude, latitude),
                               c(center_lon, center_lat))) / 1000)

  studies <- get_studies(conn)

  if (!is.null(study_names)) {
    # Explicit trial list: select by name, keep even if the location lacks
    # coordinates (left_join), regardless of type/year.
    out <- studies |>
      filter(studyName %in% study_names) |>
      left_join(locs, by = "locationDbId") |>
      arrange(year, studyName)
    missing <- setdiff(study_names, out$studyName)
    if (length(missing) > 0) {
      warning("STUDY_NAMES not found on the server: ",
              paste(missing, collapse = ", "))
    }
  } else {
    # Geographic search: trials of any allowable type within the radius / years.
    out <- studies |>
      filter(studyType %in% study_types, year %in% years) |>
      inner_join(filter(locs, distance_km <= radius_km), by = "locationDbId") |>
      arrange(distance_km, year)
  }

  write_rds(out, cache)
  out
}
