# Tier 1 (offline): trial selection logic via a fake BrAPI connection.

make_trials_conn <- function() {
  locs <- list(
    loc_rec("1", "Ithaca",  -76.50, 42.44),   # center (dist ~0)
    loc_rec("2", "FarAway", -100.0, 35.00),   # far (> radius)
    loc_rec("3", "NearNY",  -76.00, 43.00))   # within radius
  studies <- list(
    study_rec("s1", "Near_pheno_2018", "phenotyping_trial", "1", 2018),
    study_rec("s2", "Far_pheno_2018",  "phenotyping_trial", "2", 2018),
    study_rec("s3", "Near_pheno_2000", "phenotyping_trial", "1", 2000),
    study_rec("s4", "Near_geno_2018",  "genotyping",        "1", 2018),
    study_rec("s5", "Near3_pheno_2019","phenotyping_trial", "3", 2019))
  fake_conn(locations = locs, studies = studies)
}

test_that("radius mode keeps in-range phenotyping trials within the years", {
  out <- find_ny_trials(make_trials_conn(), refresh = TRUE, training_trials = NULL)
  expect_setequal(out$studyDbId, c("s1", "s5"))        # not far/old/wrong-type
  expect_true(all(out$distance_km <= 500))
  expect_true(all(c("locationName", "distance_km", "role") %in% names(out)))
  expect_true(all(out$role == "training"))
})

test_that("TRAINING_TRIALS mode selects exactly those trials regardless of type/year", {
  expect_warning(
    out <- find_ny_trials(make_trials_conn(), refresh = TRUE,
                          training_trials = c("Near_geno_2018",   # wrong type
                                              "Near_pheno_2000",  # excluded year
                                              "NOPE")),           # bogus
    "not found on the server: NOPE")
  expect_setequal(out$studyName, c("Near_geno_2018", "Near_pheno_2000"))
  expect_true(all(out$role == "training"))
})

test_that("TRAINING_TRIALS keeps a trial far outside the radius, and one with no coordinates", {
  # The radius/season/type search is bypassed ENTIRELY by an explicit list -- the
  # location join only attaches metadata. `distance_km` is therefore informational
  # here and may exceed radius_km or be NA; it is not a filter.
  locs <- list(loc_rec("1", "Ithaca",   -76.50, 42.44),
               loc_rec("2", "FarAway", -100.00, 35.00),   # ~2200 km, 4x the radius
               loc_rec("3", "NoCoords",      NA,    NA))
  studies <- list(study_rec("s1", "Far_pheno_2018", "phenotyping_trial", "2", 2018),
                  study_rec("s2", "NoCoord_trial",  "phenotyping_trial", "3", 2024))
  conn <- fake_conn(locations = locs, studies = studies)

  out <- find_ny_trials(conn, refresh = TRUE, radius_km = 500,
                        training_trials = c("Far_pheno_2018", "NoCoord_trial"))
  expect_setequal(out$studyName, c("Far_pheno_2018", "NoCoord_trial"))
  expect_gt(out$distance_km[out$studyName == "Far_pheno_2018"], 500)   # kept anyway
  expect_true(is.na(out$distance_km[out$studyName == "NoCoord_trial"]))
  expect_true(all(out$role == "training"))
})

test_that("STUDY_TYPES admits multiple types together", {
  out <- find_ny_trials(make_trials_conn(), refresh = TRUE, training_trials = NULL,
                        study_types = c("phenotyping_trial", "genotyping"))
  expect_setequal(out$studyDbId, c("s1", "s4", "s5"))  # s4 now admitted
})

test_that("test_trials are tagged 'test'; training wins on overlap", {
  out <- find_ny_trials(make_trials_conn(), refresh = TRUE,
                        training_trials = "Near_pheno_2018",
                        test_trials = c("Near3_pheno_2019", "Near_pheno_2018"))
  roles <- setNames(out$role, out$studyDbId)
  expect_equal(unname(roles["s1"]), "training")   # overlap -> training
  expect_equal(unname(roles["s5"]), "test")
  expect_equal(sum(out$studyDbId == "s1"), 1L)    # not duplicated as test
})
