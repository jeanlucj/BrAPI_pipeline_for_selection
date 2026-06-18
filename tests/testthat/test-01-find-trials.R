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
  out <- find_ny_trials(make_trials_conn(), refresh = TRUE)
  expect_setequal(out$studyDbId, c("s1", "s5"))        # not far/old/wrong-type
  expect_true(all(out$distance_km <= 500))
  expect_true(all(c("locationName", "distance_km") %in% names(out)))
})

test_that("STUDY_NAMES mode selects exactly those trials regardless of type/year", {
  expect_warning(
    out <- find_ny_trials(make_trials_conn(), refresh = TRUE,
                          study_names = c("Near_geno_2018",   # wrong type
                                          "Near_pheno_2000",  # excluded year
                                          "NOPE")),           # bogus
    "not found on the server: NOPE")
  expect_setequal(out$studyName, c("Near_geno_2018", "Near_pheno_2000"))
})

test_that("STUDY_TYPES admits multiple types together", {
  out <- find_ny_trials(make_trials_conn(), refresh = TRUE,
                        study_types = c("phenotyping_trial", "genotyping"))
  expect_setequal(out$studyDbId, c("s1", "s4", "s5"))  # s4 now admitted
})
