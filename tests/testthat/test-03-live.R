# Tier 3 (LIVE): hits the real T3/Oat server and downloads data. Skipped unless
# RUN_LIVE_TESTS is set (e.g. RUN_LIVE_TESTS=true Rscript tests/run_tests.R).
# These are slow and network-dependent; they validate the real BrAPI plumbing.

skip_if_no_live <- function() {
  if (!nzchar(Sys.getenv("RUN_LIVE_TESTS"))) {
    skip("live tests disabled (set RUN_LIVE_TESTS=true to enable)")
  }
  # T3/Oat requires a login, so without credentials these would all error rather
  # than test anything. Skip cleanly instead.
  .load_project_renviron()
  if (!nzchar(Sys.getenv("T3_USERNAME")) || !nzchar(Sys.getenv("T3_PASSWORD"))) {
    skip("live tests need T3_USERNAME / T3_PASSWORD (see .Renviron.example)")
  }
}

test_that("connect_t3 returns a usable connection", {
  skip_if_no_live()
  conn <- connect_t3()
  expect_true(is.function(conn$get) || inherits(conn, "R6"))
})

test_that("find_ny_trials returns NY-region trials (radius mode)", {
  skip_if_no_live()
  conn <- connect_t3()
  trials <- find_ny_trials(conn, refresh = TRUE)
  expect_gt(nrow(trials), 0)
  expect_true(all(c("studyDbId","studyName","distance_km") %in% names(trials)))
  expect_true(all(trials$distance_km <= RADIUS_KM, na.rm = TRUE))
})

test_that("get_phenotypes pulls observations + design for a real study", {
  skip_if_no_live()
  conn <- connect_t3()
  trials <- find_ny_trials(conn, refresh = TRUE)
  ph <- get_phenotypes(conn, head(trials$studyDbId, 3),
                       trait_names = character(0), refresh = TRUE)
  expect_true(nrow(ph$pheno) > 0)
  expect_true(all(c("germplasmDbId","germplasmName") %in% names(ph$accessions)))
})

test_that("find_and_get_genotypes builds a GRM from a single small protocol", {
  skip_if_no_live()
  conn <- connect_t3()
  trials <- find_ny_trials(conn, refresh = TRUE)
  ph <- get_phenotypes(conn, head(trials$studyDbId, 6),
                       trait_names = character(0), refresh = TRUE)
  geno <- find_and_get_genotypes(conn, ph$accessions,
                                 protocol_id = "66",      # Oat 3K array (small)
                                 pedigree_dir = NULL, refresh = TRUE)
  expect_true(isSymmetric(round(geno$G, 8)))
  expect_gt(nrow(geno$G), 0)
})
