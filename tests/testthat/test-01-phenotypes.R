# Tier 1 (offline): phenotype + design extraction via a fake connection.

make_pheno_conn <- function() {
  units <- list(
    unit_rec("ou1", "g1", "LINE1", rep = 1, block = 1, x = 1, y = 1,
             traits = list("Grain yield - g/m2|CO_350:0000260", "Plant height - cm|CO_350:0000017"),
             values = list(500, 95)),
    unit_rec("ou2", "g2", "LINE2", rep = 1, block = 2, x = 2, y = 1,
             traits = list("Grain yield - g/m2|CO_350:0000260", "Plant height - cm|CO_350:0000017"),
             values = list(420, 88)),
    unit_rec("ou3", "g1", "LINE1", rep = 2, block = 1, x = 1, y = 2,
             traits = list("Grain yield - g/m2|CO_350:0000260"),
             values = list(480)),
    unit_rec("ou4", "g2", "LINE2", rep = 2, block = 2, x = 2, y = 2,
             traits = list("Grain yield - g/m2|CO_350:0000260"),
             values = list(NA)))                       # missing value -> dropped
  fake_conn(units_by_study = list(S1 = units))
}

test_that("get_phenotypes filters traits, parses design, and lists accessions", {
  ph <- get_phenotypes(make_pheno_conn(), "S1",
                       trait_names = "Grain yield - g/m2|CO_350:0000260",
                       refresh = TRUE)
  expect_setequal(ph$pheno$trait, "Grain yield - g/m2|CO_350:0000260")  # exact match; height out
  expect_false(anyNA(ph$pheno$value))                              # NA obs dropped
  expect_equal(nrow(ph$pheno), 3)                                  # ou1,ou2,ou3 (ou4 NA dropped)
  expect_false("ou4" %in% ph$design$obsUnitDbId)                   # NA-only unit absent
  expect_true(all(c("rep","block","row","col","entryType") %in% names(ph$design)))
  expect_equal(ph$design$rep[ph$design$obsUnitDbId == "ou3"], "2")
  expect_setequal(ph$accessions$germplasmName, c("LINE1", "LINE2"))
})

test_that("empty trait_names keeps all traits", {
  ph <- get_phenotypes(make_pheno_conn(), "S1",
                       trait_names = character(0), refresh = TRUE)
  expect_true(any(stringr::str_detect(ph$pheno$trait, "height")))
})

test_that("split_by_role partitions phenotypes/accessions by trial role", {
  pheno <- list(pheno = tibble::tibble(
    studyDbId     = c("S1", "S1", "S2", "S3"),
    germplasmDbId = c("1", "2", "3", "2"),
    germplasmName = c("A", "B", "C", "B"),
    trait = "yield", value = c(1, 2, 3, 4)))
  trials <- tibble::tibble(studyDbId = c("S1", "S2", "S3"),
                           role = c("training", "test", "test"))
  s <- split_by_role(pheno, trials)
  expect_setequal(s$train_pheno$studyDbId, "S1")          # only training trial
  expect_setequal(s$train_acc$germplasmName, c("A", "B"))
  expect_setequal(s$test_acc$germplasmName, c("C", "B"))  # from S2 + S3
})

# --- per-study caching --------------------------------------------------------
# The point: adding a trial must cost one study's download, and changing the trait
# list must cost none, because the per-study files are stored unfiltered by trait.

# Each test uses its OWN study ids: the per-study cache lives in one tempdir for the
# whole session, so ids shared with another test would arrive already cached.
counting_conn <- function(studies) {
  hits <- new.env(parent = emptyenv()); hits$n <- character(0)
  units <- list(
    unit_rec("ou1", "g1", "LINE1", rep = 1, block = 1, x = 1, y = 1,
             traits = list("T1", "T2"), values = list(10, 20)),
    unit_rec("ou2", "g2", "LINE2", rep = 2, block = 1, x = 2, y = 1,
             traits = list("T1", "T2"), values = list(11, 21)))
  conn <- fake_conn(units_by_study = setNames(rep(list(units), length(studies)), studies))
  base_get <- conn$get
  conn$get <- function(call, query = NULL, ...) {
    if (grepl("/observationunits", call)) hits$n <- c(hits$n, as.character(query$studyDbId))
    base_get(call, query = query, ...)
  }
  list(conn = conn, downloaded = function() hits$n)
}

test_that("adding a study downloads only that study", {
  cc <- counting_conn(c("A1", "A2", "A3"))
  get_phenotypes(cc$conn, c("A1", "A2"), trait_names = "T1")
  expect_setequal(cc$downloaded(), c("A1", "A2"))

  get_phenotypes(cc$conn, c("A1", "A2", "A3"), trait_names = "T1")
  expect_equal(setdiff(cc$downloaded(), c("A1", "A2")), "A3")   # only the new one
})

test_that("changing the trait list re-filters without downloading anything", {
  cc <- counting_conn(c("B1", "B2"))
  ph1 <- get_phenotypes(cc$conn, c("B1", "B2"), trait_names = "T1")
  n_after_first <- length(cc$downloaded())

  ph2 <- get_phenotypes(cc$conn, c("B1", "B2"), trait_names = "T2")
  expect_length(cc$downloaded(), n_after_first)                 # no new downloads
  expect_setequal(ph1$pheno$trait, "T1")                        # and the filter changed
  expect_setequal(ph2$pheno$trait, "T2")
})

test_that("an unchanged request is served from the assembled cache", {
  cc <- counting_conn(c("C1", "C2"))
  get_phenotypes(cc$conn, c("C1", "C2"), trait_names = "T1")
  n <- length(cc$downloaded())
  get_phenotypes(cc$conn, c("C1", "C2"), trait_names = "T1")    # same request
  expect_length(cc$downloaded(), n)
})
