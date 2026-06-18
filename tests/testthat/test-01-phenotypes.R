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
                       trait_patterns = c("yield"), refresh = TRUE)
  expect_true(all(stringr::str_detect(ph$pheno$trait, "yield")))   # height filtered out
  expect_false(anyNA(ph$pheno$value))                              # NA obs dropped
  expect_equal(nrow(ph$pheno), 3)                                  # ou1,ou2,ou3 (ou4 NA dropped)
  expect_false("ou4" %in% ph$design$obsUnitDbId)                   # NA-only unit absent
  expect_true(all(c("rep","block","row","col","entryType") %in% names(ph$design)))
  expect_equal(ph$design$rep[ph$design$obsUnitDbId == "ou3"], "2")
  expect_setequal(ph$accessions$germplasmName, c("LINE1", "LINE2"))
})

test_that("empty trait_patterns keeps all traits", {
  ph <- get_phenotypes(make_pheno_conn(), "S1",
                       trait_patterns = character(0), refresh = TRUE)
  expect_true(any(stringr::str_detect(ph$pheno$trait, "height")))
})
