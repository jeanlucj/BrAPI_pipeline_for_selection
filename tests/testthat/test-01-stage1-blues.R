# Tier 1 (offline): Stage-1 per-trial BLUEs.

test_that("stage1_blues recovers genotype effects in a replicated trial", {
  set.seed(10)
  g_eff   <- setNames(seq(10, 0, length.out = 8), paste0("g", 1:8))
  rep_eff <- c(`1` = 0, `2` = 0.6, `3` = -0.5)
  blk_eff <- c(`1` = 0.4, `2` = -0.4)
  # 8 genotypes x 3 reps, 2 blocks (each 4 genotypes), with genuine rep/block
  # variance so the mixed model is well posed (no singular-fit warnings).
  d <- tidyr::crossing(germplasmName = names(g_eff), rep = 1:3) |>
    dplyr::mutate(
      studyDbId = "S1", trait = "yield",
      block = (as.integer(factor(germplasmName)) - 1L) %/% 4L + 1L,
      row = NA, col = NA,
      value = g_eff[germplasmName] + rep_eff[as.character(rep)] +
        blk_eff[as.character(block)] + rnorm(dplyr::n(), 0, 0.3))

  b <- stage1_blues(d)
  expect_true(all(c("trait","studyDbId","genotype","BLUE","SE","weight") %in% names(b)))
  expect_setequal(b$genotype, names(g_eff))
  expect_true(all(is.finite(b$weight) & b$weight > 0))
  # BLUEs track the true effects
  ord <- b$BLUE[match(names(g_eff), b$genotype)]
  expect_gt(cor(ord, g_eff), 0.9)
})

test_that("stage1_blues falls back gracefully for an unreplicated trial", {
  d <- tibble::tibble(
    studyDbId = "S2", trait = "yield",
    germplasmName = paste0("g", 1:5), rep = 1, block = 1,
    row = NA, col = NA, value = c(5, 4, 3, 2, 1))
  b <- stage1_blues(d)
  expect_equal(nrow(b), 5)
  expect_true(all(is.finite(b$weight) & b$weight > 0))   # default/median weight
})
