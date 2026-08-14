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

test_that("stage1_blues drops a study x trait with a constant response", {
  # All values identical -> no estimable genotypic signal -> would fit perfectly
  # (SE = 0, weight = Inf). Must be skipped, not emit infinite-weight BLUEs.
  d <- tidyr::crossing(germplasmName = paste0("g", 1:6), rep = 1:2) |>
    dplyr::mutate(studyDbId = "S1", trait = "survival", block = rep,
                  row = NA, col = NA, value = 100)
  b <- expect_message(stage1_blues(d), "constant")
  expect_equal(nrow(b), 0)
})

test_that("stage1_blues floors SE so duplicated-record trials get finite weights", {
  # Each genotype's records are identical duplicates but genotypes DIFFER: between-
  # genotype signal exists, yet residual var -> 0 so the raw SE -> 0 (a perfect fit).
  # The floor must keep weights finite and bounded while preserving the signal.
  base <- setNames(seq(2, 12, length.out = 6), paste0("g", 1:6))
  d <- tidyr::crossing(germplasmName = names(base), dup = 1:3) |>   # 3 identical rows
    dplyr::mutate(studyDbId = "S1", trait = "yield", rep = 1, block = 1,
                  row = NA, col = NA, value = base[germplasmName])
  b <- suppressWarnings(stage1_blues(d))                  # perfect-fit warning expected
  expect_equal(nrow(b), 6)
  expect_true(all(is.finite(b$weight) & b$weight > 0))
  expect_true(all(b$SE >= SE_FLOOR_FRAC * sd(d$value) - 1e-9))   # floor respected
  expect_gt(cor(b$BLUE[match(names(base), b$genotype)], base), 0.99)  # signal preserved
})

test_that("stage1_blues drops a redundant design factor when rep == block", {
  # rep and block are the SAME partition (labels coincide). Entering both aliases
  # their variance components and lme4 warns "nearly unidentifiable"; the guard must
  # keep only one term so the fit stays identifiable (no such warning).
  set.seed(1)
  base    <- setNames(seq(10, 0, length.out = 8), paste0("g", 1:8))
  rep_eff <- c(`1` = 0, `2` = 0.8, `3` = -0.7)
  d <- tidyr::crossing(germplasmName = names(base), rep = 1:3) |>
    dplyr::mutate(studyDbId = "S1", trait = "yield", block = rep,   # block == rep
                  row = NA, col = NA,
                  value = base[germplasmName] + rep_eff[as.character(rep)] +
                    rnorm(dplyr::n(), 0, 0.3))
  b <- expect_no_warning(stage1_blues(d))
  expect_setequal(b$genotype, names(base))
  expect_true(all(is.finite(b$weight) & b$weight > 0))
  expect_gt(cor(b$BLUE[match(names(base), b$genotype)], base), 0.9)
})
