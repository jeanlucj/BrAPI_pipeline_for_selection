# Tier 1 (offline): GRM utilities and the Wishart-EM covariance combiner.

test_that(".Gmatrix is square, symmetric, with a positive diagonal", {
  M <- make_dosage(15, 40)
  G <- .Gmatrix(M)
  expect_equal(dim(G), c(15, 15))
  expect_true(isSymmetric(unname(G)))
  expect_true(all(diag(G) > 0))
  expect_equal(rownames(G), rownames(M))
})

test_that("std_grm sets the mean diagonal to 1", {
  G <- .Gmatrix(make_dosage(12, 30))
  expect_equal(mean(diag(std_grm(G))), 1)
})

test_that("EMCovarianceCombiner stitches two overlapping GRMs", {
  g1 <- std_grm(.Gmatrix(make_dosage(8, 60, seed = 1, prefix = "")))
  g2 <- std_grm(.Gmatrix(make_dosage(8, 60, seed = 2, prefix = "")))
  rownames(g1) <- colnames(g1) <- c("A","B","C","D","E","F","G","H")  # 8
  rownames(g2) <- colnames(g2) <- c("E","F","G","H","I","J","K","L")  # overlap E..H

  combined <- unique(c(colnames(g1), colnames(g2)))                   # 12 accessions
  vi <- list(match(colnames(g1), combined), match(colnames(g2), combined))
  res <- EMCovarianceCombiner(list(g1, g2), vi, degrees_freedom = c(8, 8))

  G <- res$psi[-1, -1]
  expect_equal(dim(G), c(12, 12))                 # union, phantom dropped
  expect_true(isSymmetric(round(G, 8)))
  expect_true(all(is.finite(G)))
  # likelihood should not decrease over EM iterations
  lp <- res$loglik_path
  expect_gte(lp[length(lp)], lp[1] - 1e-6)
})

test_that("EMCovarianceCombiner validates its inputs", {
  expect_error(EMCovarianceCombiner(list(), list()), "Empty")
  expect_error(EMCovarianceCombiner(list(diag(2)), list()), "Empty|Mismatch")
})

test_that(".effective_n is bounded, scale-invariant, and tracks the spectrum", {
  expect_equal(.effective_n(diag(5)), 5)                 # k equal eigenvalues -> k
  rank1 <- tcrossprod(1:6); diag(rank1) <- diag(rank1) + 1e-9
  expect_lt(.effective_n(rank1), 1.05)                   # one dominant eigenvalue -> ~1
  G <- crossprod(matrix(rnorm(60), 6, 10))
  expect_equal(.effective_n(G), .effective_n(7 * G))     # scale-invariant
  expect_lte(.effective_n(G), nrow(G))                   # bounded by rank
})
