# Tier 1 (offline, opt-in): verify our reconstructed Wishart-EM combiner
# (EMCovarianceCombiner) against the reference implementation CovCombR::CovComb
# (Akdemir). Skipped unless CovCombR is installed
# (remotes::install_github("cran/CovCombR")).
#
# CovCombR is the slow-but-authoritative implementation (it can compute a standard
# error for every cell); we reconstructed a trimmed point-estimate-only version, so
# this guards that the reconstruction stays faithful.

test_that("EMCovarianceCombiner matches CovCombR on consistent overlapping partials", {
  skip_if_not_installed("CovCombR")

  set.seed(1)
  k <- 6; nm <- paste0("v", 1:k)
  M <- matrix(rnorm(k * k), k)
  S <- crossprod(M) + diag(k); dimnames(S) <- list(nm, nm)
  S <- S / mean(diag(S))
  A <- S[1:4, 1:4]          # vars v1..v4
  B <- S[3:6, 3:6]          # vars v3..v6  (overlap on v3, v4)

  combined <- union(colnames(A), colnames(B))
  vi <- list(match(colnames(A), combined), match(colnames(B), combined))
  invisible(capture.output(
    ours <- EMCovarianceCombiner(list(A, B), vi, degrees_freedom = c(1000, 1000),
                                 max_iter = 500, tol = 1e-8)$psi[-1, -1]))
  dimnames(ours) <- list(combined, combined)
  ref <- CovCombR::CovComb(list(A, B), nu = 1000, w = 1, lambda = 1, maxiter = 500)

  # Compare on a common scale (both standardized to mean diagonal 1), aligned by name.
  std <- function(G) (G / mean(diag(G)))[nm, nm]
  expect_lt(max(abs(std(ours) - std(ref))), 0.02)

  # Both must impute ~0 for the never-co-observed, unbridged block (v1:v2 x v5:v6):
  # related only through v3/v4, the EM (correctly) shrinks the unseen block toward 0.
  expect_lt(max(abs(std(ours)[1:2, 5:6])), 0.2)
})
