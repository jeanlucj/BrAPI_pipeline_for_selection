# Shared genomic-relationship-matrix helpers, sourced by step 4 (build per-source
# GRMs to EM-combine) and step 6 (build a kernel from a marker matrix when no
# combined GRM is supplied).

# VanRaden genomic relationship matrix from a 0/1/2 dosage matrix.
.Gmatrix <- function(M) {
  p <- colMeans(M) / 2
  W <- sweep(M, 2, 2 * p, "-")
  G <- tcrossprod(W) / (2 * sum(p * (1 - p)))
  G + diag(1e-6, nrow(G))   # nudge to positive-definite
}

# Standardize a relationship matrix so its mean diagonal is 1 (the form the EM
# combiner expects for its partial covariances).
std_grm <- function(G) G / mean(diag(G))

# Effective number of independent samples behind a relationship matrix, via
# Galwey (2009)'s "effective number of independent tests": (sum sqrt(lambda))^2 /
# sum(lambda) over the positive eigenvalues. Scale-invariant (robust to std_grm
# scaling), bounded [1, rank]. Used to set each marker GRM's EM degrees of freedom
# (markers are in LD, so a raw count overstates the independent information).
.effective_n <- function(G) {
  lambda <- eigen(G, symmetric = TRUE, only.values = TRUE)$values
  lambda <- lambda[lambda > 1e-8]                    # drop ~zero / tiny-negative
  if (length(lambda) < 1) return(1)
  (sum(sqrt(lambda)))^2 / sum(lambda)
}

# Mean imputation: fill each missing call with its marker (column) mean.
.mean_impute <- function(M) {
  cm <- colMeans(M, na.rm = TRUE)
  na <- which(is.na(M), arr.ind = TRUE)
  if (nrow(na) > 0) M[na] <- cm[na[, "col"]]
  M
}

#################################################
## This function uses glmnet to impute NA values
#################################################

# NOTES (original design, jeanlucj)
# You enter the function with a matrix, individuals in rows, markers in columns.
# I have only tested it with markers coded as -1, 0, 1 for AA, AB, and BB, but
# other codings (e.g. the pipeline's 0/1/2) work too.
# Returns a matrix of the same dimensions, but with no missing data. Imputed
# values are real numbers (not integers). This may be problematic for downstream
# mapping software (it is fine for a VanRaden GRM).
# glmnet would look at ~100 lambda penalties by default; using 10 ~doubles speed
# for a couple percent accuracy cost. 5-fold (vs 10-fold) cv ~doubles speed again.
# The biggest speed win is using only the top nPred (60) most-correlated markers
# as predictors rather than all of them.
#
# ROBUSTNESS / DETERMINISM (added when wiring this into the pipeline; see
# code/impute_diagnostic*.R for the analysis that motivated it):
#  - Each marker is imputed independently inside tryCatch: a column whose
#    cv.glmnet fails (commonly "y is constant" when a near-monomorphic marker's
#    CV training fold is all one genotype) falls back to its mean instead of
#    aborting the whole matrix. The previous `try(sapply(...))` aborted on the
#    FIRST such column, silently returning a purely mean-imputed matrix.
#  - CV folds are assigned deterministically (foldid = 1:nfolds recycled) so the
#    imputed matrix -- and hence the GRM -- is reproducible across runs.
#  - Predictor selection (the top-nPred most-associated markers) uses |Pearson
#    correlation| computed in ONE BLAS pass on the standardized mean-imputed matrix
#    (`crossprod`), rather than R's per-column cov(use="pairwise.complete.obs").
#    The original pairwise cov is O(n*m) PER imputed column -> O(n*m^2) overall and
#    dominated runtime on dense panels (e.g. ~40 min for a 17k-marker protocol);
#    the standardized matrix is built once and each column's correlations are a
#    single matrix-vector product. (Correlation rather than covariance also
#    normalizes out marker-variance differences, which is the more sensible basis
#    for "most-associated"; values are mean-imputed, not pairwise-complete.)
# This is meant to run on a reasonably large panel (see GRM_PANEL_MIN): with only
# a handful of individuals there is too little inter-individual signal for the
# regressions to beat the mean.
.impute_glmnet <- function(matNA, nfolds = 5, n_pred = 60, label = NULL) {
  cvLambda <- exp(-(2:11))
  matNoNA  <- .mean_impute(matNA)                       # baseline / per-column fallback
  nPred    <- min(n_pred, round(ncol(matNA) * 0.5))
  sumIsNA  <- apply(matNA, 2, function(v) sum(is.na(v)))
  imputeOrder <- order(sumIsNA)
  imputeOrder <- imputeOrder[sumIsNA[imputeOrder] > 0]  # skip complete columns

  # Standardize the mean-imputed matrix once so predictor correlations are a fast
  # BLAS crossprod below. Constant columns (e.g. monomorphic markers) are zeroed so
  # they contribute zero correlation (never selected) instead of NaN.
  n   <- nrow(matNoNA)
  ctr <- colMeans(matNoNA)
  Zc  <- sweep(matNoNA, 2, ctr, "-")
  sdv <- sqrt(colSums(Zc^2) / max(n - 1, 1))
  Zsc <- sweep(Zc, 2, ifelse(sdv > 0, sdv, 1), "/")
  Zsc[, sdv == 0] <- 0

  # One cv.glmnet per incomplete marker: the slowest loop in the pipeline (tens of
  # minutes on a dense panel), so it gets a bar. Ticks are batched -- redrawing on
  # every one of ~17k iterations would itself cost measurable time.
  n_glm <- 0L; n_fallback <- 0L
  pb <- pb_start(length(imputeOrder),
                 paste0(if (!is.null(label)) paste0(label, ": ") else "",
                        "imputing markers"))
  on.exit(pb_done(pb), add = TRUE)
  tick_every <- max(1L, length(imputeOrder) %/% 200L)
  done <- 0L
  for (k in imputeOrder) {
    done <- done + 1L
    if (done %% tick_every == 0L) pb_tick(pb, inc = tick_every)
    isNA <- is.na(matNA[, k])
    # Monomorphic among the observed: impute the sole observed value.
    if (sd(matNA[, k], na.rm = TRUE) == 0) {
      matNoNA[isNA, k] <- matNA[which(!isNA)[1], k]
      next
    }
    n_obs    <- sum(!isNA)
    foldid   <- rep_len(seq_len(min(nfolds, n_obs)), n_obs)   # deterministic folds
    varRange <- range(matNA[, k], na.rm = TRUE)               # clamp to observed range
    ok <- tryCatch({
      corMrk  <- abs(as.numeric(crossprod(Zsc[, k], Zsc)) / max(n - 1, 1))
      predMrk <- setdiff(order(corMrk, decreasing = TRUE)[1:nPred], k)
      cv   <- glmnet::cv.glmnet(x = matNoNA[!isNA, predMrk], y = matNA[!isNA, k],
                                foldid = foldid, lambda = cvLambda)
      pred <- predict(cv, s = "lambda.min", newx = matNoNA[isNA, predMrk, drop = FALSE])
      pred[pred < varRange[1]] <- varRange[1]
      pred[pred > varRange[2]] <- varRange[2]
      matNoNA[isNA, k] <- pred
      TRUE
    }, error = function(e) FALSE)
    if (ok) n_glm <- n_glm + 1L else n_fallback <- n_fallback + 1L
  }
  if (n_fallback > 0) {
    message(sprintf("  glmnet impute: %d markers fit, %d fell back to mean ",
                    n_glm, n_fallback),
            "(near-monomorphic / degenerate CV).")
  }
  matNoNA
}
