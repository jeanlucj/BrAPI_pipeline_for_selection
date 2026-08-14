# Part 3: the glmnet-vs-mean A/B done with a ROBUST glmnet imputation.
#
# Why this exists: the pasted .impute_glmnet wraps its whole per-column sapply in
# try(silent = TRUE). Near-monomorphic markers (a single minority genotype among
# ~17 observed) make cv.glmnet throw "y is constant" when a CV training fold is
# all one value. The FIRST such column aborts the entire sapply, so the function
# silently returns a purely mean-imputed matrix -> glmnet never actually runs on
# these partials. This version retries each column independently and falls back to
# the column mean only for the columns that genuinely fail, so we can measure what
# glmnet imputation would actually do to the partial-GRM diagonals.

library(tidyverse)
here::i_am("code/impute_diagnostic_robust.R")
source(here::here("code", "config.R"))
source(here::here("code", "grm_utils.R"))

set.seed(1)

# Robust per-column glmnet imputation (same model as .impute_glmnet, but a failing
# column falls back to its mean instead of aborting the whole matrix). Reports how
# many columns actually got a glmnet fit vs fell back.
impute_glmnet_robust <- function(matNA) {
  cvLambda <- exp(-(2:11))
  matNoNA <- apply(matNA, 2, function(v){ v[is.na(v)] <- mean(v, na.rm = TRUE); v })
  nPred <- min(60, round(ncol(matNA) * 0.5))
  sumIsNA <- apply(matNA, 2, function(v) sum(is.na(v)))
  order_k <- order(sumIsNA); order_k <- order_k[sumIsNA[order_k] > 0]
  n_glm <- 0L; n_fallback <- 0L
  for (k in order_k) {
    isNA <- is.na(matNA[, k])
    if (sd(matNA[, k], na.rm = TRUE) == 0) {
      matNoNA[isNA, k] <- matNA[which(!isNA)[1], k]; next
    }
    varRange <- range(matNA[, k], na.rm = TRUE)
    ok <- tryCatch({
      corMrk  <- abs(cov(matNA[, k], matNA, use = "pairwise.complete.obs"))
      predMrk <- setdiff(order(corMrk, decreasing = TRUE)[1:nPred], k)
      cv   <- glmnet::cv.glmnet(x = matNoNA[!isNA, predMrk], y = matNA[!isNA, k],
                                nfolds = 5, lambda = cvLambda)
      pred <- predict(cv, s = "lambda.min", newx = matNoNA[isNA, predMrk, drop = FALSE])
      pred[pred < varRange[1]] <- varRange[1]; pred[pred > varRange[2]] <- varRange[2]
      matNoNA[isNA, k] <- pred
      TRUE
    }, error = function(e) FALSE)
    if (ok) n_glm <- n_glm + 1L else n_fallback <- n_fallback + 1L
  }
  attr(matNoNA, "n_glm") <- n_glm
  attr(matNoNA, "n_fallback") <- n_fallback
  matNoNA
}

.mean_impute <- function(M) {
  cm <- colMeans(M, na.rm = TRUE); na <- which(is.na(M), arr.ind = TRUE)
  if (nrow(na) > 0) M[na] <- cm[na[, "col"]]; M
}

dosage_by_proto <- read_rds(cache_path("impute_diag_dosage.rds"))

rows <- list(); diags <- list()
for (pid in names(dosage_by_proto)) {
  M <- dosage_by_proto[[pid]]
  d_mean <- diag(std_grm(.Gmatrix(.mean_impute(M))))
  t0 <- Sys.time(); Mg <- impute_glmnet_robust(M)
  dt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  d_glm <- diag(std_grm(.Gmatrix(Mg)))
  message(sprintf("proto %s: %d acc x %d mrk | glmnet fit %d cols, fell back %d | %.1fs",
                  pid, nrow(M), ncol(M), attr(Mg, "n_glm"), attr(Mg, "n_fallback"), dt))
  rows[[pid]] <- tibble(
    proto = pid, n = nrow(M),
    cols_glmnet = attr(Mg, "n_glm"), cols_fallback = attr(Mg, "n_fallback"),
    min_mean = min(d_mean),  min_glm = min(d_glm),
    n_lt0.8_mean = sum(d_mean < 0.8), n_lt0.8_glm = sum(d_glm < 0.8),
    n_lt0.5_mean = sum(d_mean < 0.5), n_lt0.5_glm = sum(d_glm < 0.5),
    mean_diag_shift = mean(d_glm - d_mean), max_abs_shift = max(abs(d_glm - d_mean)))
  diags[[pid]] <- tibble(proto = pid, accession = rownames(M),
                         mean = d_mean, glmnet = d_glm)
}

tab <- bind_rows(rows)
cat("\n=== Robust glmnet vs mean: partial-GRM diagonals ===\n")
print(as.data.frame(tab), digits = 3)
write_rds(list(summary = tab, diags = bind_rows(diags)),
          cache_path("impute_diag_robust.rds"))
