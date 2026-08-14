# Sanity check: compare our EMCovarianceCombiner to CovCombR::CovComb on a SMALL
# real subset of the project's partial GRMs (from data/df_grid/partials.rds), to
# confirm the reconstruction agrees on realistic (not just synthetic) inputs.
#
# Kept small on purpose: CovCombR is slow on large inputs. We restrict each partial
# to a handful of shared accessions so both combiners run quickly.
#
# Usage: source(here::here("code", "covcomb_check.R"))

library(tidyverse)
here::i_am("code/covcomb_check.R")
source(here::here("code", "em_covariance_combiner.R"))
stopifnot(requireNamespace("CovCombR", quietly = TRUE))

p <- read_rds(here::here("data", "df_grid", "partials.rds"))
grms <- p$proto_grms

# Pick the ~40 accessions covered by the most partials, then take the 2-3 partials
# that overlap them best, each restricted to those names (keeps it small + connected).
allnames <- table(unlist(lapply(grms, rownames)))
hot <- names(sort(allnames, decreasing = TRUE))[seq_len(min(40, length(allnames)))]
sub <- grms |>
  map(~ { keep <- intersect(rownames(.x), hot); if (length(keep) >= 5) .x[keep, keep] else NULL }) |>
  compact()
sub <- sub[order(map_int(sub, nrow), decreasing = TRUE)][seq_len(min(3, length(sub)))]
cat("using", length(sub), "partials of sizes", paste(map_int(sub, nrow), collapse = ", "), "\n")

combined <- reduce(map(sub, colnames), union)
vi <- map(sub, ~ match(colnames(.x), combined))
invisible(capture.output(
  ours <- EMCovarianceCombiner(sub, vi, degrees_freedom = rep(1000, length(sub)),
                               max_iter = 500, tol = 1e-8)$psi[-1, -1]))
dimnames(ours) <- list(combined, combined)
ref <- CovCombR::CovComb(sub, nu = 1000, w = 1, lambda = 1, maxiter = 500)

std <- function(G) (G / mean(diag(G)))[combined, combined]
o <- std(ours); r <- std(ref)
cat("max abs diff (std to mean-diag 1):", round(max(abs(o - r)), 4), "\n")
cat("ours      : min eig", round(min(eigen(o, only.values = TRUE)$values), 4),
    " diag[min", round(min(diag(o)), 3), "max", round(max(diag(o)), 3), "]\n")
cat("CovCombR  : min eig", round(min(eigen(r, only.values = TRUE)$values), 4),
    " diag[min", round(min(diag(r)), 3), "max", round(max(diag(r)), 3), "]\n")
cat("(CovCombR applies nearPD each iteration; ours does not.)\n")
