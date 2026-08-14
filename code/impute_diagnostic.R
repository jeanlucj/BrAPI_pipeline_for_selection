# Diagnostic: does glmnet imputation reduce the frequency of low GRM diagonal
# elements relative to mean imputation in the per-protocol (partial) GRMs?
#
# Hypothesis: mean imputation pulls a heavily-missing accession toward the
# population mean, zeroing its centered marker row (W) and so deflating its
# VanRaden self-relationship (the GRM diagonal). glmnet imputation borrows
# information from correlated markers, which should preserve more of that
# accession's individuality and lift its diagonal.
#
# This is a CONTROLLED A/B: for each protocol we reconstruct the QC'd dosage
# matrix WITH its NAs intact (over the final prediction targets), then impute it
# two ways and build std_grm(.Gmatrix(.)) from each. The marker set, accession
# set and NA pattern are identical across methods, so any diagonal shift is due
# to imputation alone.
#
# Offline: reads cached VCFs under data/vcf_cache and the cached genotypes.rds;
# does not contact the server.

library(tidyverse)
here::i_am("code/impute_diagnostic.R")
source(here::here("code", "04_find_genotyping.R"))   # also brings grm_utils (.impute_glmnet)

set.seed(1)

# --- QC without imputation: replicate .qc_markers up to (not including) the
# mean-impute step, so we keep the NAs for the A/B. -----------------------------
.qc_markers_keepNA <- function(M, max_missing = MAX_MISSING, min_maf = MIN_MAF) {
  M <- M[, colMeans(is.na(M)) <= max_missing, drop = FALSE]
  M <- M[rowMeans(is.na(M)) <= max_missing, , drop = FALSE]
  af  <- colMeans(M, na.rm = TRUE) / 2
  maf <- pmin(af, 1 - af)
  M[, !is.na(maf) & maf >= min_maf, drop = FALSE]
}

.mean_impute <- function(M) {
  col_means <- colMeans(M, na.rm = TRUE)
  na_idx <- which(is.na(M), arr.ind = TRUE)
  if (nrow(na_idx) > 0) M[na_idx] <- col_means[na_idx[, "col"]]
  M
}

diag_summary <- function(d) {
  c(n = length(d), min = min(d), q05 = unname(quantile(d, .05)),
    median = median(d), mean = mean(d),
    frac_lt_0.5 = mean(d < 0.5), frac_lt_0.8 = mean(d < 0.8))
}

# --- Rebuild each protocol's QC'd dosage (NAs intact) over the targets ---------
g <- read_rds(cache_path("genotypes.rds"))
keep_samples <- rownames(g$G)                          # final prediction targets
proto_ids    <- g$protocol_ids
message("Targets: ", length(keep_samples), " accessions; protocols: ",
        paste(proto_ids, collapse = ", "))

vdir <- cache_path("vcf_cache")
files_by_proto <- set_names(lapply(proto_ids, function(pid) {
  list.files(vdir, pattern = sprintf("^proto%s_proj.*\\.vcf$", pid), full.names = TRUE)
}), proto_ids)

dosage_by_proto <- list()
for (pid in proto_ids) {
  mats <- list()
  for (f in files_by_proto[[pid]]) {
    m <- tryCatch(.vcf_to_dosage(f, keep_samples = keep_samples),
                  error = function(e) NULL)
    if (!is.null(m) && nrow(m) > 0 && ncol(m) > 0) mats[[length(mats) + 1L]] <- m
  }
  if (!length(mats)) { message("proto ", pid, ": no samples; skip"); next }
  M <- .qc_markers_keepNA(.merge_dosage(mats))
  if (nrow(M) < 2 || ncol(M) < 2) { message("proto ", pid, ": too small after QC; skip"); next }
  dosage_by_proto[[pid]] <- M
  message(sprintf("proto %s: %d acc x %d markers, %.1f%% missing",
                  pid, nrow(M), ncol(M), 100 * mean(is.na(M))))
}

write_rds(dosage_by_proto, cache_path("impute_diag_dosage.rds"))

# --- Mean-impute diagnostics (fast) --------------------------------------------
cat("\n=== Mean-imputation partial-GRM diagonals ===\n")
mean_tab <- map_dfr(names(dosage_by_proto), function(pid) {
  G <- std_grm(.Gmatrix(.mean_impute(dosage_by_proto[[pid]])))
  tibble(proto = pid, !!!as.list(diag_summary(diag(G))))
})
print(as.data.frame(mean_tab), digits = 3)
write_rds(mean_tab, cache_path("impute_diag_mean.rds"))
