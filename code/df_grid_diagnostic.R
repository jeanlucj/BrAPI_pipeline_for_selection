# Diagnostic: how the EM-combined relationship matrix G changes with the EM
# degrees-of-freedom parameters (PEDIGREE_DF, GRM_DF_MEAN). Motivated by GBLUPs
# that are not shrunk relative to the phenotypes and a combined G with diagonals
# dipping well below 1.
#
# The expensive build (download + per-protocol GRMs + pedigree partials) is run
# ONCE via .genotype_partials(); each df setting then only re-runs the cheap
# df-weighted EM combine via .combine_to_G(). Outputs go to data/df_grid/:
#   - G_ped<P>_grm<M>.rds : the combined prediction GRM for each setting
#   - df_grid_summary.csv : diagonal / off-diagonal / eigenvalue diagnostics
#
# Usage (interactive console, project root):
#   source(here::here("code", "df_grid_diagnostic.R"))

library(tidyverse)
here::i_am("code/df_grid_diagnostic.R")
purrr::walk(
  c("config.R", "01_connect.R", "02_find_trials.R", "03_get_phenotypes.R",
    "04_find_genotyping.R"),
  ~ source(here::here("code", .x))
)

out_dir <- cache_path("df_grid")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# --- Build the df-independent partials once (cached) --------------------------
parts_cache <- file.path(out_dir, "partials.rds")
if (file.exists(parts_cache)) {
  message("Reusing cached partials (", parts_cache, ").")
  parts <- read_rds(parts_cache)
} else {
  conn   <- connect_t3()
  trials <- find_ny_trials(conn)
  pheno  <- get_phenotypes(conn, trials$studyDbId)
  sets   <- split_by_role(pheno, trials)
  parts  <- .genotype_partials(conn, sets$train_acc, sets$test_acc, TEST_ACCESSIONS)
  write_rds(parts, parts_cache)
}
message(length(parts$proto_grms), " marker GRM(s) + ", length(parts$ped_partials),
        " pedigree partial(s); ", length(parts$target_names), " prediction targets.")

# --- df grid (PEDIGREE_DF, GRM_DF_MEAN); GRM_DF_STDEV left at the config value --
grid <- tibble::tribble(
  ~pedigree_df, ~grm_df_mean,
  30,  60,    # current
  60,  60,
  100, 200,
  200, 200,
  200, 400,
  400, 400)

# Diagnostics for one combined matrix.
.summarize_G <- function(G) {
  d  <- diag(G)
  od <- G[upper.tri(G)]
  ev <- eigen(G, symmetric = TRUE, only.values = TRUE)$values
  pos <- ev[ev > 1e-10]
  # how many pairs are "more related to each other than to the smaller self"?
  patho <- sum(abs(G[upper.tri(G)]) >
                 pmin(d[row(G)[upper.tri(G)]], d[col(G)[upper.tri(G)]]))
  tibble(
    n = nrow(G),
    diag_min = min(d), diag_med = median(d), diag_mean = mean(d), diag_max = max(d),
    n_diag_lt_0.5 = sum(d < 0.5),
    offdiag_min = min(od), offdiag_med = median(od), offdiag_max = max(od),
    eig_min = min(ev), n_neg_eig = sum(ev < -1e-8),
    cond = if (length(pos)) max(pos) / min(pos) else NA_real_,
    n_offdiag_gt_diag = patho)
}

summ <- pmap_dfr(grid, function(pedigree_df, grm_df_mean) {
  message("Combining: pedigree_df = ", pedigree_df, ", grm_df_mean = ", grm_df_mean)
  G <- .combine_to_G(parts, grm_df_mean = grm_df_mean,
                     grm_df_stdev = GRM_DF_STDEV, pedigree_df = pedigree_df)$G
  write_rds(G, file.path(out_dir, sprintf("G_ped%d_grm%d.rds", pedigree_df, grm_df_mean)))
  bind_cols(tibble(pedigree_df = pedigree_df, grm_df_mean = grm_df_mean,
                   grm_df_stdev = GRM_DF_STDEV),
            .summarize_G(G))
})

readr::write_csv(summ, file.path(out_dir, "df_grid_summary.csv"))
message("\nWrote ", nrow(summ), " settings to ", out_dir)
print(as.data.frame(summ), digits = 4)
invisible(summ)
