# Part 4: does estimating allele frequencies + imputation + the GRM on the FULL
# genotyped panel (then subsetting to our targets) stabilize the target diagonals,
# vs the pipeline's subset-FIRST approach?
#
# Subset-first (the pipeline): .vcf_to_dosage reads only keep_samples columns, so
# MAF, mean-impute column means and VanRaden allele frequencies are all estimated
# from a tiny panel (here 5-28). Full-panel: read ALL samples, QC + .Gmatrix on
# the whole panel, std_grm, THEN pull out the diagonals of our targets. The latter
# defines each target's self-relationship against a broad, stably-estimated
# reference rather than its few panel-mates.
#
# Offline: cached VCFs only. Runs protocols 70 and 47.

library(tidyverse)
here::i_am("code/impute_diagnostic_fullpanel.R")
source(here::here("code", "04_find_genotyping.R"))   # .vcf_to_dosage/.merge_dosage/.qc_markers

vdir  <- cache_path("vcf_cache")
files <- list(
  "70" = file.path(vdir, "proto70_proj6920.vcf"),
  "47" = c(file.path(vdir, "proto47_proj4351.vcf"),
           file.path(vdir, "proto47_proj4352.vcf")))

# Targets per protocol = rownames of the subset-first dosage from impute_diagnostic.R
dosage_sub <- read_rds(cache_path("impute_diag_dosage.rds"))
# Subset-first diagonals (mean-impute) computed earlier
sub_diags  <- read_rds(cache_path("impute_diag_robust.rds"))$diags

summarise_diag <- function(d, label) {
  tibble(view = label, n = length(d), min = min(d), q05 = quantile(d, .05),
         median = median(d), mean = mean(d),
         n_lt0.8 = sum(d < 0.8), n_lt0.5 = sum(d < 0.5))
}

report <- list(); paired <- list()
for (pid in names(files)) {
  targets <- rownames(dosage_sub[[pid]])

  # FULL panel: read every sample, merge projects, QC + GRM on the whole panel.
  mats <- lapply(files[[pid]], function(f) .vcf_to_dosage(f, keep_samples = NULL))
  mats <- mats[!map_lgl(mats, is.null)]
  M_full <- .qc_markers(.merge_dosage(mats))           # full-panel QC + mean-impute
  G_full <- std_grm(.Gmatrix(M_full))
  message(sprintf("proto %s: full panel %d acc x %d markers", pid,
                  nrow(M_full), ncol(M_full)))

  d_full_all <- diag(G_full)
  hit <- intersect(targets, rownames(G_full))
  d_full_tgt <- d_full_all[hit]

  d_sub <- sub_diags |> filter(proto == pid) |>
    { \(x) set_names(x$mean, x$accession) }()
  d_sub_tgt <- d_sub[hit]

  report[[pid]] <- bind_rows(
    summarise_diag(d_full_all,  sprintf("p%s full-panel (all %d)", pid, length(d_full_all))),
    summarise_diag(d_full_tgt,  sprintf("p%s full-panel (targets)", pid)),
    summarise_diag(d_sub_tgt,   sprintf("p%s subset-first (targets)", pid)))

  paired[[pid]] <- tibble(proto = pid, accession = hit,
                          subset_first = d_sub_tgt, full_panel = d_full_tgt) |>
    arrange(subset_first)
}

cat("\n=== Target diagonals: subset-first vs full-panel ===\n")
print(as.data.frame(bind_rows(report)), digits = 3)

cat("\n=== Per-target paired diagonals (lowest subset-first first) ===\n")
walk(paired, ~ { cat("\n-- proto", .x$proto[1], "--\n")
                 print(as.data.frame(select(.x, -proto)), digits = 3, row.names = FALSE) })

write_rds(list(report = bind_rows(report), paired = bind_rows(paired)),
          cache_path("impute_diag_fullpanel.rds"))
