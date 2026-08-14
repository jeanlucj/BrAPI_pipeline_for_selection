# Diagnostic: where do training accessions get lost between genotyping COVERAGE
# (counted by germplasmDbId, server-side) and MEMBERSHIP in the combined GRM
# (keyed by germplasmName, after VCF download + marker QC)?
#
# For each used protocol it reports the funnel
#   n_covered (dbId)  ->  n_target_in_vcf (name)  ->  n_target_in_grm (post-QC)
# and characterizes the name mismatch (composite "NAME|ID" samples, and how many
# more targets would match if VCF samples were split on "|"). Writes
#   data/match_diagnostic/funnel.csv
#   data/match_diagnostic/unmatched_targets.csv   (targets in no VCF and no pedigree)
#
# Reuses the partials saved by code/df_grid_diagnostic.R (data/df_grid/partials.rds);
# run that first if it is absent.
#
# Usage: source(here::here("code", "match_diagnostic.R"))

library(tidyverse)
here::i_am("code/match_diagnostic.R")
suppressMessages(library(httr))
purrr::walk(c("config.R", "grm_utils.R", "04_find_genotyping.R"),
            ~ source(here::here("code", .x)))

parts_cache <- cache_path("df_grid", "partials.rds")
if (!file.exists(parts_cache))
  stop("Run code/df_grid_diagnostic.R first to build ", parts_cache)
p  <- read_rds(parts_cache)
tn <- p$target_names
cat(length(tn), "target (training) accessions.\n")

out_dir <- cache_path("match_diagnostic")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

vcf_dir   <- cache_path("vcf_cache")
first_fld <- function(x) sub("\\|.*$", "", x)   # text before the first "|"

# Per-protocol funnel.
funnel <- map_dfr(p$used, function(pid) {
  files <- list.files(vcf_dir, pattern = sprintf("^proto%s_", pid), full.names = TRUE)
  samp  <- unique(unlist(lapply(files, .read_vcf_samples)))
  cov   <- p$protocols$n_covered[p$protocols$dbId == pid]
  grm   <- p$proto_grms[[pid]]
  in_vcf_exact <- sum(tn %in% samp)
  in_vcf_split <- sum(tn %in% unique(first_fld(samp)))   # if we split VCF names on "|"
  tibble(
    protocol      = pid,
    name          = p$protocols$name[p$protocols$dbId == pid],
    n_covered_dbid = if (length(cov)) cov else NA_integer_,
    n_target_in_vcf_exact = in_vcf_exact,
    n_target_in_vcf_split = in_vcf_split,
    recoverable_by_split  = in_vcf_split - in_vcf_exact,
    n_target_in_grm = sum(tn %in% rownames(grm)),
    vcf_samples     = length(samp),
    vcf_composite   = sum(grepl("|", samp, fixed = TRUE)))
})
funnel <- arrange(funnel, desc(n_covered_dbid))

# Overall: which targets reach the GRM at all (any marker partial OR pedigree)?
marker_names <- unique(unlist(lapply(p$proto_grms, rownames)))
ped_names    <- unique(unlist(lapply(p$ped_partials, rownames)))
in_any_grm   <- tn %in% union(marker_names, ped_names)

# Of the unmatched, how many would a "|"-split recover across all protocols?
all_split <- unique(first_fld(unique(unlist(
  lapply(p$used, function(pid)
    unlist(lapply(list.files(vcf_dir, sprintf("^proto%s_", pid), full.names = TRUE),
                  .read_vcf_samples)))))))
unmatched <- tn[!in_any_grm]
recoverable_split <- sum(unmatched %in% all_split)

unmatched_tbl <- tibble(germplasmName = unmatched,
                        recoverable_by_split = unmatched %in% all_split)
readr::write_csv(funnel, file.path(out_dir, "funnel.csv"))
readr::write_csv(unmatched_tbl, file.path(out_dir, "unmatched_targets.csv"))

cat("\n=== per-protocol funnel (covered dbId -> in VCF -> in GRM) ===\n")
print(as.data.frame(funnel), row.names = FALSE)
cat(sprintf("\nTargets in SOME partial (marker or pedigree): %d / %d  (so %d injected).\n",
            sum(in_any_grm), length(tn), sum(!in_any_grm)))
cat(sprintf("Of the %d unmatched, %d would match if VCF sample names were split on '|'.\n",
            length(unmatched), recoverable_split))
cat("Wrote ", file.path(out_dir, "funnel.csv"), " and unmatched_targets.csv\n")
invisible(funnel)
