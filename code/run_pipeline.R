# Interactive driver for the whole selection pipeline (steps 1..7), mirroring
# analysis/brapi_selection_pipeline.Rmd but meant to be run in a live R console.
#
# Why this exists: wflow_build() renders the Rmd in a non-interactive callr
# subprocess with knitr capturing each chunk's output, so the cli/txt progress
# bars in get_phenotypes(), find_and_get_genotypes() and the EM combiner never
# get a live terminal to draw to. Sourcing this in an interactive session does,
# so you can watch the slow download/combine steps progress.
#
# Usage (interactive console, project root):
#   source(here::here("code", "run_pipeline.R"))
# Every step caches to data/*.rds, so re-sourcing is fast and returns instantly
# (no bars). To force fresh runs and see the bars, set this in config.R:
#   PIPELINE_REFRESH <- TRUE
# All result objects (trials, pheno, geno, blues, gebv, cv, selected) are left
# in the global environment for inspection.

library(tidyverse)
here::i_am("code/run_pipeline.R")

purrr::walk(
  c("config.R", "01_connect.R", "02_find_trials.R", "03_get_phenotypes.R",
    "04_find_genotyping.R", "05_stage1_blues.R", "06_stage2_genomic_prediction.R",
    "07_select.R"),
  ~ source(here::here("code", .x))
)

# Honor an optional PIPELINE_REFRESH toggle without requiring it to exist.
refresh <- exists("PIPELINE_REFRESH") && isTRUE(PIPELINE_REFRESH)
if (refresh) message("PIPELINE_REFRESH = TRUE: re-running every cached step.")

conn <- connect_t3()

message("\n== 1. Find NY-region trials ==")
trials <- find_ny_trials(conn, refresh = refresh)
cat(sum(trials$role == "training"), "training +", sum(trials$role == "test"),
    "test trials across", n_distinct(trials$locationName), "locations\n")

message("\n== 2. Download phenotypes (progress by study) ==")
pheno <- get_phenotypes(conn, trials$studyDbId, refresh = refresh)
sets  <- split_by_role(pheno, trials)   # train_pheno, train_acc, test_acc
cat("Observations:", nrow(pheno$pheno),
    "| training accessions:", nrow(sets$train_acc),
    "| test-trial accessions:", nrow(sets$test_acc), "\n")

message("\n== 3. Genotyping coverage + combined GRM (download progress, EM bar) ==")
geno <- find_and_get_genotypes(conn, sets$train_acc, sets$test_acc,
                               TEST_ACCESSIONS, refresh = refresh)
cat("Protocols combined:", paste(geno$protocol_ids, collapse = ", "),
    "| prediction GRM:", nrow(geno$G), "x", ncol(geno$G), "accessions\n")

message("\n== 4. Stage-1 BLUEs (training trials only) ==")
blues <- stage1_blues(sets$train_pheno)
cat(nrow(blues), "BLUEs across", n_distinct(blues$studyDbId), "studies\n")

message("\n== 5. Stage-2 genomic prediction ==")
gebv <- stage2_gblup(blues, geno, refresh = refresh)
cat(nrow(gebv), "GEBVs across", n_distinct(gebv$trait), "traits\n")

message("\n== 5b. Cross-validated accuracy ==")
cv <- map_dfr(unique(gebv$trait), ~ cv_accuracy(blues, geno, .x))
print(cv)

message("\n== 6. Select parents ==")
selected <- select_parents(gebv, geno = geno, n_select = 20)
cat(sum(selected$selected), "selected; full ranking written to ",
    output_path("selected_parents.csv"), "\n", sep = "")

invisible(list(trials = trials, pheno = pheno, sets = sets, geno = geno,
               blues = blues, gebv = gebv, cv = cv, selected = selected))
