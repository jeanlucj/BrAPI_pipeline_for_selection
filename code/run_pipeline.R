# Interactive driver for the whole selection pipeline (steps 1..7), mirroring
# analysis/brapi_selection_pipeline.Rmd but meant to be run in a live R console.
#
# Why this exists: wflow_build() renders the Rmd in a non-interactive callr
# subprocess with knitr capturing each chunk's output, so the progress bars and
# status lines (see code/progress.R) never get a live terminal to draw to.
# Sourcing this in an interactive session does, so you can watch the slow
# download/fit steps progress and see how long each step took.
#
# Usage (interactive console, project root):
#   source(here::here("code", "run_pipeline.R"))
# Every step caches to data/*.rds, so re-sourcing is fast: each cached step says
# so and returns immediately. To force fresh runs, set this in config.R:
#   PIPELINE_REFRESH <- TRUE
# Reporting follows SHOW_PROGRESS (config.R, default interactive()); override with
# options(brapi.progress = TRUE/FALSE).
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

# Step banners are numbered to match the step FILES (01_connect .. 07_select),
# with cross-validation as 6b since it is part of the Stage-2 file.
timings_reset()

h <- step_start("1 Connect")
conn <- connect_t3()
step_done(h, DB_NAME)

h <- step_start("2 Trials")
trials <- find_ny_trials(conn, refresh = refresh)
step_done(h, sum(trials$role == "training"), " training + ",
          sum(trials$role == "test"), " test trials across ",
          n_distinct(trials$locationName), " locations")

h <- step_start("3 Phenotypes")
pheno <- get_phenotypes(conn, trials$studyDbId, refresh = refresh)
sets  <- split_by_role(pheno, trials)   # train_pheno, train_acc, test_acc
step_done(h, nrow(pheno$pheno), " observations | ", nrow(sets$train_acc),
          " training + ", nrow(sets$test_acc), " test-trial accessions")

h <- step_start("4 Genotyping")
geno <- find_and_get_genotypes(conn, sets$train_acc, sets$test_acc,
                               TEST_ACCESSIONS, refresh = refresh)
step_done(h, "protocols ", paste(geno$protocol_ids, collapse = ", "),
          " | prediction GRM ", nrow(geno$G), " x ", ncol(geno$G))

h <- step_start("5 Stage-1 BLUEs")
blues <- stage1_blues(sets$train_pheno)
step_done(h, nrow(blues), " BLUEs across ", n_distinct(blues$studyDbId), " studies")

h <- step_start("6 Stage-2 prediction")
gebv <- stage2_gblup(blues, geno, refresh = refresh)
step_done(h, nrow(gebv), " GEBVs across ", n_distinct(gebv$trait), " traits")

h <- step_start("6b Cross-validation")
cv <- map_dfr(unique(gebv$trait), ~ cv_accuracy(blues, geno, .x),
              .progress = pb_wrap("Cross-validation: traits"))
saveRDS(cv, output_path("trait_cvs.rds"))
step_done(h, nrow(cv), " traits cross-validated")
print(cv)

h <- step_start("7 Selection")
selected <- select_parents(gebv, cv = cv)
step_done(h, nrow(selected$parents), " accessions -> ",
          output_path("breeders_output.csv"))

print_timings()

invisible(list(trials = trials, pheno = pheno, sets = sets, geno = geno,
               blues = blues, gebv = gebv, cv = cv, selected = selected))
