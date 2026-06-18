# Step 5 (Stage 1 of the two-stage analysis): per-trial adjusted means (BLUEs).
#
# For each study x trait we fit a mixed model with genotype as a FIXED effect
# (so we get BLUEs) and the field-design factors (rep, block) as RANDOM effects.
# This adjusts genotype means for spatial/design structure, including augmented
# designs where unreplicated entries are adjusted via replicated check blocks.
# Each BLUE is carried to Stage 2 with a weight = 1 / SE^2.

library(tidyverse)
here::i_am("code/05_stage1_blues.R")
source(here::here("code", "config.R"))

# Fit one study x trait and return a tibble: genotype, BLUE, SE.
.fit_one <- function(d) {
  d <- d |>
    mutate(genotype = factor(germplasmName),
           rep = factor(rep), block = factor(block))

  # Build the random-effects part from whatever design replication exists.
  rand <- character()
  if (nlevels(droplevels(d$rep))   > 1) rand <- c(rand, "(1 | rep)")
  if (nlevels(droplevels(d$block)) > 1) rand <- c(rand, "(1 | block)")

  est <- NULL
  if (length(rand) > 0) {
    form <- as.formula(paste("value ~ 0 + genotype +", paste(rand, collapse = " + ")))
    fit  <- tryCatch(
      suppressMessages(lme4::lmer(form, data = d,
                                  control = lme4::lmerControl(check.conv.singular = "ignore"))),
      error = function(e) NULL)
    if (!is.null(fit)) {
      co  <- summary(fit)$coefficients
      est <- tibble(term = rownames(co), BLUE = co[, "Estimate"], SE = co[, "Std. Error"])
    }
  }
  # Fall back to a fixed-only model (replicated) or plot means (unreplicated).
  if (is.null(est)) {
    fit <- stats::lm(value ~ 0 + genotype, data = d)
    co  <- summary(fit)$coefficients
    est <- tibble(term = rownames(co), BLUE = co[, "Estimate"],
                  SE = if ("Std. Error" %in% colnames(co)) co[, "Std. Error"] else NA_real_)
  }

  est |>
    mutate(genotype = str_remove(term, "^genotype")) |>
    select(genotype, BLUE, SE)
}

#' Stage-1 BLUEs for every study x trait in a long phenotype table.
#'
#' @param pheno long table with columns studyDbId, germplasmName, rep, block,
#'   trait, value (the `pheno` element from get_phenotypes()).
#' @return tibble: trait, studyDbId, genotype, BLUE, SE, weight (= 1/SE^2).
#'   Studies/traits that cannot be fit are skipped with a message.
stage1_blues <- function(pheno) {
  pheno |>
    filter(!is.na(value), !is.na(germplasmName)) |>
    group_by(trait, studyDbId) |>
    group_modify(function(d, key) {
      if (n_distinct(d$germplasmName) < 2) {
        message(sprintf("Skipping %s / study %s: <2 genotypes.", key$trait, key$studyDbId))
        return(tibble())
      }
      out <- tryCatch(.fit_one(d), error = function(e) {
        message(sprintf("Skipping %s / study %s: %s", key$trait, key$studyDbId,
                        conditionMessage(e)))
        tibble()
      })
      out
    }) |>
    ungroup() |>
    # weight by inverse squared SE; default weight 1 where SE is unavailable
    mutate(weight = if_else(!is.na(SE) & SE > 0, 1 / SE^2, NA_real_)) |>
    mutate(weight = if_else(is.na(weight), stats::median(weight, na.rm = TRUE), weight)) |>
    mutate(weight = if_else(is.na(weight) | !is.finite(weight), 1, weight))
}
