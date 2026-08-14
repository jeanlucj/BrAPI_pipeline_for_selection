# Step 5 (Stage 1 of the two-stage analysis): per-trial adjusted means (BLUEs).
#
# For each study x trait we fit a mixed model with genotype as a FIXED effect
# (so we get BLUEs) and the field-design factors (rep, block) as RANDOM effects.
# This adjusts genotype means for spatial/design structure, including augmented
# designs where unreplicated entries are adjusted via replicated check blocks.
# Each BLUE is carried to Stage 2 with a weight = 1 / SE^2.
#
# Three guards keep that weight well-behaved: (1) a study x trait with a constant
# response is dropped (no estimable genotypic signal -- it would fit perfectly and
# get an infinite weight); (2) rep and block are not both entered when they are the
# same partition (aliased variance components); (3) each SE is floored at
# SE_FLOOR_FRAC of the response SD so a near-perfect fit cannot blow up 1/SE^2.

library(tidyverse)
here::i_am("code/05_stage1_blues.R")
source(here::here("code", "config.R"))

# Fit one study x trait and return a tibble: genotype, BLUE, SE.
.fit_one <- function(d) {
  d <- d |>
    mutate(genotype = factor(germplasmName),
           rep = factor(rep), block = factor(block))

  # Build the random-effects part from whatever design replication exists. Drop a
  # redundant design factor: when rep and block induce the same partition (or one is
  # nested-identical in the other -- e.g. a trial whose rep and block labels coincide),
  # entering both as crossed random intercepts aliases their variance components, which
  # lme4 flags as "nearly unidentifiable". They are redundant exactly when one factor
  # determines the other (the observed rep x block combinations number no more than the
  # levels of one factor); keep only the finer factor in that case.
  nr <- nlevels(droplevels(d$rep)); nb <- nlevels(droplevels(d$block))
  use_rep <- nr > 1; use_block <- nb > 1
  if (use_rep && use_block) {
    npair <- dplyr::n_distinct(paste(d$rep, d$block))
    if (npair == nr || npair == nb) {                 # one factor determines the other
      if (nr >= nb) use_block <- FALSE else use_rep <- FALSE
    }
  }
  rand <- c(if (use_rep) "(1 | rep)", if (use_block) "(1 | block)")

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

  # Floor each SE at SE_FLOOR_FRAC of the response SD so a (near-)perfect fit cannot
  # drive 1/SE^2 toward Inf and let one trial dominate Stage 2. Well-estimated SEs sit
  # far above this floor and are untouched; NA SEs are left for the median fallback in
  # stage1_blues(). (A fully constant trial is dropped before reaching here, so the
  # SD -- hence the floor -- is positive.)
  se_floor <- SE_FLOOR_FRAC * stats::sd(d$value)
  est |>
    mutate(genotype = str_remove(term, "^genotype"),
           SE = pmax(SE, se_floor)) |>
    select(genotype, BLUE, SE)
}

#' Stage-1 BLUEs for every study x trait in a long phenotype table.
#'
#' @param pheno long table with columns studyDbId, germplasmName, rep, block,
#'   trait, value (the `pheno` element from get_phenotypes()).
#' @return tibble: trait, studyDbId, genotype, BLUE, SE, weight (= 1/SE^2).
#'   Studies/traits that cannot be fit, or whose response is constant (no estimable
#'   genotypic signal), are skipped with a message. SE is floored at SE_FLOOR_FRAC
#'   of the response SD before weighting.
stage1_blues <- function(pheno) {
  # Typed empty result for skipped study x trait fits, so the BLUE/SE schema survives
  # row-binding even when EVERY group is skipped (else the weight step below can't
  # find the SE column).
  empty <- tibble(genotype = character(), BLUE = double(), SE = double())
  dat <- filter(pheno, !is.na(value), !is.na(germplasmName))
  # One model fit per trait x study cell; the total is knowable up front, so the
  # bar shows how far through the grid we are.
  n_fits <- nrow(distinct(dat, trait, studyDbId))
  say("Stage 1: fitting ", n_fits, " trait x study BLUE model(s) ...")
  pb <- pb_start(n_fits, "Stage 1: trait x study fits")
  on.exit(pb_done(pb), add = TRUE)
  dat |>
    group_by(trait, studyDbId) |>
    group_modify(function(d, key) {
      pb_tick(pb)
      if (n_distinct(d$germplasmName) < 2) {
        message(sprintf("Skipping %s / study %s: <2 genotypes.", key$trait, key$studyDbId))
        return(empty)
      }
      if (n_distinct(d$value) < 2) {
        message(sprintf("Skipping %s / study %s: response is constant (no estimable genotypic signal).",
                        key$trait, key$studyDbId))
        return(empty)
      }
      tryCatch(.fit_one(d), error = function(e) {
        message(sprintf("Skipping %s / study %s: %s", key$trait, key$studyDbId,
                        conditionMessage(e)))
        empty
      })
    }) |>
    ungroup() |>
    # weight by inverse squared SE; default weight 1 where SE is unavailable
    mutate(weight = if_else(!is.na(SE) & SE > 0, 1 / SE^2, NA_real_)) |>
    mutate(weight = if_else(is.na(weight), stats::median(weight, na.rm = TRUE), weight)) |>
    mutate(weight = if_else(is.na(weight) | !is.finite(weight), 1, weight))
}
