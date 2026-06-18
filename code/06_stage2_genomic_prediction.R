# Step 6 (Stage 2 of the two-stage analysis): genomic prediction with BGLR.
#
# Inputs are the Stage-1 BLUEs (one value per genotype x environment, with a
# weight = 1/SE^2) and the genotyping object from step 4 (a combined genomic
# relationship matrix `G`, and optionally a raw marker matrix). We:
#   1. remove environment main effects (environment as a fixed effect) by
#      weighted-centering BLUEs within environment;
#   2. collapse to one weighted-mean value per genotype;
#   3. fit a genomic-prediction model with BGLR over ALL genotyped accessions,
#      leaving unphenotyped candidates as NA so they receive predicted GEBVs.
#
# Default model is RKHS on the combined relationship matrix. Marker-effect models
# (BRR / BayesB) need a single raw marker matrix, which only exists for the
# single-protocol case; otherwise the pipeline falls back to RKHS on G.

library(tidyverse)
here::i_am("code/06_stage2_genomic_prediction.R")
source(here::here("code", "config.R"))
source(here::here("code", "grm_utils.R"))   # .Gmatrix(), std_grm()

# Relationship kernel for prediction: the combined GRM from step 4 when present,
# else built from the single-protocol marker matrix.
.geno_kernel <- function(geno) {
  if (!is.null(geno$G)) return(geno$G)
  if (!is.null(geno$markers)) return(.Gmatrix(geno$markers))
  stop("genotyping object has neither a combined GRM (G) nor a marker matrix.")
}

# Collapse Stage-1 BLUEs (one trait) to a weighted genotype mean after removing
# environment main effects. Returns tibble(genotype, y, w).
.genotype_means <- function(b) {
  b |>
    group_by(studyDbId) |>
    mutate(env_mean = stats::weighted.mean(BLUE, weight)) |>
    ungroup() |>
    mutate(y_adj = BLUE - env_mean) |>
    group_by(genotype) |>
    summarise(y = stats::weighted.mean(y_adj, weight),
              w = sum(weight), .groups = "drop")
}

# Fit BGLR for a single trait and return GEBVs for every candidate (rownames of
# the kernel K). For marker-effect models a raw marker matrix is required.
.predict_one <- function(b, K, markers, model, nIter, burnIn, seed) {
  cand <- rownames(K)
  gm   <- .genotype_means(filter(b, genotype %in% cand))

  y <- rep(NA_real_, length(cand)); names(y) <- cand
  w <- rep(1,         length(cand)); names(w) <- cand
  y[gm$genotype] <- gm$y
  w[gm$genotype] <- gm$w

  set.seed(seed)
  saveAt <- file.path(tempdir(), paste0("bglr_", model, "_"))

  if (model == "RKHS") {
    ETA <- list(G = list(K = K[cand, cand], model = "RKHS"))
    fm  <- BGLR::BGLR(y = y, ETA = ETA, weights = w,
                      nIter = nIter, burnIn = burnIn, verbose = FALSE, saveAt = saveAt)
    gv  <- as.vector(fm$ETA$G$u)
  } else {
    Xc  <- scale(markers[cand, , drop = FALSE], center = TRUE, scale = FALSE)
    ETA <- list(M = list(X = Xc, model = model))   # "BRR", "BayesB", "BL", ...
    fm  <- BGLR::BGLR(y = y, ETA = ETA, weights = w,
                      nIter = nIter, burnIn = burnIn, verbose = FALSE, saveAt = saveAt)
    gv  <- as.vector(Xc %*% fm$ETA$M$b)
  }

  tibble(genotype = cand, GEBV = gv, phenotyped = cand %in% gm$genotype)
}

# Resolve the model, downgrading marker-effect models to RKHS when no single raw
# marker matrix is available (e.g. an EM-combined multi-platform GRM).
.resolve_model <- function(model, geno) {
  if (model != "RKHS" && is.null(geno$markers)) {
    message("No single marker matrix (combined GRM in use); using RKHS instead of ", model, ".")
    return("RKHS")
  }
  model
}

#' Stage-2 genomic prediction across one or more traits.
#'
#' @param blues  Stage-1 output (trait, studyDbId, genotype, BLUE, weight).
#' @param geno   step-4 list with $G (combined GRM) and/or $markers.
#' @param traits which traits to predict (default: all present in `blues`).
#' @return tibble: trait, genotype, GEBV, phenotyped. Cached to data/gebv.rds.
stage2_gblup <- function(blues, geno,
                         traits = NULL, model = BGLR_MODEL,
                         nIter = BGLR_NITER, burnIn = BGLR_BURNIN,
                         seed = SEED, refresh = FALSE) {
  cache <- cache_path("gebv.rds")
  if (!refresh && file.exists(cache)) return(read_rds(cache))

  K     <- .geno_kernel(geno)
  model <- .resolve_model(model, geno)
  if (is.null(traits)) traits <- unique(blues$trait)

  out <- map_dfr(traits, function(tr) {
    b <- filter(blues, trait == tr)
    .predict_one(b, K, geno$markers, model, nIter, burnIn, seed) |>
      mutate(trait = tr, .before = 1)
  })

  write_rds(out, cache)
  out
}

#' k-fold cross-validation accuracy (correlation of predicted vs observed
#' genotype means) for one trait -- a quick sanity check for the report.
cv_accuracy <- function(blues, geno, trait, k = 5, model = BGLR_MODEL,
                        nIter = 6000, burnIn = 1000, seed = SEED) {
  K     <- .geno_kernel(geno)
  model <- .resolve_model(model, geno)
  b     <- filter(blues, trait == !!trait, genotype %in% rownames(K))
  gm    <- .genotype_means(b)
  set.seed(seed)
  folds <- sample(rep_len(1:k, nrow(gm)))
  preds <- map_dfr(1:k, function(f) {
    b_tr <- filter(b, genotype %in% gm$genotype[folds != f])
    p <- .predict_one(b_tr, K, geno$markers, model, nIter, burnIn, seed)
    filter(p, genotype %in% gm$genotype[folds == f]) |>
      select(genotype, GEBV)
  })
  obs <- gm |> select(genotype, y)
  ev  <- inner_join(preds, obs, by = "genotype")
  tibble(trait = trait, k = k, n = nrow(ev),
         accuracy = stats::cor(ev$GEBV, ev$y, use = "complete.obs"))
}
