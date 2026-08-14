# Step 6 (Stage 2 of the two-stage analysis): kernel GBLUP genomic prediction.
#
# Inputs are the Stage-1 BLUEs (one value per genotype x environment, with a
# weight = 1/SE^2) and the genotyping object from step 4 (a combined genomic
# relationship matrix `G`, or a marker matrix from which a kernel is built). We:
#   1. remove environment main effects (environment as a fixed effect) by
#      weighted-centering BLUEs within environment;
#   2. collapse to one weighted-mean value per genotype;
#   3. fit a kernel GBLUP over ALL genotyped accessions, leaving unphenotyped
#      candidates as NA so they receive predicted GEBVs.
#
# Prediction is ALWAYS kernel-based on the relationship matrix (no marker-effect
# Bayesian-alphabet models), so the single-protocol and EM-combined cases are
# handled identically. MIXED_MODEL_ENGINE picks the solver: "bglr" (Bayesian kernel
# GBLUP, fit as BRR ridge on the relationship eigen-factor -- NOT model="RKHS",
# which mishandles weights) or "sommer" (REML GBLUP).

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

# Fit BGLR RKHS (kernel GBLUP) for a single trait and return GEBVs for every
# candidate (rownames of the kernel relmat). Always kernel-based -- we deliberately
# do NOT offer marker-effect models (BRR/BayesB on a raw marker matrix), so the
# single-protocol case is treated identically to the EM-combined case.
.predict_one <- function(b, relmat, nIter, burnIn, seed) {
  cand <- rownames(relmat)
  gm   <- .genotype_means(filter(b, genotype %in% cand))

  y <- rep(NA_real_, length(cand)); names(y) <- cand
  w <- rep(1,        length(cand)); names(w) <- cand
  y[gm$genotype] <- gm$y
  w[gm$genotype] <- gm$w

  # Kernel GBLUP via ridge regression (BRR) on the eigen-factor L (relmat = L L'),
  # NOT BGLR's model = "RKHS". With weights, RKHS gives WRONG breeding values:
  # its eigendecomposition shortcut assumes iid residuals, so heteroscedastic
  # weights make its back-transformed u inflate beyond the data. BRR on L fits the
  # identical model (u = L b, so u ~ N(0, sigma2*relmat)) with the standard weighted
  # machinery and reproduces the analytic BLUP exactly (see code/bglr_rkhs_vs_brr.R).
  e <- eigen(relmat, symmetric = TRUE)
  L <- e$vectors %*% diag(sqrt(pmax(e$values, 0)))
  # BGLR sets Var(e_i) proportional to 1/weights_i^2 -- its `weights` are inverse
  # SDs, NOT inverse variances. Our w = sum(1/SE^2) is a precision, so pass sqrt(w).
  set.seed(seed)
  fm  <- BGLR::BGLR(y = y, weights = sqrt(w),
                    ETA = list(G = list(X = L, model = "BRR")),
                    nIter = nIter, burnIn = burnIn, verbose = FALSE,
                    saveAt = file.path(tempdir(), "bglr_GBLUP_"))
  gv  <- as.vector(L %*% fm$ETA$G$b)
  tibble(genotype = cand, GEBV = gv, phenotyped = cand %in% gm$genotype)
}

# Align a sommer BLUP object (vector / 1-col matrix / list-of-these) to `cand`,
# stripping the "genotype" variable-name prefix sommer prepends to level names.
.sommer_blups <- function(obj, cand) {
  if (is.list(obj) && !is.data.frame(obj)) obj <- obj[[1]]   # first random term
  if (is.matrix(obj)) { nm <- rownames(obj); v <- obj[, 1] } else { nm <- names(obj); v <- as.numeric(obj) }
  if (is.null(nm)) stop("sommer BLUPs have no level names to align to candidates.")
  nm <- sub("^genotype", "", nm)
  setNames(v, nm)[cand]
}

# Alternative engine: GBLUP via sommer (REML) instead of BGLR (Bayesian MCMC), as
# a cross-check on the GEBVs. Same contract as .predict_one: returns
# tibble(genotype, GEBV, phenotyped) over rownames(relmat). Always a relationship-kernel
# GBLUP (no marker-effect variant). Uses sommer's current API (mmes +
# vsm(ism(), Gu=relmat); sommer >= 4.3, tested on 4.4.5); BLUPs pulled via randef().
.predict_one_sommer <- function(b, relmat, seed) {
  if (!requireNamespace("sommer", quietly = TRUE)) {
    stop("engine = \"sommer\" needs the 'sommer' package: install.packages(\"sommer\").")
  }
  cand <- rownames(relmat)
  gm   <- .genotype_means(filter(b, genotype %in% cand))

  # Fit on the PHENOTYPED genotypes only, but keep ALL candidates as factor levels
  # and the full relmat as Gu, so the relationship kernel predicts the unphenotyped
  # (held-out, in CV) candidates. We must NOT feed sommer NA-y rows for them: mmes
  # excludes NA-response rows internally, which desyncs the full-length W from the
  # reduced data and aborts with a matrix-dimension mismatch. Supplying only the
  # phenotyped rows yields the kernel-only BLUP for the rest (no median-imputation
  # bias), matching the BGLR path which predicts NA rows from the kernel directly.
  dat <- data.frame(genotype = factor(gm$genotype, levels = cand),
                    y = gm$y, w = gm$w)

  set.seed(seed)
  args <- list(fixed  = y ~ 1,
               random = ~ sommer::vsm(sommer::ism(genotype), Gu = relmat),
               rcov   = ~ units, data = dat, verbose = FALSE)
  # Weight by the Stage-1 precisions to match BGLR (Var(e_i) = sigma2_e / w_i).
  # mmes forms R = sigma2_e * W^{-1}, so W = diag(w) gives exactly that. Use the
  # nrow form of diag() so a single phenotyped genotype yields a 1x1 W, not an
  # identity matrix sized by the scalar weight.
  args$W <- diag(dat$w, nrow = nrow(dat))
  if ("dateWarning" %in% names(formals(sommer::mmes))) args$dateWarning <- FALSE
  fm <- do.call(sommer::mmes, args)
  gv <- .sommer_blups(sommer::randef(fm), cand)
  tibble(genotype = cand, GEBV = unname(gv), phenotyped = cand %in% gm$genotype)
}

# Predict one trait with the chosen engine: "bglr" (Bayesian RKHS) or "sommer"
# (REML GBLUP). Both are kernel GBLUP on `relmat`; same returned shape.
.predict_engine <- function(b, relmat, engine, nIter, burnIn, seed) {
  if (engine == "sommer") .predict_one_sommer(b, relmat, seed)
  else                    .predict_one(b, relmat, nIter, burnIn, seed)
}

#' Stage-2 genomic prediction across one or more traits.
#'
#' @param blues  Stage-1 output (trait, studyDbId, genotype, BLUE, weight).
#' @param geno   step-4 list with $G (combined GRM) and/or $markers (used only to
#'   build a kernel when $G is absent; prediction is always kernel GBLUP).
#' @param traits which traits to predict (default: all present in `blues`).
#' @param engine "bglr" (Bayesian kernel GBLUP, BRR on the relationship
#'   eigen-factor) or "sommer" (REML GBLUP); defaults to `MIXED_MODEL_ENGINE`.
#' @return tibble: trait, genotype, GEBV, phenotyped. Cached to data/gebv.rds
#'   (engine "bglr") or data/gebv_sommer.rds (engine "sommer"); the cache is reused
#'   only when it still covers the requested traits over the current kernel's
#'   candidate set, else it is regenerated regardless of `refresh`.
stage2_gblup <- function(blues, geno, traits = NULL, engine = MIXED_MODEL_ENGINE,
                         nIter = BGLR_NITER, burnIn = BGLR_BURNIN,
                         seed = SEED, refresh = FALSE) {
  engine <- match.arg(engine, c("bglr", "sommer"))
  cache  <- cache_path(if (engine == "sommer") "gebv_sommer.rds" else "gebv.rds")
  relmat <- .geno_kernel(geno)
  if (is.null(traits)) traits <- unique(blues$trait)

  # NB this step deliberately does NOT use the request-keyed cache (code/cache.R) that
  # steps 2/3/4 share. Its rule is superset-tolerant on purpose -- a cache covering all
  # nine traits still serves a request for one of them -- and strict key equality would
  # throw that away and refit. Validating against the CONTENT is also stronger here:
  # it compares the cached genotypes to the kernel actually in hand, so a rebuilt G
  # invalidates the GEBVs even if every parameter looks identical.
  #
  # The cache exists to skip recomputation, NOT to serve stale predictions: reuse it
  # only when it still matches the current request -- same candidate set as the
  # relationship kernel AND covering every requested trait. A resized/regenerated G
  # or a newly requested trait regenerates even with refresh = FALSE, so switching
  # what we predict (e.g. bglr -> sommer over a rebuilt G, or adding traits) yields
  # fresh GEBVs without forcing the expensive upstream downloads to refresh too.
  # (Each engine has its own cache file, so the engine itself is part of the key.)
  if (!refresh && file.exists(cache)) {
    cached <- read_rds(cache)
    if (setequal(unique(cached$genotype), rownames(relmat)) &&
        all(traits %in% cached$trait)) {
      note_cache(cache)
      return(filter(cached, trait %in% traits))
    }
  }

  # One full model fit (a 12k-iteration MCMC under "bglr") per trait: bar over the
  # traits, plus each trait's own elapsed time as it lands, since a single fit is
  # otherwise an opaque multi-minute wait.
  # Note the kernel size once rather than per fit: every fit eigendecomposes the
  # full n x n kernel (06:.predict_one), which is most of the per-trait wait on a
  # large G, and saying so once explains it without a line per model.
  say("Stage 2: genomic prediction (", engine, ") for ", length(traits),
      " trait(s); each fit eigendecomposes the ", nrow(relmat), " x ",
      ncol(relmat), " kernel ...")
  pb <- pb_start(length(traits), paste0("Stage 2: traits (", engine, ")"))
  on.exit(pb_done(pb), add = TRUE)
  out <- map_dfr(traits, function(tr) {
    t0  <- Sys.time()
    res <- .predict_engine(filter(blues, trait == tr), relmat, engine,
                           nIter, burnIn, seed)
    pb_tick(pb)
    note(tr, ": ", .fmt_dur(as.numeric(difftime(Sys.time(), t0, units = "secs"))))
    mutate(res, trait = tr, .before = 1)
  })

  write_rds(out, cache)
  out
}

#' k-fold cross-validation accuracy (correlation of predicted vs observed
#' genotype means) for one trait -- a quick sanity check for the report.
cv_accuracy <- function(blues, geno, trait, k = 5, engine = MIXED_MODEL_ENGINE,
                        nIter = 6000, burnIn = 1000, seed = SEED) {
  engine <- match.arg(engine, c("bglr", "sommer"))
  relmat <- .geno_kernel(geno)
  b     <- filter(blues, trait == !!trait, genotype %in% rownames(relmat))
  gm    <- .genotype_means(b)
  # Need >= 2 genotyped means for a meaningful CV correlation; bail early (no
  # model fit) for traits with too little overlap with the kernel.
  if (nrow(gm) < 2) return(tibble(trait = trait, k = k, n = 0L, accuracy = NA_real_))
  set.seed(seed)
  folds <- sample(rep_len(1:k, nrow(gm)))
  # One full model fit per fold. The label carries the trait because the callers
  # loop over traits outside this function, so successive bars must be tellable apart.
  say("Cross-validating ", trait, ": ", k, " folds x one full ", engine, " fit ...")
  pb <- pb_start(k, paste0("CV folds: ", trait))
  on.exit(pb_done(pb), add = TRUE)
  preds <- map_dfr(1:k, function(f) {
    pb_tick(pb)
    b_tr <- filter(b, genotype %in% gm$genotype[folds != f])
    p <- .predict_engine(b_tr, relmat, engine, nIter, burnIn, seed)
    filter(p, genotype %in% gm$genotype[folds == f]) |>
      select(genotype, GEBV)
  })
  obs <- gm |> select(genotype, y)
  ev  <- inner_join(preds, obs, by = "genotype") |>
    filter(!is.na(GEBV), !is.na(y))
  # cor() errors on "no complete element pairs" and is undefined for <2 points;
  # a degenerate trait (e.g. too few genotypes for k folds) shouldn't abort the
  # whole report, so report NA accuracy for it instead.
  acc <- if (nrow(ev) >= 2) stats::cor(ev$GEBV, ev$y) else NA_real_
  tibble(trait = trait, k = k, n = nrow(ev), accuracy = acc)
}
