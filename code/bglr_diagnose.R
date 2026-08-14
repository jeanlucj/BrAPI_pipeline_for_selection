# Why does BGLR RKHS return inflated GEBVs while sommer is exact? Experiments:
#   (B) kernel conditioning of relmat (the injected/disconnected rows make it
#       near-singular -- BGLR RKHS drops/!handles small eigenvalues).
#   (U) BGLR UNWEIGHTED: does dropping the 1/SE^2 weights make u sane / match its
#       own analytic BLUP? (weights are the prime suspect: sigma2_e was 10.5 > var(y)).
#   (C) convergence: BGLR weighted at two seeds + a 5x-longer chain -- is u stable?
#   (S) standardized kernel (mean diag = 1, the form BGLR RKHS expects): does u
#       then match the analytic BLUP?
# A correct GBLUP engine reproduces the analytic BLUP for its own variance
# components, so each fit is scored by cor/slope/range vs hand_blup(its theta).

suppressMessages({library(tidyverse)})
here::i_am("code/bglr_diagnose.R")
source(here::here("code", "config.R")); source(here::here("code", "grm_utils.R"))

ec <- read_rds(cache_path("engine_compare.rds"))
y <- ec$y; w <- ec$w
G <- read_rds(cache_path("genotypes.rds"))$G
G <- G[names(y), names(y)]
n <- length(y)

# Analytic BLUP for fixed variance components; R = sigma2_e * diag(1/wt) (weighted)
# or sigma2_e * I (unweighted, wt = 1).
hand_blup <- function(K, s2g, s2e, wt) {
  Vi  <- solve(s2g * K + s2e * diag(1 / wt))
  one <- rep(1, nrow(K))
  mu  <- as.numeric((one %*% Vi %*% y) / (one %*% Vi %*% one))
  as.vector(s2g * K %*% Vi %*% (y - mu))
}
score <- function(engine_u, hand_u, tag) {
  cat(sprintf("  %-22s cor=%.3f slope=%.2f range[%.2f,%.2f] (hand range[%.2f,%.2f])\n",
              tag, cor(engine_u, hand_u), unname(coef(lm(engine_u ~ hand_u))[2]),
              min(engine_u), max(engine_u), min(hand_u), max(hand_u)))
}
fit_bglr <- function(K, wt = NULL, nIter = BGLR_NITER, burnIn = BGLR_BURNIN, seed = SEED) {
  set.seed(seed)
  fm <- BGLR::BGLR(y = y, weights = wt,
                   ETA = list(G = list(K = K, model = "RKHS")),
                   nIter = nIter, burnIn = burnIn, verbose = FALSE,
                   saveAt = file.path(tempdir(), "diag_"))
  list(u = setNames(as.vector(fm$ETA$G$u), names(y)),
       s2g = fm$ETA$G$varU, s2e = fm$varE)
}

cat("=== (B) Kernel conditioning of relmat ===\n")
ev <- eigen(G, symmetric = TRUE, only.values = TRUE)$values
cat(sprintf("n=%d  mean diag=%.3f  diag range[%.3f,%.3f]\n", n, mean(diag(G)),
            min(diag(G)), max(diag(G))))
cat(sprintf("eigenvalues: max=%.3f min=%.4f  #(<1e-8)=%d  #(<0)=%d  rank~%d  cond=%.1e\n",
            max(ev), min(ev), sum(ev < 1e-8), sum(ev < 0), sum(ev > 1e-8),
            max(ev) / max(min(ev[ev > 1e-8]), 1e-12)))

cat("\n=== (U) BGLR UNWEIGHTED vs analytic BLUP (R = s2e*I) ===\n")
bu <- fit_bglr(G, wt = NULL)
cat(sprintf("  varU=%.3f varE=%.3f lambda=%.2f\n", bu$s2g, bu$s2e, bu$s2e / bu$s2g))
score(bu$u, hand_blup(G, bu$s2g, bu$s2e, rep(1, n)), "BGLR-unwt vs hand")
cat(sprintf("  cor(BGLR-unwt u, sommer GEBV)=%.3f\n", cor(bu$u, ec$gebv_s)))

cat("\n=== (C) Convergence: weighted, two seeds + long chain ===\n")
bw1 <- fit_bglr(G, wt = w, seed = 1234)
bw2 <- fit_bglr(G, wt = w, seed = 99)
bwl <- fit_bglr(G, wt = w, seed = 1234, nIter = 60000, burnIn = 15000)
cat(sprintf("  seed1234 vs seed99: cor(u)=%.3f  | sd1=%.2f sd2=%.2f\n",
            cor(bw1$u, bw2$u), sd(bw1$u), sd(bw2$u)))
cat(sprintf("  short vs long chain: cor(u)=%.3f | varE short=%.2f long=%.2f | u-range long[%.1f,%.1f]\n",
            cor(bw1$u, bwl$u), bw1$s2e, bwl$s2e, min(bwl$u), max(bwl$u)))
score(bwl$u, hand_blup(G, bwl$s2g, bwl$s2e, w), "BGLR-wt-long vs hand")

cat("\n=== (S) Standardized kernel (mean diag = 1), weighted ===\n")
Ks <- std_grm(G)
bs <- fit_bglr(Ks, wt = w)
cat(sprintf("  varU=%.3f varE=%.3f lambda=%.2f\n", bs$s2g, bs$s2e, bs$s2e / bs$s2g))
score(bs$u, hand_blup(Ks, bs$s2g, bs$s2e, w), "BGLR-std-K vs hand")
cat(sprintf("  cor(BGLR-stdK u, sommer GEBV)=%.3f\n", cor(bs$u, ec$gebv_s)))
