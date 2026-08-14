# Is BGLR's weighted RKHS itself the problem, or is the BGLR-vs-sommer gap just
# Bayesian-posterior-mean vs REML-plug-in?
#
# Fit the SAME variance model two ways under weights = sqrt(w):
#   - RKHS on G            (eigendecomposition shortcut; assumes iid residuals)
#   - BRR ridge on L, G=LL' (u = L b, Var(u)=sigma2_b*G; standard weighted machinery)
# Both are u ~ N(0, sigma2*G). Compare each to the analytic plug-in BLUP for its own
# variance components and to sommer (exact, ground truth). A valid shrinkage BLUP
# never exceeds the data spread; we flag ranges that do.

suppressMessages(library(tidyverse))
here::i_am("code/bglr_rkhs_vs_brr.R")
source(here::here("code", "config.R")); source(here::here("code", "grm_utils.R"))

ec <- read_rds(cache_path("engine_compare.rds"))
y <- ec$y; w <- ec$w; gebv_s <- ec$gebv_s              # sommer = ground truth (weighted)
G <- read_rds(cache_path("genotypes.rds"))$G; G <- G[names(y), names(y)]
n <- length(y)
cat(sprintf("data y range [%.2f, %.2f] (spread %.2f)\n", min(y), max(y), diff(range(y))))

# eigen factor: G = L L'
e <- eigen(G, symmetric = TRUE); L <- e$vectors %*% diag(sqrt(pmax(e$values, 0)))

# analytic plug-in GBLUP: u = s2g G (s2g G + s2e diag(1/w))^-1 (y - mu_GLS)
hand <- function(s2g, s2e) {
  Vi  <- solve(s2g * G + s2e * diag(1 / w)); one <- rep(1, n)
  mu  <- as.numeric((one %*% Vi %*% y) / (one %*% Vi %*% one))
  as.vector(s2g * G %*% Vi %*% (y - mu))
}
rep_ <- function(a, b) sprintf("cor=%.3f slope=%.2f", cor(a, b), unname(coef(lm(b ~ a))[2]))
rng  <- function(u) sprintf("[%.2f,%.2f] sd=%.2f%s", min(u), max(u), sd(u),
                            if (diff(range(u)) > diff(range(y))) "  *EXCEEDS DATA*" else "")

set.seed(SEED)
fm_rkhs <- BGLR::BGLR(y = y, weights = sqrt(w), ETA = list(G = list(K = G, model = "RKHS")),
                      nIter = BGLR_NITER, burnIn = BGLR_BURNIN, verbose = FALSE,
                      saveAt = file.path(tempdir(), "rk_"))
u_rkhs <- as.vector(fm_rkhs$ETA$G$u); s2g_rk <- fm_rkhs$ETA$G$varU; s2e_rk <- fm_rkhs$varE

set.seed(SEED)
fm_brr <- BGLR::BGLR(y = y, weights = sqrt(w), ETA = list(M = list(X = L, model = "BRR")),
                     nIter = BGLR_NITER, burnIn = BGLR_BURNIN, verbose = FALSE,
                     saveAt = file.path(tempdir(), "br_"))
u_brr <- as.vector(L %*% fm_brr$ETA$M$b); s2g_br <- fm_brr$ETA$M$varB; s2e_br <- fm_brr$varE

cat("\n=== weighted (sqrt(w)) RKHS vs BRR (same u ~ N(0, s2*G) model) ===\n")
cat(sprintf("RKHS : varU=%.3f varE=%.3f lambda=%.1f | u %s\n", s2g_rk, s2e_rk, s2e_rk/s2g_rk, rng(u_rkhs)))
cat(sprintf("BRR  : varB=%.3f varE=%.3f lambda=%.1f | u %s\n", s2g_br, s2e_br, s2e_br/s2g_br, rng(u_brr)))
cat(sprintf("sommer (ground truth)                  | u %s\n", rng(gebv_s)))
cat("\n  RKHS u vs analytic(RKHS varcomps):  ", rep_(hand(s2g_rk, s2e_rk), u_rkhs), "\n")
cat("  BRR  u vs analytic(BRR  varcomps):  ", rep_(hand(s2g_br, s2e_br), u_brr), "\n")
cat("  RKHS u vs sommer:                   ", rep_(gebv_s, u_rkhs), "\n")
cat("  BRR  u vs sommer:                   ", rep_(gebv_s, u_brr), "\n")
