# Experiments to explain the BGLR vs sommer GEBV difference (Freeze damage).
#   (1) Fit both engines (now BOTH weighted by 1/SE^2) and report variance
#       components + the implied genetic ratio sigma2_g/(sigma2_g+sigma2_e).
#   (3) Hand-compute the GBLUP BLUP with FIXED variance components and compare to
#       each engine's GEBVs: if each engine reproduces the hand solution for its
#       OWN (sigma2_g, sigma2_e), the engines are both correct and the difference
#       is purely variance estimation; if an engine does NOT reproduce its own
#       hand solution, that engine has a parameterization/scale problem.
# Offline (cached phenotypes/trials/genotypes); refits both models.

suppressMessages({library(tidyverse); library(httr)})
here::i_am("code/engine_compare.R")
source(here::here("code", "04_find_genotyping.R"))
source(here::here("code", "03_get_phenotypes.R"))
source(here::here("code", "05_stage1_blues.R"))
source(here::here("code", "06_stage2_genomic_prediction.R"))

trait <- "Freeze damage severity - 0-9 Rating|CO_350:0005001"
trials <- read_rds(cache_path("ny_trials.rds")); pheno <- read_rds(cache_path("phenotypes.rds"))
geno   <- read_rds(cache_path("genotypes.rds")); sets <- split_by_role(pheno, trials)
blues  <- stage1_blues(sets$train_pheno)

G    <- geno$G
cand <- rownames(G)
gm   <- .genotype_means(filter(blues, trait == !!trait, genotype %in% cand))
y <- setNames(rep(NA_real_, length(cand)), cand); y[gm$genotype] <- gm$y
w <- setNames(rep(1,        length(cand)), cand); w[gm$genotype] <- gm$w
pheno_idx <- !is.na(y)
cat(sprintf("Freeze: %d candidates, %d phenotyped. y sd=%.3f\n",
            length(cand), sum(pheno_idx), sd(y[pheno_idx])))

## ---- Fit BGLR RKHS (weighted) ----------------------------------------------
set.seed(SEED)
fm_b <- BGLR::BGLR(y = y, weights = w,
                   ETA = list(G = list(K = G, model = "RKHS")),
                   nIter = BGLR_NITER, burnIn = BGLR_BURNIN, verbose = FALSE,
                   saveAt = file.path(tempdir(), "cmp_bglr_"))
s2g_b <- fm_b$ETA$G$varU; s2e_b <- fm_b$varE
gebv_b <- setNames(as.vector(fm_b$ETA$G$u), cand)

## ---- Fit sommer GBLUP (weighted) -------------------------------------------
gv_s_tib <- .predict_one_sommer(filter(blues, trait == !!trait), G, SEED)
gebv_s <- setNames(gv_s_tib$GEBV, gv_s_tib$genotype)[cand]
# refit once more to grab varcomps (predict_one_sommer doesn't return the fit)
dat <- data.frame(genotype = factor(cand, levels = cand), y = y, w = w)
set.seed(SEED)
fm_s <- sommer::mmes(y ~ 1, random = ~ sommer::vsm(sommer::ism(genotype), Gu = G),
                     rcov = ~ units, data = dat, W = diag(dat$w),
                     verbose = FALSE, dateWarning = FALSE)
vc <- summary(fm_s)$varcomp
s2g_s <- vc[grep("genotype", rownames(vc))[1], "VarComp"]
s2e_s <- vc[grep("units",    rownames(vc))[1], "VarComp"]

cat("\n=== (1) Variance components (both weighted) ===\n")
ratio <- function(g, e) g / (g + e)
cat(sprintf("BGLR  : sigma2_g=%.3f  sigma2_e=%.3f  g/(g+e)=%.3f  lambda(e/g)=%.2f\n",
            s2g_b, s2e_b, ratio(s2g_b, s2e_b), s2e_b/s2g_b))
cat(sprintf("sommer: sigma2_g=%.3f  sigma2_e=%.3f  g/(g+e)=%.3f  lambda(e/g)=%.2f\n",
            s2g_s, s2e_s, ratio(s2g_s, s2e_s), s2e_s/s2g_s))

## ---- (3) Hand GBLUP with fixed variance components -------------------------
# All targets phenotyped here, so Z = I over `cand`; BLUP:
#   V = s2g*G + s2e*diag(1/w);  mu = (1' Vinv y)/(1' Vinv 1);  u = s2g*G*Vinv*(y-mu)
hand_blup <- function(s2g, s2e) {
  Vi  <- solve(s2g * G + s2e * diag(1 / w))
  one <- rep(1, length(cand))
  mu  <- as.numeric((one %*% Vi %*% y) / (one %*% Vi %*% one))
  as.vector(s2g * G %*% Vi %*% (y - mu))
}
u_hand_b <- hand_blup(s2g_b, s2e_b)
u_hand_s <- hand_blup(s2g_s, s2e_s)

cmp <- function(a, b) c(cor = cor(a, b), slope = unname(coef(lm(b ~ a))[2]),
                        max_abs_diff = max(abs(a - b)), sd_engine = sd(b), sd_hand = sd(a))
cat("\n=== (3) Engine GEBV vs hand-GBLUP with that engine's OWN variance components ===\n")
cat("BGLR   vs hand(BGLR theta):  ");  print(round(cmp(u_hand_b, gebv_b), 3))
cat("sommer vs hand(sommer theta):"); print(round(cmp(u_hand_s, gebv_s), 3))
cat("\nCross-check: hand(sommer theta) vs hand(BGLR theta) cor=",
    round(cor(u_hand_s, u_hand_b), 3),
    " | engine-vs-engine cor=", round(cor(gebv_b, gebv_s), 3), "\n")
cat(sprintf("\nRanges: BGLR [%.1f,%.1f]  sommer [%.2f,%.2f]  hand(BGLR) [%.1f,%.1f]  hand(sommer) [%.2f,%.2f]\n",
            min(gebv_b), max(gebv_b), min(gebv_s), max(gebv_s),
            min(u_hand_b), max(u_hand_b), min(u_hand_s), max(u_hand_s)))

write_rds(list(s2g_b=s2g_b, s2e_b=s2e_b, s2g_s=s2g_s, s2e_s=s2e_s,
               gebv_b=gebv_b, gebv_s=gebv_s, u_hand_b=u_hand_b, u_hand_s=u_hand_s,
               y=y, w=w), cache_path("engine_compare.rds"))
