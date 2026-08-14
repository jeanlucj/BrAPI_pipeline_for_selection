## Possible issue: `model = "RKHS"` with observation `weights` returns inflated random effects

Dear BGLR Gustavo and Paulino,
[I had Claude write this up.  Forgive awkwardness.]
First, thank you for BGLR — we use it heavily and appreciate the work. While doing
two-stage genomic prediction (stage‑1 BLUEs fed into a stage‑2 GBLUP, with
observation weights proportional to the stage‑1 precisions), we ran into behavior
that we believe may be a bug, and we wanted to share a minimal reproducible example
and our analysis in case it is useful. We may well be misunderstanding something —
if so, we would be grateful to be corrected.

**Summary.** When `model = "RKHS"` is combined with observation `weights`, the
returned random effects (`$ETA[[k]]$u`) appear to be systematically **inflated** —
in our example ~6× — and exceed the true signal, which a shrinkage predictor should
not do. The variance components are estimated correctly, and the *unweighted* RKHS
fit is correct. Fitting the **mathematically identical model** as a ridge regression
(`model = "BRR"`) on the eigen‑factor `L` (where `K = L Lᵀ`), with the same weights,
gives the correct, weight‑consistent estimates. This makes us think the weights are
not fully propagated into the RKHS random‑effect computation, even though they are
used correctly for variance‑component estimation.

Environment: R 4.5.2, BGLR 1.1.4.

## Minimal reproducible example

```r
library(BGLR)

set.seed(1)
n <- 200

## A positive-definite relationship matrix (mean diagonal 1) and its eigen-factor
## L, so that K = L %*% t(L).
X   <- matrix(rnorm(n * 400), n)
K   <- tcrossprod(scale(X)) / ncol(X);  K <- K / mean(diag(K))
EVD <- eigen(K, symmetric = TRUE)
L   <- EVD$vectors %*% diag(sqrt(pmax(EVD$values, 0)))

## Simulate g ~ N(0, K) with HETEROSCEDASTIC residuals. The BGLR help states the
## residual variance of observation i is proportional to 1/weights_i^2, so we
## generate residuals to match that convention.
sg2 <- 1; se2 <- 1
g   <- as.vector(L %*% rnorm(n, sd = sqrt(sg2)))
w   <- runif(n, 1, 10)                       # weights, a 10x spread
y   <- 2 + g + rnorm(n, sd = sqrt(se2) / w)  # Var(e_i) = se2 / w_i^2

ni <- 12000; bi <- 2000

## (1) RKHS on K, weighted
f_rkhs <- BGLR(y, ETA = list(list(K = K, model = "RKHS")),
               weights = w, nIter = ni, burnIn = bi, verbose = FALSE)
u_rkhs <- as.vector(f_rkhs$ETA[[1]]$u)

## (2) The SAME model as a ridge regression on L: u = L b, so u ~ N(0, var * K).
##     Identical prior on the random effect, but the standard (non-eigendecomposition)
##     machinery. Same weights.
f_brr <- BGLR(y, ETA = list(list(X = L, model = "BRR")),
              weights = w, nIter = ni, burnIn = bi, verbose = FALSE)
u_brr <- as.vector(L %*% f_brr$ETA[[1]]$b)

## (3) UNWEIGHTED RKHS, as a control.
f_unwt <- BGLR(y, ETA = list(list(K = K, model = "RKHS")),
               nIter = ni, burnIn = bi, verbose = FALSE)
u_unwt <- as.vector(f_unwt$ETA[[1]]$u)

## Analytic GBLUP (BLUP at fixed variance components) for reference.
gblup <- function(s2g, s2e, R) {
  Vi <- solve(s2g * K + R); mu <- sum(Vi %*% y) / sum(Vi)
  as.vector(s2g * K %*% Vi %*% (y - mu))
}
a_rkhs <- gblup(f_rkhs$ETA[[1]]$varU, f_rkhs$varE, f_rkhs$varE * diag(1 / w^2))
a_brr  <- gblup(f_brr$ETA[[1]]$varB,  f_brr$varE,  f_brr$varE  * diag(1 / w^2))
a_unwt <- gblup(f_unwt$ETA[[1]]$varU, f_unwt$varE, f_unwt$varE * diag(n))

## Variance components are identical for RKHS and BRR:
c(varU_RKHS = f_rkhs$ETA[[1]]$varU, varE_RKHS = f_rkhs$varE,
  varB_BRR  = f_brr$ETA[[1]]$varB,  varE_BRR  = f_brr$varE)

## A BLUP is shrunken and should not exceed the true signal -- but weighted RKHS does:
rbind(true_g    = range(g),
      RKHS_wt   = range(u_rkhs),
      BRR_wt    = range(u_brr),
      RKHS_unwt = range(u_unwt))

## Agreement with the analytic BLUP at each fit's OWN variance components:
c(cor = cor(u_rkhs, a_rkhs), slope = coef(lm(u_rkhs ~ a_rkhs))[2])  # RKHS  weighted
c(cor = cor(u_brr,  a_brr),  slope = coef(lm(u_brr  ~ a_brr))[2])   # BRR   weighted
c(cor = cor(u_unwt, a_unwt), slope = coef(lm(u_unwt ~ a_unwt))[2])  # RKHS  unweighted
```

## Observed output (seed 1; MCMC numbers vary slightly, the pattern does not)

```
var comps  RKHS: varU=0.67 varE=5.16 | BRR: varB=0.67 varE=5.21

range of estimates           low      high
  true_g                   -2.89      2.28
  RKHS  (weighted)        -14.95     13.81     <- ~6x the true signal
  BRR   (weighted)         -1.78      2.09
  RKHS  (unweighted)       -2.13      1.94

agreement with analytic BLUP (cor, slope):
  RKHS  weighted     cor = 0.948   slope = 6.04     <- inflated, ~ lambda = varE/varU
  BRR   weighted     cor = 1.000   slope = 1.00     <- exact
  RKHS  unweighted   cor = 1.000   slope = 1.00     <- exact
```

## Why this looks like a problem

- **RKHS and BRR fit the identical model** `u ~ N(0, σ²·K)` (since `K = L Lᵀ`, a
  ridge on `L` puts exactly this prior on `u = L b`). They estimate the **same
  variance components**. They should therefore return the same random effects.
- **They do, when unweighted, and they do not, when weighted.** The weighted RKHS
  estimates are ~6× larger than both the true breeding values and the BRR estimates,
  and they exceed the range of the true signal — which a shrinkage estimator cannot
  legitimately do.
- The weighted RKHS estimates also fail to match the **analytic BLUP computed at the
  fit's own posterior‑mean variance components** (slope ≈ 6), whereas weighted BRR
  and unweighted RKHS match it exactly (slope = 1). The inflation factor tracks
  `λ = varE/varU`, which grows as the weights spread the residual variances.

Because the **variance‑component estimation is correct** and only the **random‑effect
output is wrong, and only under weights**, our (external, unverified) hypothesis is
that the weights are correctly used in the likelihood but are not fully accounted for
when the RKHS path forms/back‑transforms `u` from the eigenvector representation
(`K = U D Uᵀ`), which is exact for homoscedastic residuals.

## Workaround (for other users who hit this)

Fit the equivalent model as a ridge regression on an eigen‑factor of the kernel,
which handles weights correctly:

```r
EVD <- eigen(K, symmetric = TRUE)
L   <- EVD$vectors %*% diag(sqrt(pmax(EVD$values, 0)))   # K = L Lᵀ
fm  <- BGLR(y, ETA = list(list(X = L, model = "BRR")), weights = w, ...)
u   <- as.vector(L %*% fm$ETA[[1]]$b)                    # u ~ N(0, var · K)
```

## Questions

1. Is this expected behavior — i.e., are observation `weights` intended to be
   unsupported (or handled differently than we assume) for `model = "RKHS"`?
2. I assume this is unintended. Do you want me to open a GitHub issue?

Thank you very much for maintaining BGLR.
