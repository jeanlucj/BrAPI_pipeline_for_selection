# Shared genomic-relationship-matrix helpers, sourced by step 4 (build per-source
# GRMs to EM-combine) and step 6 (build a kernel from a marker matrix when no
# combined GRM is supplied).

# VanRaden genomic relationship matrix from a 0/1/2 dosage matrix.
.Gmatrix <- function(M) {
  p <- colMeans(M) / 2
  W <- sweep(M, 2, 2 * p, "-")
  G <- tcrossprod(W) / (2 * sum(p * (1 - p)))
  G + diag(1e-6, nrow(G))   # nudge to positive-definite
}

# Standardize a relationship matrix so its mean diagonal is 1 (the form the EM
# combiner expects for its partial covariances).
std_grm <- function(G) G / mean(diag(G))

# Effective number of independent samples behind a relationship matrix, via
# Galwey (2009)'s "effective number of independent tests": (sum sqrt(lambda))^2 /
# sum(lambda) over the positive eigenvalues. Scale-invariant (robust to std_grm
# scaling), bounded [1, rank]. Used to set each marker GRM's EM degrees of freedom
# (markers are in LD, so a raw count overstates the independent information).
.effective_n <- function(G) {
  lambda <- eigen(G, symmetric = TRUE, only.values = TRUE)$values
  lambda <- lambda[lambda > 1e-8]                    # drop ~zero / tiny-negative
  if (length(lambda) < 1) return(1)
  (sum(sqrt(lambda)))^2 / sum(lambda)
}
