# Tier 1 (offline): Stage-2 helpers (no BGLR) and parent selection.

test_that(".genotype_means removes the environment main effect", {
  base <- c(g1 = 3, g2 = 1, g3 = 2)
  b <- tibble::tibble(
    trait = "yield",
    studyDbId = rep(c("S1", "S2"), each = 3),
    genotype  = rep(names(base), 2),
    BLUE      = c(base, base + 50),     # S2 shifted by a +50 environment effect
    weight    = 1)
  gm <- .genotype_means(b)
  expect_setequal(gm$genotype, names(base))
  expect_equal(gm$w[gm$genotype == "g1"], 2)            # summed weights
  expect_equal(gm$genotype[which.max(gm$y)], "g1")      # ranking follows base effect
  expect_equal(mean(gm$y), 0, tolerance = 1e-8)         # env-centered
})

test_that(".geno_kernel prefers G, falls back to markers, else errors", {
  G <- diag(3); dimnames(G) <- list(c("a","b","c"), c("a","b","c"))
  expect_identical(.geno_kernel(list(G = G)), G)
  M <- make_dosage(4, 20)
  expect_equal(dim(.geno_kernel(list(G = NULL, markers = M))), c(4, 4))
  expect_error(.geno_kernel(list(G = NULL, markers = NULL)), "neither")
})

test_that(".resolve_model downgrades marker-effect models without a marker matrix", {
  expect_equal(.resolve_model("BRR", list(markers = NULL)), "RKHS")
  expect_equal(.resolve_model("BRR", list(markers = make_dosage(2, 2))), "BRR")
  expect_equal(.resolve_model("RKHS", list(markers = NULL)), "RKHS")
})

test_that("select_parents ranks, flags redundancy, and writes a CSV", {
  gebv <- tibble::tibble(
    trait = "yield",
    genotype = c("A","B","C","D"),
    GEBV = c(3, 2.9, 1, 0),               # A and B both high
    phenotyped = c(TRUE, TRUE, FALSE, FALSE))
  G <- diag(4); dimnames(G) <- list(c("A","B","C","D"), c("A","B","C","D"))
  G["A","B"] <- G["B","A"] <- 0.9         # B closely related to higher-ranked A

  out <- file.path(tempdir(), "sel.csv")
  res <- select_parents(gebv, geno = list(G = G), n_select = 2,
                        relatedness_threshold = 0.5, outfile = out)
  expect_equal(res$genotype[1], "A")                    # ranked by index
  expect_true(all(res$rank == seq_len(nrow(res))))
  expect_equal(sum(res$selected), 2)
  expect_equal(res$redundant_with[res$genotype == "B"], "A")
  expect_true(is.na(res$redundant_with[res$genotype == "A"]))
  expect_true(file.exists(out))
})
