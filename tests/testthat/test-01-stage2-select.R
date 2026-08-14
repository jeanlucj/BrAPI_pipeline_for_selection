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

test_that("select_parents writes the two-block breeders_output.csv (un-standardized)", {
  tn <- c("yield", "height")                            # full "trait names"
  gebv <- tidyr::expand_grid(genotype = c("A", "B", "C"), trait = tn) |>
    dplyr::mutate(GEBV = c(10, 1,  6, 2,  2, 3),         # A: 10/1, B: 6/2, C: 2/3
                  phenotyped = genotype %in% c("A", "B"))
  weights <- setNames(c(2, -1), tn)                     # index = 2*yield - 1*height
  shorts  <- setNames(c("Yld", "Ht"), tn)
  cv <- tibble::tibble(trait = tn, accuracy = c(0.5, 0.3))

  out <- file.path(tempdir(), "breeders.csv")
  res <- select_parents(gebv, cv = cv, trait_weights = weights,
                        trait_names = tn, trait_short = shorts, outfile = out)

  # Index is the RAW weighted sum (no standardization); sorted descending.
  expect_equal(res$parents$accession, c("A", "B", "C"))     # 19 > 10 > 1
  expect_equal(res$parents$Index, c(2*10 - 1, 2*6 - 2, 2*2 - 3))
  expect_equal(res$parents$In_Training, c(1L, 1L, 0L))       # C not phenotyped
  expect_setequal(names(res$parents), c("accession", "In_Training", "Yld", "Ht", "Index"))

  # n_select truncates to the top n by index.
  expect_equal(nrow(select_parents(gebv, cv = cv, trait_weights = weights,
                                   trait_names = tn, trait_short = shorts,
                                   n_select = 2, outfile = tempfile())$parents), 2)

  # CSV shape: header + spacer column + block 2.
  raw <- readLines(out)
  expect_equal(raw[1], "accession,In_Training,Yld,Ht,Index,,Trait,Weight,CV_accuracy")
  expect_match(raw[2], ",,Yld,2,0.5$")                       # blank spacer, block-2 row 1
  expect_true(file.exists(out))
})

test_that("select_parents with NULL weights drops the Index and Weight columns", {
  tn <- c("yield", "height")
  gebv <- tidyr::expand_grid(genotype = c("A", "B", "C"), trait = tn) |>
    dplyr::mutate(GEBV = c(10, 1,  6, 2,  2, 3),
                  phenotyped = genotype %in% c("A", "B"))
  shorts <- setNames(c("Yld", "Ht"), tn)
  cv <- tibble::tibble(trait = tn, accuracy = c(0.5, 0.3))

  out <- file.path(tempdir(), "breeders_noweights.csv")
  res <- select_parents(gebv, cv = cv, trait_weights = NULL,
                        trait_names = tn, trait_short = shorts, outfile = out)

  # No selection index: Index column gone, rows ordered by accession, GEBVs present.
  expect_setequal(names(res$parents), c("accession", "In_Training", "Yld", "Ht"))
  expect_equal(res$parents$accession, c("A", "B", "C"))      # sorted by accession
  expect_false("Weight" %in% names(res$traits))              # block 2 drops Weight
  expect_equal(readLines(out)[1],
               "accession,In_Training,Yld,Ht,,Trait,CV_accuracy")
})
