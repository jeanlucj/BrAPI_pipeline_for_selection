# Tier 1 (offline): step-4 helpers exercised with synthetic VCFs / matrices.

test_that(".vcf_to_dosage round-trips genotypes (incl. NA) and unique-ifies IDs", {
  D <- make_dosage(6, 8); D[1, 1] <- NA               # one missing call
  vcf <- write_test_vcf(tempfile(fileext = ".vcf"), D)
  M <- .vcf_to_dosage(vcf, target_density = NULL)
  expect_equal(nrow(M), 6); expect_equal(ncol(M), 8)
  expect_false(anyDuplicated(colnames(M)) > 0)        # IDs made unique
  expect_equal(unname(M[rownames(D), ][, order(colnames(M))][, 1]),
               unname(D[, 1]))                          # values preserved (col 1)
  expect_true(is.na(M["L1", ][1]) || any(is.na(M)))    # NA carried through
})

test_that(".vcf_to_dosage can restrict to keep_samples", {
  D <- make_dosage(6, 5)
  vcf <- write_test_vcf(tempfile(fileext = ".vcf"), D)
  M <- .vcf_to_dosage(vcf, keep_samples = c("L2", "L4"), target_density = NULL)
  expect_setequal(rownames(M), c("L2", "L4"))
})

test_that("marker thinning counts and thins to ~target density", {
  D <- make_dosage(4, 100)
  vcf <- write_test_vcf(tempfile(fileext = ".vcf"), D)
  expect_equal(.count_markers(vcf), 100)
  th <- .thin_to_target(vcf, target = 25)              # floor(100/25)=4 -> ~25 kept
  expect_false(is.null(th$tmp))
  expect_lte(.count_markers(th$tmp), 30)
  expect_gte(.count_markers(th$tmp), 20)
  expect_true(is.null(.thin_to_target(vcf, target = 100)$tmp))  # <2x target: no thin
})

test_that(".read_vcf_samples reads the #CHROM header", {
  D <- make_dosage(3, 2)
  vcf <- write_test_vcf(tempfile(fileext = ".vcf"), D)
  expect_equal(.read_vcf_samples(vcf), c("L1", "L2", "L3"))
})

test_that(".merge_dosage unions markers and de-duplicates accessions", {
  a <- make_dosage(3, 3, seed = 1, prefix = "X"); colnames(a) <- c("m1","m2","m3")
  b <- make_dosage(3, 3, seed = 2, prefix = "Y"); colnames(b) <- c("m3","m4","m5")
  M <- .merge_dosage(list(a, b))
  expect_setequal(colnames(M), c("m1","m2","m3","m4","m5"))
  expect_equal(nrow(M), 6)
  expect_true(is.na(M["X1", "m4"]))                    # X rows lack m4/m5
})

test_that(".qc_markers drops minority-annotation markers before accessions", {
  # 10 accessions on markers m1..m5, 4 on disjoint m6..m10 -> every accession ~50% NA
  M <- matrix(NA_real_, 14, 10,
              dimnames = list(paste0("a", 1:14), paste0("m", 1:10)))
  set.seed(7)
  M[1:10, 1:5]  <- sample(0:2, 50, replace = TRUE)
  M[11:14, 6:10] <- sample(0:2, 20, replace = TRUE)
  Q <- .qc_markers(M, max_missing = 0.5, min_maf = 0.0)
  expect_true(all(colnames(Q) %in% paste0("m", 1:5)))  # majority markers kept
  expect_false(any(paste0("a", 11:14) %in% rownames(Q)))# minority accessions dropped
  expect_false(anyNA(Q))                                # imputed
})

test_that(".read_pedigree_group reconstructs a dense symmetric matrix", {
  csv <- tempfile(fileext = ".csv")
  readr::write_csv(tibble::tibble(
    germplasmDbId_i = c(101, 102, 101),
    germplasmDbId_j = c(101, 102, 102),
    relationship    = c(1,   1,   0.5)), csv)
  A <- .read_pedigree_group(csv)
  expect_equal(dim(A), c(2, 2))
  expect_true(isSymmetric(unname(A)))
  expect_equal(A["101", "102"], 0.5)
  expect_equal(diag(A), c(`101` = 1, `102` = 1))
})

test_that(".pedigree_partials only returns overlapping groups, standardized", {
  dir <- file.path(tempdir(), "ped"); dir.create(dir, showWarnings = FALSE)
  readr::write_csv(tibble::tibble(
    germplasmDbId_i = c(101, 102), germplasmDbId_j = c(101, 102),
    relationship = c(1, 1)),
    file.path(dir, "T3_Oat_group2.csv"))
  hit  <- .pedigree_partials(dir, setNames(c("LineA","LineB"), c("101","102")))
  miss <- .pedigree_partials(dir, setNames("Z", "999"))
  expect_length(hit, 1)
  expect_equal(mean(diag(hit[[1]])), 1)
  expect_setequal(rownames(hit[[1]]), c("LineA", "LineB"))
  expect_length(miss, 0)
  expect_length(.pedigree_partials(NULL, setNames("Z","1")), 0)  # disabled
})

test_that(".coverage_table parses BrAPI content; .safe_coverage swallows errors", {
  content <- list(results = list(
    counts = list(ranked_genotyping_protocols = list("8", "61"),
                  accessions_by_genotyping_protocol = list(`8` = 30, `61` = 25),
                  accessions_total = 40),
    lookups = list(genotyping_protocols = list(`8` = "GBS POGI", `61` = "Diversity"))))
  ct <- .coverage_table(content, "protocol")
  expect_equal(nrow(ct), 2)
  expect_equal(ct$dbId, c("8", "61"))
  expect_equal(ct$n_covered, c(30L, 25L))
  expect_equal(ct$n_total, c(40L, 40L))
  expect_equal(nrow(.safe_coverage(function(...) stop("boom"), "x", "protocol")), 0)
})
