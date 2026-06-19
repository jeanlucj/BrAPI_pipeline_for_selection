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

test_that(".read_pedigree_group densifies only the kept-dbid intersection", {
  csv <- tempfile(fileext = ".csv")
  readr::write_csv(tibble::tibble(
    germplasmDbId_i = c(101, 102, 103, 101, 102),
    germplasmDbId_j = c(101, 102, 103, 102, 103),
    relationship    = c(1,   1,   1,   0.5, 0.25)), csv)
  A <- .read_pedigree_group(csv, keep_dbids = c("101", "102"))
  expect_setequal(rownames(A), c("101", "102"))   # 103 never materialized
  expect_equal(A["101", "102"], 0.5)
  expect_equal(diag(A), c(`101` = 1, `102` = 1))
})

test_that(".pedigree_partials standardizes groups and restricts to keep_names", {
  dir <- file.path(tempfile(), "T3_Oat"); dir.create(dir, recursive = TRUE)
  csv <- file.path(dir, "T3_Oat_group2.csv")
  readr::write_csv(tibble::tibble(
    germplasmDbId_i = c(101, 102, 103), germplasmDbId_j = c(101, 102, 103),
    relationship = c(1, 1, 1)), csv)
  members <- setNames(list(c("101", "102", "103")), csv)
  d2n <- setNames(c("LineA", "LineB", "Bridge"), c("101", "102", "103"))

  hit <- .pedigree_partials(members, d2n, keep_names = c("LineA", "LineB", "Bridge"))
  expect_length(hit, 1)
  expect_equal(mean(diag(hit[[1]])), 1)
  expect_setequal(rownames(hit[[1]]), c("LineA", "LineB", "Bridge"))

  # keep_names drops the non-kept accession (103/Bridge) before densifying.
  sub <- .pedigree_partials(members, d2n, keep_names = c("LineA", "LineB"))
  expect_setequal(rownames(sub[[1]]), c("LineA", "LineB"))

  expect_length(.pedigree_partials(members, d2n, keep_names = "Z"), 0)  # no overlap
  expect_length(.pedigree_partials(list(),  d2n, keep_names = "LineA"), 0)  # disabled
})

test_that(".germplasm_name_map loads dbid->name from the companion cache", {
  parent <- tempfile(); dir.create(parent)
  dir <- file.path(parent, "T3_Oat"); dir.create(dir)
  saveRDS(list(list(germplasmDbId = 7, germplasmName = "LineG"),
               list(germplasmDbId = 8, germplasmName = "LineH")),
          file.path(parent, "germplasm_cache_T3_Oat.rds"))
  m <- .germplasm_name_map(dir)
  expect_equal(m[c("7", "8")], c(`7` = "LineG", `8` = "LineH"))
  expect_length(.germplasm_name_map(file.path(parent, "Absent")), 0)  # no cache
  expect_length(.germplasm_name_map(NULL), 0)                          # disabled
})

test_that(".resolve_names_to_dbids maps names via the cache and warns on misses", {
  parent <- tempfile(); dir.create(parent)
  dir <- file.path(parent, "T3_Oat"); dir.create(dir)
  saveRDS(list(list(germplasmDbId = 7, germplasmName = "LineG"),
               list(germplasmDbId = 8, germplasmName = "LineH")),
          file.path(parent, "germplasm_cache_T3_Oat.rds"))
  expect_message(
    res <- .resolve_names_to_dbids(c("LineG", "Ghost"), dir),
    "could not be resolved")
  expect_equal(res$germplasmDbId[res$germplasmName == "LineG"], "7")
  expect_false("Ghost" %in% res$germplasmName)                 # unresolved dropped
})

test_that(".center_dfs re-centers on mean_df and caps the spread", {
  # observed SD (5) below the cap (15): spread preserved, mean re-centered.
  out <- .center_dfs(c(20, 25, 30), mean_df = 60, sd_df = 15)
  expect_equal(mean(out), 60)
  expect_equal(sd(out), sd(c(20, 25, 30)))
  expect_equal(order(out), order(c(20, 25, 30)))             # ordering preserved

  # observed SD large: spread capped at sd_df.
  capped <- .center_dfs(c(10, 50, 400), mean_df = 60, sd_df = 15)
  expect_equal(mean(capped), 60)
  expect_equal(sd(capped), 15)

  expect_equal(.center_dfs(42, 60, 15), 60)                  # single GRM -> mean_df
  expect_equal(.center_dfs(c(7, 7, 7), 60, 15), rep(60, 3))  # zero SD -> mean_df
  expect_true(all(.center_dfs(c(1, 2, 1000), 20, 50) > 0))   # floored positive
})

test_that(".pedigree_group_members reads the groups rds, falls back to the CSV", {
  parent <- tempfile(); dir.create(parent)
  dir <- file.path(parent, "T3_Oat"); dir.create(dir)
  csv <- file.path(dir, "T3_Oat_group1.csv")
  readr::write_csv(tibble::tibble(germplasmDbId_i = c(7, 8), germplasmDbId_j = c(7, 8),
                                  relationship = c(1, 1)), csv)
  rds <- file.path(dir, "T3_Oat_pedigree_groups.rds")
  saveRDS(list(group1 = c("7", "8", "9")), rds)        # 9 not in CSV -> proves source
  m <- .pedigree_group_members(dir)
  expect_named(m, csv)
  expect_setequal(m[[csv]], c("7", "8", "9"))           # from rds

  file.remove(rds)
  expect_setequal(.pedigree_group_members(dir)[[csv]], c("7", "8"))  # CSV fallback
  expect_length(.pedigree_group_members(NULL), 0)       # disabled
})

test_that("malformed pedigree companions fail loudly (not silently)", {
  # group CSV missing a required column
  bad_csv <- tempfile(fileext = ".csv")
  readr::write_csv(tibble::tibble(germplasmDbId_i = "1", relationship = 1), bad_csv)
  expect_error(.read_pedigree_group(bad_csv), "missing column")

  # germplasm cache that isn't a list of BrAPI records
  parent <- tempfile(); dir.create(parent)
  dir <- file.path(parent, "T3_Oat"); dir.create(dir)
  saveRDS(data.frame(germplasmDbId = 1, germplasmName = "X"),
          file.path(parent, "germplasm_cache_T3_Oat.rds"))
  expect_error(.germplasm_name_map(dir), "germplasm records")

  # groups rds that isn't a NAMED list
  dirp <- file.path(tempfile(), "T3_Oat"); dir.create(dirp, recursive = TRUE)
  saveRDS(list(c("1", "2"), c("3")),
          file.path(dirp, "T3_Oat_pedigree_groups.rds"))
  expect_error(.pedigree_group_members(dirp), "named list")
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
