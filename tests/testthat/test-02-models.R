# Tier 2 (offline, heavier): real BGLR genomic prediction and a fully mocked
# multi-protocol EM-combine integration. No network. A few seconds to run.

test_that("stage2_gblup (real BGLR/RKHS) predicts all candidates incl. unphenotyped", {
  M <- make_dosage(30, 200, seed = 3)
  G <- .Gmatrix(M)
  geno <- list(G = G, markers = NULL)

  set.seed(3)
  u <- as.numeric(scale(M %*% rnorm(ncol(M))))          # polygenic genetic values
  names(u) <- rownames(M)
  train <- rownames(M)[1:20]                             # 10 left unphenotyped
  blues <- tibble::tibble(trait = "yield", studyDbId = "S1", genotype = train,
                          BLUE = u[train] + rnorm(20, 0, 0.3), weight = 1)

  gebv <- stage2_gblup(blues, geno, engine = "bglr", nIter = 800, burnIn = 150,
                       refresh = TRUE)
  expect_equal(nrow(gebv), 30)                           # every candidate scored
  expect_equal(sum(gebv$phenotyped), 20)
  expect_true(all(is.finite(gebv$GEBV)))
  expect_gt(cor(gebv$GEBV, u[gebv$genotype]), 0.3)       # tracks true values
})

test_that("cv_accuracy returns a finite correlation", {
  M <- make_dosage(25, 150, seed = 5)
  geno <- list(G = .Gmatrix(M), markers = NULL)
  set.seed(5)
  u <- as.numeric(scale(M %*% rnorm(ncol(M)))); names(u) <- rownames(M)
  blues <- tibble::tibble(trait = "yield", studyDbId = "S1",
                          genotype = rownames(M), BLUE = u + rnorm(25, 0, 0.3), weight = 1)
  cv <- cv_accuracy(blues, geno, "yield", k = 3, engine = "bglr",
                    nIter = 600, burnIn = 100)
  expect_equal(nrow(cv), 1)
  expect_true(is.finite(cv$accuracy) && abs(cv$accuracy) <= 1)
})

test_that("cv_accuracy returns NA (not an error) for a degenerate trait", {
  M <- make_dosage(25, 150, seed = 5)
  geno <- list(G = .Gmatrix(M), markers = NULL)
  # A trait whose genotypes don't overlap the kernel -> no complete pairs.
  blues <- tibble::tibble(trait = "yield", studyDbId = "S1",
                          genotype = paste0("ghost", 1:25), BLUE = rnorm(25), weight = 1)
  cv <- expect_no_error(cv_accuracy(blues, geno, "yield", k = 3))
  expect_equal(nrow(cv), 1)
  expect_true(is.na(cv$accuracy))
})

# ---- mocked multi-protocol find_and_get_genotypes (EM combine path) ----------

make_geno_conn <- function() {
  # Two protocols sharing samples N5..N7 + bridges B1,B2 (genotyped in both).
  p1 <- make_dosage(9, 40, seed = 11, prefix = "")
  rownames(p1) <- c(paste0("N", 1:7), "B1", "B2")
  p2 <- make_dosage(8, 40, seed = 12, prefix = "")
  rownames(p2) <- c(paste0("N", 5:10), "B1", "B2")
  dosages <- list(`101` = p1, `102` = p2)

  cov_content <- function(type) {
    n <- sprintf
    list(results = list(
      counts = setNames(list(list("501","502"),
                             setNames(list(7L, 6L), c("501","502")), 10L),
                        c(n("ranked_genotyping_%ss", type),
                          n("accessions_by_genotyping_%s", type),
                          "accessions_total")),
      lookups = setNames(list(setNames(list("ProtoA","ProtoB"), c("501","502"))),
                         n("genotyping_%ss", type))))
  }
  list(
    filter_geno_protocols = function(accessions) list(content = cov_content("protocol")),
    filter_geno_projects  = function(accessions) list(content = cov_content("project")),
    vcf_archived_list = function(genotyping_protocol_id) {
      proj <- if (genotyping_protocol_id == 501) 101 else 102
      data.frame(protocol_id = genotyping_protocol_id,
                 protocol_name = "P", project_id = proj,
                 project_name = paste0("proj", proj),
                 file_name = paste0("f", proj, ".vcf"),
                 stringsAsFactors = FALSE)
    },
    vcf_archived = function(output, genotyping_project_id, file_name) {
      write_test_vcf(output, dosages[[as.character(genotyping_project_id)]])
      invisible(NULL)
    })
}

test_that("find_and_get_genotypes EM-combines protocols then subsets to targets", {
  acc <- tibble::tibble(germplasmDbId = as.character(1:10),
                        germplasmName = paste0("N", 1:10))
  geno <- find_and_get_genotypes(make_geno_conn(), acc,
                                 protocol_id = NULL, pedigree_dir = NULL,
                                 refresh = TRUE)
  expect_setequal(geno$protocol_ids, c("501", "502"))
  expect_null(geno$markers)                              # multi-protocol -> no raw markers
  expect_true(isSymmetric(round(geno$G, 8)))
  # G is subset to the training targets; bridges B1/B2 informed the combine but
  # are NOT prediction targets, so they are dropped.
  expect_setequal(rownames(geno$G), paste0("N", 1:10))
  expect_false(any(c("B1", "B2") %in% rownames(geno$G)))
  expect_true(all(is.finite(geno$G)))
})

test_that("find_and_get_genotypes accepts explicit EM degrees-of-freedom args", {
  acc <- tibble::tibble(germplasmDbId = as.character(1:10),
                        germplasmName = paste0("N", 1:10))
  geno <- find_and_get_genotypes(make_geno_conn(), acc, protocol_id = NULL,
                                 pedigree_dir = NULL, grm_df_mean = 200,
                                 grm_df_stdev = 15, pedigree_df = 100,
                                 refresh = TRUE)
  expect_true(isSymmetric(round(geno$G, 8)))
  expect_setequal(rownames(geno$G), paste0("N", 1:10))
})

test_that("find_and_get_genotypes predicts a TEST accession bridged via pedigree", {
  # "PB" is genotyped in one protocol and sits in a pedigree group with N1/N2.
  # Passed as a TEST accession (by name), it should be resolved, combined, and
  # retained in the final (subset) G. Distinct protocol/project ids avoid colliding
  # with make_geno_conn()'s on-disk VCF cache.
  p1 <- make_dosage(10, 40, seed = 11)
  rownames(p1) <- c(paste0("N", 1:7), "B1", "B2", "PB")
  p2 <- make_dosage(8, 40, seed = 12)
  rownames(p2) <- c(paste0("N", 5:10), "B1", "B2")
  dosages <- list(`301` = p1, `302` = p2)

  cov_content <- function(type) {
    n <- sprintf
    list(results = list(
      counts = setNames(list(list("601", "602"),
                             setNames(list(7L, 6L), c("601", "602")), 10L),
                        c(n("ranked_genotyping_%ss", type),
                          n("accessions_by_genotyping_%s", type), "accessions_total")),
      lookups = setNames(list(setNames(list("ProtoA", "ProtoB"), c("601", "602"))),
                         n("genotyping_%ss", type))))
  }
  conn <- list(
    filter_geno_protocols = function(accessions) list(content = cov_content("protocol")),
    filter_geno_projects  = function(accessions) list(content = cov_content("project")),
    vcf_archived_list = function(genotyping_protocol_id) {
      proj <- if (genotyping_protocol_id == 601) 301 else 302
      data.frame(protocol_id = genotyping_protocol_id, protocol_name = "P",
                 project_id = proj, project_name = paste0("proj", proj),
                 file_name = paste0("f", proj, ".vcf"), stringsAsFactors = FALSE)
    },
    vcf_archived = function(output, genotyping_project_id, file_name) {
      write_test_vcf(output, dosages[[as.character(genotyping_project_id)]]); invisible(NULL)
    })

  parent <- tempfile(); dir.create(parent)
  dir <- file.path(parent, "T3_Oat"); dir.create(dir)
  readr::write_csv(tibble::tibble(
    germplasmDbId_i = c(201, 1, 2, 201, 201, 1),
    germplasmDbId_j = c(201, 1, 2, 1,   2,   2),
    relationship    = c(1,   1, 1, 0.5, 0.5, 0.5)),
    file.path(dir, "T3_Oat_group2.csv"))
  saveRDS(list(group2 = c("201", "1", "2")),
          file.path(dir, "T3_Oat_pedigree_groups.rds"))
  saveRDS(list(list(germplasmDbId = 201, germplasmName = "PB"),
               list(germplasmDbId = 1,   germplasmName = "N1"),
               list(germplasmDbId = 2,   germplasmName = "N2")),
          file.path(parent, "germplasm_cache_T3_Oat.rds"))

  acc <- tibble::tibble(germplasmDbId = as.character(1:10),
                        germplasmName = paste0("N", 1:10))
  geno <- find_and_get_genotypes(conn, acc, test_names = "PB",
                                 protocol_id = NULL, pedigree_dir = dir,
                                 refresh = TRUE)
  expect_true("PB" %in% rownames(geno$G))      # test accession, retained via combine
  expect_false(any(c("B1", "B2") %in% rownames(geno$G)))  # bridges still dropped
  expect_true(all(is.finite(geno$G)))
})

test_that("find_and_get_genotypes single protocol keeps a raw marker matrix", {
  acc <- tibble::tibble(germplasmDbId = as.character(1:7),
                        germplasmName = paste0("N", 1:7))   # all in protocol 501
  geno <- find_and_get_genotypes(make_geno_conn(), acc,
                                 protocol_id = "501", pedigree_dir = NULL,
                                 refresh = TRUE)
  expect_equal(geno$protocol_ids, "501")
  expect_false(is.null(geno$markers))                    # single protocol, no injection
  expect_setequal(rownames(geno$markers), rownames(geno$G))
  expect_true(isSymmetric(round(geno$G, 8)))
})

test_that("find_and_get_genotypes injects ungenotyped training, drops ungenotyped test", {
  train <- tibble::tibble(germplasmDbId = as.character(c(1:10, 99)),
                          germplasmName = c(paste0("N", 1:10), "GHOST"))  # GHOST: no geno
  test  <- tibble::tibble(germplasmDbId = "77", germplasmName = "ABSENT") # not genotyped
  geno <- find_and_get_genotypes(make_geno_conn(), train, test_accessions = test,
                                 protocol_id = NULL, pedigree_dir = NULL,
                                 refresh = TRUE)
  expect_true("GHOST" %in% rownames(geno$G))             # training -> force-injected
  expect_false("ABSENT" %in% rownames(geno$G))           # test, ungenotyped -> dropped
  expect_false(any(c("B1", "B2") %in% rownames(geno$G))) # bridges dropped
  # injected row: diagonal = mean diagonal of the genotyped block, off-diagonals 0
  others <- setdiff(rownames(geno$G), "GHOST")
  expect_equal(geno$G["GHOST", "GHOST"], mean(diag(geno$G[others, others])))
  expect_true(all(geno$G["GHOST", others] == 0))
  expect_null(geno$markers)                              # injection -> markers dropped
})
