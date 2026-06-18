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

  gebv <- stage2_gblup(blues, geno, nIter = 800, burnIn = 150, refresh = TRUE)
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
  cv <- cv_accuracy(blues, geno, "yield", k = 3, nIter = 600, burnIn = 100)
  expect_equal(nrow(cv), 1)
  expect_true(is.finite(cv$accuracy) && abs(cv$accuracy) <= 1)
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

test_that("find_and_get_genotypes EM-combines two protocols into one GRM", {
  acc <- tibble::tibble(germplasmDbId = as.character(1:10),
                        germplasmName = paste0("N", 1:10))
  geno <- find_and_get_genotypes(make_geno_conn(), acc,
                                 protocol_id = NULL, pedigree_dir = NULL,
                                 refresh = TRUE)
  expect_setequal(geno$protocol_ids, c("501", "502"))
  expect_null(geno$markers)                              # multi-protocol -> no raw markers
  expect_true(isSymmetric(round(geno$G, 8)))
  # combined set = our accessions present + bridges B1,B2 = N1..N10 + B1,B2
  expect_setequal(rownames(geno$G), c(paste0("N", 1:10), "B1", "B2"))
  expect_true(all(is.finite(geno$G)))
})

test_that("find_and_get_genotypes single protocol keeps a raw marker matrix", {
  acc <- tibble::tibble(germplasmDbId = as.character(1:10),
                        germplasmName = paste0("N", 1:10))
  geno <- find_and_get_genotypes(make_geno_conn(), acc,
                                 protocol_id = "501", pedigree_dir = NULL,
                                 refresh = TRUE)
  expect_equal(geno$protocol_ids, "501")
  expect_false(is.null(geno$markers))                    # single protocol -> markers kept
  expect_true(isSymmetric(round(geno$G, 8)))
})
