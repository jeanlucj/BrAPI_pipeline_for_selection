# Tier 2 (offline): step 4's intermediate caches -- per-protocol partial GRMs, the
# pedigree partials, and the EM-combined G before subsetting.
#
# These exist because rebuilding G re-parses gigabytes of VCFs, re-runs the
# per-marker imputation and re-derives every partial -- minutes apiece -- even when
# the only thing that changed was an EM weight. So the tests care about two things:
# that the expensive work is genuinely skipped, and that skipping it changes nothing
# about the answer.
#
# NB "delete the VCFs and see if it still works" is NOT a valid check: the flow needs
# the VCF headers for bridge detection regardless, so the cache does not make the
# files disposable. What it must avoid is re-PARSING them, so that is what is counted.

pc_conn <- function() {
  p1 <- make_dosage(9, 40, seed = 11, prefix = "")
  rownames(p1) <- c(paste0("N", 1:7), "B1", "B2")
  p2 <- make_dosage(8, 40, seed = 12, prefix = "")
  rownames(p2) <- c(paste0("N", 5:10), "B1", "B2")
  dosages <- list(`101` = p1, `102` = p2)

  cov_content <- function(type) {
    n <- sprintf
    list(results = list(
      counts = setNames(list(list("501", "502"),
                             setNames(list(7L, 6L), c("501", "502")), 10L),
                        c(n("ranked_genotyping_%ss", type),
                          n("accessions_by_genotyping_%s", type),
                          "accessions_total")),
      lookups = setNames(list(setNames(list("ProtoA", "ProtoB"), c("501", "502"))),
                         n("genotyping_%ss", type))))
  }
  list(
    filter_geno_protocols = function(accessions) list(content = cov_content("protocol")),
    filter_geno_projects  = function(accessions) list(content = cov_content("project")),
    vcf_archived_list = function(genotyping_protocol_id) {
      proj <- if (genotyping_protocol_id == 501) 101 else 102
      data.frame(protocol_id = genotyping_protocol_id, protocol_name = "P",
                 project_id = proj, project_name = paste0("proj", proj),
                 file_name = paste0("f", proj, ".vcf"), stringsAsFactors = FALSE)
    },
    vcf_archived = function(output, genotyping_project_id, file_name) {
      write_test_vcf(output, dosages[[as.character(genotyping_project_id)]])
    })
}

pc_acc <- function(n = 10) tibble::tibble(germplasmDbId = as.character(seq_len(n)),
                                          germplasmName = paste0("N", seq_len(n)))

# Count VCF->dosage parses: the expensive work the partial cache exists to skip.
count_parses <- function(env = parent.frame()) {
  hits <- new.env(parent = emptyenv()); hits$n <- 0L
  orig <- .vcf_to_dosage
  .vcf_to_dosage <<- function(...) { hits$n <- hits$n + 1L; orig(...) }
  withr::defer(.vcf_to_dosage <<- orig, envir = env)
  function() hits$n
}

# Bust only the top-level genotypes.rds, so step 4 re-derives G but may reuse
# everything beneath it.
bust_top <- function() {
  unlink(c(cache_path("genotypes.rds"), cache_key_file(cache_path("genotypes.rds"))))
}

test_that("re-deriving G reuses the partials: no VCF is parsed again, and G is identical", {
  parses <- count_parses()
  g1 <- find_and_get_genotypes(pc_conn(), pc_acc(), protocol_id = NULL,
                               pedigree_dir = NULL, refresh = TRUE)
  after_first <- parses()
  expect_gt(after_first, 0)                       # it really did parse the first time

  bust_top()
  g2 <- find_and_get_genotypes(pc_conn(), pc_acc(), protocol_id = NULL,
                               pedigree_dir = NULL)
  expect_equal(parses(), after_first)             # nothing re-parsed
  expect_equal(g2$G, g1$G)                        # identical, not merely similar
  expect_setequal(g2$protocol_ids, g1$protocol_ids)
})

test_that("changing only the EM weights reuses the partials and redoes the combine", {
  find_and_get_genotypes(pc_conn(), pc_acc(), protocol_id = NULL,
                         pedigree_dir = NULL, refresh = TRUE)
  parts_before <- file.mtime(list.files(cache_path("partials"), full.names = TRUE))
  expect_gt(length(parts_before), 0)

  parses <- count_parses()
  bust_top()
  g <- find_and_get_genotypes(pc_conn(), pc_acc(), protocol_id = NULL,
                              pedigree_dir = NULL, grm_df_mean = 120)
  expect_equal(parses(), 0L)                      # partials untouched ...
  expect_equal(file.mtime(list.files(cache_path("partials"), full.names = TRUE)),
               parts_before)
  expect_true(isSymmetric(round(g$G, 8)))         # ... but the combine still ran
})

test_that("the combined-G cache is reused when nothing upstream changed", {
  find_and_get_genotypes(pc_conn(), pc_acc(), protocol_id = NULL,
                         pedigree_dir = NULL, refresh = TRUE)
  combined_mtime <- file.mtime(cache_path("combined_G.rds"))
  bust_top()
  find_and_get_genotypes(pc_conn(), pc_acc(), protocol_id = NULL, pedigree_dir = NULL)
  expect_equal(file.mtime(cache_path("combined_G.rds")), combined_mtime)
})

test_that("refresh = TRUE rebuilds every layer", {
  find_and_get_genotypes(pc_conn(), pc_acc(), protocol_id = NULL,
                         pedigree_dir = NULL, refresh = TRUE)
  parses <- count_parses()
  find_and_get_genotypes(pc_conn(), pc_acc(), protocol_id = NULL,
                         pedigree_dir = NULL, refresh = TRUE)
  expect_gt(parses(), 0L)
})

test_that("the protocol key ignores accessions the protocol does not carry", {
  # .protocol_grm() builds its panel from `intersect(keep_samples, all_samples)`, so
  # requesting an accession this protocol never genotyped must not invalidate it --
  # otherwise every protocol rebuilds whenever any accession is added anywhere.
  k1 <- .protocol_key("501", character(0), c("N1", "N2"), c("N1", "N2", "N3"),
                      character(0))
  k2 <- .protocol_key("501", character(0), c("N1", "N2", "GHOST"),
                      c("N1", "N2", "N3"), character(0))
  expect_length(.key_diff(k1, k2), 0)

  # ... but an accession it DOES carry changes the panel, so it must invalidate.
  k3 <- .protocol_key("501", character(0), c("N1", "N2", "N3"),
                      c("N1", "N2", "N3"), character(0))
  expect_gt(length(.key_diff(k1, k3)), 0)
})

test_that("the file stamp ignores mtime but notices a changed size", {
  # mtime changes for reasons unrelated to content (a copy, an rsync, a restore),
  # and rebuilding every partial for that would cost hours.
  f <- tempfile(fileext = ".vcf")
  writeLines("x", f); s1 <- .file_stamp(f)
  Sys.setFileTime(f, Sys.time() + 3600)
  expect_equal(.file_stamp(f), s1)
  writeLines(c("x", "y"), f)
  expect_false(identical(.file_stamp(f), s1))
})

test_that("markers still come back for the single-protocol case from cache", {
  acc <- tibble::tibble(germplasmDbId = as.character(1:7),
                        germplasmName = paste0("N", 1:7))
  g1 <- find_and_get_genotypes(pc_conn(), acc, protocol_id = "501",
                               pedigree_dir = NULL, refresh = TRUE)
  expect_false(is.null(g1$markers))

  parses <- count_parses()
  bust_top()
  g2 <- find_and_get_genotypes(pc_conn(), acc, protocol_id = "501",
                               pedigree_dir = NULL)
  expect_equal(parses(), 0L)
  expect_equal(g2$markers, g1$markers)            # dosage survived in the partial cache
})
