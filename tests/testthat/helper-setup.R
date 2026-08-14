# Shared setup + synthetic-data builders for the test suite.
# testthat sources helper*.R into the test environment before running tests.

suppressMessages(library(tidyverse))

# Never narrate the test suite. SHOW_PROGRESS already defaults to interactive()
# (FALSE under Rscript), but the suite is also run from RStudio, so pin it.
options(brapi.progress = FALSE)

# Source the whole pipeline (offline: defines functions, does not connect).
for (f in c("config.R", "grm_utils.R", "em_covariance_combiner.R",
            "01_connect.R", "02_find_trials.R", "03_get_phenotypes.R",
            "04_find_genotyping.R", "05_stage1_blues.R",
            "06_stage2_genomic_prediction.R", "07_select.R")) {
  source(here::here("code", f))
}

# Redirect caches/outputs to a throwaway dir so tests never touch data/ or
# output/ in the repo. The pipeline functions are source()'d into the global
# environment (so is cache_path/output_path from config.R); we must override
# them THERE with `<<-`, not in this helper's environment, or the functions
# would keep resolving the real data/ paths.
.TEST_TMP <- file.path(tempdir(), "brapi_tests")
dir.create(.TEST_TMP, showWarnings = FALSE, recursive = TRUE)
cache_path  <<- function(...) { p <- file.path(.TEST_TMP, ...); dir.create(dirname(p), showWarnings = FALSE, recursive = TRUE); p }
output_path <<- function(...) { p <- file.path(.TEST_TMP, ...); dir.create(dirname(p), showWarnings = FALSE, recursive = TRUE); p }

# --- synthetic-data builders -------------------------------------------------

# Write a minimal vcfR-readable VCF from a dosage matrix (accessions x markers,
# values 0/1/2/NA). All IDs are "." so .vcf_to_dosage exercises its unique-ID fix.
write_test_vcf <- function(path, dosage) {
  gt_of <- c("0" = "0/0", "1" = "0/1", "2" = "1/1")
  hdr <- c("##fileformat=VCFv4.2",
           "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
           paste(c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT",
                   rownames(dosage)), collapse = "\t"))
  body <- vapply(seq_len(ncol(dosage)), function(j) {
    calls <- ifelse(is.na(dosage[, j]), "./.", gt_of[as.character(dosage[, j])])
    paste(c("1", j * 100L, ".", "A", "G", ".", "PASS", ".", "GT", calls), collapse = "\t")
  }, character(1))
  writeLines(c(hdr, body), path)
  path
}

# A reproducible accessions x markers dosage matrix.
make_dosage <- function(n_acc = 20, n_mar = 50, seed = 1, prefix = "L") {
  set.seed(seed)
  M <- matrix(sample(0:2, n_acc * n_mar, replace = TRUE), n_acc, n_mar)
  dimnames(M) <- list(paste0(prefix, seq_len(n_acc)),
                      paste0("m", seq_len(n_mar)))
  M
}

# A RELATED population of inbred lines: `n_fam` biparental crosses among `n_founders`
# founders, each cross giving doubled-haploid progeny (dosages 0/2, as oat lines are
# essentially homozygous). Lines from one cross share ~half their genome; lines from
# different crosses are related through shared founders.
#
# Why not make_dosage(): it draws every line's markers independently, so the lines are
# mutually UNRELATED, every GRM off-diagonal is ~0, and no line can be predicted from
# any other. Genomic prediction of an unphenotyped line is then impossible BY
# CONSTRUCTION -- held-out accuracy is ~0 no matter how correct the code is. Relatedness
# is the thing being exploited, so the fixture has to contain some.
make_related_dosage <- function(n_acc = 200, n_mrk = 2000, n_founders = 20,
                                n_fam = 20, seed = 1, prefix = "L") {
  set.seed(seed)
  p <- stats::runif(n_mrk, 0.1, 0.9)                       # founder allele frequencies
  H <- matrix(stats::rbinom(n_founders * n_mrk, 1, rep(p, each = n_founders)),
              n_founders, n_mrk)                            # founder haplotypes
  parents <- t(replicate(n_fam, sample.int(n_founders, 2)))
  fam_of  <- rep_len(seq_len(n_fam), n_acc)
  D <- vapply(seq_len(n_acc), function(i) {
    pa <- parents[fam_of[i], ]
    from_a <- stats::rbinom(n_mrk, 1, 0.5) == 1             # which parent each locus came from
    2L * ifelse(from_a, H[pa[1], ], H[pa[2], ])             # doubled haploid: 0 or 2
  }, numeric(n_mrk)) |> t()
  dimnames(D) <- list(paste0(prefix, seq_len(n_acc)), paste0("m", seq_len(n_mrk)))
  attr(D, "family") <- setNames(fam_of, rownames(D))
  D
}

# Multi-trial phenotypes with a REAL genetic basis, so that the marker-based GRM is
# the correct kernel for the trait and Stage 2 has actual signal to find. Without
# this, phenotypes are independent of the markers, REML drives the genetic variance
# component to the boundary, and sommer returns BLUPs of exactly 0 while BGLR (whose
# prior keeps sigma2_g > 0) returns something non-zero -- which looks like an engine
# bug and is not one.
#
# Generative model (every marker is causal; effects additive):
#   g_i    = sum_j W_ij beta_j        W = column-centered dosage, beta_j ~ N(0, 1),
#                                     then rescaled so var(g) = h2
#   y_ijkl = mu + trial_j + rep_k + block_l + g_i + e,   var(e) = 1 - h2
# over the related population from make_related_dosage().
#
# `prop_unphenotyped` of the lines appear in the marker matrix (hence in G) but in no
# trial, so Stage 2 has genuinely unphenotyped candidates to predict -- and their true
# genetic values are known, which makes an honest accuracy check possible.
#
# @return list(D, beta, g, ph, phenotyped, unphenotyped, h2): `D` the dosage matrix
#   (build G from THIS, or the kernel and the phenotypes are decoupled), `g` the named
#   true genetic values, `ph` the long phenotype table stage1_blues() expects.
simulate_trials <- function(n_acc = 200, n_mrk = 2000, h2 = 0.5, n_trials = 2,
                            n_rep = 2, n_block = 4, prop_unphenotyped = 0.2,
                            n_founders = 20, n_fam = 20,
                            trait = "Yield", mu = 100, seed = 1) {
  D <- make_related_dosage(n_acc, n_mrk, n_founders, n_fam, seed = seed)
  set.seed(seed + 1L)

  # True breeding values: every marker contributes. Centering matters -- it is what
  # makes VanRaden's G the covariance of these g's.
  W    <- sweep(D, 2, colMeans(D), "-")
  beta <- stats::rnorm(n_mrk)
  g    <- as.vector(W %*% beta)
  g    <- (g - mean(g)) / stats::sd(g) * sqrt(h2)  # var(g) = h2 exactly
  names(g) <- rownames(D)

  lines <- rownames(D)
  n_out <- round(prop_unphenotyped * n_acc)
  unphenotyped <- if (n_out > 0) sort(sample(lines, n_out)) else character(0)
  phenotyped   <- setdiff(lines, unphenotyped)

  trial_eff <- stats::rnorm(n_trials, 0, 1)
  rep_eff   <- stats::rnorm(n_rep, 0, 0.3)
  block_eff <- stats::rnorm(n_block, 0, 0.3)

  ph <- tidyr::expand_grid(studyDbId = paste0("s", seq_len(n_trials)),
                           rep = seq_len(n_rep),
                           germplasmName = phenotyped) |>
    # Blocks are RE-RANDOMIZED within each trial x rep, as in a real incomplete-block
    # design. Giving a line the same block everywhere would confound block with
    # genotype and make lme4's block variance unidentifiable (degenerate Hessian).
    dplyr::group_by(studyDbId, rep) |>
    dplyr::mutate(block = sample(rep_len(seq_len(n_block), dplyr::n()))) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      trait = trait,
      value = mu + trial_eff[as.integer(sub("^s", "", studyDbId))] +
        rep_eff[rep] + block_eff[block] + g[germplasmName] +
        stats::rnorm(dplyr::n(), 0, sqrt(1 - h2))) |>
    dplyr::select(studyDbId, germplasmName, rep, block, trait, value)

  list(D = D, beta = beta, g = g, ph = ph,
       phenotyped = phenotyped, unphenotyped = unphenotyped, h2 = h2)
}

# Fake BrAPI connection whose $get(call, ...) returns canned locations/studies/
# observationunits. `locations`/`studies`/`units_by_study` are plain lists shaped
# like BrAPI $data elements.
fake_conn <- function(locations = list(), studies = list(), units_by_study = list()) {
  list(get = function(call, query = NULL, page = NULL, pageSize = NULL) {
    if (grepl("/locations", call)) return(list(data = locations))
    if (grepl("/studies",   call)) return(list(data = studies))
    if (grepl("/observationunits", call)) {
      sid <- as.character(query$studyDbId)
      return(list(data = units_by_study[[sid]] %||% list()))
    }
    stop("fake_conn: unexpected call ", call)
  })
}

loc_rec <- function(id, name, lon, lat, country = "USA") {
  list(locationDbId = id, locationName = name, countryCode = country,
       coordinates = list(geometry = list(coordinates = list(lon, lat, 100))))
}
study_rec <- function(id, name, type, loc, year) {
  list(studyDbId = id, studyName = name, studyType = type,
       locationDbId = loc, seasons = list(as.character(year)), trialName = name)
}
# One observation unit with embedded observations for the given traits/values.
unit_rec <- function(ou, gdbid, gname, rep, block, x, y, traits, values) {
  list(
    germplasmDbId = gdbid, germplasmName = gname, observationUnitDbId = ou,
    observationUnitPosition = list(
      positionCoordinateX = x, positionCoordinateY = y, entryType = "test",
      observationLevelRelationships = list(
        list(levelName = "rep",   levelCode = as.character(rep)),
        list(levelName = "block", levelCode = as.character(block)),
        list(levelName = "plot",  levelCode = ou))),
    observations = Map(function(t, v)
      list(observationVariableName = t, value = as.character(v)), traits, values))
}
