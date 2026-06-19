# Step 4: find the genotyping protocols/projects that cover our accessions and
# build a combined genomic relationship matrix (GRM) for genomic prediction.
#
# Coverage: Breedbase ranks genotyping protocols AND projects by how many of a
# supplied accession list they genotype (conn$filter_geno_protocols / _projects).
#
# Markers: T3's BrAPI genotyping endpoints are not implemented and on-the-fly VCF
# generation (conn$vcf) is slow/unreliable, so we download pre-generated ARCHIVED
# VCFs (conn$vcf_archived). Large GBS files are thinned genome-wide to ~
# TARGET_DENSITY markers and read with only the needed sample columns to bound
# memory.
#
# Combining platforms: with GENO_PROTOCOL_ID = NULL we build one standardized GRM
# per covering protocol and EM-combine them (Wishart-EM, em_covariance_combiner.R)
# into a single GRM, using accessions genotyped on >= 2 platforms as bridges so
# the EM can impute cross-platform relationships. A precomputed pedigree
# relationship matrix (PEDIGREE_DIR) can be added as an extra partial to stitch
# disjoint marker GRMs. A single integer GENO_PROTOCOL_ID uses just that protocol
# (and also returns its raw marker matrix for marker-effect models).

library(tidyverse)
here::i_am("code/04_find_genotyping.R")
source(here::here("code", "01_connect.R"))
source(here::here("code", "grm_utils.R"))
source(here::here("code", "em_covariance_combiner.R"))

# --- Coverage ----------------------------------------------------------------

# Turn a filter_geno_* result's $content into a tidy coverage tibble.
# `type` is "protocol" or "project". Returns 0-row tibble if nothing matched
# (the BrAPI.R wrapper itself errors on empty results, so callers guard it).
.coverage_table <- function(content, type) {
  res     <- content$results
  ranked  <- res$counts[[sprintf("ranked_genotyping_%ss", type)]]
  lookups <- res$lookups[[sprintf("genotyping_%ss", type)]]
  counts  <- res$counts[[sprintf("accessions_by_genotyping_%s", type)]]
  total   <- res$counts$accessions_total %||% NA_integer_
  if (is.null(ranked) || length(ranked) == 0) {
    return(tibble(rank = integer(), dbId = character(), name = character(),
                  n_covered = integer(), n_total = integer()))
  }
  ids <- map_chr(ranked, as.character)
  tibble(
    rank      = seq_along(ids),
    dbId      = ids,
    name      = map_chr(ids, ~ lookups[[.x]] %||% NA_character_),
    n_covered = map_int(ids, ~ as.integer(counts[[.x]] %||% 0L)),
    n_total   = as.integer(total)
  )
}

# Safely run a coverage ranking, returning a 0-row tibble instead of erroring
# when no accession is genotyped.
.safe_coverage <- function(fn, accessions, type) {
  out <- tryCatch(fn(accessions), error = function(e) NULL)
  if (is.null(out)) return(.coverage_table(list(results = list()), type))
  .coverage_table(out$content, type)
}

# --- VCF reading / thinning --------------------------------------------------

# Sample names from a VCF's #CHROM header (no genotype data read).
.read_vcf_samples <- function(path) {
  con <- if (grepl("\\.gz$", path)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(con))
  repeat {
    line <- readLines(con, n = 1)
    if (!length(line)) stop("No #CHROM header found in ", basename(path))
    if (startsWith(line, "#CHROM")) return(strsplit(line, "\t")[[1]][-(1:9)])
  }
}

# Count data (variant) lines, streaming so memory is bounded.
.count_markers <- function(path, block = 100000L) {
  con <- if (grepl("\\.gz$", path)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(con))
  n <- 0L
  repeat {
    lines <- readLines(con, n = block)
    if (!length(lines)) break
    n <- n + sum(!startsWith(lines, "#"))
  }
  n
}

# Stream a VCF and write every n-th marker (variant line) to `out`, preserving
# all header lines. Ported from build_grm_for_cv00.R: thinning genome-wide bounds
# memory for projects with >1M markers. Returns markers written.
.thin_vcf <- function(path, n, out, block = 100000L) {
  con_in  <- if (grepl("\\.gz$", path)) gzfile(path, "rt") else file(path, "rt")
  con_out <- file(out, "wt")
  on.exit({ close(con_in); close(con_out) })
  data_seen <- 0L; kept <- 0L
  repeat {
    lines <- readLines(con_in, n = block)
    if (!length(lines)) break
    is_hdr <- startsWith(lines, "#")
    global_data <- data_seen + cumsum(!is_hdr)            # 1-based data-line index
    keep <- is_hdr | (!is_hdr & ((global_data - 1L) %% n == 0L))
    writeLines(lines[keep], con_out)
    data_seen <- data_seen + sum(!is_hdr)
    kept <- kept + sum(keep & !is_hdr)
  }
  kept
}

# Thin a VCF to ~target markers (keep every n-th, n = floor(n_markers/target)).
# Returns a path to read (a temp file if thinned, else the original) and the temp
# path to clean up (or NULL).
.thin_to_target <- function(path, target) {
  if (is.null(target) || target <= 0) return(list(path = path, tmp = NULL))
  nm <- .count_markers(path)
  n  <- max(1L, as.integer(floor(nm / target)))
  if (n <= 1L) return(list(path = path, tmp = NULL))
  tmp  <- tempfile(fileext = ".vcf")
  kept <- .thin_vcf(path, n, tmp)
  message(sprintf("  thinning %s: %d -> %d markers (every %dth)",
                  basename(path), nm, kept, n))
  list(path = tmp, tmp = tmp)
}

# VCF -> accession x marker dosage matrix (count of ALT allele), optionally
# thinned to target_density and restricted to keep_samples (by germplasmName) to
# bound memory. Returns NULL if none of keep_samples are present.
.vcf_to_dosage <- function(vcf_path, keep_samples = NULL,
                           target_density = TARGET_DENSITY) {
  th <- .thin_to_target(vcf_path, target_density)
  on.exit(if (!is.null(th$tmp)) unlink(th$tmp))

  cols <- NULL
  if (!is.null(keep_samples)) {
    samples <- .read_vcf_samples(vcf_path)   # original header valid after thinning
    keep_cols <- 9 + which(samples %in% keep_samples)
    if (!length(keep_cols)) return(NULL)
    cols <- c(1:9, keep_cols)
  }
  vcf <- if (is.null(cols)) vcfR::read.vcfR(th$path, verbose = FALSE)
         else                vcfR::read.vcfR(th$path, cols = cols, verbose = FALSE)
  if (is.null(vcf@gt) || ncol(vcf@gt) < 2) return(NULL)

  # vcfR uses the ID column for marker names and requires it unique; many T3 VCFs
  # leave ID as "." or repeat it, so build unique CHROM_POS ids.
  ids <- vcf@fix[, "ID"]
  if (any(is.na(ids)) || anyDuplicated(ids)) {
    vcf@fix[, "ID"] <- make.unique(paste(vcf@fix[, "CHROM"], vcf@fix[, "POS"], sep = "_"))
  }
  gt <- vcfR::extract.gt(vcf, element = "GT")            # markers x samples
  g <- matrix(NA_real_, nrow(gt), ncol(gt), dimnames = dimnames(gt))
  g[gt %in% c("0/0", "0|0")] <- 0
  g[gt %in% c("0/1", "1/0", "0|1", "1|0")] <- 1
  g[gt %in% c("1/1", "1|1")] <- 2
  t(g)                                                   # accessions x markers
}

# Merge dosage matrices (same protocol = same markers) by union of markers and
# row-binding accessions; de-duplicate accessions (keep first occurrence).
.merge_dosage <- function(mats) {
  if (length(mats) == 1) return(mats[[1]])
  all_markers <- unique(unlist(lapply(mats, colnames)))
  aligned <- lapply(mats, function(m) {
    full <- matrix(NA_real_, nrow(m), length(all_markers),
                   dimnames = list(rownames(m), all_markers))
    full[, colnames(m)] <- m
    full
  })
  M <- do.call(rbind, aligned)
  M[!duplicated(rownames(M)), , drop = FALSE]
}

# Marker QC: drop high-missing markers/accessions, MAF filter, mean-impute.
# Markers are filtered BEFORE accessions: when projects use different genome
# annotations (different marker names), merging leaves every accession ~50%
# missing, so dropping the minority-annotation markers first is what makes the
# accessions look complete again.
.qc_markers <- function(M, max_missing = MAX_MISSING, min_maf = MIN_MAF) {
  M <- M[, colMeans(is.na(M)) <= max_missing, drop = FALSE]
  M <- M[rowMeans(is.na(M)) <= max_missing, , drop = FALSE]
  af  <- colMeans(M, na.rm = TRUE) / 2
  maf <- pmin(af, 1 - af)
  M <- M[, !is.na(maf) & maf >= min_maf, drop = FALSE]
  col_means <- colMeans(M, na.rm = TRUE)
  na_idx <- which(is.na(M), arr.ind = TRUE)
  if (nrow(na_idx) > 0) M[na_idx] <- col_means[na_idx[, "col"]]
  M
}

# --- Marker download ---------------------------------------------------------

# Download a protocol's archived project VCF(s) to the on-disk cache (so re-runs
# skip the download). Returns the cached file paths; per-file failures are warned
# and skipped (so a multi-GB panel can't abort an all-protocols run).
.download_protocol_files <- function(conn, protocol_id) {
  arch <- tryCatch(
    conn$vcf_archived_list(genotyping_protocol_id = as.integer(protocol_id)),
    error = function(e) NULL)
  if (!is.data.frame(arch) || nrow(arch) == 0) return(character(0))

  vdir <- cache_path("vcf_cache")
  dir.create(vdir, showWarnings = FALSE, recursive = TRUE)
  paths <- character(0)
  for (i in seq_len(nrow(arch))) {
    dest <- file.path(vdir, paste0("proto", protocol_id, "_proj",
                                   arch$project_id[i], ".vcf"))
    if (!file.exists(dest) || file.size(dest) == 0) {
      message("  downloading protocol ", protocol_id, " project ",
              arch$project_id[i], " (", i, "/", nrow(arch),
              "); large GBS files may take a while...")
      ok <- tryCatch({
        conn$vcf_archived(output = dest,
                          genotyping_project_id = as.integer(arch$project_id[i]),
                          file_name = arch$file_name[i])
        TRUE
      }, error = function(e) FALSE)
      if (!ok || !file.exists(dest) || file.size(dest) == 0) {
        warning("Skipping project ", arch$project_id[i], " of protocol ",
                protocol_id, " (download failed; may be too large).")
        unlink(dest)
        next
      }
    }
    paths <- c(paths, dest)
  }
  paths
}

# Build one standardized GRM (+ its dosage) for a protocol from its cached VCFs.
.protocol_grm <- function(files, keep_samples) {
  mats <- list()
  for (f in files) {
    m <- tryCatch(.vcf_to_dosage(f, keep_samples = keep_samples),
                  error = function(e) { warning("parse failed: ", basename(f),
                                                " (", conditionMessage(e), ")"); NULL })
    if (!is.null(m) && nrow(m) > 0 && ncol(m) > 0) mats[[length(mats) + 1L]] <- m
  }
  if (!length(mats)) return(NULL)
  dosage <- .qc_markers(.merge_dosage(mats))
  if (nrow(dosage) < 2 || ncol(dosage) < 2) return(NULL)
  grm <- std_grm(.Gmatrix(dosage))
  # m_eff = effective # of independent samples (Galwey), used to set this GRM's EM
  # degrees of freedom relative to the other partials.
  list(dosage = dosage, grm = grm, m_eff = .effective_n(grm))
}

# --- Pedigree partials -------------------------------------------------------

# Files in pedigree_dir come from the sibling BrAPI_pedigree_relmat project (see
# CLAUDE.md, Cross-repo dependencies). A *missing* companion degrades gracefully;
# a present-but-malformed one is a contract violation -- fail loudly with an
# actionable message rather than silently building a wrong/empty GRM.
.pedigree_contract_error <- function(file, detail) {
  stop("Unexpected shape in BrAPI_pedigree_relmat output '", file, "': ", detail,
       "\n  The sibling project's output contract may have changed; see CLAUDE.md ",
       "(Cross-repo dependencies) and update code/04_find_genotyping.R.",
       call. = FALSE)
}

# Assert a pedigree group CSV has the expected triplet columns.
.require_pedigree_cols <- function(d, csv) {
  req <- c("germplasmDbId_i", "germplasmDbId_j", "relationship")
  if (!all(req %in% names(d))) {
    .pedigree_contract_error(csv, paste0(
      "missing column(s) ", paste(setdiff(req, names(d)), collapse = ", "),
      "; expected ", paste(req, collapse = ", "), "."))
  }
  invisible(d)
}

# Reconstruct a dense relationship matrix from a sparse-triplet group CSV
# (germplasmDbId_i, germplasmDbId_j, relationship; upper triangle, symmetric).
# When keep_dbids is supplied, the triplets are filtered to pairs where BOTH
# endpoints are relevant *before* the dense matrix is allocated, so we never
# materialize the full group (some groups are large). Returns a matrix over only
# the relevant germplasmDbIds present in the CSV (possibly 0x0).
.read_pedigree_group <- function(csv, keep_dbids = NULL) {
  # Read every column as character so ids are preserved exactly as written (a
  # double parse would turn round-number dbIds into "1e+06" and break matching);
  # relationship is coerced back to numeric after the contract check.
  d <- readr::read_csv(csv, show_col_types = FALSE,
                       col_types = readr::cols(.default = readr::col_character()))
  .require_pedigree_cols(d, csv)
  d <- mutate(d, relationship = as.numeric(relationship))
  if (!is.null(keep_dbids)) {
    d <- filter(d, germplasmDbId_i %in% keep_dbids,
                   germplasmDbId_j %in% keep_dbids)
  }
  ids <- union(d$germplasmDbId_i, d$germplasmDbId_j)
  A <- matrix(0, length(ids), length(ids), dimnames = list(ids, ids))
  i <- match(d$germplasmDbId_i, ids)
  j <- match(d$germplasmDbId_j, ids)
  A[cbind(i, j)] <- d$relationship
  A[cbind(j, i)] <- d$relationship
  A
}

# germplasmDbId -> germplasmName for ALL pedigree accessions, from the companion
# germplasm cache the pedigree repo writes next to pedigree_dir
# (output/germplasm_cache_<id>.rds, a raw BrAPI germplasm dump). This is what lets
# a non-phenotyped accession be recognized as a pedigree bridge: pedigree CSVs are
# keyed by germplasmDbId but VCF samples (and keep_samples) are germplasmName.
# Returns a named character vector (names = dbId); empty if the cache is absent,
# in which case only our own accessions get pedigree treatment (legacy behavior).
.germplasm_name_map <- function(pedigree_dir) {
  if (is.null(pedigree_dir)) return(character(0))
  cache <- file.path(dirname(pedigree_dir),
                     paste0("germplasm_cache_", basename(pedigree_dir), ".rds"))
  if (!file.exists(cache)) {
    message("Germplasm cache not found (", cache,
            "); pedigree bridges limited to phenotyped accessions.")
    return(character(0))
  }
  recs <- read_rds(cache)
  if (!is.list(recs) || is.data.frame(recs) || !length(recs) ||
      !is.list(recs[[1]]) ||
      !all(c("germplasmDbId", "germplasmName") %in% names(recs[[1]]))) {
    .pedigree_contract_error(cache, paste0(
      "not a non-empty list of BrAPI germplasm records with ",
      "$germplasmDbId / $germplasmName."))
  }
  ids  <- map_chr(recs, ~ as.character(.x$germplasmDbId %||% NA))
  nms  <- map_chr(recs, ~ .x$germplasmName %||% NA_character_)
  ok   <- !is.na(ids) & !is.na(nms)
  set_names(nms[ok], ids[ok])
}

# Resolve germplasmName -> germplasmDbId via the pedigree germplasm cache (first
# occurrence wins for non-unique names). Returns a tibble(germplasmDbId,
# germplasmName) for the names that resolve; warns about any that don't (those can
# still be predicted if they appear as a VCF sample in a selected protocol).
.resolve_names_to_dbids <- function(names, pedigree_dir) {
  names <- unique(names[!is.na(names)])
  if (!length(names)) return(tibble(germplasmDbId = character(), germplasmName = character()))
  d2n <- .germplasm_name_map(pedigree_dir)          # dbId -> name
  n2d <- d2n[!duplicated(d2n)]                       # keep first dbId per name
  n2d <- set_names(names(n2d), unname(n2d))          # name -> dbId
  hit <- names[names %in% names(n2d)]
  miss <- setdiff(names, hit)
  if (length(miss)) {
    message(length(miss), " TEST_ACCESSIONS could not be resolved to a germplasmDbId ",
            "(predicted only if present in a selected protocol/pedigree by name): ",
            paste(utils::head(miss, 10), collapse = ", "),
            if (length(miss) > 10) ", ..." else "")
  }
  tibble(germplasmDbId = unname(n2d[hit]), germplasmName = hit)
}

# Per-group accession membership (germplasmDbId vectors) for the CSV-backed
# pedigree groups in pedigree_dir, needed to detect bridges BEFORE densifying.
# Reads the companion `<id>_pedigree_groups.rds` (cheap) rather than the group
# CSVs, which can be >1GB; falls back to the CSV's id columns if the rds lacks a
# group. Returns a named list: csv path -> character dbId vector.
.pedigree_group_members <- function(pedigree_dir) {
  if (is.null(pedigree_dir) || !dir.exists(pedigree_dir)) {
    if (!is.null(pedigree_dir)) {
      message("Pedigree dir not found (", pedigree_dir, "); skipping pedigree stitch.")
    }
    return(list())
  }
  csvs <- list.files(pedigree_dir, pattern = "_group\\d+\\.csv$", full.names = TRUE)
  groups_rds <- file.path(pedigree_dir,
                          paste0(basename(pedigree_dir), "_pedigree_groups.rds"))
  groups <- list()
  if (file.exists(groups_rds)) {
    groups <- read_rds(groups_rds)
    if (!is.list(groups) || is.null(names(groups))) {
      .pedigree_contract_error(groups_rds, paste0(
        "not a named list of germplasmDbId vectors (elements named 'group<N>')."))
    }
  }
  set_names(
    map(csvs, function(csv) {
      key <- sub("^.*_(group\\d+)\\.csv$", "\\1", basename(csv))
      ids <- groups[[key]]
      if (is.null(ids)) {                              # rds missing this group
        d <- readr::read_csv(csv, show_col_types = FALSE,
                             col_types = readr::cols(.default = readr::col_character()))
        .require_pedigree_cols(d, csv)
        ids <- union(d$germplasmDbId_i, d$germplasmDbId_j)
      }
      as.character(ids)
    }),
    csvs)
}

# Build EM-ready standardized pedigree partials, one per CSV-backed group,
# restricted to keep_names (our accessions + bridges) and relabeled by
# germplasmName via dbid_to_name (the full map from .germplasm_name_map overlaid
# with our own). group_members is the .pedigree_group_members() list. Groups are
# densified over only the kept germplasmDbIds (see .read_pedigree_group), so the
# 1GB group CSV never becomes a full dense matrix.
.pedigree_partials <- function(group_members, dbid_to_name, keep_names) {
  partials <- list()
  for (csv in names(group_members)) {
    dbids <- group_members[[csv]]
    nm    <- dbid_to_name[dbids]
    keep_dbids <- dbids[!is.na(nm) & nm %in% keep_names]
    if (length(keep_dbids) < 2) next
    A <- .read_pedigree_group(csv, keep_dbids)        # densified over overlap only
    if (nrow(A) < 2) next
    # Relabel dbId -> name; drop accessions with no name or a duplicate name
    # (T3 germplasmNames are not guaranteed unique -- first occurrence wins).
    nm2  <- unname(dbid_to_name[rownames(A)])
    keep <- !is.na(nm2) & !duplicated(nm2)
    A <- A[keep, keep, drop = FALSE]
    if (nrow(A) < 2) next
    rownames(A) <- colnames(A) <- nm2[keep]
    partials[[length(partials) + 1L]] <- std_grm(A)
  }
  if (length(partials)) {
    message("Pedigree: ", length(partials),
            " group matrix(es) overlap our accessions/bridges.")
  } else {
    message("Pedigree: no precomputed group matrices overlap our accessions/bridges.")
  }
  partials
}

# Turn per-GRM effective-sample-size values into EM degrees of freedom: preserve
# their relative ordering but re-center on `mean_df` and cap the spread at
# `sd_df` (the lesser of `sd_df` and the observed SD). This keeps the absolute
# marker-vs-pedigree weighting under explicit control rather than letting raw
# eigenvalue spectra dictate it. Floored at 1 so every df stays positive.
.center_dfs <- function(m_eff, mean_df, sd_df) {
  if (length(m_eff) <= 1) return(rep(mean_df, length(m_eff)))
  obs_sd <- stats::sd(m_eff)
  if (!is.finite(obs_sd) || obs_sd == 0) return(rep(mean_df, length(m_eff)))
  z <- (m_eff - mean(m_eff)) / obs_sd
  pmax(mean_df + z * min(sd_df, obs_sd), 1)
}

# --- Driver ------------------------------------------------------------------

#' Identify covering genotyping protocols/projects and build the prediction GRM.
#'
#' @param conn             BrAPI connection.
#' @param train_accessions tibble (germplasmDbId, germplasmName) for the training
#'   accessions. Always present in the returned `G` (injected if they have neither
#'   genotypes nor pedigree).
#' @param test_accessions  tibble (germplasmDbId, germplasmName) of extra accessions
#'   to predict (e.g. from TEST_TRIALS); predicted only if they end up in the GRM.
#' @param test_names       germplasmName vector of extra accessions to predict
#'   (TEST_ACCESSIONS); their dbIds are resolved via the pedigree germplasm cache so
#'   their protocols are also covered.
#' @param protocol_id NULL = combine ALL covering protocols (EM); an id = that
#'   protocol only.
#' @return list(protocols, projects, protocol_ids, G, markers). `G` is the combined
#'   GRM by germplasmName, SUBSET to the prediction targets (training + predictable
#'   test); `markers` is the raw dosage matrix only for the single-protocol/no-
#'   pedigree/no-injection case (else NULL). Cached to data/genotypes.rds.
find_and_get_genotypes <- function(conn, train_accessions,
                                   test_accessions = NULL,
                                   test_names = TEST_ACCESSIONS,
                                   protocol_id = GENO_PROTOCOL_ID,
                                   pedigree_dir = PEDIGREE_DIR,
                                   refresh = FALSE) {
  cache <- cache_path("genotypes.rds")
  if (!refresh && file.exists(cache)) return(read_rds(cache))

  clean <- function(a) {
    if (is.null(a) || !nrow(a)) return(tibble(germplasmDbId = character(),
                                              germplasmName = character()))
    a |> transmute(germplasmDbId = as.character(germplasmDbId),
                   germplasmName = as.character(germplasmName)) |>
      filter(!is.na(germplasmDbId)) |> distinct()
  }
  train <- clean(train_accessions)
  test  <- clean(test_accessions)
  resolved <- .resolve_names_to_dbids(test_names, pedigree_dir)   # TEST_ACCESSIONS -> dbId

  # Coverage / dbid_to_name span training + test (protocols cover both); the
  # prediction target names are everyone we ultimately want in G; training names
  # are force-included even if ungenotyped.
  acc <- bind_rows(train, test, resolved) |> distinct()
  dbids        <- unique(acc$germplasmDbId)
  dbid_to_name <- setNames(acc$germplasmName, acc$germplasmDbId)
  target_names <- unique(c(train$germplasmName, test$germplasmName,
                           resolved$germplasmName, test_names))
  target_names <- target_names[!is.na(target_names)]
  force_names  <- unique(train$germplasmName)
  our_names    <- target_names                        # build keep-set (+ bridges below)

  protocols <- .safe_coverage(conn$filter_geno_protocols, dbids, "protocol")
  projects  <- .safe_coverage(conn$filter_geno_projects,  dbids, "project")
  if (nrow(protocols) == 0) {
    stop("No genotyping protocol covers any of the supplied accessions.")
  }

  proto_ids <- if (is.null(protocol_id)) protocols$dbId[protocols$n_covered > 0]
               else as.character(protocol_id)

  # Download archived VCFs per protocol (cached); keep protocols that yielded files.
  files_by_proto <- set_names(
    map(proto_ids, ~ .download_protocol_files(conn, .x),
        .progress = "Genotypes: downloading VCFs by protocol"),
    proto_ids)
  files_by_proto <- files_by_proto[lengths(files_by_proto) > 0]
  if (!length(files_by_proto)) {
    stop("Could not download any archived VCF for the selected protocol(s). ",
         "Set GENO_PROTOCOL_ID to a smaller protocol (e.g. an Oat 3K array).")
  }

  # Pedigree accessions, discovered up front (by germplasmName) so they can act as
  # bridges: dbid_to_name is widened with the pedigree germplasm cache, and each
  # CSV-backed group's membership is read cheaply (no densifying).
  ped_d2n <- .germplasm_name_map(pedigree_dir)
  ped_d2n[names(dbid_to_name)] <- dbid_to_name        # our own names take precedence
  ped_members <- .pedigree_group_members(pedigree_dir)
  ped_group_names <- map(ped_members, function(ids) {
    nm <- ped_d2n[ids]; unique(nm[!is.na(nm)])
  })

  # Bridges: accessions (by name) appearing in >= 2 partials, counting each
  # genotyping protocol AND each pedigree group. This keeps an accession seen in
  # only one protocol when it also sits in a pedigree group, so it can bridge that
  # protocol's marker GRM to the pedigree relationship matrix.
  samples_by_proto <- map(files_by_proto, ~ unique(unlist(lapply(.x, .read_vcf_samples))))
  occ <- table(unlist(c(samples_by_proto, ped_group_names), use.names = FALSE))
  bridges <- names(occ)[occ >= 2]
  keep_samples <- union(our_names, bridges)

  # One standardized GRM per protocol.
  proto_grms <- list(); proto_dosage <- list(); m_effs <- numeric(0); used <- character(0)
  for (pid in names(files_by_proto)) {
    message("Building GRM for protocol ", pid)
    res <- .protocol_grm(files_by_proto[[pid]], keep_samples)
    if (is.null(res)) { message("  no usable markers; skipping"); next }
    proto_grms[[pid]]   <- res$grm
    proto_dosage[[pid]] <- res$dosage
    m_effs <- c(m_effs, res$m_eff)
    used   <- c(used, pid)
  }
  if (!length(proto_grms)) stop("No protocol produced a usable marker GRM.")

  # Pedigree partials (precomputed group matrices overlapping our accessions or
  # bridges), labeled by germplasmName via the widened map.
  ped_partials <- .pedigree_partials(ped_members, ped_d2n, keep_samples)

  # EM degrees of freedom: marker GRMs from their effective sample sizes
  # (re-centered on GRM_DF_MEAN, spread capped at GRM_DF_STDEV); pedigree fixed.
  grm_dfs  <- .center_dfs(m_effs, GRM_DF_MEAN, GRM_DF_STDEV)
  partials <- c(unname(proto_grms), ped_partials)
  df_all   <- c(grm_dfs, rep(PEDIGREE_DF, length(ped_partials)))

  if (length(partials) == 1L) {
    G <- partials[[1]]
  } else {
    combined_names <- unique(unlist(lapply(partials, colnames)))
    var_indices <- lapply(partials, function(g) match(colnames(g), combined_names))
    res <- EMCovarianceCombiner(partial_covs = partials, var_indices = var_indices,
                                degrees_freedom = df_all)
    G <- res$psi[-1, -1]                      # drop the phantom variable
    G <- (G + t(G)) / 2                        # clean tiny EM floating-point asymmetry
    rownames(G) <- colnames(G) <- combined_names
  }

  # Restrict to the prediction targets: bridges/extra accessions informed the EM
  # combine but are dropped from the matrix we predict on. Then force-include any
  # training accession that has neither genotypes nor pedigree, with diagonal =
  # mean diagonal of the kept matrix and zero off-diagonals (a prior-only line).
  keep <- intersect(target_names, rownames(G))
  G <- G[keep, keep, drop = FALSE]
  inject <- setdiff(force_names, rownames(G))
  if (length(inject)) {
    diag_val  <- if (nrow(G)) mean(diag(G)) else 1
    all_names <- c(rownames(G), inject)
    G2 <- matrix(0, length(all_names), length(all_names),
                 dimnames = list(all_names, all_names))
    if (nrow(G)) G2[rownames(G), colnames(G)] <- G   # existing block; cross terms 0
    G2[cbind(inject, inject)] <- diag_val            # injected self-relationships
    G <- G2
    message("Injected ", length(inject), " ungenotyped/unpedigreed training ",
            "accession(s) into G (diag = ", round(diag_val, 3), ", off-diag = 0).")
  }

  # Raw marker matrix only for the single-protocol, no-pedigree, no-injection case
  # (injected accessions have no marker row), subset to the retained accessions so
  # marker-effect BGLR models stay consistent with G.
  markers <- NULL
  if (length(used) == 1L && length(ped_partials) == 0L && !length(inject)) {
    m <- proto_dosage[[used]]
    markers <- m[rownames(m) %in% rownames(G), , drop = FALSE]
  }

  out <- list(protocols    = protocols,
              projects     = projects,
              protocol_ids = used,
              G            = G,
              markers      = markers)
  write_rds(out, cache)
  out
}
