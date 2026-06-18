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
  list(dosage = dosage, grm = std_grm(.Gmatrix(dosage)), n = nrow(dosage))
}

# --- Pedigree partials -------------------------------------------------------

# Reconstruct a dense relationship matrix from a sparse-triplet group CSV
# (germplasmDbId_i, germplasmDbId_j, relationship; upper triangle, symmetric).
.read_pedigree_group <- function(csv) {
  d <- readr::read_csv(csv, show_col_types = FALSE)
  ids <- as.character(union(d$germplasmDbId_i, d$germplasmDbId_j))
  A <- matrix(0, length(ids), length(ids), dimnames = list(ids, ids))
  i <- match(as.character(d$germplasmDbId_i), ids)
  j <- match(as.character(d$germplasmDbId_j), ids)
  A[cbind(i, j)] <- d$relationship
  A[cbind(j, i)] <- d$relationship
  A
}

# Load precomputed pedigree group matrices that overlap our accessions, relabel
# them by germplasmName, and return them as EM-ready standardized partials.
# dbid_to_name maps germplasmDbId -> germplasmName for our accessions.
.pedigree_partials <- function(pedigree_dir, dbid_to_name) {
  if (is.null(pedigree_dir) || !dir.exists(pedigree_dir)) {
    if (!is.null(pedigree_dir)) {
      message("Pedigree dir not found (", pedigree_dir, "); skipping pedigree stitch.")
    }
    return(list())
  }
  csvs <- list.files(pedigree_dir, pattern = "_group\\d+\\.csv$", full.names = TRUE)
  partials <- list()
  for (csv in csvs) {
    A <- .read_pedigree_group(csv)
    keep <- rownames(A) %in% names(dbid_to_name)      # overlap with our accessions
    if (sum(keep) < 2) next
    A <- A[keep, keep, drop = FALSE]
    rownames(A) <- colnames(A) <- unname(dbid_to_name[rownames(A)])
    partials[[length(partials) + 1L]] <- std_grm(A)
  }
  if (length(partials)) {
    message("Pedigree: ", length(partials),
            " group matrix(es) overlap our accessions.")
  } else {
    message("Pedigree: no precomputed group matrices overlap our accessions.")
  }
  partials
}

# --- Driver ------------------------------------------------------------------

#' Identify covering genotyping protocols/projects and build a combined GRM.
#'
#' @param conn       BrAPI connection.
#' @param accessions tibble with germplasmDbId and germplasmName (e.g. the
#'   `accessions` element from get_phenotypes()).
#' @param protocol_id NULL = combine ALL covering protocols (EM); an id = that
#'   protocol only.
#' @return list(protocols, projects, protocol_ids, G, markers). `G` is the
#'   combined GRM (accession x accession, by germplasmName); `markers` is the raw
#'   dosage matrix only for the single-protocol/no-pedigree case (else NULL).
#'   Cached to data/genotypes.rds.
find_and_get_genotypes <- function(conn, accessions,
                                   protocol_id = GENO_PROTOCOL_ID,
                                   pedigree_dir = PEDIGREE_DIR,
                                   refresh = FALSE) {
  cache <- cache_path("genotypes.rds")
  if (!refresh && file.exists(cache)) return(read_rds(cache))

  acc <- accessions |>
    transmute(germplasmDbId = as.character(germplasmDbId),
              germplasmName = as.character(germplasmName)) |>
    filter(!is.na(germplasmDbId)) |>
    distinct()
  dbids       <- unique(acc$germplasmDbId)
  our_names   <- unique(acc$germplasmName)
  dbid_to_name <- setNames(acc$germplasmName, acc$germplasmDbId)

  protocols <- .safe_coverage(conn$filter_geno_protocols, dbids, "protocol")
  projects  <- .safe_coverage(conn$filter_geno_projects,  dbids, "project")
  if (nrow(protocols) == 0) {
    stop("No genotyping protocol covers any of the supplied accessions.")
  }

  proto_ids <- if (is.null(protocol_id)) protocols$dbId[protocols$n_covered > 0]
               else as.character(protocol_id)

  # Download archived VCFs per protocol (cached); keep protocols that yielded files.
  files_by_proto <- set_names(map(proto_ids, ~ .download_protocol_files(conn, .x)), proto_ids)
  files_by_proto <- files_by_proto[lengths(files_by_proto) > 0]
  if (!length(files_by_proto)) {
    stop("Could not download any archived VCF for the selected protocol(s). ",
         "Set GENO_PROTOCOL_ID to a smaller protocol (e.g. an Oat 3K array).")
  }

  # Bridges: accessions (by name) genotyped in >= 2 of the chosen protocols.
  samples_by_proto <- map(files_by_proto, ~ unique(unlist(lapply(.x, .read_vcf_samples))))
  all_samples <- unlist(samples_by_proto, use.names = FALSE)
  bridges <- unique(all_samples[duplicated(all_samples)])
  keep_samples <- union(our_names, bridges)

  # One standardized GRM per protocol.
  proto_grms <- list(); proto_dosage <- list(); dfs <- integer(0); used <- character(0)
  for (pid in names(files_by_proto)) {
    message("Building GRM for protocol ", pid)
    res <- .protocol_grm(files_by_proto[[pid]], keep_samples)
    if (is.null(res)) { message("  no usable markers; skipping"); next }
    proto_grms[[pid]]   <- res$grm
    proto_dosage[[pid]] <- res$dosage
    dfs  <- c(dfs, res$n)
    used <- c(used, pid)
  }
  if (!length(proto_grms)) stop("No protocol produced a usable marker GRM.")

  # Pedigree partials (precomputed group matrices overlapping our accessions).
  ped_partials <- .pedigree_partials(pedigree_dir, dbid_to_name)

  partials <- c(unname(proto_grms), ped_partials)
  df_all   <- c(dfs, rep(PEDIGREE_DF, length(ped_partials)))

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

  # Raw marker matrix only when a single protocol and no pedigree partial was
  # combined (so marker-effect BGLR models remain possible).
  markers <- if (length(used) == 1L && length(ped_partials) == 0L) proto_dosage[[used]] else NULL

  out <- list(protocols    = protocols,
              projects     = projects,
              protocol_ids = used,
              G            = G,
              markers      = markers)
  write_rds(out, cache)
  out
}
