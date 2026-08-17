# Evaluation tooling: step through the pipeline one module at a time.
#
# This file is NOT part of the pipeline. Nothing in code/ sources it -- it is loaded
# by hand when you sit down to evaluate the code (see EVALUATION.md §3). That is
# deliberate: arm_evaluation() calls debug(), and a function left armed will block
# the next background run or wflow_build() on the debugger prompt.
#
#   arm_evaluation("stage2")   # debug() every function in the group
#   ... run something, step with n / s / c / Q ...
#   disarm_evaluation()        # ALWAYS, before any unattended run
#
# peek(x) prints a one-line health summary of whatever is flowing between steps and
# returns it unchanged, so it can be dropped inline.

library(tidyverse)
here::i_am("code/evaluation.R")

# --- the groups --------------------------------------------------------------

# Function names, grouped the way you actually evaluate them: one group per
# conceptual module, not one per file (step 4's file alone spans marker parsing,
# GRM building, the pedigree stitch and the combine, which are four separate
# things to reason about). Names that do not exist in the session are skipped, so
# optional pieces (canonicalize_to_primary comes from T3_brapi_helpers, sommer
# helpers need the package) never break arming.
EVAL_GROUPS <- list(
  # --- offline ---------------------------------------------------------------
  config     = c("config_lines", "config_traits", ".config_path", ".config_raw"),
  progress   = c("step_start", "step_done", "say", "note", "note_cache",
                 "pb_start", "pb_tick", "pb_done", "pb_wrap", "time_it",
                 "print_timings"),
  grm        = c(".Gmatrix", "std_grm", ".effective_n", ".mean_impute",
                 ".impute_glmnet"),
  combine    = c("EMCovarianceCombiner", "initialize_psi",
                 "compute_conditional_expectation", "compute_log_likelihood",
                 ".center_dfs"),
  markers    = c(".vcf_to_dosage", ".merge_dosage", ".qc_markers", ".thin_vcf",
                 ".thin_to_target", ".count_markers", ".read_vcf_samples"),
  stage1     = c("stage1_blues", ".fit_one"),
  stage2     = c("stage2_gblup", ".predict_engine", ".predict_one",
                 ".predict_one_sommer", ".geno_kernel", ".genotype_means",
                 ".sommer_blups"),
  cv         = c("cv_accuracy"),
  select     = c("select_parents", ".block_matrix"),

  # --- online (needs a live, logged-in connection) ---------------------------
  connect    = c("connect_t3", "t3_login", ".load_project_renviron"),
  trials     = c("find_ny_trials", "get_locations", "get_studies", "cache_key"),
  phenotypes = c("get_phenotypes", "split_by_role", ".level_code"),
  coverage   = c(".safe_coverage", ".coverage_table", "build_alias_lookup",
                 "canonicalize_to_primary", ".synonyms_available"),
  pedigree   = c(".germplasm_name_map", ".resolve_names_to_dbids",
                 ".pedigree_group_members", ".read_pedigree_group",
                 ".pedigree_partials", ".require_pedigree_cols",
                 ".pedigree_contract_error"),
  genotyping = c("find_and_get_genotypes", ".genotype_partials", ".combine_to_G",
                 ".protocol_grm", ".download_protocol_files")
)

# Aggregate: everything step 4 touches. Handy once the individual groups are clean
# and you want to watch one protocol go all the way through.
EVAL_GROUPS$step4 <- unique(c(EVAL_GROUPS$markers, EVAL_GROUPS$grm,
                              EVAL_GROUPS$pedigree, EVAL_GROUPS$combine,
                              EVAL_GROUPS$genotyping))

# Which groups need a live server -- drives eval_groups()' Off/online column.
EVAL_ONLINE <- c("connect", "trials", "phenotypes", "coverage", "pedigree",
                 "genotyping", "step4")

# Fast -> slow. Offline modules first so a cheap bug surfaces before an expensive
# one; within the live layer, the cheapest server work first (see EVALUATION.md §5).
EVAL_ORDER <- c("config", "progress", "grm", "combine", "markers",
                "stage1", "stage2", "cv", "select",
                "connect", "trials", "phenotypes", "coverage", "pedigree",
                "genotyping")

# --- arming ------------------------------------------------------------------

.eval_members <- function(groups) {
  unknown <- setdiff(groups, names(EVAL_GROUPS))
  if (length(unknown)) {
    stop("unknown evaluation group(s): ", paste(unknown, collapse = ", "),
         "\n  known: ", paste(names(EVAL_GROUPS), collapse = ", "), call. = FALSE)
  }
  unique(unlist(EVAL_GROUPS[groups], use.names = FALSE))
}

#' debug() every function in one or more groups.
#'
#' debug(), not debugonce(): the flag persists until undebug()'d, and debugonce()
#' can be neither cancelled nor detected with isdebugged(). The exists()/isdebugged()
#' guards make this idempotent, and silently skip functions that are not loaded.
#'
#' @param groups Group names (see eval_groups()).
#' @return The number of functions newly armed, invisibly.
arm_evaluation <- function(groups) {
  fns <- .eval_members(groups)
  armed <- 0L
  for (fn in fns) {
    if (exists(fn, mode = "function") && !isdebugged(get(fn))) {
      debug(get(fn)); armed <- armed + 1L
    }
  }
  message("armed ", armed, " function(s) in group(s): ", paste(groups, collapse = ", "),
          "\n  debugger keys:  n = next line   s = step into   c = finish this call   Q = quit",
          "\n  disarm_evaluation() when done -- ALWAYS before an unattended run.")
  invisible(armed)
}

#' undebug() every function in the given groups (default: all of them).
#' @return The number of functions disarmed, invisibly.
disarm_evaluation <- function(groups = names(EVAL_GROUPS)) {
  fns <- .eval_members(groups)
  off <- 0L
  for (fn in fns) {
    if (exists(fn, mode = "function") && isdebugged(get(fn))) {
      undebug(get(fn)); off <- off + 1L
    }
  }
  message("disarmed ", off, " function(s).")
  invisible(off)
}

#' Which functions are currently armed (across every group).
armed_functions <- function() {
  fns <- .eval_members(names(EVAL_GROUPS))
  fns[map_lgl(fns, ~ exists(.x, mode = "function") && isdebugged(get(.x)))]
}

#' Print the menu of groups, in evaluation order.
eval_groups <- function() {
  ordered <- c(EVAL_ORDER, setdiff(names(EVAL_GROUPS), EVAL_ORDER))
  out <- tibble(
    group   = ordered,
    n_fns   = map_int(ordered, ~ length(EVAL_GROUPS[[.x]])),
    loaded  = map_int(ordered, ~ sum(map_lgl(EVAL_GROUPS[[.x]],
                                             \(f) exists(f, mode = "function")))),
    layer   = if_else(ordered %in% EVAL_ONLINE, "online", "offline"))
  print(out, n = nrow(out))
  arm <- armed_functions()
  if (length(arm)) {
    message("NOTE: ", length(arm), " function(s) currently armed: ",
            paste(utils::head(arm, 6), collapse = ", "),
            if (length(arm) > 6) ", ..." else "")
  }
  invisible(out)
}

# --- peek: one-line health summaries -----------------------------------------

.pk <- function(...) cat(paste0(...), "\n", sep = "")
.rng <- function(v) {
  v <- v[is.finite(v)]
  if (!length(v)) return("none finite")
  sprintf("%.4g .. %.4g", min(v), max(v))
}

#' One-line health summary of an object flowing between pipeline steps.
#'
#' Recognizes the five shapes this pipeline moves around and calls out the
#' failure signatures that otherwise pass silently. Returns `x` unchanged
#' (invisibly), so it can be dropped inline: `blues <- peek(stage1_blues(ph))`.
#'
#' @param x          pheno tibble / dosage matrix / relationship matrix / BLUEs /
#'   GEBV tibble.
#' @param accessions For a dosage matrix: the names you EXPECT to match, so the
#'   rowname overlap (the synonym/name-mismatch signal) can be reported.
peek <- function(x, accessions = NULL) {
  if (is.matrix(x))                                    .peek_matrix(x, accessions)
  else if (is.data.frame(x))                           .peek_df(x)
  else if (is.list(x) && all(map_lgl(x, is.matrix)))   walk2(x, names(x) %||% seq_along(x),
                                                             ~ .peek_matrix(.x, accessions, .y))
  else if (is.numeric(x))                              .peek_numeric(x)
  else .pk("peek: ", class(x)[1], " of length ", length(x))
  invisible(x)
}

.peek_numeric <- function(v) {
  .pk("numeric: n=", length(v), " distinct=", n_distinct(v),
      " NA=", sum(is.na(v)), " range=", .rng(v))
  if (n_distinct(v[is.finite(v)]) == 1) .pk("  !! every finite value identical -- degenerate")
}

# A dosage matrix (accessions x markers, 0/1/2) and a relationship matrix are told
# apart by squareness + matching dimnames, since they need different diagnostics.
.peek_matrix <- function(m, accessions = NULL, label = NULL) {
  is_rel <- nrow(m) == ncol(m) && !is.null(rownames(m)) &&
    identical(rownames(m), colnames(m))
  tag <- if (is.null(label)) "" else paste0("[", label, "] ")
  if (is_rel) {
    d  <- diag(m); off <- m[upper.tri(m)]
    sym <- isTRUE(all.equal(m, t(m), tolerance = 1e-8))
    .pk(tag, "relationship matrix: ", nrow(m), " x ", ncol(m),
        " symmetric=", sym, " diag mean=", sprintf("%.3f", mean(d)),
        " diag=", .rng(d), " offdiag=", .rng(off))
    if (!sym) .pk("  !! NOT symmetric -- the EM combine or a relabel went wrong")
    if (anyNA(m)) .pk("  !! ", sum(is.na(m)), " NA cells")
    # Injected training accessions: diagonal at the mean, off-diagonals exactly 0.
    inj <- which(rowSums(abs(m)) - abs(d) == 0)
    if (length(inj)) {
      .pk("  ", length(inj), " row(s) with ALL off-diagonals exactly 0 (injected, ",
          "unrelated to everyone): ", paste(utils::head(rownames(m)[inj], 4), collapse = ", "))
    }
    if (abs(mean(d) - 1) < 1e-8) {
      .pk("  note: mean diagonal is exactly 1 -- std_grm()'d, or a pedigree/RKHS ",
          "kernel. A RAW VanRaden GRM on inbred oats should be > 1 (1 + F).")
    }
  } else {
    vals <- sort(unique(as.vector(m[!is.na(m)])))
    .pk(tag, "dosage matrix: ", nrow(m), " accessions x ", ncol(m), " markers  missing=",
        sprintf("%.1f%%", 100 * mean(is.na(m))),
        "  values={", paste(utils::head(vals, 6), collapse = ","), "}")
    if (length(vals) && (min(vals) < 0 || max(vals) > 2)) {
      .pk("  !! values outside {0,1,2} -- .Gmatrix()'s p = colMeans(M)/2 assumes 0/1/2")
    }
    if (!is.null(accessions)) {
      ov <- length(intersect(rownames(m), accessions))
      .pk("  rowname overlap with ", length(accessions), " expected: ", ov)
      if (ov == 0 && nrow(m) > 0) {
        .pk("  !! ZERO overlap on a non-empty matrix -- synonym / name-mismatch ",
            "(VCF samples under a preliminary name); see USE_SYNONYMS")
      }
    }
  }
}

.peek_df <- function(d) {
  nm <- names(d)
  if (all(c("trait", "genotype", "GEBV") %in% nm))          .peek_gebv(d)
  else if (all(c("genotype", "BLUE", "SE") %in% nm))        .peek_blues(d)
  else if (all(c("trait", "value", "studyDbId") %in% nm))   .peek_pheno(d)
  else .pk("tibble: ", nrow(d), " x ", ncol(d), " [", paste(nm, collapse = ", "), "]")
  .na_report(d)
}

# An all-NA column is the canonical silent failure: rep/block all NA means the
# Stage-1 model quietly degrades to plot means.
.na_report <- function(d) {
  na <- map_int(d, ~ sum(is.na(.x)))
  bad <- names(na)[na == nrow(d) & nrow(d) > 0]
  hit <- na[na > 0]
  if (length(hit)) {
    .pk("  NA per column: ",
        paste(sprintf("%s=%d", names(hit), hit), collapse = "  "))
  }
  if (length(bad)) .pk("  !! ALL NA: ", paste(bad, collapse = ", "))
}

.peek_pheno <- function(d) {
  .pk("phenotypes: ", nrow(d), " obs  studies=", n_distinct(d$studyDbId),
      " traits=", n_distinct(d$trait),
      " accessions=", n_distinct(d$germplasmName),
      "  value=", .rng(d$value))
}

.peek_blues <- function(d) {
  .pk("BLUEs: ", nrow(d), " rows  traits=", n_distinct(d$trait),
      " studies=", n_distinct(d$studyDbId), " genotypes=", n_distinct(d$genotype),
      "  BLUE=", .rng(d$BLUE), "  SE=", .rng(d$SE))
  if ("weight" %in% names(d)) {
    w <- d$weight[is.finite(d$weight)]
    .pk("  weight (1/SE^2)=", .rng(d$weight),
        if (length(w) && max(w) / max(min(w), 1e-12) > 1e6)
          "   !! spans >1e6 -- one trial can dominate Stage 2 (see SE_FLOOR_FRAC)" else "")
    if (any(!is.finite(d$weight))) .pk("  !! non-finite weights")
  }
}

.peek_gebv <- function(d) {
  n_pred <- if ("phenotyped" %in% names(d)) sum(!d$phenotyped) else NA_integer_
  .pk("GEBVs: ", nrow(d), " rows  traits=", n_distinct(d$trait),
      " genotypes=", n_distinct(d$genotype),
      "  predicted-only=", n_pred, " of ", n_distinct(d$genotype))
  if (!is.na(n_pred) && n_pred == 0) {
    .pk("  !! every candidate was phenotyped -- nothing NEW is being predicted")
  }
  d |>
    group_by(trait) |>
    summarise(range = .rng(GEBV), sd = sd(GEBV, na.rm = TRUE), .groups = "drop") |>
    pwalk(function(trait, range, sd) {
      .pk("  ", trait, ": ", range, if (isTRUE(sd < 1e-10)) "   !! constant GEBVs" else "")
    })
}
