# Step 7: rank candidates and write the breeders' output spreadsheet.
#
# The selection index is a raw, UN-standardized weighted sum of per-trait GEBVs:
# GEBVs are already appropriately shrunk given their uncertainty, so standardizing
# them across traits would distort that and is not done here. The weights
# (TRAIT_WEIGHTS) carry both direction (sign) and relative scale/importance.
#
# Output: output/breeders_output.csv -- an intentionally spreadsheet-ish file with
# two side-by-side blocks separated by one blank column.
#   Block 1 (one row per predicted accession, sorted by Index, descending):
#     accession | In_Training (1/0) | <one GEBV column per trait> | Index
#   Block 2 (one row per trait):
#     Trait (short name) | Weight | CV_accuracy (cross-validated prediction accuracy)
# When TRAIT_WEIGHTS is NULL there is no index: block 1 drops the Index column (rows
# sorted by accession) and block 2 drops the Weight column.

library(tidyverse)
here::i_am("code/07_select.R")
source(here::here("code", "config.R"))   # TRAIT_* defaults, output_path()

# Coerce a data frame to a character matrix with its column names as the first
# row; NA -> "". Used to lay the two blocks side by side.
.block_matrix <- function(df) {
  body <- vapply(df, function(col) {
    col <- as.character(col); col[is.na(col)] <- ""; col
  }, character(nrow(df)))
  if (nrow(df) == 1L) body <- matrix(body, nrow = 1L)   # vapply drops 1-row to a vector
  rbind(colnames(df), body)
}

#' Rank candidates and write the breeders' output CSV.
#'
#' @param gebv          stage2_gblup() output (trait, genotype, GEBV, phenotyped).
#' @param cv            cv_accuracy() rows (trait, accuracy, ...); fills block 2's
#'                      CV_accuracy (NA where a trait is absent). Optional.
#' @param n_select      NULL = report every predicted accession; an integer keeps the
#'                      top n by selection index.
#' @param trait_weights named numeric (full trait name -> index weight); NULL skips
#'   the selection index entirely (unranked per-trait GEBVs, no Weight column).
#' @param trait_names   trait order (full names); GEBV columns and trait rows follow it.
#' @param trait_short   optional named vector (full name -> short label) for column
#'                      headers; full names are used where a label is missing.
#' @param outfile       CSV path.
#' @return invisibly list(parents = <block 1 tibble>, traits = <block 2 tibble>);
#'   also writes `outfile`.
select_parents <- function(gebv, cv = NULL, n_select = NULL,
                           trait_weights = TRAIT_WEIGHTS,
                           trait_names = TRAIT_NAMES,
                           trait_short = TRAIT_SHORT_NAMES,
                           outfile = output_path("breeders_output.csv")) {
  if (is.null(trait_names)) trait_names <- unique(gebv$trait)
  traits <- intersect(trait_names, unique(gebv$trait))   # config order, present only
  if (!length(traits)) stop("None of TRAIT_NAMES are present in gebv$trait.")
  short_of <- function(t) {
    s <- if (!is.null(trait_short)) trait_short[[t]] else NULL
    if (is.null(s) || is.na(s)) t else s
  }
  shorts <- unname(vapply(traits, short_of, character(1)))   # unnamed: used as col names

  g <- filter(gebv, trait %in% traits)
  # With no trait weights there is no selection index: the output is just the
  # per-trait GEBVs, unranked (e.g. an all-traits prediction run where selection
  # weights have not been chosen yet). The Index column and block-2 Weight column
  # are both omitted in that case.
  has_weights <- !is.null(trait_weights)

  # Per-accession: training flag + (when weighted) raw weighted-sum selection index.
  idx <- g |>
    group_by(genotype) |>
    summarise(In_Training = as.integer(any(phenotyped)),
              Index = if (has_weights) sum(trait_weights[trait] * GEBV) else NA_real_,
              .groups = "drop")

  # Wide GEBVs, columns in trait order, relabeled to short names.
  wide <- g |>
    select(genotype, trait, GEBV) |>
    pivot_wider(names_from = trait, values_from = GEBV) |>
    select(genotype, all_of(traits))
  names(wide) <- c("genotype", shorts)

  block1 <- left_join(idx, wide, by = "genotype")
  if (has_weights) {
    block1 <- block1 |>
      arrange(desc(Index)) |>
      select(accession = genotype, In_Training, all_of(shorts), Index)
  } else {
    block1 <- block1 |>
      arrange(genotype) |>
      select(accession = genotype, In_Training, all_of(shorts))
  }
  if (!is.null(n_select)) block1 <- head(block1, n_select)
  block1 <- mutate(block1, across(all_of(c(shorts, if (has_weights) "Index")),
                                  ~ round(.x, 4)))

  # Block 2: per-trait weight (omitted when unweighted) and cross-validated accuracy.
  cv_acc <- if (!is.null(cv)) setNames(cv$accuracy, cv$trait) else NULL
  block2 <- tibble(Trait = unname(shorts))
  if (has_weights) block2$Weight <- unname(trait_weights[traits])
  block2$CV_accuracy <- if (is.null(cv_acc)) NA_real_ else round(unname(cv_acc[traits]), 4)

  # Lay the blocks side by side with one blank spacer column, padding the shorter
  # block with blank rows so the two heights match.
  b1 <- .block_matrix(block1)
  b2 <- .block_matrix(block2)
  H  <- max(nrow(b1), nrow(b2))
  pad <- function(m) if (nrow(m) < H) rbind(m, matrix("", H - nrow(m), ncol(m))) else m
  full <- cbind(pad(b1), matrix("", H, 1), pad(b2))
  readr::write_csv(as.data.frame(full, stringsAsFactors = FALSE), outfile,
                   col_names = FALSE)

  invisible(list(parents = block1, traits = block2))
}
