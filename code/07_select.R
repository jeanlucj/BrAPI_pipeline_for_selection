# Step 7: rank candidates and select the best accessions for crossing.
#
# Builds a selection index from the Stage-2 GEBVs (standardized per trait, then
# combined with per-trait weights), ranks all genotyped candidates, and flags
# top picks that are highly related to a better-ranked pick (so the parent set
# is not dominated by near-duplicates).

library(tidyverse)
here::i_am("code/07_select.R")
source(here::here("code", "grm_utils.R"))   # .Gmatrix() for the markers fallback

#' Select parents from Stage-2 GEBVs.
#'
#' @param gebv          stage2_gblup() output (trait, genotype, GEBV, phenotyped).
#' @param geno          step-4 list with $G (combined GRM) and/or $markers, used
#'                      for the relatedness/redundancy flag (optional).
#' @param n_select      number of parents to flag as selected.
#' @param trait_weights named numeric vector (trait -> weight; +1 higher-is-better,
#'                       -1 lower-is-better). Default: +1 for every trait.
#' @param relatedness_threshold  relationship above which two picks are "redundant".
#' @return ranked tibble (also written to output/selected_parents.csv).
select_parents <- function(gebv, geno = NULL, n_select = 20,
                           trait_weights = NULL, relatedness_threshold = 0.5,
                           outfile = output_path("selected_parents.csv")) {
  traits <- unique(gebv$trait)
  if (is.null(trait_weights)) trait_weights <- setNames(rep(1, length(traits)), traits)

  # per-trait standardized GEBV, weighted, summed into an index
  scored <- gebv |>
    group_by(trait) |>
    mutate(z = as.numeric(scale(GEBV))) |>
    ungroup() |>
    mutate(wz = z * trait_weights[trait])

  index <- scored |>
    group_by(genotype) |>
    summarise(index = sum(wz, na.rm = TRUE),
              phenotyped = any(phenotyped), .groups = "drop")

  # wide per-trait GEBVs for the report
  wide <- gebv |>
    select(trait, genotype, GEBV) |>
    pivot_wider(names_from = trait, values_from = GEBV,
                names_prefix = "GEBV_")

  ranked <- index |>
    left_join(wide, by = "genotype") |>
    arrange(desc(index)) |>
    mutate(rank = row_number(), selected = rank <= n_select, .before = 1)

  # flag redundancy among selected picks using the relationship matrix: the
  # combined GRM from step 4 if present, else one built from the marker matrix.
  ranked$redundant_with <- NA_character_
  Grel <- if (!is.null(geno$G)) geno$G
          else if (!is.null(geno$markers)) .Gmatrix(geno$markers) else NULL
  if (!is.null(Grel)) {
    sel <- ranked |> filter(selected, genotype %in% rownames(Grel))
    if (nrow(sel) > 1) {
      G <- Grel[sel$genotype, sel$genotype]
      for (i in 2:nrow(sel)) {
        better <- sel$genotype[seq_len(i - 1)]
        close  <- better[G[sel$genotype[i], better] > relatedness_threshold]
        if (length(close) > 0) {
          ranked$redundant_with[ranked$genotype == sel$genotype[i]] <-
            str_c(close, collapse = "; ")
        }
      }
    }
  }

  readr::write_csv(ranked, outfile)
  ranked
}
