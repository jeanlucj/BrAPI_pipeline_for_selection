# Diagnostic for the Freeze-damage GEBVs: are the extreme-GEBV accessions the ones
# that actually have marker/pedigree data, the ones with low GRM diagonals, or the
# phenotyped ones? Offline: reads cached genotypes.rds + gebv.rds.

library(tidyverse)
here::i_am("code/gebv_diagnostic.R")
source(here::here("code", "config.R"))

geno  <- read_rds(cache_path("genotypes.rds"))
gebv  <- read_rds(cache_path("gebv.rds"))
G     <- geno$G
trait <- "Freeze damage severity - 0-9 Rating|CO_350:0005001"

# Per-accession: GRM diagonal, and "informed" = has any nonzero off-diagonal
# (injected/mean-only accessions have all off-diagonals == 0 by construction).
d   <- diag(G)
off <- rowSums(abs(G)) - abs(d)
info <- tibble(genotype = rownames(G), diag = d, offsum = off,
               informed = off > 1e-8)

fz <- gebv |> filter(trait == !!trait) |> select(genotype, GEBV, phenotyped)
df <- info |> left_join(fz, by = "genotype")

cat("G:", nrow(G), "accessions |", sum(info$informed), "informed (nonzero off-diag),",
    sum(!info$informed), "injected/mean-only\n")
cat("Freeze GEBV range:", round(range(df$GEBV, na.rm = TRUE), 2), "\n\n")

cat("=== GEBV by informed status ===\n")
df |> group_by(informed) |>
  summarise(n = n(), gebv_min = min(GEBV), gebv_max = max(GEBV),
            gebv_sd = sd(GEBV), abs_gebv_mean = mean(abs(GEBV)), .groups = "drop") |>
  as.data.frame() |> print(digits = 3)

cat("\n=== GEBV by phenotyped status ===\n")
df |> group_by(phenotyped) |>
  summarise(n = n(), gebv_min = min(GEBV), gebv_max = max(GEBV),
            abs_gebv_mean = mean(abs(GEBV)), .groups = "drop") |>
  as.data.frame() |> print(digits = 3)

cat("\n=== Top 12 |GEBV|: informed? phenotyped? diagonal? ===\n")
df |> arrange(desc(abs(GEBV))) |> head(12) |>
  mutate(GEBV = round(GEBV, 2), diag = round(diag, 3)) |>
  select(genotype, GEBV, informed, phenotyped, diag) |> as.data.frame() |> print()

cat("\n=== Among INFORMED accessions: does a low diagonal go with extreme GEBV? ===\n")
inf <- df |> filter(informed)
cat("cor(diag, GEBV)   =", round(cor(inf$diag, inf$GEBV), 3), "\n")
cat("cor(diag, |GEBV|) =", round(cor(inf$diag, abs(inf$GEBV)), 3), "\n")
cat("\nLowest-diagonal informed accessions and their GEBV:\n")
inf |> arrange(diag) |> head(10) |>
  mutate(GEBV = round(GEBV, 2), diag = round(diag, 3)) |>
  select(genotype, diag, GEBV, phenotyped) |> as.data.frame() |> print()
