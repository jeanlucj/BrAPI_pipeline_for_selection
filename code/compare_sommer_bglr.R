# Cross-check the BGLR (Bayesian) Stage-2 GEBVs against a sommer (REML GBLUP) fit
# on the same relationship kernel, for Freeze damage. Offline: uses cached
# phenotypes / trials / genotypes; reconstructs Stage-1 BLUEs locally. Writes a
# scatterplot to output/sommer_vs_bglr_freeze.png.

library(tidyverse)
here::i_am("code/compare_sommer_bglr.R")
for (f in c("config.R", "02_find_trials.R", "03_get_phenotypes.R",
            "grm_utils.R", "05_stage1_blues.R",
            "06_stage2_genomic_prediction.R")) source(here::here("code", f))

trait <- "Freeze damage severity - 0-9 Rating|CO_350:0005001"

trials <- read_rds(cache_path("ny_trials.rds"))
pheno  <- read_rds(cache_path("phenotypes.rds"))
geno   <- read_rds(cache_path("genotypes.rds"))
sets   <- split_by_role(pheno, trials)
blues  <- stage1_blues(sets$train_pheno)

# BGLR GEBVs from the cached pipeline run; sommer GEBVs fit now.
# "informed" = accession has genomic/pedigree connections (nonzero off-diagonal in
# G); the rest were injected with a prior-only diagonal (off-diag 0).
G <- geno$G
informed <- tibble(genotype = rownames(G),
                   informed = (rowSums(abs(G)) - abs(diag(G))) > 1e-8)

bglr_g <- read_rds(cache_path("gebv.rds")) |>
  filter(trait == !!trait) |> select(genotype, GEBV_bglr = GEBV)
somm_g <- stage2_gblup(blues, geno, traits = trait, engine = "sommer",
                       refresh = TRUE) |>
  select(genotype, GEBV_sommer = GEBV)

cmp <- bglr_g |> inner_join(somm_g, by = "genotype") |>
  left_join(informed, by = "genotype")
r  <- cor(cmp$GEBV_bglr, cmp$GEBV_sommer)
rs <- cor(cmp$GEBV_bglr, cmp$GEBV_sommer, method = "spearman")
cat(sprintf("Freeze damage: %d accessions | Pearson r = %.3f | Spearman = %.3f\n",
            nrow(cmp), r, rs))
cat(sprintf("BGLR range  [%.1f, %.1f] (sd %.2f); sommer range [%.2f, %.2f] (sd %.3f)\n",
            min(cmp$GEBV_bglr), max(cmp$GEBV_bglr), sd(cmp$GEBV_bglr),
            min(cmp$GEBV_sommer), max(cmp$GEBV_sommer), sd(cmp$GEBV_sommer)))

cat("\nGEBV sd by informed status (markers/pedigree vs injected prior-only):\n")
cmp |> group_by(informed) |>
  summarise(n = n(), sd_bglr = sd(GEBV_bglr), sd_sommer = sd(GEBV_sommer),
            .groups = "drop") |> as.data.frame() |> print(digits = 3)

# Free scales: the ~20x scale gap makes a fixed 1:1 plot uninformative.
p <- ggplot(cmp, aes(GEBV_bglr, GEBV_sommer, color = informed)) +
  geom_vline(xintercept = 0, colour = "grey85") +
  geom_hline(yintercept = 0, colour = "grey85") +
  geom_point(alpha = 0.8) +
  labs(title = "Freeze damage GEBVs: sommer (REML) vs BGLR (Bayesian)",
       subtitle = sprintf("n = %d | Pearson r = %.3f, Spearman = %.3f | note differing axis scales",
                          nrow(cmp), r, rs),
       x = "BGLR GEBV", y = "sommer GEBV",
       color = "has markers/pedigree") +
  theme_bw()
ggsave(output_path("sommer_vs_bglr_freeze.png"), p, width = 6.5, height = 5.5, dpi = 120)
cat("\nWrote ", output_path("sommer_vs_bglr_freeze.png"), "\n", sep = "")
write_rds(cmp, cache_path("sommer_vs_bglr_freeze.rds"))
