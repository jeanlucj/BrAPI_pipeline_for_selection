# Why do so few training accessions end up genomically connected in G?
# Classify each training accession by whether it (a) is genotyped on a USED
# protocol (its name appears in that protocol's cached VCF), (b) sits in a pedigree
# group, and (c) ended up "informed" in G (nonzero off-diagonal) vs injected.
# Offline: cached VCFs + pedigree files only. Also reports covering protocols that
# were ranked but NOT used (e.g. VCF download skipped), since accessions covered
# only by those are silently disconnected.

library(tidyverse)
here::i_am("code/coverage_diagnostic.R")
source(here::here("code", "04_find_genotyping.R"))   # .read_vcf_samples, pedigree helpers
source(here::here("code", "02_find_trials.R"))
source(here::here("code", "03_get_phenotypes.R"))

trials <- read_rds(cache_path("ny_trials.rds"))
pheno  <- read_rds(cache_path("phenotypes.rds"))
geno   <- read_rds(cache_path("genotypes.rds"))
sets   <- split_by_role(pheno, trials)

train <- sets$train_acc |> distinct(germplasmDbId, germplasmName)
train_names <- train$germplasmName
G <- geno$G

# (a) genotyped on a USED protocol = name appears in that protocol's cached VCF.
vdir <- cache_path("vcf_cache")
vcf_names <- unique(unlist(lapply(geno$protocol_ids, function(pid) {
  fs <- list.files(vdir, pattern = sprintf("^proto%s_proj.*\\.vcf$", pid), full.names = TRUE)
  unlist(lapply(fs, .read_vcf_samples))
})))

# (b) pedigree group membership (by name). Widen the cache dbid->name map with our
# own training names so our accessions resolve even if absent from the cache.
ped_d2n <- .germplasm_name_map(PEDIGREE_DIR)
ped_d2n[train$germplasmDbId] <- train$germplasmName
ped_members <- .pedigree_group_members(PEDIGREE_DIR)
ped_names <- unique(unlist(lapply(ped_members, function(ids) {
  nm <- ped_d2n[ids]; nm[!is.na(nm)]
})))

# (c) informed in G = nonzero off-diagonal.
informed_names <- rownames(G)[(rowSums(abs(G)) - abs(diag(G))) > 1e-8]

cls <- train |> transmute(
  germplasmName,
  has_markers  = germplasmName %in% vcf_names,
  in_pedigree  = germplasmName %in% ped_names,
  in_G         = germplasmName %in% rownames(G),
  informed     = germplasmName %in% informed_names)

cat(sprintf("Training accessions: %d (distinct names: %d)\n",
            nrow(train), n_distinct(train_names)))
cat(sprintf("In G: %d | informed (nonzero off-diag): %d | injected: %d\n",
            sum(cls$in_G), sum(cls$informed), sum(cls$in_G & !cls$informed)))

cat("\n=== Coverage of training accessions ===\n")
# (compute the union/neither BEFORE summarise so the reducing sums don't shadow
# the has_markers / in_pedigree columns mid-expression)
cls2 <- cls |> mutate(markers_or_ped = has_markers | in_pedigree,
                      neither = !has_markers & !in_pedigree)
cls2 |> summarise(
  n              = n(),
  has_markers    = sum(has_markers),
  in_pedigree    = sum(in_pedigree),
  markers_or_ped = sum(markers_or_ped),
  neither        = sum(neither)) |>
  as.data.frame() |> print()

cat("\n=== Cross-tab: connected-source vs informed-in-G ===\n")
cls |> count(has_markers, in_pedigree, informed) |> as.data.frame() |> print()

cat("\n=== LEAK: genotyped and/or pedigreed but NOT informed in G ===\n")
leak <- cls |> filter((has_markers | in_pedigree) & !informed)
cat(nrow(leak), "accessions have a marker/pedigree source but are disconnected in G:\n")
print(as.data.frame(head(leak, 20)))

cat("\n=== Covering protocols RANKED but NOT used (e.g. VCF skipped) ===\n")
used <- as.character(geno$protocol_ids)
geno$protocols |>
  mutate(used = dbId %in% used) |>
  filter(n_covered > 0) |>
  select(rank, dbId, name, n_covered, n_total, used) |>
  as.data.frame() |> print()

write_rds(cls, cache_path("coverage_diag.rds"))
