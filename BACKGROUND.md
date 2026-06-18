# BACKGROUND

The theory behind the methods and the reasoning behind the design decisions. For
usage see [README.md](README.md); for structure see [DESIGN.md](DESIGN.md).

## The breeding question

We want to choose **parents for crossing** that will improve a trait (e.g. grain
yield) in the **New York** environment. The best available evidence is the union of
many public oat field trials plus marker data on their accessions. Genomic
prediction lets us rank not only phenotyped lines but also their **genotyped but
unphenotyped relatives**, which is exactly the candidate pool a breeder draws
crossing parents from.

## Data source: BrAPI / T3/Oat

T3/Oat is a Breedbase instance exposing the **Breeding API (BrAPI)**. We use the
`TriticeaeToolbox/BrAPI.R` wrapper. Decisions and quirks discovered against the live
server:

- **Anonymous access** is sufficient for the phenotype and (archived) genotype data
  used here; login is wired as a fallback only.
- **`/studies` ignores a server-side `locationDbId` filter**, so all studies are
  pulled once and filtered client-side.
- **Genotyping BrAPI endpoints are not implemented** on T3/Oat (`/allelematrix`
  returns an error page), so marker data must come from Breedbase VCF downloads.
- On-the-fly VCF generation (`conn$vcf`) is **slow and stalls** for large protocols,
  so we prefer **pre-generated archived VCFs** (`conn$vcf_archived`).
- Archived GBS matrices can be **multiple GB**; the largest (diversity panels) may
  fail mid-transfer. Hence marker **thinning** and **skip-on-failure**.

## Two-stage genomic prediction

We use the standard two-stage approach (rather than a single huge model) because it
is robust and modular:

**Stage 1 — per-trial BLUEs (`lme4`).** Within each trial × trait we fit genotype as
a *fixed* effect and the field design (replicate, block) as *random* effects, giving
Best Linear Unbiased Estimates adjusted for local spatial/design structure. This
also handles augmented designs (unreplicated entries adjusted via replicated check
blocks). Each BLUE carries a weight = 1/SE², propagating its precision to Stage 2;
trials without replication fall back to a fixed-only model or plot means.

**Stage 2 — genomic prediction (`BGLR`).** Environment main effects are removed
(environment as a fixed effect, via weighted centering), BLUEs are weight-averaged to
one value per genotype, and a genomic model is fit over **all** genotyped accessions
— phenotyped lines train the model, unphenotyped candidates receive predicted GEBVs.

### Why BGLR, and RKHS by default
`BGLR` is a well-established Bayesian whole-genome-regression package. The default
model is **RKHS** (reproducing-kernel Hilbert space) using a genomic **relationship
matrix** as the kernel. RKHS is the natural choice here because the multi-platform
combine (below) produces a *relationship matrix*, not a single aligned marker
matrix — there is no common marker set to feed a marker-effect model. When a single
protocol is used a raw marker matrix is available, so marker-effect models
(BRR/BayesB) remain selectable; otherwise the pipeline transparently uses RKHS on the
combined GRM.

## Combining multiple genotyping platforms (EM covariance combine)

Different genotyping protocols (GBS runs, SNP arrays) genotype overlapping but
distinct sets of accessions on **disjoint marker sets**, so their genomic
relationship matrices cannot simply be averaged. We combine them with a **Wishart-EM
covariance combiner** (Akdemir), reused from the sibling `T3Predictathon2026`
project:

- Each protocol contributes a standardized GRM (mean diagonal = 1) as a *partial
  covariance* over the accessions it genotyped.
- Accessions genotyped on **≥ 2 platforms act as bridges**: the EM uses them to
  *impute* relationships between accessions that were never co-genotyped, treating
  the unobserved blocks as missing data in a multivariate Wishart model.
- The result is one combined GRM spanning the union of accessions. Pairs with no
  shared platform and no bridge stay near the prior — a faithful reflection of
  genuine lack of connectivity rather than a fabricated relationship.

`GENO_PROTOCOL_ID = NULL` triggers this multi-platform combine across all covering
protocols; a single id skips it (one GRM, faster, with a usable marker matrix).

### Marker thinning to a target density
GBS files can carry >1M markers and exceed memory. `thin_vcf()` streams a VCF and
keeps every n-th marker genome-wide, with `n = floor(n_markers / TARGET_DENSITY)`
(default target 10000). Relationship estimates are stable well below the full marker
count, and per-protocol standardization absorbs differing marker counts across
platforms, so thinning costs little. Thinning only engages when a file has ≥ 2×
the target.

### Marker QC ordering
Within a protocol, project files can use **different genome annotations** (different
marker names), so merging leaves every accession ~50% missing. QC therefore drops
high-missing **markers before accessions** — removing the minority-annotation markers
first makes the accessions look complete again instead of discarding all of them.

## Pedigree stitch

A **pedigree (numerator) relationship matrix** can stitch together marker GRMs that
remain disconnected — pedigree links span accessions that no shared platform covers.
The pipeline reads **precomputed** pedigree relationship matrices from the sibling
`BrAPI_pedigree_relmat` project (`PEDIGREE_DIR`) and adds any matrix overlapping our
accessions to the EM combine as an extra partial covariance, weighted by
`PEDIGREE_DF`.

**Current limitation for T3/Oat:** that project only precomputed matrices for its
*small* pedigree groups; the large connected component (~27,800 accessions, where
the NY breeding lines mostly live) exceeded its size cap and was skipped. So the
pedigree stitch rarely engages for Oat today. The mechanism is correct and fully
exercised where matrices exist (e.g. wheat); to enable it for Oat, rerun the pedigree
project with a higher `max_relmat_size`.

## Selection

Stage-2 GEBVs are standardized per trait and combined into a weighted **selection
index** (default equal weights, `+1` = higher-is-better). Candidates are ranked, and
top picks that are highly related (via the combined GRM) to a better-ranked pick are
flagged `redundant_with` — so a breeder can avoid choosing near-identical parents and
keep diversity in the crossing block.

## Targeting "New York"

The environment is targeted at the **trial-selection** stage: trials are chosen by
proximity to a New York center point (or by an explicit `STUDY_NAMES` list). Because
Stage 2 removes environment main effects and predicts a genomic value across those
NY-region environments, the resulting GEBVs reflect performance in that target
environment rather than a global mean.
