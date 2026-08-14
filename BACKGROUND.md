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

- **Authentication is required.** T3/Oat no longer serves this data anonymously, and
  an unauthenticated call returns an *empty* result rather than a 401, so a missing
  login looks exactly like a missing dataset. `connect_t3()` therefore logs in on
  every run (`T3_USERNAME`/`T3_PASSWORD` from a gitignored `.Renviron`) and fails
  loudly when credentials are absent or rejected — the latter also needs an explicit
  check, since `conn$login()` merely warns on a wrong password.
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

Because the weight is 1/SE², a trial that fits *perfectly* (residual variance → 0)
would hand Stage 2 a near-infinite weight and dominate it — exactly what a trait with
no within-genotype variation does (e.g. winter survival scored 100 on every plot, or
duplicated records). Three guards prevent this: (1) a trial × trait whose response is
**constant** is dropped (no estimable genotypic signal — it carries no information,
not infinite precision); (2) `replicate` and `block` are not both entered as random
effects when they describe the **same partition** (their labels coincide), which would
alias the two variance components; (3) every SE is **floored** at `SE_FLOOR_FRAC` of
the trial × trait response SD, so a near-perfect fit cannot blow up 1/SE². Well-
estimated BLUEs sit far above the floor and are untouched.

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
- Accessions appearing in **≥ 2 partials act as bridges** — counting each
  genotyping platform *and* each pedigree group (below). The EM uses them to
  *impute* relationships between accessions that were never co-genotyped, treating
  the unobserved blocks as missing data in a multivariate Wishart model. Because a
  pedigree group counts as a partial, an accession seen on only **one** platform is
  still retained when it also sits in a pedigree group — it bridges that platform's
  marker GRM to the pedigree matrix, even if it was never phenotyped.
- The result is one combined GRM spanning the union of accessions. Pairs with no
  shared platform and no bridge stay near the prior — a faithful reflection of
  genuine lack of connectivity rather than a fabricated relationship.

### Degrees of freedom = how much each partial is trusted
Each partial enters the EM with a **degrees-of-freedom** value νᵢ that acts as its
relative weight in the M-step (Ψ = Σνᵢ·E[…] / Σνᵢ). In the Wishart model νᵢ is the
number of independent samples behind the matrix (Akdemir et al. 2020, *Front. Plant
Sci.* 11:947; 2023, *Axioms* 12:161): a single common ν leaves the combined estimate
unchanged (it only scales the standard errors), but **differing** νᵢ reweight the
estimate. A marker GRM is estimated over *markers*, not individuals, so its df should
reflect the number of markers — but markers are in strong LD, so the raw count badly
overstates the independent information. We therefore set each marker GRM's df from its
**effective number of independent samples** (Galwey 2009's `(Σ√λ)²/Σλ` over the GRM
eigenvalues, scale-invariant and bounded by the rank).

That measure is used only for the GRMs' *relative* ordering: the values are then
**re-centered on `GRM_DF_MEAN` and their spread capped at `GRM_DF_STDEV`** (the lesser
of `GRM_DF_STDEV` and the observed SD). Pedigree partials keep a fixed `PEDIGREE_DF`,
so `GRM_DF_MEAN` vs `PEDIGREE_DF` is the explicit knob for how much more marker
matrices are trusted than pedigree (default 60 vs 30 ≈ 2×). This deliberately bounds
an artifact of the eigenvalue measure — more *diverse* panels have flatter spectra and
thus a higher effective N — so panel diversity nudges the relative weighting but cannot
dominate it.

`GENO_PROTOCOL_ID = NULL` triggers this multi-platform combine across all covering
protocols; a single id skips it (one GRM, faster, with a usable marker matrix).

### Marker thinning to a target density
GBS files can carry >1M markers and exceed memory. `thin_vcf()` streams a VCF and
keeps every n-th marker genome-wide, with `n = floor(n_markers / TARGET_DENSITY)`
(default target 10000). Relationship estimates are stable well below the full marker
count, and per-protocol standardization absorbs differing marker counts across
platforms, so thinning costs little. Thinning only engages when a file has ≥ 2×
the target.

### Marker keys and QC ordering
A marker is keyed by **canonical CHROM_POS** (a leading `chr` stripped), built the same
way for every VCF, because the `ID` column is inconsistent across T3 projects (some
carry SNP names, others `.`); keying on `ID` made the *same* physical marker look
distinct between projects, so merging a multi-project protocol exploded into ~2× markers
that were ~50% missing and gutted accessions in QC. With CHROM_POS keys the shared array
markers align across a protocol's projects.

Even so, project files can still use genuinely **different panels** (e.g. an array vs a
GBS run), so merging leaves real missingness. QC therefore drops high-missing
**markers before accessions** — removing the minority-panel markers first makes the
accessions look complete again instead of discarding all of them.

### Estimating each partial GRM on a full panel
A VanRaden GRM centers and scales markers by **allele frequencies estimated from the
accessions in the matrix**, and the diagonal is each accession's self-relationship
*relative to that panel's* frequencies. If we restrict to only our handful of relevant
accessions (targets + bridges) before computing frequencies, imputing, and building the
GRM, those quantities are estimated from a few individuals and `std_grm` then pins their
mean diagonal to 1 — so a diagonal becomes "distance from these few panel-mates" rather
than a stable property. In practice this produced wildly low diagonals for some
accessions, and the apparent ranking even **inverted** versus a broad-panel estimate (an
accession that looked under-related on a tiny panel was actually among the *most*
divergent on the full panel). The same instability corrupts the off-diagonals that feed
the EM combine. So each protocol's imputation + GRM are estimated on a **panel** of all
its relevant accessions plus non-relevant fillers up to `max(GRM_PANEL_MIN, n_relevant)`,
then the GRM is subset back to the relevant accessions. (Diagnostic: `code/impute_diagnostic*.R`.)

### Imputing the residual missing calls
Missing genotype calls are filled before the GRM is built (`GRM_IMPUTE`). The default
`.impute_glmnet` regresses each incomplete marker on its most-correlated markers via
cross-validated elastic net, borrowing information across markers — which is only
worthwhile on a reasonably large panel (hence the panel above). It is **robust**: a
marker whose CV fails (commonly a near-monomorphic marker whose CV training fold is all
one genotype) falls back to its column mean rather than aborting, and CV folds are fixed
so the imputed matrix is reproducible. The simpler `.mean_impute` fills each missing call
with the marker mean; mean-filling pulls heavily-missing accessions toward the population
mean, but at this panel size missingness turned out **not** to be the driver of the low
diagonals — the panel-size effect above was.

## Pedigree stitch

A **pedigree (numerator) relationship matrix** can stitch together marker GRMs that
remain disconnected — pedigree links span accessions that no shared platform covers.
The pipeline reads **precomputed** pedigree relationship matrices from the sibling
`BrAPI_pedigree_relmat` project (`PEDIGREE_DIR`) and adds each group overlapping our
accessions **or bridges** to the EM combine as an extra partial covariance, weighted
by `PEDIGREE_DF`.

**Two keying spaces meet here.** Pedigree group CSVs are keyed by `germplasmDbId`,
but marker GRMs (and `keep_samples`) are keyed by `germplasmName`. To decide whether
a *non-phenotyped* accession bridges a platform to a pedigree group, the pipeline
needs names for accessions beyond our own. It gets them from the companion
`germplasm_cache_<id>.rds` (a raw BrAPI germplasm dump the pedigree project writes
next to `PEDIGREE_DIR`); if that cache is absent it degrades gracefully to
phenotyped-only pedigree treatment. Group **membership** is read from the cheap
`<id>_pedigree_groups.rds` so bridge detection never densifies the group, and the
group matrix itself is reconstructed over only the kept (our + bridge) accessions —
so the large Oat component (~27,800 accessions, a 1 GB triplet CSV) never becomes a
full dense matrix.

## Selection

Stage-2 GEBVs are combined into a weighted **selection index**,
`Index = Σ_trait TRAIT_WEIGHTS × GEBV`. The GEBVs are **not** standardized: a GEBV is
already shrunk in proportion to its uncertainty, so standardizing per trait would
inflate poorly-estimated values and distort the index. The weights therefore carry both
direction (sign) and the relative scale/importance of each trait.

Results are written to `output/breeders_output.csv`, an intentionally spreadsheet-ish
file with two side-by-side blocks (separated by a blank column): block 1 lists every
predicted accession (sorted by index) with a training-membership flag, its per-trait
GEBVs, and its index; block 2 lists, per trait, the index weight and the
cross-validated prediction accuracy (`cv_accuracy`).

## Targeting "New York"

The environment is targeted at the **trial-selection** stage: training trials are
chosen by proximity to a New York center point (or by an explicit `TRAINING_TRIALS`
list). Because Stage 2 removes environment main effects and predicts a genomic value
across those NY-region environments, the resulting GEBVs reflect performance in that
target environment rather than a global mean.

## Who gets predicted

Prediction is scoped explicitly so the relationship matrix stays small and relevant.
The target set is the union of the training-trial accessions (always predicted),
`TEST_TRIALS` accessions, and `TEST_ACCESSIONS`. The pipeline builds the full
EM-combined relationship matrix (using bridges and pedigree to estimate cross-platform
relatedness) and then **subsets it to the targets** — bridge/extra accessions inform
the combine but are not carried into prediction. Test accessions are predicted only if
they land in the matrix (via genotype or pedigree); training accessions are always
predicted, and one with neither genotype nor pedigree is **injected** with a
prior-only row (diagonal = mean diagonal, zero off-diagonals), so it still receives a
GEBV shrunk toward the mean. `TEST_TRIALS` phenotypes are held out of training (a
trial listed in both `TRAINING_TRIALS` and `TEST_TRIALS` is used for training).
