# Evaluation guide — stepping through the selection pipeline

This is the document for **evaluating and validating** the pipeline from the R console
— for me now and in six months. `README.md` is for running it; `DESIGN.md` is the
static architecture; `BACKGROUND.md` is why it is built this way; `code/README.md` is
the per-function reference. This file is the **runbook**: it lists the modules, gives a
fast→slow order to check them in, and walks each one step by step — one console line at
a time — saying what each step should return, what to eyeball, and the subtle failure
signatures to hunt for.

The pipeline is long (seven steps, two keying spaces, four data sources) and most of it
runs unattended, so the goal is that **working through this document is faster and more
thorough than re-reading all the code**, and that it maximizes the chance of catching
the bugs that make the pipeline *look* like it works.

The objective, to keep in view the whole time: produce a **ranked list of oat
accessions to use as parents for crossing in New York**, whose GEBVs are trustworthy
enough to act on. The failure mode being hunted is not a crash — it is the run that
finishes cleanly while silently **dropping** accessions (a test accession absent from
`G` is dropped without a word), **mis-scaling** them (weighted RKHS inflates breeding
values ~λ-fold), **mis-keying** them (coverage is counted by `germplasmDbId`, markers
match by `germplasmName`), or estimating a relationship matrix off the wrong panel.
Every one of those produces a plausible-looking parent list.

------------------------------------------------------------------------

## 1. The modules, in one screen

Twelve files in `code/`, one driver, one report. The **group** column is the name you
pass to `arm_evaluation()` (§2) to `debug()` that module.

| File (`code/…`) | Owns | Off/online | Group |
|----|----|----|----|
| `config.R` | every parameter; reads the long selection lists from `data/config/*.txt` | offline | `config` |
| `progress.R` | console status lines, bars, per-step timings | offline | `progress` |
| `grm_utils.R` | VanRaden GRM, standardization, Galwey effective *n*, marker imputation | offline | `grm` |
| `em_covariance_combiner.R` | the Wishart-EM combiner (verbatim copy from `T3Predictathon2026`) | offline | `combine` |
| `04_find_genotyping.R` (VCF half) | VCF sample/marker reading, thinning, dosage merge, marker QC | offline | `markers` |
| `05_stage1_blues.R` | per-trial BLUEs (lme4) + `1/SE²` weights | offline | `stage1` |
| `06_stage2_genomic_prediction.R` | kernel GBLUP → GEBVs (BGLR or sommer) | offline\* | `stage2` |
| `06_…` (`cv_accuracy`) | k-fold cross-validated accuracy | offline\* | `cv` |
| `07_select.R` | selection index → `breeders_output.csv` | offline | `select` |
| `01_connect.R` | BrAPI connection + login | **online** | `connect` |
| `02_find_trials.R` | locations, studies, training/test roles | **online** | `trials` |
| `03_get_phenotypes.R` | observations + field design; `split_by_role()` | **online** | `phenotypes` |
| `synonyms.R` + coverage helpers | alias→primary map; protocol/project coverage tables | **online** | `coverage` |
| `04_…` (pedigree half) | the sibling repo's group matrices + its output contract | **online** | `pedigree` |
| `04_…` (orchestration) | download → per-protocol GRM → EM combine → subset/inject | **online** | `genotyping` |

\* Stages 2 and CV are offline *once* `blues` and `geno` exist — they fit models, they
do not call the server. `arm_evaluation("step4")` is an aggregate over `markers`,
`grm`, `pedigree`, `combine` and `genotyping` for when you want to watch one protocol
go all the way through.

**The one structural thing to hold in your head:** there are **two keying spaces**.
Genotyping *coverage* is resolved server-side by `germplasmDbId`; everything
marker-, GRM- and BLUE-side is keyed by `germplasmName`. Step 04 holds the only bridge
between them (`dbid_to_name`, widened by the pedigree germplasm cache). Most of the
silent-loss bugs in §6 live at that seam.

------------------------------------------------------------------------

## 2. How to evaluate: the plan

Three principles, in priority order:

1. **Offline before online.** Prove the math on synthetic data (and the test suite, §7)
   before spending a single download. Levels L1–L8 need no server and no `data/`.
2. **Cached before cold.** Steps 2/3/4/6 cache to `data/*.rds`; the live levels run
   against what you already have and are near-instant. Exactly **one** level (L14)
   deliberately goes cold, and it is scoped to the small Oat 3K protocol.
3. **One group armed at a time.** `arm_evaluation("stage2")`, run it, step through
   *only* Stage 2, `disarm_evaluation()`. Arming everything makes the debugger unusable.

### The three tools

| Tool | Use it to | When |
|----|----|----|
| `arm_evaluation(group)` / `disarm_evaluation()` (`code/evaluation.R`) | `debug()` one module and watch it fill in line by line | stepping through logic |
| `peek(x)` (`code/evaluation.R`) | print a one-line health summary of the object flowing between steps — shape, NA counts, diagonal, weight range, name overlap — and pass it through unchanged | inspecting an intermediate |
| the `code/*diagnostic*.R` / `*_compare.R` scripts (§9) | re-derive a result **independently** of the pipeline | confirming a real result is real |

`eval_groups()` prints the whole menu. The debugger keys, once: **`n`** next line,
**`s`** step into a call, **`c`** finish this function, **`Q`** quit.
**Always `disarm_evaluation()` before an unattended run or `wflow_build()`** —
otherwise every call blocks on the debugger prompt. `armed_functions()` tells you
what is still armed.

Progress bars and the debugger fight over the console: disarm before you watch a bar,
arm when you want to step.

------------------------------------------------------------------------

## 3. Bootstrap (paste once per session)

``` r
library(tidyverse)
purrr::walk(c("config.R", "01_connect.R", "02_find_trials.R", "03_get_phenotypes.R",
              "04_find_genotyping.R", "05_stage1_blues.R",
              "06_stage2_genomic_prediction.R", "07_select.R"),
            ~ source(here::here("code", .x)))
source(here::here("code", "evaluation.R"))   # NOT sourced by the pipeline, on purpose
eval_groups()                                 # confirm the tooling loaded
options(brapi.progress = TRUE)                # narrate steps even if not interactive
```

> Source the step files by name, **not** with a `list.files("code")` glob:
> `run_pipeline.R` runs the entire pipeline on source, and the diagnostics in §9 run
> their own experiments.

For the **offline** levels (L1–L8) also load the synthetic-data builders the test suite
already provides — no new fixtures, no network:

``` r
source(here::here("tests", "testthat", "helper-setup.R"))   # make_dosage, write_test_vcf, fake_conn
```

> **Caveat** `helper-setup.R` redirects `cache_path()` and
> `output_path()` to a tempdir so tests never touch `data/`. That is exactly what you
> want offline — but it means an online level run afterwards would read and write the
> *wrong* cache. Before the online levels, restore them with
> `source(here::here("code", "config.R"))` (or start a fresh session).

For the **online** levels (L9 on), open one connection and reuse it. T3/Oat requires a
login, so `.Renviron` must hold `T3_USERNAME` / `T3_PASSWORD` (copy `.Renviron.example`):

``` r
conn <- connect_t3()        # reads .Renviron itself; no restart needed
nchar(conn$auth_token)      # -> a non-zero token length
```

> **Auth symptoms.** No credentials → an error of class `t3_missing_credentials`.
> Credentials the server rejects → `t3_bad_credentials`, preceded by BrAPI's own
> `Incorrect Password` warning. Anything else that comes back *empty* is a data
> question, not an auth question — an unauthenticated T3/Oat call returns an empty
> result rather than a 401, which is precisely why the login is mandatory.

------------------------------------------------------------------------

## 4. Level-by-level walkthrough

Each level: **arm** the group, run the numbered console lines (each assigns a
Global-Environment variable), and check the four annotations — **→ returns** (the shape
you should see), **🔍 eyeball** (what healthy looks like), **🚩 red flag** (the subtle
failure to hunt), and the italic *Assumption:* the code is making at that step (test
**it**, not just the output). `disarm_evaluation()` when you finish a level.

The numbers quoted below are from the synthetic fixtures at the seeds given, so they
are reproducible: if yours differ, either the fixture or the code moved.

### L1 — `config` (offline, instant)

``` r
arm_evaluation("config")
tr <- config_traits("Trait_Sel2026.txt")       # step: watch the column split
tt <- config_lines("Trial_Sel2026_Intersect20NoYld.txt")
ac <- config_lines("Acc_Sel2026.txt")
disarm_evaluation()
str(tr); length(tt); length(ac)
try(config_lines("no_such_file.txt"))          # MUST error
```

- **→ returns** `tr` is `list(names=, weights=, short=)` — with the current one-column
  trait file, `weights` and `short` are `NULL`; `tt` is 9 study names; `ac` is 816
  germplasm names.
- **🔍 eyeball** `TRAIT_WEIGHTS`/`TRAIT_SHORT_NAMES`, when present, are **named vectors
  keyed by `TRAIT_NAMES`** — `identical(names(tr$weights), tr$names)`. Names with a pipe
  (`APPLER|CIAV2680`, `Grain yield - g/m2|CO_350:0000260`) come through **verbatim**.
- **🚩 red flag** a missing file returning `NULL` instead of erroring. `TRAINING_TRIALS
  <- NULL` means "fall back to the geographic radius search" — a typo'd filename would
  silently run a completely different analysis. Also: a half-filled weight column must
  error, not silently zero part of the selection index.
- *Assumption:* trait names are matched **exactly**, not as substrings, against
  `observationVariableName`. A trait in the file that does not exist on the server is
  not an error anywhere — it simply yields no observations (check at L11).

### L2 — `progress` (offline, instant)

``` r
options(brapi.progress = FALSE)
x <- time_it("nothing", 6 * 7)        # silent, returns 42
options(brapi.progress = TRUE)
timings_reset(); h <- step_start("demo"); Sys.sleep(0.2); step_done(h, "ok"); print_timings()
```

- **→ returns** `x == 42` either way; `pipeline_timings()` gains one row per `step_done()`.
- **🔍 eyeball** with reporting off, nothing prints and nothing changes — the helpers are
  dropped inline all through the pipeline, so they must be transparent.
- **🚩 red flag** a bar that redraws every iteration of a long loop (it costs real time);
  a bar that keeps drawing after its function returned (the `.envir` contract).
- *Assumption:* `SHOW_PROGRESS <- interactive()` — so a batch `Rscript` run, the test
  suite and `wflow_build()` are silent by design. If you expected narration and got
  none, that is why.

### L3 — `grm` (offline; the relationship-matrix crux)

``` r
arm_evaluation("grm")
D  <- make_dosage(n_acc = 30, n_mar = 200, seed = 1)   # 0/1/2, accessions x markers
G  <- .Gmatrix(D)                                      # raw VanRaden
Gs <- std_grm(G)                                       # mean diagonal -> 1
.effective_n(Gs)
M  <- D; M[sample(length(M), 200)] <- NA
i1 <- .impute_glmnet(M); i2 <- .impute_glmnet(M)
disarm_evaluation()
peek(Gs); identical(i1, i2); anyNA(i1)
```

- **→ returns** `G`, `Gs` are 30×30 symmetric; `.effective_n(Gs)` ≈ **27.9** of 30;
  `i1` has no `NA` left and `identical(i1, i2)` is `TRUE`.
- **🔍 eyeball** the **raw** diagonal mean is **1.308** here and `std_grm()` takes it to
  exactly 1. On real inbred oat data the raw VanRaden diagonal is `1 + F` and should sit
  comfortably **above 1** — a raw diagonal near 1.0 on real data is suspicious, not
  reassuring. To see that directly, build a GRM from the *inbred* population L6 uses:
  `mean(diag(.Gmatrix(simulate_trials(seed = 1)$D)))` is **exactly 2.0**, i.e. `F = 1`
  for fully homozygous lines.
- **🚩 red flag** `.impute_glmnet` not reproducible across calls (it fixes its CV folds
  deterministically on purpose); a diagonal mean of exactly 1 where you expected a *raw*
  GRM (something already standardized it); dosages outside `{0,1,2}` — `.Gmatrix()`'s
  `p <- colMeans(M)/2` assumes that coding and produces a meaningless `p` otherwise.
- *Assumption:* `.impute_glmnet` retries each marker independently and falls back to the
  column mean on failure, reporting the fallback count. A large fallback count means the
  panel is too small or too monomorphic for imputation to be earning its keep — compare
  against `GRM_IMPUTE = "mean"` (§9, `impute_diagnostic.R`).

### L4 — `combine` (offline; the EM weighting)

``` r
arm_evaluation("combine")
A  <- std_grm(.Gmatrix(make_dosage(20, 150, seed = 2)))
B  <- std_grm(.Gmatrix(make_dosage(20, 150, seed = 3)))
nm <- union(colnames(A), colnames(B))
vi <- list(match(colnames(A), nm), match(colnames(B), nm))
res <- EMCovarianceCombiner(partial_covs = list(A, B), var_indices = vi,
                            degrees_freedom = c(60, 30))
disarm_evaluation()
dim(res$psi); length(nm)                       # -> 21 21   and   20
.center_dfs(c(20, 25, 30), mean_df = 60, sd_df = 15)   # -> 55 60 65
```

- **→ returns** `res$psi` is **one row/col larger** than `length(nm)` — the leading
  **phantom variable** the combiner adds. `.combine_to_G()` drops it with `psi[-1, -1]`.
- **🔍 eyeball** after dropping the phantom row and symmetrizing, `all.equal(G, t(G))`.
  `.center_dfs` preserves the *ordering* of the effective sample sizes while re-centering
  them on `GRM_DF_MEAN` and capping their spread at `GRM_DF_STDEV`.
- **🚩 red flag** a combined matrix one row/col too big (the phantom variable survived —
  every accession would be shifted by one); the EM bar running to its full 100 iterations
  without the `Converged at iteration N` message (it converges in single digits on
  well-conditioned partials).
- *Assumption:* the degrees of freedom are the **weights**. `GRM_DF_MEAN` (60) vs
  `PEDIGREE_DF` (30) is the statement "trust markers about twice as much as pedigree".
  That ratio is a *choice*, not a measurement — `df_grid_diagnostic.R` (§9) shows what it
  does to `G`.
- *Assumption:* partials are stitched only through accessions they **share**. Two
  genuinely disjoint partials produce zero cross-blocks, which is correct and not a bug —
  but it means those accessions cannot inform each other's GEBVs at all.

### L5 — `markers` (offline; the keying crux)

``` r
arm_evaluation("markers")
D   <- make_dosage(12, 40, seed = 4)
v   <- write_test_vcf(tempfile(fileext = ".vcf"), D)   # writes every ID as "."
.read_vcf_samples(v); .count_markers(v)
dos <- .vcf_to_dosage(v)
M   <- dos; M[, 1:5] <- NA; M[1, 6:20] <- NA
q   <- .qc_markers(M, max_missing = 0.5, min_maf = 0)
disarm_evaluation()
head(colnames(dos), 3)                  # -> "1_100" "1_200" "1_300"
peek(dos, accessions = rownames(D))     # overlap 12 of 12
peek(dos, accessions = c("NOPE1"))      # overlap 0 -> the mismatch alarm
```

- **→ returns** `dos` is 12 accessions × 40 markers of `{0,1,2}`; markers are named
  **`CHROM_POS`** (`1_100`, `1_200`, …), *not* from the VCF `ID` column — the fixture
  writes every `ID` as `"."` precisely to exercise that. QC here keeps 35 of 40 markers
  and all 12 accessions.
- **🔍 eyeball** `peek(dos, accessions = )` reports the rowname overlap. Twelve of twelve
  is health; **zero overlap on a non-empty matrix is the synonym/name-mismatch
  signature** — the VCF carries our accessions under a preliminary name later demoted to
  a synonym (see `USE_SYNONYMS`, L12).
- **🚩 red flag** markers keyed by the VCF `ID` column: T3 projects key inconsistently
  (SNP names in one, `.` in another), which makes *identical* markers look distinct and
  silently breaks the multi-project merge — the merged matrix then looks ~50% missing.
  Also: QC dropping **accessions before markers**. The order is deliberate — a protocol's
  projects can use different genome annotations, leaving every accession ~50% missing
  until the minority-annotation markers go first.
- *Assumption:* large VCFs are thinned genome-wide to `TARGET_DENSITY` (10 000) by
  keeping every *n*-th variant line. That is uniform in *position*, not in linkage — fine
  for a GRM, wrong if you ever want specific markers.

### L6 — `stage1` (offline; where the weights are born, and where truth is defined)

The fixture from here on is a **simulated population with a known answer**, because
L7's whole job is to decide whether Stage 2 is right, and you cannot do that without
knowing what right is. `simulate_trials()` (`tests/testthat/helper-setup.R`) builds:

- **a related population** — `n_fam` biparental crosses among `n_founders` founders,
  doubled-haploid progeny, so dosages are 0/2 (oat lines are essentially homozygous) and
  lines from one cross share about half their genome;
- **a genetic value per line**, `g = Wβ` with **every marker causal** (`β ~ N(0,1)`, `W`
  the centred dosage), rescaled so `var(g) = h²`;
- **phenotypes** `y = μ + trial + rep + block + g + e`, `var(e) = 1 − h²`, in a
  randomized incomplete-block design;
- **a held-out fraction** (`prop_unphenotyped`) present in the markers but in **no
  trial** — genuinely unphenotyped candidates whose true `g` you nonetheless know.

``` r
arm_evaluation("stage1")
sim <- simulate_trials(n_acc = 200, n_mrk = 2000, h2 = 0.5, n_trials = 2, seed = 1)
bl  <- stage1_blues(sim$ph)
disarm_evaluation()
peek(sim$ph); peek(bl)
length(sim$unphenotyped); var(sim$g)
cor(bl$BLUE, sim$g[bl$genotype])                          # pooled: LOW, and that is correct
bl |> dplyr::group_by(studyDbId) |>
  dplyr::summarise(cor = cor(BLUE, sim$g[genotype]))      # within trial: ~0.81, ~0.83
```

- **→ returns** `sim$ph` is 640 rows (160 phenotyped lines × 2 trials × 2 reps), 40 lines
  held out, `var(sim$g)` exactly `h²` = 0.5. `bl` is 320 rows (160 genotypes × 2 studies)
  with `SE` ≈ 0.66–0.68 and `weight = 1/SE²` ≈ 2.15–2.26.
- **🔍 eyeball** the **within-trial** correlation between BLUE and the true `g` is
  ~0.81–0.83, which is what theory says it should be: with `h² = 0.5` and 2 reps,
  `sqrt(0.5 / (0.5 + 0.5/2))` ≈ 0.82. If Stage 1 is not recovering `g` this well, nothing
  in L7 is interpretable and the problem is *here*.
- **🔍 eyeball** the **pooled** correlation is much lower (~0.41) — and that is not a bug.
  Pooling BLUEs across trials mixes in the trial means, which have nothing to do with `g`.
  It is precisely why Stage 2 removes environment main effects before averaging
  (`.genotype_means()`; on this fixture that lifts the correlation to ~0.90).
- **🔍 eyeball** `peek(bl)` prints the **weight range**. Healthy weights span perhaps an
  order of magnitude. Watch for the trait×study cells that were *skipped* (`<2
  genotypes`, constant response, or a fit error) — those messages are the only sign that
  part of your data never entered the model.
- **🚩 red flag** lme4 warning `degenerate Hessian with 1 negative eigenvalue`. In a
  simulated design that usually means a design factor is **confounded with genotype** —
  e.g. every line sitting in the same block in every rep, so the block variance is not
  identifiable. `simulate_trials()` re-randomizes blocks within each trial × rep to
  avoid exactly that; a real trial with that structure gives the same warning and the
  same unusable block term.
- *Assumption:* every marker is causal and all effects are additive, so a VanRaden GRM is
  **exactly** the right kernel for this trait. Real traits are not like this — the fixture
  is deliberately the easy case, so that a failure in L7 is unambiguous.
- **🚩 red flag** a weight range spanning **>1e6**. A trial that fits (near-)perfectly —
  duplicated records, a degenerate design — drives its residual variance toward 0, every
  SE toward 0 and `1/SE²` toward `Inf`, letting that one trial dominate Stage 2 entirely.
  `SE_FLOOR_FRAC` (0.01 of the trial×trait response SD) exists to cap exactly this; set
  it to 0 and re-run to see whether it is currently doing any work.
- **🚩 red flag** `rep`/`block` **all `NA`** in the input (check with `peek(pheno$pheno)`
  at L11): the model silently degrades to plot means and the SEs become meaningless.
- *Assumption:* Stage 1 is fit **per trial**, so the BLUEs carry no across-trial
  adjustment — the two stages are joined only through `weight`. A genotype in one trial
  only is estimated from that trial alone.

### L7 — `stage2` + `cv` (offline; the scaling crux)

`G` is built from **`sim$D`** — the same markers that generated the phenotypes. That is
the whole point: build it from any other `make_dosage()` call and the kernel and the
trait are decoupled, so nothing below means anything.

``` r
arm_evaluation("stage2")
G    <- std_grm(.Gmatrix(sim$D))
gebv <- stage2_gblup(bl, list(G = G), nIter = 3000, burnIn = 1000, refresh = TRUE)
gs   <- stage2_gblup(bl, list(G = G), engine = "sommer", refresh = TRUE)
disarm_evaluation()
peek(gebv); peek(gs)
cor(gebv$GEBV, gs$GEBV)                                     # engines agree?      ~0.998
cor(gebv$GEBV, sim$g[gebv$genotype])                        # vs TRUTH, all       ~0.82
ho <- gebv$genotype %in% sim$unphenotyped                   # the honest number
cor(gebv$GEBV[ho], sim$g[gebv$genotype[ho]])                #                     ~0.49
cv <- cv_accuracy(bl, list(G = G), "Yield", nIter = 600, burnIn = 100); cv   #    ~0.52
```

- **→ returns** one row per candidate per trait (`trait, genotype, GEBV, phenotyped`);
  `cv` is a one-row tibble `trait, k, n, accuracy`. Both engines return **non-zero**
  GEBVs, correlated with each other at ~0.998.
- **🔍 eyeball** `peek(gebv)` reports **predicted-only** = 40 of 200 — the held-out lines,
  present in `G` but in no trial. That is the case that matters in production, and it is
  the only one where accuracy is not inflated by the line's own phenotype.
- **🔍 eyeball** **four numbers that should tell one story**: engine-vs-engine (~0.998),
  GEBV-vs-truth over all lines (~0.82, flattered by the phenotyped 80%), GEBV-vs-truth on
  the **held-out** lines (~0.49), and `cv_accuracy()` (~0.52). The last two should land in
  the same neighbourhood — that is the check that tells you whether `cv_accuracy()`, the
  only one of the four available on real data, can be trusted.
- **🚩 red flag — the one that motivated this fixture.** **sommer GEBVs all exactly 0.**
  REML has put σ²_g on the boundary: the kernel explains none of the phenotypic variance,
  so the BLUPs collapse to zero. On *this* fixture that is a failure. On signal-free data
  it is the **correct** answer — and BGLR cannot reproduce it, because its prior keeps
  σ²_g strictly positive, so it always returns something non-zero. A BGLR-vs-sommer
  "disagreement" therefore means nothing until you have confirmed σ²_g > 0: check the
  variance components with `code/engine_compare.R` (§9) before concluding anything.
- **🚩 red flag** held-out accuracy near zero **while** the phenotyped lines look fine.
  Either the kernel and the phenotypes are decoupled (a different `make_dosage()` call),
  or the population has no structure to exploit: markers drawn independently per line
  leave every GRM off-diagonal ≈ 0, and an unrelated line is unpredictable **by
  construction**, however correct the code is. `simulate_trials()` builds a
  founder/family population precisely so this number can be non-zero — check
  `range(G[upper.tri(G)])` (here ≈ −0.21 … 0.64) before blaming Stage 2.
- **🚩 red flag — the big one.** Stage 2 fits kernel GBLUP as **BRR on the relationship
  eigen-factor**, deliberately *not* BGLR's `model = "RKHS"`. RKHS's eigendecomposition
  shortcut assumes iid residuals, so with our `1/SE²` weights it inflates GEBVs ~λ-fold.
  Unweighted RKHS is fine; weighted is not. If GEBVs come back on roughly the phenotype
  scale instead of shrunk, suspect this first and run `code/bglr_rkhs_vs_brr.R` (§9).
- **🚩 red flag** BGLR's `weights` are **inverse SDs** (`Var ∝ 1/weights²`), not inverse
  variances. Our `w` is a precision, so the code passes `sqrt(w)`. Passing `w` directly
  is a silent mis-weighting, not an error.
- **🚩 red flag** a stale cache serving the wrong candidates. `stage2_gblup` reuses
  `data/gebv.rds` **only** when its genotype set equals `rownames(relmat)` *and* it
  covers every requested trait — resize `G` or add a trait and it must regenerate even
  with `refresh = FALSE`. Verify that: `stage2_gblup(bl, list(G = G[1:50, 1:50]))` must
  come back with 50 genotypes, not the cached 200.
- *Assumption:* prediction is **always** kernel-based; there is no marker-effect
  (BRR/BayesB-on-markers) path, so single-protocol and EM-combined runs are handled
  identically. `geno$markers` exists only to build a kernel when `G` is absent.
- *Assumption (performance, not correctness):* `eigen(relmat)` is recomputed inside
  **every** trait fit and **every** CV fold — 9 traits × 5 folds is 45+ decompositions of
  the same matrix. Expect the wait; it is not a hang.

The independent oracles for this level are `code/engine_compare.R` (variance components
plus a hand-computed GBLUP at fixed variances — whichever engine fails to reproduce the
hand solution is the broken one), `code/compare_sommer_bglr.R` (the same comparison on
the real cached data) and `code/bglr_rkhs_vs_brr.R` (RKHS vs BRR under weights). See §9.

### L8 — `select` (offline, instant)

``` r
arm_evaluation("select")
sel <- select_parents(gebv, cv = cv, trait_weights = NULL, trait_names = NULL,
                      trait_short = NULL, outfile = file.path(tempdir(), "bo.csv"))
disarm_evaluation()
names(sel$parents); nrow(sel$parents)
```

- **→ returns** with `trait_weights = NULL`: `accession`, `In_Training`, and one
  **un-standardized** GEBV column per trait — and **no `Index` column** (rows sorted by
  accession). With weights set, an `Index` column appears and rows sort by it.
- **🔍 eyeball** the GEBV columns are on each trait's own scale — the index is
  `Σ weight × GEBV` over *un-standardized* values, so the weights carry the unit
  conversion as well as the priority. A weight vector chosen as if the GEBVs were
  standardized will silently rank on whichever trait has the largest numeric spread.
- **🚩 red flag** `In_Training` all 1 (nothing new is being proposed); an `Index` column
  present when you expected the plain all-traits dump, or vice versa.
- *Assumption:* block 2 of the CSV reports each trait's **cross-validated accuracy** —
  which is a property of the training data, not of any individual accession. A
  high-index accession backed by a trait with near-zero CV accuracy is not a
  recommendation.

> **Everything below is online.** Restore the real cache paths first if you loaded
> `helper-setup.R` (`source(here::here("code","config.R"))`), then open `conn` (§3) and
> follow §5: cheapest work first.

### L9 — `connect` (online; one round-trip)

``` r
arm_evaluation("connect")
conn <- connect_t3()          # step into t3_login() to watch the token arrive
disarm_evaluation()
nchar(conn$auth_token)
```

- **→ returns** an R6 BrAPI connection with a non-empty `auth_token`.
- **🔍 eyeball** the token is non-empty. That is the *only* proof the login worked.
- **🚩 red flag** `conn$login()` returns normally on a **rejected password** — it merely
  warns and leaves `auth_token` `NULL`. And called with empty arguments it *prompts*
  (`readline`/`askpass`), which hangs any non-interactive run. `t3_login()` guards both;
  if you ever call `conn$login()` by hand, you own those two traps.
- *Assumption:* `.Renviron` is read by `connect_t3()` itself, so it works from any
  working directory and picks up edits without restarting R.

### L10 — `trials` (online; two paged GETs, then client-side filtering)

``` r
arm_evaluation("trials")
trials <- find_ny_trials(conn)             # cached: watch for the "using cached" line
disarm_evaluation()
dplyr::count(trials, role)
dplyr::filter(trials, is.na(distance_km))  # explicit trials with no coordinates
```

- **→ returns** one row per study with `studyDbId, studyName, studyType, locationDbId,
  year, distance_km, role`. With `TRAINING_TRIALS` set from `data/config/`, expect
  exactly those 9 names tagged `training`.
- **🔍 eyeball** the first line names the **mode**: `Trial selection: 9 named
  TRAINING_TRIALS -- the radius / season / study-type search is bypassed`, or
  `Trial selection: radius search -- within 500 km of (...)`. Then every configured trial
  name should appear. A name not found on the server is a **warning**, not an error —
  read the console, because the run continues with fewer training trials than you asked for.
- **🔍 eyeball — not a red flag.** Both modes fetch `/locations` and compute
  `distance_km`, so in explicit mode you will still see `Fetching locations …` and a
  `distance_km` column whose values may exceed `RADIUS_KM` or be `NA`. That join supplies
  **metadata only** (`locationName`, the report's map); it filters nothing. Confirm with
  `dplyr::filter(trials, distance_km > RADIUS_KM)` — rows there are correct in explicit
  mode and impossible in radius mode.
- **🔍 eyeball** the cache line. Either `using cached data/ny_trials.rds …` — in which
  case it was built for **this** request — or a `built for a different request —
  rebuilding` message naming what changed. Both are informative; silence would not be.
  Verify the invalidation works while you are here:
  ```r
  find_ny_trials(conn, training_trials = head(TRAINING_TRIALS, 2))   # must rebuild
  find_ny_trials(conn)                                              # must rebuild back
  find_ny_trials(conn)                                              # must now hit cache
  ```
- **🚩 red flag** a step returning trials you did not configure. That was the pre-keying
  failure: the cache filename was the whole key, so `TRAINING_TRIALS` could change from
  3 names to 9 and `find_ny_trials()` would keep handing back the old three forever. If
  you ever see a count that disagrees with `length(TRAINING_TRIALS)`, check the sidecar
  (`readRDS(cache_key_file(cache_path("ny_trials.rds")))`) before anything else.
- **🚩 red flag** `/studies` **ignores a server-side `locationDbId` filter**, so every
  study is pulled and filtered client-side. If you ever "optimize" that into a server
  query, the filter will appear to work and quietly return everything.
- *Assumption:* a trial named in **both** `TRAINING_TRIALS` and `TEST_TRIALS` counts as
  **training** — its phenotypes are used. If you meant to hold it out, it is not held out.
- *Assumption:* with `TRAINING_TRIALS` non-`NULL`, the radius/year/type search is bypassed
  entirely — `CENTER_LAT`/`RADIUS_KM`/`YEARS`/`STUDY_TYPES` do nothing.

### L11 — `phenotypes` (online; cached, else ~30 min)

``` r
arm_evaluation("phenotypes")
pheno <- get_phenotypes(conn, trials$studyDbId)
sets  <- split_by_role(pheno, trials)
disarm_evaluation()
peek(pheno$pheno)
dplyr::count(pheno$pheno, trait)
length(intersect(sets$train_pheno$studyDbId,
                 dplyr::filter(trials, role == "test")$studyDbId))   # MUST be 0
```

- **→ returns** `pheno` = `list(pheno, design, accessions)`; `sets` =
  `train_pheno / train_acc / test_acc`.
- **🔍 eyeball** the cache lines: studies already in `data/pheno_cache/` are reused and
  only genuinely new ones download, so a re-run after adding one trial should report
  one download, not the whole set. Changing `TRAIT_NAMES` must download **nothing** —
  the per-study files are stored unfiltered and the trait filter is applied at read time.
- **🔍 eyeball** `peek(pheno$pheno)` — `value` finite, and **`rep`/`block` not all
  `NA`** (they can legitimately be absent, but all-`NA` means Stage 1 silently degrades
  to plot means: the canonical subtle bug). `count(trait)` should show only your
  configured traits, and *which* of them actually have data.
- **🚩 red flag** the last line returning anything but **0** — test-trial phenotypes
  leaking into training. Also: a configured trait that yields zero rows is invisible
  otherwise, since trait matching is exact and a near-miss simply matches nothing.
- *Assumption:* `germplasmDbId` and `germplasmName` are both carried here, and this is
  where the dbId↔name map that step 04 depends on originates.

### L12 — `coverage` (online; the dbId↔name seam)

``` r
arm_evaluation("coverage")
dbids <- unique(c(sets$train_acc$germplasmDbId, sets$test_acc$germplasmDbId))
protocols <- .safe_coverage(conn$filter_geno_protocols, dbids, "protocol")
lk <- build_alias_lookup(conn, sets$train_acc$germplasmName)
disarm_evaluation()
head(protocols, 10)
length(lk); sum(names(lk) != unname(lk))      # aliases, and how many are real synonyms
```

- **→ returns** `protocols` ranked by `n_covered`; `lk` an alias→primary lookup (empty
  when `USE_SYNONYMS` is off or `T3_brapi_helpers` is unavailable).
- **🔍 eyeball** compare `protocols$n_covered` against how many of those accessions
  actually turn up as VCF samples at L14. **They will not match**, and that gap is
  expected: the server's coverage count for the Oat 3K array far exceeds what the
  archived VCF actually contains.
- **🚩 red flag** treating coverage as membership. Coverage is counted server-side by
  **`germplasmDbId`**; membership in `G` is decided by **`germplasmName`** matching a VCF
  sample after download and QC. `match_diagnostic.R` (§9) reports that funnel per
  protocol — run it whenever "why is this accession not in `G`?" comes up.
- *Assumption:* synonym canonicalization is the correct guard but, empirically, recovers
  **≈0** extra matches for the current T3/Oat targets — most disconnected accessions are
  genuinely absent from the archived VCFs. If you see a large synonym gain, something
  changed upstream and is worth understanding.

### L13 — `pedigree` (online-ish; reads the sibling repo's output)

``` r
arm_evaluation("pedigree")
d2n     <- .germplasm_name_map(PEDIGREE_DIR)
members <- .pedigree_group_members(PEDIGREE_DIR)
res     <- .resolve_names_to_dbids(head(TEST_ACCESSIONS, 20), PEDIGREE_DIR)
disarm_evaluation()
length(d2n); length(members); nrow(res)
```

- **→ returns** `d2n` a dbId→name character vector; `members` a list of dbId vectors, one
  per `<id>_group<N>.csv`; `res` a tibble of the names that resolved.
- **🔍 eyeball** unresolved `TEST_ACCESSIONS` are reported by count and sample — they can
  still be predicted if they appear as a VCF sample by name, so this is information, not
  failure.
- **🚩 red flag** this step depends on the **output contract** of
  `../BrAPI_pedigree_relmat` (group CSVs, `<id>_pedigree_groups.rds`,
  `germplasm_cache_<id>.rds` one level above `PEDIGREE_DIR`). Missing companions degrade
  *gracefully and silently* to phenotyped-only pedigree; a present-but-malformed one
  fails loudly via `.pedigree_contract_error()`. Silence here means "no pedigree", which
  is easy to mistake for "pedigree fine".
- *Assumption:* a pedigree group counts as a partial for **bridge detection**, so an
  accession genotyped on only one platform is kept when it also sits in a pedigree group.
  That is what stitches otherwise-disjoint marker GRMs together.

### L14 — `genotyping`, the one cold path (online; minutes, not hours)

Scoped deliberately: a single **small** protocol, no pedigree, forced refresh. This is
the level where you watch download → thin → QC → impute → GRM actually happen.

``` r
arm_evaluation("genotyping")
geno <- find_and_get_genotypes(conn, sets$train_acc, sets$test_acc, TEST_ACCESSIONS,
                               protocol_id = "66",       # Oat 3K array (small)
                               pedigree_dir = NULL, refresh = TRUE)
disarm_evaluation()
peek(geno$G)
setdiff(sets$train_acc$germplasmName, rownames(geno$G))   # MUST be empty
length(setdiff(TEST_ACCESSIONS, rownames(geno$G)))        # silently dropped test names
```

- **→ returns** `geno$G` subset to the **prediction targets** (training ∪ predictable
  test), plus `protocols`, `projects`, `protocol_ids`, and `markers` only in the
  single-protocol / no-pedigree / no-injection case.
- **🔍 eyeball** `peek(geno$G)` — square, symmetric, and it names any row whose
  off-diagonals are **exactly 0**: those are **injected** training accessions with no
  genotype and no pedigree (diagonal = mean diagonal). They are predicted, but from
  nothing; their GEBV is the population mean. Count them.
- **🚩 red flag** the second-to-last line returning names: training accessions are
  supposed to be **force-kept**. The last line is the mirror image and is *expected* to
  be non-zero — test accessions absent from `G` are **dropped without a word**, so this
  is the only place that loss is visible. Write the number down.
- **🚩 red flag** each protocol's imputation and VanRaden GRM are estimated on a **panel**
  (relevant accessions + fillers up to `max(GRM_PANEL_MIN, n_relevant)`), then subset back.
  Estimating on only the relevant handful makes the diagonals an artifact of who shares
  the panel. Watch the `panel: N relevant + M filler` line; if `M` is 0 on a small target
  set, the diagonals are not trustworthy (`impute_diagnostic_fullpanel.R`, §9).
- *Assumption:* `GENO_PROTOCOL_ID = NULL` (the pipeline default) means "combine **all**
  covering protocols via the EM combiner". Pinning `"66"` as above is an evaluation
  convenience — it is a *different* `G` from the one a real run builds.

### L15 — `flow` end to end (watch the funnel)

``` r
disarm_evaluation()                 # never run this armed
source(here::here("code", "run_pipeline.R"))
```

- **→ returns** `trials, pheno, sets, geno, blues, gebv, cv, selected` in the global
  environment, `output/breeders_output.csv` written, and a per-step timing table.
- **🔍 eyeball** the **funnel**, at each stage, and reconcile it top to bottom:

  | Count | From |
  |---|---|
  | accessions phenotyped | `nrow(sets$train_acc)` + `nrow(sets$test_acc)` |
  | accessions covered (server, by dbId) | `geno$protocols$n_covered` |
  | accessions in `G` | `nrow(geno$G)` |
  | of those, informed vs injected | `peek(geno$G)`'s zero-off-diagonal report |
  | accessions predicted | `dplyr::n_distinct(gebv$genotype)` |
  | rows in the output | `nrow(selected$parents)` |

  Every drop between two lines should have an explanation you can name.
- **🚩 red flag** cross-validated accuracy that is high for a trait with very few
  phenotyped genotypes; every GEBV identical within a trait; a timing table where a step
  you expected to be cached took minutes.
- *Assumption:* the report (`analysis/brapi_selection_pipeline.Rmd`) runs the same
  functions with the same caches — so if the Rmd and the driver disagree, it is a cache
  or config difference, not a code difference.

------------------------------------------------------------------------

## 5. Fast → slow within the live layer

| Cache | Written by | Cost to rebuild | Keyed on (a change here rebuilds it) |
|----|----|----|----|
| `data/ny_trials.rds` | step 02 | seconds–a minute (2 paged GETs) | `TRAINING_TRIALS`, `TEST_TRIALS`; + radius/year/type **only** in radius mode |
| `data/pheno_cache/<sid>.rds` | step 03 | one study's download | nothing — one file per study, **unfiltered by trait**, reused forever |
| `data/phenotypes.rds` | step 03 | seconds from the per-study files | the study set, `TRAIT_NAMES` |
| `data/synonym_map.rds` | `synonyms.R` | one `/search/germplasm` call | the accession name set, `USE_SYNONYMS` |
| `data/vcf_cache/*.vcf` | step 04 | minutes to **hours** (GB-scale GBS files) | nothing — keyed by protocol/project, reused forever |
| `data/genotypes.rds` | step 04 | minutes (cached VCFs) to hours (cold) | the accessions, `GENO_PROTOCOL_ID`, `PEDIGREE_DIR`, EM df, marker QC/impute, `SEED` |
| `data/gebv.rds` / `gebv_sommer.rds` | step 06 | minutes per trait | validated by **content**: `G`'s candidate set + the requested traits |

Each step records its request in a sidecar `data/<name>.key.rds` and rebuilds when it
differs — so **you do not delete caches by hand**; change the config and the affected
steps rebuild themselves, saying what changed. Two things worth knowing while
evaluating: **deleting one cache forces every downstream step**, and a cache with **no
sidecar** (built by hand, or before this mechanism) is treated as stale and rebuilt.

------------------------------------------------------------------------

## 6. Subtle-bug catalogue

| Signature | What it means | Expose it with |
|----|----|----|
| GEBVs on roughly the phenotype scale, not shrunk | weighted **RKHS** instead of BRR on the eigen-factor — inflated ~λ-fold | `code/bglr_rkhs_vs_brr.R`; L7 |
| GEBVs off by a smooth factor vs sommer | BGLR `weights` are inverse **SDs**; passing precision `w` instead of `sqrt(w)` | `code/engine_compare.R`, `compare_sommer_bglr.R` |
| **sommer GEBVs all exactly 0** while BGLR is non-zero | REML put σ²_g on the boundary — the kernel explains nothing. Correct behaviour on signal-free data; BGLR's prior cannot return exactly 0, so the engines *always* look like they disagree there | L7 on `simulate_trials()`; variance components via `engine_compare.R` |
| held-out accuracy ≈ 0 while phenotyped lines look fine | no relatedness to exploit (GRM off-diagonals ≈ 0), or the kernel and phenotypes are decoupled | `range(G[upper.tri(G)])`; L7 |
| non-empty dosage matrix, **zero rowname overlap** | synonym / name mismatch (VCF under a preliminary name) | `peek(dos, accessions=)`; `match_diagnostic.R` |
| coverage count ≫ accessions in `G` | the dbId-vs-name seam, or the archived VCF is thinner than the server claims | `match_diagnostic.R`, `coverage_diagnostic.R` |
| every accession ~50% missing after merge | markers keyed by VCF `ID` instead of canonical `CHROM_POS`; or QC dropped accessions before markers | L5 |
| GRM diagonals dipping below 1 / unstable | imputation + allele frequencies estimated on too few relevant accessions | `impute_diagnostic*.R`; the `panel:` line at L14 |
| rows in `G` with all off-diagonals exactly 0 | **injected** training accessions — predicted from nothing, GEBV = population mean | `peek(geno$G)` |
| `TEST_ACCESSIONS` missing from the output | absent from `G`; dropped **silently** by design | the `setdiff` at L14 |
| BLUE weights spanning >1e6 | a degenerate trial driving `1/SE²` → ∞ and dominating Stage 2 | `peek(blues)`; `SE_FLOOR_FRAC` |
| `rep`/`block` all `NA` | Stage 1 silently degraded to plot means | `peek(pheno$pheno)` |
| combined `G` one row/col too big | the EM phantom variable was not dropped | L4 |
| pedigree silently absent | a missing sibling-repo companion file degrades gracefully | L13 |
| a step returns results that do not match the current config | a cache served for a *different* request — the failure request-keying exists to prevent | the sidecar: `readRDS(cache_key_file(cache_path("<name>.rds")))`; L10 |
| CV accuracy high on a trait with few phenotyped genotypes | too little data for the correlation to mean anything | `cv` table at L15 |

The meta-rule, worth more than any single check: **trust agreement between two
independent derivations; distrust a single number.**

------------------------------------------------------------------------

## 7. The test suite (offline, deterministic)

``` bash
Rscript tests/run_tests.R                       # tiers 1 + 2 (offline)
RUN_LIVE_TESTS=true Rscript tests/run_tests.R   # also tier 3 (needs .Renviron creds)
```

| Command | Expected |
|----|----|
| `tests/run_tests.R` | `FAIL 0`, `WARN 0`, `SKIP 4`, `PASS 291` — the 4 skips are the live tier |
| with `RUN_LIVE_TESTS=true` and credentials | the same, with the 4 live tests running instead of skipping |
| with `RUN_LIVE_TESTS=true` and **no** credentials | still 4 skips, not 4 errors |

If a count drifts after a change, that is the regression signal — reconcile it before
moving on. Tier 1 covers the pure helpers on synthetic inputs, tier 2 runs the real
models and a fully mocked `find_and_get_genotypes`, tier 3 hits the live server.
`test-01-evaluation.R` additionally asserts that **every function named in
`EVAL_GROUPS` exists** — that is what keeps this document honest as the code moves.

------------------------------------------------------------------------

## 8. What the code touches

| Path | Role | Tracked? |
|----|----|----|
| `.Renviron` | `T3_USERNAME` / `T3_PASSWORD` | **no** (`.Renviron.example` is) |
| `data/config/*.txt` | the selection lists that define a run | **yes** — the one tracked part of `data/` |
| `data/*.rds`, `data/vcf_cache/` | every step's cache | no |
| `output/breeders_output.csv`, `trait_cvs.rds` | results | no |
| `../BrAPI_pedigree_relmat/output/<id>/` | pedigree group matrices + companions | external contract (L13) |
| `../T3_brapi_helpers/` | synonym helpers | external, optional |

------------------------------------------------------------------------

## 9. Validation tooling — the independent oracles

`peek()` and the debugger tell you what the pipeline *did*. These scripts tell you
whether it was *right*, by re-deriving the answer another way. Run them with
`source(here::here("code", "<script>"))`; all but the noted ones work offline from the
caches.

| Script (`code/…`) | What it independently re-derives | Read a failure as |
|----|----|----|
| `bglr_rkhs_vs_brr.R` | the same variance model fitted as RKHS-on-`G` vs BRR-on-eigen-factor, under weights | if they disagree, the weighted-RKHS inflation is live again |
| `bglr_diagnose.R` | kernel conditioning + unweighted-vs-weighted BGLR | near-singular `G` (injected/disconnected rows) rather than a coding bug |
| `engine_compare.R` | variance components + a **hand-computed** GBLUP at fixed variances | whichever engine fails to reproduce the hand solution is the broken one |
| `compare_sommer_bglr.R` | the whole Stage 2 under sommer REML instead of BGLR | a systematic scale gap = weighting; scatter = MCMC noise |
| `covcomb_check.R` | our EM combiner vs `CovCombR::CovComb` on a **real** partial subset | a combiner bug, not a data problem |
| `df_grid_diagnostic.R` | how `G` responds to `PEDIGREE_DF` × `GRM_DF_MEAN` | GEBVs that are unshrunk / diagonals below 1 are a *weighting choice*, not a bug |
| `impute_diagnostic.R` | glmnet vs mean imputation, effect on low GRM diagonals | mean imputation deflating heavily-missing accessions' diagonals |
| `impute_diagnostic_fullpanel.R` | full-panel vs subset-first frequency/imputation estimation | diagonals that are an artifact of who shared the panel |
| `impute_diagnostic_robust.R` | the robust per-marker retry vs the all-or-nothing original | one bad marker aborting a whole protocol's imputation |
| `coverage_diagnostic.R` | per training accession: genotyped? in a pedigree group? informed in `G` or injected? | the real reason `G` is sparser than coverage suggested |
| `match_diagnostic.R` | the funnel from coverage (dbId) → VCF sample (name) → GRM membership | exactly where accessions are lost across the two keying spaces |
| `gebv_diagnostic.R` | whether extreme GEBVs belong to the data-rich, the low-diagonal, or the phenotyped | extreme GEBVs concentrated on low-diagonal accessions = a kernel problem |

------------------------------------------------------------------------

## 10. Extending

1. Add the function to its step file, keeping the existing conventions (tidyverse,
   `|>`, `here::here()`, `pkg::fn()`).
2. Extend `tests/testthat/` — an **oracle** test (compare against an independently
   computed answer) rather than a snapshot.
3. Update the living docs: `README.md` (usage), `DESIGN.md` (structure),
   `BACKGROUND.md` (why), `code/README.md` (the function table).
4. If it is slow, give it a bar (`code/progress.R`, and see `CLAUDE.md` for the two cli
   gotchas).
5. If it is worth stepping through, **add its name to the right group in
   `code/evaluation.R::EVAL_GROUPS`** so `arm_evaluation()` reaches it — and add a level
   here. `test-01-evaluation.R` will fail if a name in `EVAL_GROUPS` stops existing, but
   nothing can tell you about a function you never added.
