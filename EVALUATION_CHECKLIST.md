# Selection-pipeline evaluation checklist

Tick these off over time. Each line names the `arm_evaluation()` group; the full
step-by-step walkthrough lives in **`EVALUATION.md`** at the matching level (L1…L15).
The order is offline → cached-online → cold on purpose: work top to bottom so a cheap
bug surfaces before an expensive one. `disarm_evaluation()` at the end of each level,
and **always before an unattended run or `wflow_build()`**.

Bootstrap once (see `EVALUATION.md` §3), then:

## Offline — no server, no `data/`, do these first

- [ ] **L1 `config`** — trait weights/short names are named vectors keyed by
  `TRAIT_NAMES`; pipe-bearing names survive verbatim; a **missing config file errors**
  rather than falling back to the radius search.
- [ ] **L2 `progress`** — silent and value-transparent with reporting off; timings
  accumulate one row per step.
- [ ] **L3 `grm`** — raw VanRaden diagonal > 1 (`1 + F`), `std_grm()` mean diagonal
  exactly 1, `.effective_n` ≈ 27.9/30 on the fixture; `.impute_glmnet` reproducible and
  leaves no `NA`.
- [ ] **L4 `combine`** — `psi` is one row/col larger than the name set (**phantom
  variable dropped**); result symmetric; converges well before iteration 100;
  `.center_dfs(c(20,25,30), 60, 15)` → `55 60 65`.
- [ ] **L5 `markers`** — markers keyed **`CHROM_POS`** (`1_100`…), not the VCF `ID`;
  QC drops markers *before* accessions; `peek(dos, accessions=)` raises the zero-overlap
  alarm on a deliberate mismatch.
- [ ] **L6 `stage1`** — on `simulate_trials()`: within-trial `cor(BLUE, true g)` ≈ 0.82
  (theory: `sqrt(h²/(h² + (1−h²)/n_rep))`); `weight = 1/SE²` spans a sane range; no lme4
  `degenerate Hessian` warning; skipped trait×study cells read and understood; check
  whether `SE_FLOOR_FRAC` is currently doing any work.
- [ ] **L7 `stage2` + `cv`** — `G` built from **`sim$D`**; both engines **non-zero** and
  agreeing (~0.998); held-out `cor(GEBV, true g)` ≈ 0.49 and `cv_accuracy()` ≈ 0.52 land
  in the same neighbourhood; GEBVs **shrunk**, not on the phenotype scale (the
  weighted-RKHS trap); `sqrt(w)` passed to BGLR; the cache **regenerates** when `G`'s
  candidate set changes. All-zero sommer output ⇒ check σ²_g before suspecting the code.
- [ ] **L8 `select`** — `TRAIT_WEIGHTS = NULL` ⇒ no `Index` column; the index is over
  **un-standardized** GEBVs, so the weights carry the unit conversion.
- [ ] **Test suite** — `Rscript tests/run_tests.R` → `FAIL 0`, `SKIP 4`, `PASS 310`.

## Online — cached first, cheapest server work first

- [ ] **L9 `connect`** — non-empty `conn$auth_token` (the only proof the login worked);
  the two `conn$login()` traps understood (silent on a bad password, prompts on empty).
- [ ] **L10 `trials`** — the trial count matches `length(TRAINING_TRIALS)`; the cache
  either reports a hit or names what changed and rebuilds; every configured trial name
  found (a miss is only a *warning*); a trial in both lists counts as **training**;
  `/studies` filtered client-side.
- [ ] **L11 `phenotypes`** — adding a trial downloads only that trial; changing
  `TRAIT_NAMES` downloads nothing; `peek(pheno$pheno)`: `value` finite, **`rep`/`block` not
  all `NA`**; every configured trait actually has data; **zero** test-trial studies in
  `train_pheno`.
- [ ] **L12 `coverage`** — server coverage (by `germplasmDbId`) reconciled against GRM
  membership (by `germplasmName`); synonym gain ≈ 0 as expected.
- [ ] **L13 `pedigree`** — the sibling repo's companions present; silence means "no
  pedigree", not "pedigree fine".
- [ ] **L14 `genotyping` (the one cold path)** — `protocol_id = "66"`, `refresh = TRUE`:
  no training accession missing from `G`; **count the injected rows** (all off-diagonals
  exactly 0) and the **silently dropped** `TEST_ACCESSIONS`; the `panel:` line shows real
  filler.
- [ ] **L15 `flow`** — `run_pipeline.R` end to end; every drop in the funnel
  (phenotyped → covered → in `G` → informed → predicted → output rows) has a named
  explanation; timing table matches expectations.

## Whole system

- [ ] **An independent oracle agrees** — at least one script from §9 re-derives a
  headline result (start with `match_diagnostic.R` for membership and
  `compare_sommer_bglr.R` for the GEBV scale).
- [ ] **Report renders** — `workflowr::wflow_build("analysis/brapi_selection_pipeline.Rmd")`
  agrees with the driver's numbers.
- [ ] **Clean shutdown** — `armed_functions()` returns `character(0)` so no background
  run will block on the debugger.

------------------------------------------------------------------------

*Notes / anomalies noticed (the "that's a little funny" observations — the point of the
exercise):*

-

-

-
