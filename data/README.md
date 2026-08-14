# Data

Cached BrAPI pulls and intermediate results (`ny_trials.rds`, `phenotypes.rds`,
`genotypes.rds`, `gebv*.rds`, `vcf_cache/`, …). Everything here is regenerable and
**gitignored** — with one exception, below.

## `config/` — the selection lists (tracked in git)

The long lists that define a run live here as plain text rather than as literal
vectors in `code/config.R`, and are read by `config_lines()` / `config_traits()`:

| File (example) | Setting |
|----------------|---------|
| `Trial_Sel2026_Intersect20NoYld.txt` | `TRAINING_TRIALS` (and/or `TEST_TRIALS`) — `studyName`s |
| `Acc_Sel2026.txt` | `TEST_ACCESSIONS` — `germplasmName`s |
| `Trait_Sel2026.txt` | `TRAIT_NAMES` + optional `TRAIT_WEIGHTS` / `TRAIT_SHORT_NAMES` |

**Format:** one value per line. Blank lines and lines starting with `#` are ignored,
surrounding whitespace is trimmed, duplicates are dropped with a message, and a
trailing newline is optional. Everything else is used **verbatim** — a
`germplasmName` such as `APPLER|CIAV2680` keeps its pipe. A file named in `config.R`
but missing here is an **error**, not an empty list.

**Trait file:** 1 to 3 **tab-separated** columns — trait name, selection-index weight,
short label:

```
Grain yield - g/m2|CO_350:0000260	0.5	Yield g/m2
Freeze damage severity - 0-9 Rating|CO_350:0005001	-5	Freeze_damage 0-9
```

A column that is used must be filled on every line; a column omitted everywhere
yields `NULL` (no selection index / full trait names). See README.md section 3.
