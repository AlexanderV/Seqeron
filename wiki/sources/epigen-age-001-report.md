---
type: source
title: "Validation report: EPIGEN-AGE-001 (epigenetic age — Horvath 2013 multi-tissue + Horvath 2018 skin&blood + Levine 2018 PhenoAge DNAm clocks, EpigeneticsAnalyzer.CalculateEpigeneticAge / CalculateSkinBloodAge / CalculatePhenoAge)"
tags: [validation, epigenetics, governance]
doc_path: docs/Validation/reports/EPIGEN-AGE-001.md
sources:
  - docs/Validation/reports/EPIGEN-AGE-001.md
source_commit: 4416f109abc62e50a964a9613c0361828869e163
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Validation report: EPIGEN-AGE-001

The two-stage **validation write-up** for test unit **EPIGEN-AGE-001** (epigenetic age — DNA-methylation
"clocks" estimating chronological age from CpG β-values), validated 2026-06-24. This is the *report*
artifact that feeds one row of the [[validation-ledger]]; it records the validator's independent
**verdict** on both the algorithm description (Stage A) and the shipped code (Stage B), and the wider
campaign is [[validation-and-testing]]. The clocks, their two-stage model, transform branches, oracles
and edge cases are synthesized in the concept [[epigenetic-age-horvath-clock]]; [[test-unit-registry]]
defines the unit. Distinct from [[epigen-age-001-evidence]] — the pre-implementation evidence artifact
sourced from `docs/Evidence/` — this page is the independent two-stage re-validation verdict.

## Verdict

**Stage A: PASS · Stage B: PASS · End state: CLEAN.** No defect found; no code changed. This
re-validation **supersedes** the 2026-06-15 report, which predated the embedding of the skin&blood and
PhenoAge clocks (commit `08b38201`): the multi-tissue clock (commit `9edac5ff`) was re-checked and the
two new clocks were independently validated for the first time. The filtered suite
(`EpigeneticsAnalyzer_CalculateEpigeneticAge_Tests`) ran **34 passed, 0 failed**; build 0 warnings /
0 errors. No code changed this session, so per protocol the full suite was not re-run (last left green
by the 2026-06-23 commits).

## Canonical methods & source under test

In `EpigeneticsAnalyzer.cs` (age-clock region) plus three coefficient-table partials:

- `CalculateEpigeneticAge(methylation)` (`:1145–1152`) — built-in **Horvath 2013 353-CpG multi-tissue**
  clock; delegates to the caller-supplied overload with the embedded table + intercept `0.695507258`.
- `CalculateEpigeneticAge(methylation, coefficients, intercept)` (`:1171–1197`) — caller-supplied path:
  null/null/empty guards, `Y = intercept + Σ coef·β` over CpGs in the table, then `HorvathAntiTransform(Y)`.
- `HorvathAntiTransform` (`:1206–1211`) — the `anti.trafo`: `Y<0 ? 21·exp(Y)−1 : 21·Y+20`
  (`HorvathAdultAge = 20.0`). Byte-for-byte the sourced transform.
- `CalculateSkinBloodAge(methylation)` (`:1228–1235`) — built-in **Horvath 2018 391-CpG skin&blood**;
  delegates to the `anti.trafo` overload with intercept **−0.447119319**.
- `CalculatePhenoAge(methylation)` / `CalculatePhenoAge(methylation, coefficients, intercept)`
  (`:1252–1303`) — **Levine 2018 513-CpG DNAm PhenoAge**: `age = intercept + Σ weight·β`, returned
  **without** any transform (intercept 60.664).
- Coefficient tables: `EpigeneticsAnalyzer.HorvathClock.cs` (353), `.SkinBloodClock.cs` (391),
  `.PhenoAgeClock.cs` (513).
- Tests: `EpigeneticsAnalyzer_CalculateEpigeneticAge_Tests` (34 tests).

## Stage A — description (algorithm faithfulness)

Confirmed against independently fetched references: **biolearn `biolearn/model.py`** (the reference
Python reimplementation), the biolearn `Horvath2.csv` (skin&blood) and `PhenoAge.csv` coefficient
tables, the aldringsvitenskap/epigeneticclock reference R (`horvath2013.R`, `AdditionalFile3.csv`), and
the three primary papers — **Horvath 2013** (*Genome Biol.* 14:R115, 353 CpGs), **Horvath et al. 2018**
(*Aging* 10(7):1758, 391 CpGs), **Levine et al. 2018** (*Aging* 10(4):573, 513 CpGs).

Formula checks (all ✓):

- **`anti.trafo`** `F⁻¹(Y)` = `21·exp(Y)−1` for `Y<0`, else `21·Y+20` (adult.age = 20) — matches
  biolearn `anti_trafo` byte-for-byte and the 2013 reference R. The **strict `<`** at the boundary puts
  `Y=0` on the linear branch → exactly 20.
- **Multi-tissue & skin&blood:** `age = anti.trafo(intercept + Σ coef·β)`; the skin&blood intercept
  −0.447119319 is the value biolearn feeds into `anti_trafo`.
- **PhenoAge:** `age = intercept + Σ weight·β`, returned **untransformed** (biolearn has no `transform`
  field for PhenoAge → identity; Levine 2018 linear predictor is already in years).
- **Non-clock CpG handling:** only CpGs present in the coefficient table contribute (inner join / matrix
  product in both reference implementations).

**Independent coefficient-table cross-check** (field-by-field numeric diff of the embedded C# tables
against the biolearn CSVs): skin&blood 391 vs 391 rows, key sets equal, **max abs coef diff 0.0**;
PhenoAge 513 vs 513 rows, key sets equal, **max abs coef diff 0.0**. The multi-tissue 353-CpG table was
verified byte-identical in the prior session and is unchanged; spot-checked probes this session all
matched, including the scientific-notation entry `cg00431549 = 8.83e-6`. Intercepts confirmed against
source: 0.695507258 / −0.447119319 / 60.664.

**Independent age cross-check** (hand-computed in Python from the validated formulas, every value matched
the spec/tests): multi-tissue empty map → `21·Y+20 = 34.605652418`; multi-tissue `cg00864867 β=1` →
68.200697468; multi-tissue two-probe negative predictor → exp branch `10.936325872311789`; skin&blood
empty map → exp branch `12.428819664840216`; skin&blood `cg12140144 β=1` → 18.309250637525345; PhenoAge
empty map → 60.664; PhenoAge `cg15611364 β=1` → 123.78815047; PhenoAge two probes β=0.5 → 70.22137867.
A **negative control** confirms PhenoAge must **not** be anti-transformed: `anti.trafo(60.664) = 1293.944`
is the value the code must avoid.

## Stage B — implementation

Every code path above was re-traced and reproduced exactly. The parameterless built-in overloads
(E6 / SB5 / PA8) assert equality with the explicit overloads called with the built-in tables/intercepts
(within 1e-12). **PA6 is a negative control** locking in "no transform": PhenoAge ≠
`HorvathAntiTransform(60.664)` (= 1293.944). Coefficient-table counts confirmed by grep (353/391/513) and
by the numeric diff above.

**Test-quality audit (hard gate):** E1/SB1/PA1 assert exact counts (353/391/513) and exact intercepts;
E2/SB2/PA2 assert exact named-probe coefficients incl. the scientific-notation entry — **sourced, not
code echoes**. Both `anti.trafo` branches, boundary `Y=0`, intercept-only, monotonicity (S3) and the
PhenoAge no-transform negative control are covered; all three error contracts (null map, null coef,
empty coef) are exercised per clock. Filtered suite green 34/0; build clean.

## Findings

- **No defect, no code change (state CLEAN).** All three embedded tables are numerically identical to the
  independent biolearn source (max abs diff 0.0 over the 391 and 513 pairs; multi-tissue byte-identical in
  the prior session); transforms and intercepts match the reference implementation exactly. PhenoAge
  correctly applies **no** anti-transform; both Horvath clocks share the same adult.age=20 `anti.trafo`.
- **Out of scope (unchanged):** `PredictImprintedGenes` is listed under the unit in the registry but is
  unrelated to the age clocks and has no retrieved authoritative basis; not modified, not tested here.
