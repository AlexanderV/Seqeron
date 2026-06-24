# Validation Report: EPIGEN-AGE-001 — Epigenetic Age Estimation (DNA-methylation clocks)

- **Validated:** 2026-06-24   **Area:** Epigenetics
- **Canonical method(s):**
  `EpigeneticsAnalyzer.CalculateEpigeneticAge(methylation)` (built-in Horvath 2013 353-CpG multi-tissue),
  `CalculateEpigeneticAge(methylation, coefficients, intercept)` (caller-supplied, anti.trafo path),
  `CalculateSkinBloodAge(methylation)` (built-in Horvath 2018 391-CpG skin&blood),
  `CalculatePhenoAge(methylation)` / `CalculatePhenoAge(methylation, coefficients, intercept)` (Levine 2018 513-CpG, no transform),
  `HorvathAntiTransform(transformedAge)`.
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS

This re-validation supersedes the 2026-06-15 report, which predated the embedding of the
skin&blood and PhenoAge clocks (commit `08b38201`). The multi-tissue clock (commit `9edac5ff`)
was re-checked; the two new clocks were independently validated for the first time.

## Stage A — Description

### Sources opened this session (independent of repo artifacts)

1. **biolearn `biolearn/model.py`** (fetched verbatim via GitHub raw) — the reference Python
   reimplementation of these clocks. Confirms verbatim:
   - `anti_trafo(x, adult_age=20)` (line 21–23): `x < 0 → (1+20)·exp(x)−1, else (1+20)·x+20`.
   - **Horvathv2 (skin&blood)** (line 437–438): `"transform": lambda sum: anti_trafo(sum - 0.447119319)`
     → uses anti.trafo (adult_age=20) with intercept **−0.447119319**.
   - **PhenoAge** (line 382–390): the model has **no `transform` field** → defaults to
     `no_transform` (identity, line 1479–1486). Intercept is the `intercept` row in `PhenoAge.csv`
     (= 60.664) consumed as coefficient×1 (line 1496). → DNAm PhenoAge = intercept + Σ weight·β, untransformed.
2. **biolearn `Horvath2.csv`** — skin&blood coefficient table; 391 CpG rows.
3. **biolearn `PhenoAge.csv`** — PhenoAge weight table; 513 CpG rows + an `intercept,60.664` row.
4. **aldringsvitenskap/epigeneticclock `horvath2013.R` / `StepwiseAnalysis.R` / `AdditionalFile3.csv`**
   (prior session) — the 2013 `trafo`/`anti.trafo` (adult.age=20), `predictedAge =
   anti.trafo(intercept + meth·coef)`, and the 353-CpG `CoefficientTraining` table (intercept 0.695507258).
5. **Horvath S (2013) Genome Biology 14:R115**, **Horvath et al. (2018) Aging 10(7):1758-1775**,
   **Levine et al. (2018) Aging 10(4):573-591** — the three primary papers (353 / 391 / 513 CpGs).

### Formula check

- **anti.trafo F⁻¹(Y):** `21·exp(Y)−1` for `Y < 0`, else `21·Y + 20` (adult.age=20). Matches
  biolearn `anti_trafo` byte-for-byte and the 2013 reference R. The strict `<` at the boundary
  places Y=0 on the linear branch → exactly 20. ✓
- **Multi-tissue & skin&blood:** age = `anti.trafo(intercept + Σ coef·β)`. The skin&blood intercept
  −0.447119319 is the value biolearn feeds into `anti_trafo`. ✓
- **PhenoAge:** age = `intercept + Σ weight·β`, returned **untransformed** (biolearn: no transform
  field; Levine 2018 Methods linear predictor in years). ✓
- **Non-clock CpG handling:** only CpGs present in the coefficient table contribute (matrix product /
  inner join in both reference implementations). ✓

### Independent coefficient-table cross-check (numeric, all rows)

A field-by-field numeric comparison of the embedded C# tables against the biolearn CSVs:

| Clock | Embedded n | Source n | Key sets equal | Max abs coef diff |
|-------|-----------|----------|----------------|-------------------|
| Skin&blood (Horvath2.csv) | 391 | 391 | yes | **0.0** |
| PhenoAge (PhenoAge.csv)   | 513 | 513 | yes | **0.0** |

The multi-tissue table (353 CpGs, intercept 0.695507258) was verified byte-identical against the
Springer supplement + GitHub mirror in the prior session and is unchanged. Spot-checked probes
this session (E2/SB2/PA2) all matched, including the scientific-notation entry cg00431549 = 8.83e-6.

Intercepts confirmed against source: multi-tissue 0.695507258; skin&blood −0.447119319 (biolearn
transform); PhenoAge 60.664 (`intercept` row of PhenoAge.csv).

### Independent age cross-check (hand-computed from the validated formulas)

| Case | Linear predictor | Branch | Expected age |
|------|------------------|--------|--------------|
| Multi-tissue, empty map | 0.695507258 | linear | 21·Y+20 = **34.605652418** |
| Multi-tissue, cg00864867 β=1 | 2.295271308 | linear | **68.200697468** |
| Multi-tissue, cg09809672+cg27544190 β=1 | −0.564936093 | exp | 21·e^Y−1 = **10.936325872311789** |
| Skin&blood, empty map | −0.447119319 | exp | 21·e^Y−1 = **12.428819664840216** |
| Skin&blood, cg12140144 β=1 | −0.083938149 | exp | **18.309250637525345** |
| PhenoAge, empty map | 60.664 | none | **60.664** |
| PhenoAge, cg15611364 β=1 | — | none | 60.664+63.12415047 = **123.78815047** |
| PhenoAge, two probes β=0.5 | — | none | **70.22137867** |
| (PhenoAge negative-control) anti.trafo(60.664) | — | linear | 1293.944 (must NOT be applied) |

All recomputed in Python from the validated formulas; every value matches the spec/tests.

**Stage A findings:** All three clocks (transforms, intercepts, coefficient tables) match the
independent biolearn reference implementation and the primary papers. PhenoAge correctly uses **no**
anti-transform; both Horvath clocks use the same adult.age=20 anti.trafo. **PASS.**

## Stage B — Implementation

### Code path reviewed
- `EpigeneticsAnalyzer.cs:1206-1211` `HorvathAntiTransform`: `Y<0 ? (1+20)·exp(Y)−1 : (1+20)·Y+20`,
  `HorvathAdultAge = 20.0`. Byte-for-byte the sourced anti.trafo. ✓
- `EpigeneticsAnalyzer.cs:1171-1197` caller-supplied multi-tissue/anti.trafo overload: null/null/empty
  guards, `Y = intercept + Σ coef·β` over CpGs in the table, then `HorvathAntiTransform(Y)`. ✓
- `EpigeneticsAnalyzer.cs:1145-1152` parameterless multi-tissue overload delegates with the built-in
  353-CpG table + 0.695507258. ✓
- `EpigeneticsAnalyzer.cs:1228-1235` `CalculateSkinBloodAge` delegates to the anti.trafo overload with
  the built-in 391-CpG table + intercept −0.447119319. ✓
- `EpigeneticsAnalyzer.cs:1252-1303` `CalculatePhenoAge` (both overloads): `age = intercept + Σ weight·β`,
  returned **without** any transform; same null/empty contract. ✓
- Coefficient tables: `EpigeneticsAnalyzer.HorvathClock.cs` (353), `.SkinBloodClock.cs` (391),
  `.PhenoAgeClock.cs` (513) — counts confirmed by grep (353/391/513) and by numeric table diff.

### Cross-verification vs code
Filtered run of `EpigeneticsAnalyzer_CalculateEpigeneticAge_Tests`: **34 passed, 0 failed**. Every
expected value is one of the externally hand-computed constants above (34.605652418, 68.200697468,
10.936325872311789, 12.428819664840216, 18.309250637525345, 60.664, 123.78815047, 70.22137867, …).

### Variant/delegate consistency
- E6 / SB5 / PA8 assert the parameterless built-in overloads equal the explicit overloads called with
  the built-in tables/intercepts (within 1e-12). ✓
- PA6 is a negative control: it asserts PhenoAge ≠ `HorvathAntiTransform(60.664)` (= 1293.944),
  locking in "no transform". ✓

### Test quality audit (HARD gate)
- Tables: E1/SB1/PA1 assert exact counts (353/391/513) and exact intercepts; E2/SB2/PA2 assert exact
  named-probe coefficients incl. the scientific-notation entry. Sourced, not code echoes. ✓
- Both anti.trafo branches, boundary Y=0, intercept-only, monotonicity (S3), and the PhenoAge
  no-transform negative control are covered; all three error contracts (null map, null coef, empty coef)
  exercised per clock. ✓
- Filtered suite green (34/0); build 0 warnings / 0 errors. No code changed this session, so the full
  suite was not re-run (protocol requires it only on code change); the 2026-06-23 commits last left it green.

### Findings / defects
- **None.** All three embedded tables are numerically identical to the independent biolearn source
  (max abs diff 0.0 over 391 and 513 pairs; multi-tissue byte-identical in the prior session). Transforms
  and intercepts match the reference implementation exactly.
- Out of scope (unchanged): `PredictImprintedGenes` is listed under the unit in the registry but is
  unrelated to the age clocks and has no retrieved authoritative basis; not modified, not tested here.

## Verdict & follow-ups
- **Stage A: PASS.** Three clocks validated against the biolearn reference implementation + primary
  papers: multi-tissue (353, intercept 0.695507258, anti.trafo), skin&blood (391, intercept −0.447119319,
  anti.trafo), PhenoAge (513, intercept 60.664, **no** transform).
- **Stage B: PASS.** Code realises each formula exactly; embedded tables numerically identical to source;
  34/34 unit tests pass; delegate overloads consistent; PhenoAge no-transform locked by a negative control.
- **End-state: CLEAN.** No defect found. No code changed.
