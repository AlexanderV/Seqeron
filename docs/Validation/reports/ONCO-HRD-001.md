# Validation Report: ONCO-HRD-001 — Homologous Recombination Deficiency (HRD) composite genomic-scar score

- **Validated:** 2026-06-16   **Area:** Oncology
- **Canonical method(s):** `OncologyAnalyzer.CalculateHRDScore(int loh, int tai, int lst)`, `OncologyAnalyzer.ClassifyHRDStatus(int score)`, `OncologyAnalyzer.DetectHRD(HrdComponents)` (plus public const `HrdHighScoreThreshold`, enum `HrdStatus`, records `HrdComponents`/`HrdResult`)
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS

## Scope

The Registry/TestSpec/Evidence for ONCO-HRD-001 deliberately scope this unit to the **composite-sum +
threshold classification**: it takes the three already-computed genomic-scar component counts
(LOH, TAI, LST) and returns their unweighted sum and the HRD-high classification against the 42 cutoff.
Per-segment derivation of the components from raw allele-specific copy-number data (Abkevich 15 Mb LOH
segmentation, Birkbak sub-telomere/centromere TAI geometry, Popova LST smoothing) is **out of scope here**
and lives in the separate ONCO-LOH-001 / ONCO-CNA-001 units. The methods `DetectLOH`,
`CalculateHrdLohScore`, and `CalculateLOHFraction` in the same file are part of the `#region Loss of
heterozygosity (ONCO-LOH-001)` block and are therefore not validated under this unit. (The checklist's
method table for ONCO-HRD-001 lists aspirational `CalculateLOHScore`/`CalculateTAIScore`/`CalculateLSTScore`
names that do not exist as methods; the actual implemented surface is the three composite methods above.)

## Stage A — Description

### Sources opened & what they confirm (retrieved this session)

1. **Telli ML et al. (2016), Clin Cancer Res 22(15):3764–3773** — primary source.
   The PubMed page (https://pubmed.ncbi.nlm.nih.gov/26957554/) was reCAPTCHA-gated on WebFetch, but the
   abstract content was retrieved this session via WebSearch (which returned the indexed abstract text):
   verbatim *"a combined homologous recombination deficiency (HRD) score, calculated as an unweighted sum
   of LOH …, TAI …, and LST … scores"* and *"HR deficiency was defined as HRD score ≥42 or BRCA1/2
   mutation."* → confirms **both** the unweighted-sum formula and the inclusive **≥ 42** cutoff.
2. **Stewart MD et al. (2022), Oncologist 27(3):167–174** — peer-reviewed review (WebFetched
   https://pmc.ncbi.nlm.nih.gov/articles/PMC8914493/): myChoice CDx *"determines HR status by … genomic
   instability … (gLOH + TAI + LST…)"* → confirms the three-way sum independently. Cutoff: *"other trials
   on niraparib and olaparib used a cutoff of 42"* (and notes a veliparib trial used 33). → corroborates 42.
3. **scarHRD reference implementation** (WebFetched https://github.com/sztup/scarHRD) — independent
   reference tool. Documented worked examples show `HRD-sum` is the literal arithmetic sum of the three
   columns (see cross-check below). → independent confirmation of the unweighted sum.

### Formula check

- `HRD = LOH + TAI + LST` (unweighted) — matches Telli 2016 verbatim, Stewart 2022 (`gLOH+TAI+LST`),
  and scarHRD's `HRD-sum` column.
- `status = HrdHigh iff score ≥ 42` (boundary inclusive) — matches Telli 2016 ("≥42").

### Edge-case semantics check

- Boundary 42 → HRD-high; 41 → HRD-negative: directly from the inclusive "≥42" definition.
- All-zero (near-diploid) → score 0 → HRD-negative: sum of zeros, below cutoff. Well-defined.
- Components are non-negative event counts (LOH regions / telomeric AIs / LSTs per Abkevich/Birkbak/Popova);
  a negative count is invalid → defined as `ArgumentOutOfRangeException`. This is a sound, sourced convention
  (counts cannot be negative); not "implementation-defined".

### Independent cross-check (numbers)

scarHRD documentation worked examples (https://github.com/sztup/scarHRD), independent of this repo:

| HRD-LOH | Telomeric AI | LST | scarHRD HRD-sum | LOH+TAI+LST |
|---------|--------------|-----|-----------------|-------------|
| 1 | 2 | 0 | 3 | 1+2+0 = 3 ✓ |
| 25 | 35 | 33 | 93 | 25+35+33 = 93 ✓ |

Boundary triples (derived from Telli 2016's sum + ≥42 cutoff, the same arithmetic, independent of code):

| LOH | TAI | LST | Sum | Status (≥42 inclusive) |
|-----|-----|-----|-----|------------------------|
| 20 | 15 | 12 | 47 | HRD-high |
| 14 | 14 | 14 | 42 | HRD-high (boundary) |
| 14 | 13 | 14 | 41 | HRD-negative |
| 5 | 4 | 3 | 12 | HRD-negative |
| 0 | 0 | 0 | 0 | HRD-negative |

### Findings / divergences

None. The description (algorithm doc, TestSpec, Evidence) matches the external primary literature and an
independent reference implementation. Invariants INV-01..INV-04 are genuine (integer addition is
commutative; cutoff inclusive). The scope simplification (components supplied as inputs) is honestly
declared and does not affect correctness of the sum/threshold.

## Stage B — Implementation

### Code path reviewed

`src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.cs:1752–1862`:
- `HrdHighScoreThreshold = 42` (line 1759).
- `CalculateHRDScore` (1806–1825): validates each component `< 0` → `ArgumentOutOfRangeException`; returns
  `loh + tai + lst`. Exactly the unweighted sum.
- `ClassifyHRDStatus` (1836–1844): validates `score < 0` → throws; returns `score >= 42 ? HrdHigh : HrdNegative`.
  Inclusive `>=` matches the source ("≥42").
- `DetectHRD` (1856–1860): `score = CalculateHRDScore(...)`; `new HrdResult(components, score, ClassifyHRDStatus(score))`.
  Composition contract holds; validation delegated to `CalculateHRDScore`.

### Formula realised correctly?

Yes — exact integer sum, no weighting, no approximation; inclusive `>=` comparison against the sourced 42.

### Cross-verification table recomputed vs code

Traced by hand against the code and confirmed by the passing tests: 20+15+12=47→HrdHigh; 5+4+3=12→HrdNegative;
14+14+14=42→HrdHigh; 14+13+14=41→HrdNegative; 0+0+0=0→HrdNegative; ClassifyHRDStatus(42)=HrdHigh,
(41)=HrdNegative, (100)=HrdHigh, (0)=HrdNegative. All match the externally-sourced values.

### Variant/delegate consistency

`DetectHRD` reuses `CalculateHRDScore` and `ClassifyHRDStatus` directly (no duplicated logic), so INV-04
holds by construction. Negative-component rejection is shared via `CalculateHRDScore`.

### Test quality audit

File `tests/Seqeron/Seqeron.Genomics.Tests/OncologyAnalyzer_CalculateHRDScore_Tests.cs` (16 `[Test]` methods):

- **Sourced, exact expectations** — every assertion uses `Is.EqualTo` with an exact value (47, 12, 42, 41,
  0) or an exact enum (`HrdHigh`/`HrdNegative`), each annotated with the Telli-2016 source. No
  `Greater`/`AtLeast`/`Contains`/range/tolerance anywhere. A deliberately-wrong implementation (e.g. weighted
  sum, or `>42` strict cutoff) would fail M1/C1 and M3/M6 respectively, so the tests are not code-echoes.
- **All public surface exercised** — `CalculateHRDScore` (M1, M2, C1, S2), `ClassifyHRDStatus` (M3, M4, M5,
  S3, S4), `DetectHRD` (M6, M7, M8, M9, S1, negative), `HrdHighScoreThreshold` constant locked to 42.
- **All Stage-A branches/edges** — sum formula, inclusive boundary 42 vs 41, well-above (100), all-zero
  near-diploid (S1/S4), commutativity (C1), and the invalid negative-input error path on all three component
  positions (S2) plus score (S3) plus the end-to-end path (`DetectHRD_NegativeComponent_Throws`).
- **Honest green** — full unfiltered suite `Failed: 0, Passed: 6632` (one pre-existing benchmark skipped,
  unrelated). No skip/ignore/weakening was applied.

### Findings / defects

None. No code change and no test change were required: the implementation realises the validated
description exactly and the tests already lock the externally-sourced values across every branch and edge case.

## Verdict & follow-ups

- **Stage A: PASS** — sum formula and ≥42 inclusive cutoff confirmed against Telli 2016 (primary),
  Stewart 2022 (review), and the scarHRD reference implementation (worked examples 3 and 93 reproduce the
  literal sum).
- **Stage B: PASS** — code computes the exact unweighted sum and inclusive `>=42` classification; tests are
  exact, sourced, and cover every method and Stage-A edge/error case.
- **Test-quality gate: PASS** — exact sourced assertions, full surface and all branches covered, no
  green-washing, full unfiltered suite green (6632/0).
- **End-state: CLEAN.** No defect found; algorithm fully functional within its declared (composite-sum +
  threshold) scope. Component-derivation from raw segments remains, by design, in ONCO-LOH-001 / ONCO-CNA-001.
- **Follow-ups:** none for ONCO-HRD-001. (Cosmetic: the checklist method table lists non-existent
  `CalculateLOHScore`/`CalculateTAIScore`/`CalculateLSTScore` names — left as-is; out of scope for this unit.)
