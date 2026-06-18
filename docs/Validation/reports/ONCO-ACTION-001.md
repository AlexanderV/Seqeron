# Validation Report: ONCO-ACTION-001 — Clinical Actionability Assessment (OncoKB Therapeutic Levels of Evidence)

- **Validated:** 2026-06-16   **Area:** Oncology
- **Canonical method(s):** `OncologyAnalyzer.ClassifyActionabilityLevel`, `AssessActionability`, `GetTherapyRecommendations`, `CompareLevels`, `IsStandardCare` (`src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.cs:6918-7180`)
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS
- **End-state:** ✅ CLEAN

## Stage A — Description

### Sources opened & what they confirm (all retrieved this session)

1. **oncokb-annotator README** (WebFetch of `raw.githubusercontent.com/oncokb/oncokb-annotator/master/README.md`; corroborated by WebSearch) — verbatim column orderings:
   - `HIGHEST_LEVEL`: **"LEVEL_R1 > LEVEL_1 > LEVEL_2 > LEVEL_3A > LEVEL_3B > LEVEL_4 > LEVEL_R2"**
   - `HIGHEST_SENSITIVE_LEVEL`: **"LEVEL_1 > LEVEL_2 > LEVEL_3A > LEVEL_3B > LEVEL_4"**
   - `HIGHEST_RESISTANCE_LEVEL`: **"LEVEL_R1 > LEVEL_R2"**
   - Definition: `HIGHEST_LEVEL` = "The highest level of evidence for therapeutic implications."
2. **OncoKB Curation SOP v3** (downloaded 17 MB PDF via `curl`, text-extracted with `pypdf`) — verbatim:
   - Level 1 = "FDA-recognized biomarkers predictive of response to an FDA-approved drug in a specified indication"
   - Level 2 = "Standard care biomarkers recommended by the NCCN or other professional guidelines predictive of response to an FDA-approved drug in a specified indication"
   - Level R1 = "Standard care biomarkers predictive of resistance to an FDA-approved drug in this indication"
   - Grouping: "The highest levels of evidence, **Levels 1 and 2**, refer to the standard implications for sensitivity to an FDA-approved drug. Additionally, **Level R1** refers to the standard implications for resistance to an FDA-approved drug. **Levels 3A, 3B and 4** refer to the investigational implications for sensitivity… **Level R2** includes investigational implications for resistance…" ⇒ standard-care set = {1, 2, R1}.
   - Refinement note: "this system was refined to deprioritize the significance of standard care biomarkers when present in indications outside of the FDA-approved/NCCN listed indication" ⇒ Level 3A (in-indication clinical evidence) ranks above Level 3B (other-indication standard care).
   - Conflicting-data note: Levels 1/2/R1 "are categorized by their inclusion in either the FDA or NCCN guidelines, and therefore conflicting data is not relevant."
3. **OncoKB FAQ (`faq.oncokb.org/llms-full.txt`)** — verbatim Level 2 and Level 3B ("Standard care or investigational biomarker predictive of response to an FDA-approved or investigational drug in another indication").
4. **Chakravarty 2017 (PMC5586540)** — confirms the system: "Potential treatment implications are stratified by the level of evidence that a specific molecular alteration is predictive of drug response…" (full level table truncated in the PMC excerpt; level prose obtained from SOP v3 + FAQ + README above).

### Formula / ordering check
The combined precedence **R1 > 1 > 2 > 3A > 3B > 4 > R2** and the two sub-axis orders match the annotator README verbatim. The standard-care grouping {1, 2, R1} matches the SOP v3 prose verbatim. 3A > 3B ordering is confirmed by the SOP refinement note. No formula divergence.

### Edge-case semantics
- No leveled association ⇒ no highest level (annotator leaves HIGHEST_LEVEL empty). Modeled as `None`/`NotActionable` — documented assumption A1, faithful to the observable behaviour.
- A variant may carry both a sensitivity and a resistance association; the two axes are reported separately (SOP v3 + README) — the combined axis interleaves them (R1 above 1; R2 below 4).

### Independent cross-check (numbers, all from README order)
| Input set | Axis | Expected | Source |
|-----------|------|----------|--------|
| {1, R1} | combined | R1 | README HIGHEST_LEVEL (R1 > 1) |
| {4, R2} | combined | Level 4 | README (4 > R2) |
| {2, 3A} | sensitive | Level 2 | README HIGHEST_SENSITIVE_LEVEL |
| {3A, 3B, 4} | sensitive | Level 3A | README |
| {R1, R2} | resistance | R1 | README HIGHEST_RESISTANCE_LEVEL |
| {} | all | None / not actionable | README (empty HIGHEST_LEVEL) |

### Findings / divergences
None. The repo's TestSpec/Evidence orderings match the independently-retrieved README; the level definitions match SOP v3 + FAQ verbatim. Stage A **PASS**.

## Stage B — Implementation

### Code path reviewed
- `OncoKbLevel` enum (`:6926`) — ascending actionability `None < R2 < Level4 < Level3B < Level3A < Level2 < Level1 < R1`; the ordinal order **exactly** encodes the combined precedence R1 > 1 > 2 > 3A > 3B > 4 > R2 with None lowest.
- `CompareLevels` (`:7080`) = `((int)a).CompareTo((int)b)` — correct because the enum ordinals are the precedence.
- `HighestLevel` (`:7091`) — linear scan, `allowed` HashSet filter for the sensitivity/resistance sub-axes; `None` start value ⇒ empty/no-match ⇒ `None`.
- `SensitivityLevels`={1,2,3A,3B,4}, `ResistanceLevels`={R1,R2}, `StandardCareLevels`={1,2,R1} — all match the sources.
- `ClassifyActionabilityLevel` / `AssessActionability` / `GetTherapyRecommendations` / `IsStandardCare` — all delegate to the above; order preserved per input; recommendations sorted descending via `CompareLevels`.

### Formula realised correctly?
Yes. The maxima on each axis are the sourced README orders; the standard-care predicate is the sourced SOP grouping. Validation (null variants, null associations at construction) matches the contract.

### Cross-verification table recomputed vs code
All six cross-check rows above are exercised by tests and pass against the actual code (full suite green).

### Variant/delegate consistency
`ClassifyActionabilityLevel` and the combined field of `AssessActionability` both call `HighestLevel(..., allowed: null)` — consistent. `GetTherapyRecommendations` sorts via the same `CompareLevels`.

### Test quality audit (HARD gate)
- **Sourced expectations:** every expected value is keyed to the annotator README order / SOP grouping, not to code output. M1 checks every adjacent pair of the full chain; M2/M3 check the R1/R2 interleaving — a wrong implementation (e.g. ordinal flip, R1 treated as low) fails these.
- **No green-washing:** all assertions use `Is.EqualTo` exact values; no skips, no widened tolerances, no weakened assertions.
- **Coverage gaps found and fixed (test-only, additive, no code change):**
  1. Public `CompareLevels` was never directly exercised → added `CompareLevels_RealisesCombinedOrderAndNoneIsLowest` (full pairwise strict-monotonic chain over the 8-value ascending order + reflexivity; signs sourced to README + "None not actionable").
  2. `GetTherapyRecommendations` sort never exercised resistance/sensitivity interleaving (only sensitivity levels) → added `GetTherapyRecommendations_InterleavesResistanceAndSensitivity` ({R2,3A,R1,1} ⇒ R1,1,3A,R2).
  3. `AssessActionability` combined axis on a {sensitivity, resistance} mix not directly asserted → added `AssessActionability_Level4AndR2_CombinedIsLevel4`.
  - Fixture 16 → 19 tests.
- **Honest green:** full unfiltered suite **Failed: 0, Passed: 6694** (was 6691); `dotnet build` 0 errors, 0 new warnings (4 pre-existing NUnit2007 warnings only in unrelated `ApproximateMatcher_EditDistance_Tests.cs`).

### Findings / defects
No code or description defect. Three additive test-coverage improvements only. Stage B **PASS**.

## Verdict & follow-ups
- **Stage A: PASS · Stage B: PASS · End-state: ✅ CLEAN.**
- Test-quality gate: **PASS** (sourced expectations, no green-washing, coverage gaps closed, honest full-suite green).
- No follow-ups.
