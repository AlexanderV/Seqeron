# Validation Report: META-TAXA-001 — Significant Taxa Detection (Mann–Whitney U / Wilcoxon rank-sum)

- **Validated:** 2026-06-15   **Area:** Metagenomics
- **Canonical method(s):** `MetagenomicsAnalyzer.MannWhitneyU(group1, group2, useContinuityCorrection)`,
  `MetagenomicsAnalyzer.FindSignificantTaxa(profiles, groups, pThreshold, useContinuityCorrection)`
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS-WITH-NOTES (two test-quality defects found and fixed this session)

> Scope note: the per-session prompt's generic hints (lineage parsing, 7-rank hierarchy, LCA,
> Greengenes/SILVA `k__;p__;c__` strings, ETE3/QIIME) describe a *different* algorithm. The repo
> defines META-TAXA-001 consistently across the checklist, TestSpec, Evidence and source as the
> **Mann–Whitney U / Wilcoxon rank-sum significant-taxa test** under the normal approximation.
> It was validated as defined. (Taxonomy-string parsing exists separately in `ParseTaxonomyString`
> used by `ClassifyReads`, which is not this unit.)

## Stage A — Description

### Sources opened & what they confirm
- **scipy 1.13.1 `scipy.stats.mannwhitneyu`** (installed locally, executed this session) — the reference
  implementation. Reproduced every documented value (see cross-checks).
- **Wikipedia "Mann–Whitney U test"** (cited in Evidence, citing Mann & Whitney 1947) — formula set.
- **Mann & Whitney (1947)**, Ann. Math. Statist. 18(1):50–60 — primary definition of U and the
  normal-approximation mean/variance.
- **A&S 7.1.26 erf** (constants in `StatisticsHelper.Erf`) — confirms the |ε| ≤ 1.5×10⁻⁷ error bound,
  hence p-values match the exact normal CDF only to ≈1×10⁻⁶ (the test tolerance `ErfPTolerance`).
- **Xia & Sun (2017)** PMC6128532 — domain justification (rank-sum used per taxon for differential abundance).

### Formula check
All match the cited sources exactly:
- `U1 = R1 − n1(n1+1)/2`; `U2 = n1·n2 − U1`; `U1 + U2 = n1·n2` (Mann & Whitney 1947) — code lines 1384–1386.
- `m_U = n1·n2/2`; `σ_U = sqrt(n1·n2·(n1+n2+1)/12)` — code lines 1388, 1390–1392.
- Tie correction `σ² = n1·n2/12 · [(n+1) − Σ(t³−t)/(n(n−1))]` with midranks — code lines 1364–1377, 1390–1392.
  This is algebraically identical to the Evidence/Wikipedia form
  `σ² = n1·n2·(n+1)/12 − n1·n2·Σ(t³−t)/(12·n·(n−1))` (verified by factoring n1·n2/12).
- Continuity correction subtracts 0.5 from |U − m_U| (SciPy default on) — code lines 1407–1408.
- Two-tailed p = `2·(1 − Φ(|z|))` clamped to [0,1] — code line 1411.

### Edge-case semantics check
- **All-tied / zero-variance** → σ → 0, z defined as 0, p = 1 (no evidence against H₀) — code lines 1396–1401.
  Sourced (degenerate σ→0). Confirmed scipy returns p≈1 in this case as well.
- **Empty / null group** → throws (ArgumentException / ArgumentNullException) — code lines 1346–1349. Correct
  (test undefined for n=0).
- **Absent taxon in a profile** → filled with abundance 0 (ASM-03, standard for abundance tables) — code line 1462–1463.

### Independent cross-check (numbers — all recomputed this session with scipy 1.13.1)
| Quantity | Sourced value | Method |
|----------|---------------|--------|
| U1 (x=[19,22,16,29,24], y=[20,11,17,12]) | 17.0 | scipy.statistic |
| U2 | 3.0 | n1·n2 − U1 |
| σ_U | 4.08248290463863 | sqrt(200/12) |
| z (no cc) | 1.7146428199482247 | (17−10)/σ |
| z (cc) | 1.5921683328090657 | (17−10−0.5)/σ |
| p (no cc) | 0.08641073297370001 | scipy asymptotic two-sided |
| p (cc) | 0.11134688653314039 | scipy asymptotic two-sided |
| Tortoise/hare U_T, U_H | 11, 25 | scipy + Mann&Whitney 1947 |
| M7 ties σ ([1,1,2] vs [1,2,2]) | 2.012461179749811 | sqrt(4.05), hand-derived tie correction |

All match the TestSpec/Evidence claims (the spec's p(cc)=…41 vs scipy's …39 differ in the last digit
only, well within the 1×10⁻⁶ erf tolerance).

### Findings / divergences
None. The description is biologically and mathematically correct and matches the reference implementation.

## Stage B — Implementation

### Code path reviewed
`src/Seqeron/Algorithms/Seqeron.Genomics.Metagenomics/MetagenomicsAnalyzer.cs:1293–1477`
(region "Significant Taxa Detection (Mann–Whitney U)"), plus `StatisticsHelper.NormalCDF`/`Erf`
(`Infrastructure/StatisticsHelper.cs:13–33`, A&S 7.1.26 with the exact sourced constants).

### Formula realised correctly?
Yes. The implementation computes the validated U, tie-corrected σ, continuity correction, and
two-tailed asymptotic p exactly. Verified by tracing and by independent scipy reproduction:

| Case | Code output (this session) | scipy 1.13.1 | Match |
|------|---------------------------|--------------|-------|
| SciPy ex. U1/U2 | 17 / 3 | 17 / 3 | exact |
| SciPy ex. z (no cc) | 1.7146428199482247 | 1.7146428199482247 | exact |
| SciPy ex. p (no cc) | — (asserted 0.0864107329737 ±1e-6) | 0.08641073297370001 | ✓ |
| SciPy ex. p (cc) | — (asserted 0.11134688653314041 ±1e-6) | 0.11134688653314039 | ✓ |
| Tortoise/hare U | 11 / 25 | 11 / 25 | exact |
| M7 ties z (no cc) | 0.7453559924999298 | matches hand σ=2.0124611797498108 | ✓ |
| M9 TaxonA p (cc) | 0.08085551657571322 | 0.08085559837005224 | ✓ (≤1e-6) |
| M9 TaxonB p (cc) | 0.6192567546069347 | 0.6192567541768621 | ✓ (≤1e-6) |
| S3 TaxonZ p (cc) | 0.22067149635775185 | 0.22067136191984682 | ✓ (≤1e-6) |
| M9b 4v4 separated p (cc) | (asserted 0.030382821976577504 ±1e-6) | 0.030382821976577504 | ✓ |

### Variant/delegate consistency
`FindSignificantTaxa` delegates per-taxon to `MannWhitneyU` (with absent→0 fill), takes the larger U,
and reuses the same p; deterministic ascending-p ordering. Consistent with the canonical core.

### Test quality audit (HARD gate)
Two Stage-B **test-quality defects** found (no code defect):

1. **M9 (`FindSignificantTaxa_SeparatedAndOverlapping_FlagsCorrectly`) was code-echoing.** It asserted
   `a.Significant == (a.PValue < 0.05)` and `b.Significant == (b.PValue < 0.05)` — tautologies that
   restate the method's own contract against itself and would pass against *any* implementation
   (e.g. one that computes the wrong p, or always returns Significant=false). The only sourced facts
   it carried were the relative ordering. **Fixed:** rewrote to assert the externally-sourced exact
   p-values (TaxonA 0.08085559837005224, TaxonB 0.6192567541768621 from scipy, within erf tolerance)
   and the correct significance flags (both **False** at p=0.05 — a 3-vs-3 comparison cannot reach 0.05).
   Also **added M9b** (`FindSignificantTaxa_FullySeparated_FlagsSignificant`, 4-vs-4) to genuinely
   exercise the Significant=**true** branch with the sourced p=0.030382821976577504; the prior test
   never reached a significant outcome despite its name.

2. **S3 (`FindSignificantTaxa_AbsentTaxon_TreatedAsZeroAbundance`) was green-washed.** It asserted only
   `p > 0 and ≤ 1` — a wide range where the exact value is known and where the absent→0 fill behaviour
   it claims to test is not actually pinned (a wrong fill would still give a p in (0,1]). **Fixed:** now
   asserts the exact sourced value p = 0.22067136191984682 (scipy `mannwhitneyu([0,0],[50,60])`), which
   is reproducible only if the absent taxon is filled with 0 — genuinely locking ASM-03.

All other tests (M1–M8, S1, S2, C1, null/empty/mismatch/single-group/invalid-label) assert exact
sourced values or true invariants and exercise every public method, both branches of the continuity
flag, the tie path, the degenerate σ→0 path, and all documented error cases.

**Honest green:** full unfiltered `dotnet test` = **6558 passed, 0 failed, 0 skipped-relevant**
(MFE benchmark intentionally skipped); `dotnet build` 0 errors; edited test file warning-free.

### Findings / defects
- Defect (test-quality, fixed): M9 code-echo; S3 range-instead-of-exact. Both rewritten to sourced
  exact values; M9b added to cover the significant branch. No implementation defect.

## Verdict & follow-ups
- **Stage A: PASS.** Description matches Mann & Whitney (1947) and scipy verbatim.
- **Stage B: PASS-WITH-NOTES.** Implementation is correct against the reference; two test-quality
  defects were found and completely fixed in-session.
- **End-state: CLEAN.** No remaining gaps; suite green at 6558.
</content>
</invoke>
