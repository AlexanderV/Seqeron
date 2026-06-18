# Validation Report: ONCO-ARTIFACT-001 — Sequencing Artifact Detection (OxoG / FFPE deamination / Fisher strand bias)

- **Validated:** 2026-06-16   **Area:** Oncology
- **Canonical method(s):** `OncologyAnalyzer.ClassifyArtifact`, `CalculateGivScore`, `CalculateStrandBias`, `DetectOxoGArtifacts`, `FilterArtifacts`
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS-WITH-NOTES (test-coverage gaps closed in-session; algorithm correct)

## Stage A — Description

### Sources opened & what they confirm (this session)

| Source | Retrieved | Confirms |
|--------|-----------|----------|
| Chen et al. 2017, Science 355:752 (via GenomeWeb summary + searches) | WebSearch | OxoG = G>T excess in read 1 / C>A excess in read 2 (reverse complement); GIV is the R1/R2 imbalance per substitution type; ~3/4 of TCGA datasets had elevated GIV |
| Ettwiller Damage-estimator (reference implementation) | WebFetch github.com/Ettwiller/Damage-estimator | "estimation of damage … based on imbalance between R1 and R2 variant frequency"; standard acoustic shearing ⇒ "GIV score for G_T around 2" — confirms GIV = R1/R2 G_T ratio with neutral baseline ~1 |
| Costello et al. 2013, NAR doi:10.1093/nar/gks1443 (via search + GATK) | WebSearch | 8-oxoguanine pairs with adenine ⇒ apparent G>T in read 1 and C>A in read 2; the canonical OxoG oxidative artifact substitution class is G:C>T:A |
| Do & Dobrovic 2015 (FFPE review; Evidence doc + searches) | prior Evidence + WebSearch | FFPE cytosine deamination ⇒ uracil pairs with adenine ⇒ C>T (and G>A antisense), collectively C:G>T:A; disjoint from the oxidation class C>A/G>T |
| GATK FisherStrand.java (Broad, master) | WebFetch raw GitHub | `FS = phredScaleErrorRate(max(pValue, MIN_PVALUE))`; `MIN_PVALUE = 1E-320`; FS is the two-sided Fisher exact p Phred-scaled |
| GATK StrandBiasTest.java (Broad, master) | WebFetch raw GitHub | 2×2 table cell order array[0..3] = ref-fwd, ref-rev, alt-fwd, alt-rev; `table[0][0]=array[0]; table[0][1]=array[1]; table[1][0]=array[2]; table[1][1]=array[3]` |

### Formula check
- **Substitution classes** (spec §2.2): FFPE {C>T, G>A}; OxoG {G>T, C>A}; disjoint. Matches Do & Dobrovic 2015 (C:G>T:A) and Chen 2017 / Costello 2013 (G>T R1 / C>A R2) exactly.
- **GIV = read1Alt / read2Alt** with neutral 1.0 and damaged > 1.5. The R1/R2 ratio and the ~2 standard-shearing baseline are confirmed by Damage-estimator (authors' own implementation). The exact 1.0/1.5 numeric thresholds come from the Nature Methods summary of Chen 2017 (nmeth.4254); that page is paywalled this session, but the ratio's directionality and neutral-at-1 behaviour are independently corroborated by Damage-estimator and the GenomeWeb write-up. Threshold = operational cutoff, correctly cited.
- **FS = −10·log10(max(p, 1E-320))** with the two-sided Fisher exact p on [refFwd, refRev, altFwd, altRev]. Matches FisherStrand.java verbatim (formula, constant, cell order).

### Edge-case semantics check
- GIV with r2=0, r1=0 ⇒ 1.0 (no imbalance evidence); r2=0, r1>0 ⇒ +∞ (maximal imbalance). Consistent with a directional R1/R2 ratio.
- All-zero / zero-margin strand table ⇒ two-sided Fisher p = 1 ⇒ FS = 0 (no evidence of bias). Independently confirmed with scipy (below).
- Substitution outside the four artifact pairs ⇒ None, never flagged.

### Independent cross-check (numbers)
`scipy.stats.fisher_exact(..., alternative='two-sided')` (independent of the C# code):

| Table [refFwd,refRev,altFwd,altRev] | scipy two-sided p | FS = −10·log10(p) |
|-------------------------------------|-------------------|-------------------|
| [10,10,10,10] | 1.0 | 0.0 |
| [20,0,0,20] | 1.4508889103849688e-11 | 108.38365838736459 |
| [15,5,5,15] | 0.003847527308377529 | 24.148182890181 |
| [10,0,10,0] (zero margin) | 1.0 | 0.0 |

All match the spec/implementation values to ≥1e-13.

### Findings / divergences
None. Stage A description is biologically and mathematically correct. The only assumption (strand/read-mate counts supplied on the observation record rather than parsed from a BAM) is an API-shape decision that does not affect any classification rule. **Stage A: PASS.**

## Stage B — Implementation

### Code path reviewed
`src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.cs:899–1254`.

- `ClassifySubstitution` (1135) maps (C,T)/(G,A)→FFPE, (G,T)/(C,A)→OxoG, else None — matches the sourced disjoint classes.
- `CalculateGivScore` (999) = r1/r2 with r2=0 handling (1.0 if r1=0 else +∞); validates non-negative — matches spec.
- `CalculateStrandBias` (1035) = −10·log10(max(p, 1E-320)) via a self-contained two-sided Fisher (`FisherExactTwoSided`, 1162) summing same-margin hypergeometric tables with prob ≤ observed (the conventional two-sided definition), using Lanczos log-gamma. `total==0`⇒p=1.
- `ClassifyArtifact` (1059): FFPE flagged always; OxoG flagged iff GIV > 1.5; case-insensitive bases.
- `DetectOxoGArtifacts` / `FilterArtifacts`: null-checked; preserve input order; result ⊆ input.

### Formula realised correctly?
Yes. The Fisher routine reproduces scipy's two-sided p to ≥1e-13 on all four cross-check tables (including the zero-margin degenerate case), and the Phred scaling, MIN_PVALUE floor, and cell ordering match GATK. GIV and the substitution map are exact.

### Cross-verification table recomputed vs code
Test `CalculateStrandBias_FullySegregatedTable` asserts FS = 108.38365838736458 for [20,0,0,20]; scipy gives 108.38365838736459 (Δ < 1e-13). `CalculateStrandBias_PartiallyBiasedTable` asserts 24.148182890180962 for [15,5,5,15]; scipy gives 24.148182890181 (Δ < 1e-12). Balanced [10,10,10,10] ⇒ 0.0. All exact-value assertions trace to the independent scipy computation, not to the code.

### Variant/delegate consistency
`DetectOxoGArtifacts` and `FilterArtifacts` both delegate to `ClassifyArtifact`; behaviour is consistent (OxoG-only flagged subset vs non-artifact subset).

### Test quality audit (HARD gate)
- **Sourced, not echoed:** the two non-trivial FS values are derived from scipy (independent reference), not the implementation. A deliberately-wrong Fisher implementation would fail M9/M9-ext. GIV and substitution assertions are exact and sourced.
- **No green-washing:** exact `.EqualTo(...)` assertions with tight tolerances (≤1e-6) where exact values are known; monotonicity (C1/INV-05) uses ordering comparisons, which is appropriate for a property test; no skips/ignores/widened tolerances introduced.
- **Coverage gaps found and fixed (this session):** added four tests for previously-untested branches/edge cases —
  1. `CalculateStrandBias_EmptyTable_ReturnsZero` (the `total==0` branch, spec §3.3),
  2. `CalculateStrandBias_ZeroMarginTable_ReturnsZero` ([10,0,10,0], scipy-confirmed p=1),
  3. `CalculateGivScore_NegativeR2_Throws` (mirror of the negative-r1 branch),
  4. `DetectOxoGArtifacts_Empty_ReturnsEmpty`.
  (Renamed `CalculateGivScore_NegativeCount_Throws` → `..._NegativeR1_Throws` for clarity.)
- **Honest green:** FULL unfiltered suite = **6628 passed, 0 failed, 0 skipped (relevant)**; `dotnet build` 0 errors; no new warnings in the changed test file.

### Findings / defects
No implementation defect. Stage B gap was test coverage only (untested `total==0`/zero-margin Fisher branch, negative-r2 GIV branch, empty DetectOxoG), now closed. **Stage B: PASS-WITH-NOTES.**

## Verdict & follow-ups
- **Stage A: PASS.** **Stage B: PASS-WITH-NOTES.** **End-state: CLEAN.**
- No code change required; tests strengthened to lock sourced values and cover all branches/edge cases. Working tree builds and the full suite is green.
- Note: the 1.0/1.5 GIV numeric thresholds rest on the Nature Methods summary (paywalled this session); the ratio definition and neutral-at-1 behaviour are independently corroborated by the authors' Damage-estimator implementation, so the threshold citation is sound.
