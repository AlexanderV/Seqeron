# Significant Taxa Detection (Mann–Whitney U / Wilcoxon Rank-Sum)

| Field | Value |
|-------|-------|
| Algorithm Group | Metagenomics |
| Test Unit ID | META-TAXA-001 |
| Related Projects | Seqeron.Genomics.Metagenomics |
| Implementation Status | Production |
| Last Reviewed | 2026-06-13 |

## 1. Overview

Significant taxa detection identifies which microbial taxa differ in abundance between two sample groups (e.g., case vs control). For each taxon, the per-sample abundances of the two groups are compared with the Mann–Whitney U (Wilcoxon rank-sum) test, the standard non-parametric approach for microbiome differential abundance [3]. The test is non-parametric: it ranks the pooled observations and makes no normality assumption, which suits compositional abundance data. P-values are computed from the asymptotic normal approximation of the U statistic with midrank tie handling and an optional continuity correction [1][2].

## 2. Scientific / Formal Basis

### 2.1 Domain Context

In two-group microbiome/metagenomic comparisons, each sample yields an abundance for every observed taxon. To ask "is taxon X differentially abundant?", the abundances in group 1 are compared with those in group 2. Because abundances are non-normal and often heavily tied/zero-inflated, the non-parametric Wilcoxon rank-sum test is widely used to find statistically significant taxa [3].

### 2.2 Core Model

Pool the `n1` observations of group 1 and `n2` observations of group 2 (`n = n1 + n2`) and assign ranks; tied values receive the average (midrank) of the positions they occupy [1]. With `R1` the rank sum of group 1:

- `U1 = R1 − n1(n1+1)/2`, and `U2 = n1·n2 − U1`, so `U1 + U2 = n1·n2` [1][2].

Under H₀ the U statistic is asymptotically normal with [1]:

- mean `m_U = n1·n2 / 2`,
- standard deviation `σ_U = sqrt( n1·n2·(n1+n2+1)/12 − n1·n2·Σ(t_k³ − t_k) / (12·n·(n−1)) )`, where `t_k` is the number of tied observations at distinct value `k` (the right-hand subtraction is the tie correction; it is 0 when there are no ties) [1].

The z-score uses the larger U and an optional continuity correction of 0.5 [2]:

- `z = (|U − m_U| − cc) / σ_U`, `cc ∈ {0, 0.5}`,

and the two-tailed p-value is `p = 2·(1 − Φ(z))` where Φ is the standard normal CDF [1].

### 2.3 Modeling Assumptions

| ID | Assumption | Consequence if Violated |
|----|------------|--------------------------|
| ASM-01 | Observations within each group are independent draws from a continuous-ish distribution | Heavy ties (e.g., many zeros) make the asymptotic approximation conservative; exact methods preferred at very small n |
| ASM-02 | Sample sizes large enough for the normal approximation | At very small n the asymptotic p-value is approximate (an exact permutation p-value would differ) [2] |
| ASM-03 | A taxon absent in a sample has abundance 0 in that sample | Mislabeling absence as missing would change rank assignment |

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | `U1 + U2 = n1·n2` | Definition of U2 = n1·n2 − U1 [1][2] |
| INV-02 | `0 ≤ U ≤ n1·n2` | Follows from INV-01 with both U ≥ 0 [1] |
| INV-03 | p-value ∈ [0, 1] | `2·(1 − Φ(z))` for z ≥ 0, clamped |
| INV-04 | Significant ⇔ `PValue < pThreshold` | Method contract [3] |
| INV-05 | All-tied groups → p = 1, not significant | σ_U → 0 ⇒ no evidence against H₀ |
| INV-06 | Swapping group1/group2 leaves p unchanged | z uses `max(U1,U2)` / `|U − m_U|`, symmetric in the groups [1] |

### 2.5 Comparison with Related Methods

| Aspect | Mann–Whitney U (this) | Welch's t-test |
|--------|-----------------------|----------------|
| Distributional assumption | None (rank-based) | Approximate normality |
| Sensitivity to ties/zeros | Handled via tie correction | Affected by variance estimate |
| Statistic | Rank-sum U | Mean difference / pooled SE |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `group1` / `group2` | `IReadOnlyList<double>` | required | Observations per group (`MannWhitneyU`) | non-null, non-empty |
| `profiles` | `IReadOnlyList<IReadOnlyDictionary<string,double>>` | required | Per-sample taxon→abundance maps (`FindSignificantTaxa`) | non-null |
| `groups` | `IReadOnlyList<int>` | required | Group label (1 or 2) per profile, index-aligned | non-null; same length as profiles; each ∈ {1,2}; both groups non-empty |
| `pThreshold` | `double` | 0.05 | Significance cutoff | typically (0,1) |
| `useContinuityCorrection` | `bool` | true | Subtract 0.5 from `|U − m_U|` (SciPy default) [2] | — |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `MannWhitneyResult.U1`, `.U2` | `double` | U statistics for group1, group2 |
| `MannWhitneyResult.Z` | `double` | Normal-approx z-score `(|U − m_U| − cc)/σ_U` |
| `MannWhitneyResult.PValue` | `double` | Two-tailed asymptotic p-value, [0,1] |
| `SignificantTaxon` | record | `(Taxon, U, Z, PValue, Significant)` per taxon, ascending p-value |

### 3.3 Preconditions and Validation

`MannWhitneyU` throws `ArgumentNullException` for a null group and `ArgumentException` for an empty group. `FindSignificantTaxa` throws `ArgumentNullException` for null `profiles`/`groups`, and `ArgumentException` for length mismatch, a label other than 1/2, or a missing group; empty `profiles` returns an empty list. Abundances are numeric; a taxon absent from a profile is treated as abundance 0.

## 4. Algorithm

### 4.1 High-Level Steps

1. Pool both samples; sort by value.
2. Assign midranks (average rank within each tie block) and accumulate `Σ(t³−t)` over tie blocks.
3. Compute `R1`, then `U1 = R1 − n1(n1+1)/2` and `U2 = n1·n2 − U1`.
4. Compute `m_U`, tie-corrected variance, and z-score from `max(U1,U2)` with optional continuity correction.
5. Two-tailed p-value `2·(1 − Φ(z))`, clamped to [0,1]; degenerate σ=0 → p=1.
6. (`FindSignificantTaxa`) Run steps 1–5 per taxon over the two label-defined groups; flag `PValue < pThreshold`; sort ascending by p-value.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

- Continuity correction constant 0.5 [2]; variance divisor 12 from `σ_U² = n1·n2·(n1+n2+1)/12` [1].
- Φ is computed via the repository `StatisticsHelper.NormalCDF`, which uses Abramowitz & Stegun erf approximation 7.1.26 [4].

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `MannWhitneyU` | O(n log n) | O(n) | dominated by the pooled sort, n = n1+n2 |
| `FindSignificantTaxa` | O(t · s log s) | O(s) | t taxa, s = samples; per-taxon sort over the group vectors |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [MetagenomicsAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Metagenomics/MetagenomicsAnalyzer.cs)

- `MetagenomicsAnalyzer.MannWhitneyU(group1, group2, useContinuityCorrection)`: core rank-sum test, returns U1/U2/z/p.
- `MetagenomicsAnalyzer.FindSignificantTaxa(profiles, groups, pThreshold, useContinuityCorrection)`: per-taxon two-group test, returns ordered `SignificantTaxon` list.

### 5.2 Current Behavior

Ties use midranks and the tie-corrected variance. The z-score is taken from the larger of U1/U2 so the result is symmetric in the two groups. The continuity correction is on by default (SciPy parity). When every observation is tied (variance ≤ 0) the test reports z=0, p=1. No multiple-testing correction is applied across taxa (caller's responsibility). No substring/pattern search is involved, so the repository suffix tree is **not applicable** to this unit.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- `U1 = R1 − n1(n1+1)/2`, `U2 = n1·n2 − U1` [1][2].
- `m_U = n1·n2/2`, tie-corrected `σ_U` [1].
- z-score with continuity correction 0.5 and two-tailed `p = 2·(1 − Φ(z))` [1][2].

**Intentionally simplified:**

- Exact (permutation) p-value: only the asymptotic normal approximation is provided; **consequence:** at very small n the p-value differs slightly from an exact test [2].
- Φ via A&S 7.1.26 (|ε| ≤ 1.5×10⁻⁷); **consequence:** p-values match an exact normal CDF to ≈1×10⁻⁶ [4].

**Not implemented:**

- Multiple-testing (FDR) correction across taxa; **users should rely on:** applying a correction (e.g., Benjamini–Hochberg) to the returned p-values externally [3].

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Empty group | `ArgumentException` | Test undefined for n=0 |
| Null group/profiles/groups | `ArgumentNullException` | Input validation |
| All observations identical | z=0, p=1, not significant | σ_U → 0 (INV-05) |
| Taxon absent in some profiles | treated as abundance 0 | ASM-03 |
| Group label ∉ {1,2} | `ArgumentException` | Only two-group comparison supported |
| Empty profiles list | empty result | No taxa to test |

### 6.2 Limitations

Two-group comparison only; asymptotic (not exact) p-values; no multiplicity correction; rank-based test ignores the compositional structure of microbiome data (a known limitation of non-parametric per-taxon tests [3]).

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
var mw = MetagenomicsAnalyzer.MannWhitneyU(
    new double[] { 19, 22, 16, 29, 24 },
    new double[] { 20, 11, 17, 12 });
// mw.U1 == 17, mw.U2 == 3, mw.PValue ≈ 0.1113 (SciPy reference)
```

**Numerical walk-through:** x=[19,22,16,29,24] (n1=5), y=[20,11,17,12] (n2=4). Pooled ranks give R1, U1 = 17, U2 = 3. m_U = 10, σ_U = sqrt(200/12) = 4.0825. With continuity: z = (17−10−0.5)/4.0825 = 1.5922 → p ≈ 0.1113; without: z = 1.7146 → p ≈ 0.0864 [2].

### 7.3 Related Tests, Evidence, or Documents

- Tests: [MetagenomicsAnalyzer_FindSignificantTaxa_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/MetagenomicsAnalyzer_FindSignificantTaxa_Tests.cs) — covers `INV-01`..`INV-06`
- Evidence: [META-TAXA-001-Evidence.md](../../../docs/Evidence/META-TAXA-001-Evidence.md)

## 8. References

1. Mann, H.B.; Whitney, D.R. 1947. On a Test of Whether one of Two Random Variables is Stochastically Larger than the Other. *Annals of Mathematical Statistics* 18(1):50–60. https://doi.org/10.1214/aoms/1177730491 (formulas accessed via https://en.wikipedia.org/wiki/Mann%E2%80%93Whitney_U_test , 2026-06-13)
2. SciPy developers. scipy.stats.mannwhitneyu — SciPy Manual. https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.mannwhitneyu.html (accessed 2026-06-13)
3. Xia, Y.; Sun, J. 2017. Hypothesis testing and statistical analysis of microbiome. *Genes & Diseases* 4(3):138–148. https://pmc.ncbi.nlm.nih.gov/articles/PMC6128532/ (accessed 2026-06-13)
4. Abramowitz, M.; Stegun, I.A. 1964. Handbook of Mathematical Functions, formula 7.1.26 (erf approximation). Documented at https://www.johndcook.com/blog/python_erf/ (accessed 2026-06-13)
