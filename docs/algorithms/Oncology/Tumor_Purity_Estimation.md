# Tumor Purity Estimation

| Field | Value |
|-------|-------|
| Algorithm Group | Oncology |
| Test Unit ID | ONCO-PURITY-001 |
| Related Projects | Seqeron.Genomics.Oncology |
| Implementation Status | Production |
| Last Reviewed | 2026-06-14 |

## 1. Overview

Tumor purity ρ (also π) is the fraction of cells in a bulk sequencing sample that are tumour cells, the remainder being normal/stromal cells. These estimators recover ρ from the variant allele frequencies (VAFs) of clonal somatic mutations by inverting the closed-form expected-VAF relation that links VAF to purity, mutation multiplicity, and local copy number [1][2]. The computation is exact (deterministic, closed form) for a given (VAF, multiplicity, copy-number) state; the canonical special case — a clonal heterozygous somatic SNV at a copy-neutral diploid locus — gives ρ = 2·VAF [1]. It should be used when somatic SNV calls (and, for non-diploid loci, allele-specific copy-number state) are available; it is not a copy-ratio or B-allele-frequency segmentation method.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

A bulk tumour sample mixes tumour cells (fraction ρ) with normal diploid cells (fraction 1−ρ). For a somatic mutation present in m copies of the tumour genome on a segment of total copy number n_tot, the alternate-allele reads come only from tumour cells, while the total reads come from both populations. The expected VAF therefore depends jointly on ρ, m, and n_tot [1][2]. ABSOLUTE [3] and FACETS [4] use the same mixture to convert allelic fractions to per-cancer-cell allele counts and to estimate purity/ploidy.

### 2.2 Core Model

CNAqc gives the expected VAF of a clonal mutation (cancer cell fraction c = 1) as [1][2]:

> v = m·π / [ 2(1−π) + π·(n_A + n_B) ]

where m is the mutation multiplicity, π the tumour purity, and n_A + n_B = n_tot the tumour total copy number; the term 2(1−π) is the contribution of the healthy diploid normal cells and π·n_tot is the tumour contribution [2]. FACETS independently confirms the normal-diploid term: it mixes a normal (1,1) genotype with the aberrant (m,p) genotype at cellular fraction Φ via m* = mΦ + (1−Φ) [4].

Inverting for purity:

> π = 2·v / [ m + v·(2 − n_tot) ]

For a clonal **heterozygous** SNV at a **copy-neutral diploid** locus (m = 1, n_tot = 2) this collapses to v = π/2, i.e. **ρ = 2·v** [1]. CNAqc's worked example states a real purity of 60% corresponds to VAF 30% (and the 55–65% purity band to the 27.5–32.5% VAF band) [1].

### 2.3 Modeling Assumptions

| ID | Assumption | Consequence if Violated |
|----|------------|--------------------------|
| ASM-01 | Mutations are clonal (cancer cell fraction c = 1) | A subclonal mutation (c<1) has depressed VAF (v ∝ c); purity is underestimated [1] |
| ASM-02 | For the VAF-only estimator: m = 1 and n_tot = 2 (copy-neutral diploid heterozygous) | On amplified/LOH segments ρ = 2·VAF is wrong; use the allele-specific `EstimatePurity` with explicit (m, n_tot) [1] |
| ASM-03 | Normal cells are diploid (2 copies) | The 2(1−π) term and the closed form no longer hold |

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | 0 ≤ purity ≤ 1 | Purity is a fraction of cells [3]; inputs yielding ρ outside [0,1] are rejected |
| INV-02 | For m=1, n_tot=2: purity = 2·VAF exactly | v = π/2 special case of §2.2 [1] |
| INV-03 | Inversion recovers the purity that generated the VAF | π = 2v/[m+v(2−n_tot)] is the exact algebraic inverse of the §2.2 relation |
| INV-04 | VAF = 0 ⇒ purity = 0; estimate is monotone non-decreasing in VAF for fixed m, n_tot | ρ = 2v closed form [1] |

### 2.5 Comparison with Related Methods

| Aspect | This estimator (VAF/closed-form) | ABSOLUTE / FACETS |
|--------|----------------------------------|-------------------|
| Input | Clonal SNV VAFs (+ optional allele-specific CN) | Genome-wide segmented copy ratios + BAF + SNVs |
| Output | Point purity from the closed-form relation | Joint purity + ploidy + absolute CN by model fitting |
| Multiplicity | Supplied per variant | Inferred jointly |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| variants (`EstimatePurityFromVAF`) | `IEnumerable<VariantObservation>` | required | Clonal heterozygous diploid somatic SNVs | non-null, non-empty; valid read counts; each VAF ≤ 0.5 |
| vaf (`EstimatePurityFromVaf`) | `double` | required | Single clonal het diploid SNV VAF | ∈ [0, 0.5] |
| variants (`EstimatePurity`) | `IEnumerable<PurityVariant>` | required | Clonal SNVs with VAF, multiplicity m, total CN n_tot | non-null, non-empty; m ≥ 1; n_tot ≥ 1; VAF ∈ [0,1] |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| (return) | `double` | Estimated tumour purity ρ ∈ [0, 1]; median of per-variant estimates for the collection overloads |

### 3.3 Preconditions and Validation

Null collections throw `ArgumentNullException`; empty collections throw `ArgumentException` (purity undefined). A VAF outside [0,1] throws `ArgumentOutOfRangeException`; for the diploid model a VAF > 0.5 (implying ρ > 1) throws `ArgumentOutOfRangeException`. For the allele-specific overload, m < 1 or n_tot < 1, or any (VAF, m, n_tot) combination yielding ρ outside [0,1] (including a non-positive denominator), throws `ArgumentOutOfRangeException`. Read counts are validated via `CalculateVAF` (alt/total) as in ONCO-VAF-001.

## 4. Algorithm

### 4.1 High-Level Steps

1. For each variant compute its VAF (read-count overload) or read its supplied VAF.
2. Map the VAF to a per-variant purity: ρ = 2·v (diploid model) or ρ = 2v/[m + v(2 − n_tot)] (allele-specific).
3. Validate ρ ∈ [0, 1].
4. Aggregate per-variant purities by their median.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

- Normal total copy number fixed at 2 (autosomal diploid normal) [2][4].
- Heterozygous-diploid factor: ρ = 2·VAF (m=1, n_tot=2) [1].
- Median chosen over mean for robustness to subclonal/outlier VAFs (aggregation policy; does not alter the single-variant formula).

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| EstimatePurity / EstimatePurityFromVAF | O(n log n) | O(n) | n = #variants; dominated by the median sort. O(1) per variant. |
| EstimatePurityFromVaf | O(1) | O(1) | single closed-form evaluation |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [OncologyAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.cs)

- `OncologyAnalyzer.EstimatePurityFromVAF(IEnumerable<VariantObservation>)`: median of ρ = 2·VAF over clonal het diploid SNVs.
- `OncologyAnalyzer.EstimatePurityFromVaf(double)`: single-VAF closed form ρ = 2·VAF.
- `OncologyAnalyzer.EstimatePurity(IEnumerable<PurityVariant>)`: median of the allele-specific inversion ρ = 2v/[m+v(2−n_tot)].

### 5.2 Current Behavior

Collection overloads aggregate per-variant purities by median (lower-mid average for even counts). The VAF-only overload fixes the copy state to copy-neutral diploid heterozygous; non-diploid states are handled only by the allele-specific overload with explicit (m, n_tot). No search/matching is involved, so the repository suffix tree is not applicable (N/A).

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- ρ = 2·VAF for a clonal heterozygous copy-neutral diploid SNV (m=1, n_tot=2) [1].
- π = 2v/[m + v(2 − n_tot)], the exact inverse of v = m·π/[2(1−π)+π·n_tot] [1][2].
- Normal contribution fixed at 2 copies weighted (1−π) [2][4].

**Intentionally simplified:**

- Aggregation uses the median of per-variant point estimates; **consequence:** no model-fit confidence interval or joint ploidy estimate (unlike ABSOLUTE/FACETS) is produced.

**Not implemented:**

- Joint purity+ploidy+absolute-CN model fitting; **users should rely on:** ABSOLUTE [3] / FACETS [4] for genome-wide joint inference.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | VAF-only estimator fixes m=1, n_tot=2 | Assumption | Wrong on amplified/LOH loci | accepted | ASM-02; use `EstimatePurity` for other states |
| 2 | Median aggregation | Assumption | Robust central estimate, not a fitted value | accepted | does not change the per-variant formula |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| VAF = 0 | purity = 0 | ρ = 2·0 = 0 (INV-04) |
| VAF = 0.5 (diploid) | purity = 1.0 | ρ = 2·0.5 = 1 |
| VAF > 0.5 (diploid model) | ArgumentOutOfRangeException | ρ > 1 impossible (INV-01) |
| Empty collection | ArgumentException | purity undefined [1] |
| null collection | ArgumentNullException | guard |
| n_tot < 1 or m < 1 (allele-specific) | ArgumentOutOfRangeException | formula domain |

### 6.2 Limitations

Purity below ~0.1 approaches sequencing noise and is reported as a small value, not validated against a detection-limit model. The VAF-only estimator assumes clonality (c=1) and copy-neutral diploid heterozygosity; subclonal or amplified-segment variants must use the allele-specific overload with the correct (m, n_tot). No confidence interval, ploidy, or whole-genome-doubling handling is provided here (see ONCO-PLOIDY-001).

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
// Clonal heterozygous somatic SNV at a copy-neutral diploid locus, VAF 0.30.
double purity = OncologyAnalyzer.EstimatePurityFromVaf(0.30); // 0.60  (ρ = 2·VAF)

// Allele-specific: a 2:1 segment (n_tot=3), m=2, VAF 2/3 ⇒ fully pure tumour.
double rho = OncologyAnalyzer.EstimatePurity(new[]
{
    new OncologyAnalyzer.PurityVariant(2.0 / 3.0, Multiplicity: 2, TumorTotalCopyNumber: 3)
}); // 1.0
```

**Numerical walk-through:** CNAqc 2:1 example at π=1 [1]: m=1 ⇒ v = 1·1/[0 + 1·3] = 1/3; m=2 ⇒ v = 2/3. Inverting m=2, v=2/3, n_tot=3: π = 2·(2/3)/[2 + (2/3)(2−3)] = (4/3)/(4/3) = 1.0.

### 7.3 Related Tests, Evidence, or Documents

- Tests: [OncologyAnalyzer_EstimatePurity_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Oncology/OncologyAnalyzer_EstimatePurity_Tests.cs) — covers INV-01..INV-04
- Evidence: [ONCO-PURITY-001-Evidence.md](../../../docs/Evidence/ONCO-PURITY-001-Evidence.md)
- Related algorithms: [Variant_Allele_Frequency](../Oncology/Variant_Allele_Frequency.md)

## 8. References

1. Antonello A, Bergamin R, Calonaci N, Househam J, Milite S, Williams MJ, Anselmi F, d'Onofrio A, Sundaram V, Sosinsky A, Cross WCH, Caravagna G. 2024. Computational validation of clonal and subclonal copy number alterations from bulk tumor sequencing using CNAqc. Genome Biology 25(1):38. https://doi.org/10.1186/s13059-024-03170-5
2. CNAqc package vignette (Caravagna lab). Quality control of allele-specific copy numbers, mutations and tumour purity. https://caravagnalab.github.io/CNAqc/articles/CNAqc.html
3. Carter SL, Cibulskis K, Helman E, McKenna A, Shen H, Zack T, Laird PW, Onofrio RC, Winckler W, Weir BA, Beroukhim R, Pellman D, Levine DA, Lander ES, Meyerson M, Getz G. 2012. Absolute quantification of somatic DNA alterations in human cancer. Nature Biotechnology 30(5):413–421. https://doi.org/10.1038/nbt.2203
4. Shen R, Seshan VE. 2016. FACETS: allele-specific copy number and clonal heterogeneity analysis tool for high-throughput DNA sequencing. Nucleic Acids Research 44(16):e131. https://doi.org/10.1093/nar/gkw520
