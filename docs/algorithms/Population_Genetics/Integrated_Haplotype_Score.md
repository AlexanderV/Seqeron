# Integrated Haplotype Score (iHS) — Selection Signature Detection

| Field | Value |
|-------|-------|
| Algorithm Group | Population Genetics |
| Test Unit ID | POP-SELECT-001 |
| Related Projects | Seqeron.Genomics.Population |
| Implementation Status | Simplified |
| Last Reviewed | 2026-06-13 |

## 1. Overview

The integrated Haplotype Score (iHS) detects recent positive selection from phased haplotype
data by contrasting the rate of haplotype-homozygosity decay around the ancestral versus the
derived allele of a focal SNP [1]. A site under a recent selective sweep carries an unusually
long, conserved haplotype around the favored allele; iHS quantifies that asymmetry. The
computation is deterministic and exact for the formal definition: it builds the Extended
Haplotype Homozygosity (EHH) decay curves [2], integrates each into iHH by the trapezoidal rule,
and reports the unstandardized score ln(iHH_A/iHH_D) [1]. A separate frequency-binned
standardization converts these into approximately standard-normal, cross-comparable scores [1],
and a genome-wide scan summarizes windows by the proportion of SNPs with |iHS| > 2 [1].

## 2. Scientific / Formal Basis

### 2.1 Domain Context

Under a recent positive sweep, the selected allele rises in frequency faster than recombination
can break down the haplotype it sits on, so haplotype homozygosity extends further than expected
under neutrality [1][2]. Comparing the favored (typically derived) allele against the alternative
(typically ancestral) allele controls for local recombination rate and demography, because both
alleles share the same genomic background [1].

### 2.2 Core Model

**Extended Haplotype Homozygosity (EHH).** For chromosomes carrying a core allele *c*, extended
to marker *x_i*, EHH is the probability that two randomly chosen core-carrying chromosomes are
identical over the interval [2][3]:

EHH_c(x_i) = Σ_{h ∈ H_c(x_i)} C(n_h, 2) / C(n_c, 2)

where n_h is the count of each distinct extended haplotype *h* and n_c the number of
core-carrying chromosomes; C(n,2) = n(n−1)/2 [3]. The rehh form
EHH = (1/(n_a(n_a−1))) Σ_k n_k(n_k−1) is algebraically identical [4].

**Integrated EHH (iHH).** iHH is the area under the EHH curve, joined by straight lines
(trapezoidal rule) against marker positions, summed in both directions away from the core SNP,
truncated where EHH first drops below 0.05 [1][3][4]:

iHH_c = Σ_i ½(EHH_c(x_{i−1}) + EHH_c(x_i))·g(x_{i−1}, x_i)

over downstream and upstream markers, where g is the distance between adjacent markers [3].

**Unstandardized iHS** [1]:

unstandardized iHS = ln(iHH_A / iHH_D)

(ancestral numerator, derived denominator). selscan adopts the reciprocal
ln(iHH_1/iHH_0) = ln(iHH_D/iHH_A) and explicitly notes the sign difference from Voight [3].

**Standardized iHS** [1]:

iHS = (ln(iHH_A/iHH_D) − E_p[ln(iHH_A/iHH_D)]) / SD_p[ln(iHH_A/iHH_D)]

where the expectation and standard deviation are taken over SNPs whose derived allele frequency
*p* falls in the same bin; the result is approximately standard normal [1].

### 2.3 Modeling Assumptions

| ID | Assumption | Consequence if Violated |
|----|------------|--------------------------|
| ASM-01 | Haplotypes are correctly phased and ancestral state is known | Mis-phasing or mis-polarization flips/biases the iHS sign and magnitude [1] |
| ASM-02 | Standardization bins contain enough SNPs at each derived-allele frequency | Sparse bins give unstable E/SD estimates and noisy standardized scores [1] |

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | EHH ∈ [0,1]; EHH = 1 for a single chromosome, 0 for all-distinct haplotypes | ratio of pair counts C(n_h,2)/C(n_c,2) [3] |
| INV-02 | iHH ≥ 0 | sum of non-negative trapezoid areas [1][3] |
| INV-03 | Balanced EHH decay ⇒ unstandardized iHS = 0 | iHH_A/iHH_D ≈ 1 ⇒ ln(1) = 0 [1] |
| INV-04 | ln(iHH_A/iHH_D) = −ln(iHH_D/iHH_A) | logarithm of reciprocal; matches Voight vs selscan sign note [1][3] |
| INV-05 | Standardized scores within a bin have mean 0 (sd 1 when >1 element) | subtract bin mean, divide by bin sd [1] |
| INV-06 | Scan ProportionExtreme = ExtremeCount / SnpCount ∈ [0,1] | count of |iHS|>2 over window size [1] |

### 2.5 Comparison with Related Methods

| Aspect | iHS | Tajima's D / Fay & Wu's H |
|--------|-----|--------------------------|
| Signal | Haplotype-length asymmetry around an allele [1] | Site-frequency-spectrum distortion [1] |
| Within vs across populations | Within one population, derived vs ancestral [1] | Within one population, allele frequencies |
| Optimal sweep stage | Incomplete/ongoing sweep (intermediate frequency) [1] | Completed sweeps / frequency skew |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `extendedHaplotypes` | `IReadOnlyList<string>` | required | Extended haplotype per core-carrying chromosome (`CalculateEhh`) | non-null |
| `haplotypes` | `IReadOnlyList<string>` | required | Phased haplotypes, one allele char per marker (`CalculateIHS`) | non-null, equal length = positions.Count, core allele ∈ {'0','1'} |
| `positions` | `IReadOnlyList<int>` | required | Marker positions (genomic or genetic units) | non-null, indexed parallel to haplotype chars |
| `coreIndex` | `int` | required | Index of focal SNP | 0 ≤ coreIndex < positions.Count |
| `scores` | `IReadOnlyList<(double Unstandardized, double DerivedAlleleFrequency)>` | required | iHS values with derived freq (`StandardizeIHS`) | non-null |
| `binCount` | `int` | 20 | Equal-width derived-frequency bins over (0,1) | ≥ 1; default width 0.05 [4] |
| `standardizedScores` | `IReadOnlyList<double>` | required | Per-SNP standardized iHS in genomic order (`ScanForSelection`) | non-null |
| `windowSize` | `int` | 50 | Consecutive SNPs per window | ≥ 1; Voight used 50 [1] |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `CalculateEhh` → | `double` | EHH ∈ [0,1] |
| `IhsResult.UnstandardizedIHS` | `double` | ln(iHH_A/iHH_D) (Voight sign) [1] |
| `IhsResult.IhhAncestral` / `IhhDerived` | `double` | Integrated EHH areas |
| `IhsResult.DerivedAlleleFrequency` | `double` | Frequency of derived core allele |
| `StandardizeIHS` → | `IReadOnlyList<double>` | Frequency-binned z-standardized iHS, input order |
| `SelectionScanWindow` | record | WindowIndex, SnpCount, ExtremeCount, ProportionExtreme (|iHS|>2) |

### 3.3 Preconditions and Validation

Indexing is 0-based. Core alleles use the polarized encoding `'0'` = ancestral, `'1'` = derived;
any other core character throws `ArgumentException`. Null inputs throw `ArgumentNullException`;
`coreIndex` out of range throws `ArgumentOutOfRangeException`; inconsistent haplotype lengths or a
monomorphic focal SNP (only one allele present) throw `ArgumentException`. `binCount` and
`windowSize` below 1 throw `ArgumentOutOfRangeException`. An empty EHH sample returns 0; a
single-chromosome sample returns EHH = 1.

## 4. Algorithm

### 4.1 High-Level Steps

1. Partition chromosomes at the focal SNP into ancestral ('0') and derived ('1') sets.
2. For each set, walk outward in both directions; at each marker compute EHH over the extended
   substring via `CalculateEhh`.
3. Accumulate the trapezoidal area between consecutive markers; stop a direction once EHH < 0.05.
4. Sum both directions to obtain iHH_A and iHH_D; return ln(iHH_A/iHH_D).
5. (Optional) Standardize scores within derived-allele-frequency bins.
6. (Optional) Scan SNP-ordered standardized scores in fixed windows, reporting the proportion
   with |iHS| > 2.

### 4.2 Decision Rules, Scoring, Reference Tables

| Constant | Value | Source |
|----------|-------|--------|
| EHH integration cutoff | 0.05 | Voight (2006) §M&M; rehh `limehh` default [1][4] |
| Extreme iHS threshold | 2.0 (|iHS| > 2) | Voight (2006) §M&M [1] |
| Default window size | 50 SNPs | Voight (2006) §M&M [1] |
| Default bin width | 0.05 (20 bins) | rehh `freqbin` convention [4] |

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `CalculateEhh` | O(n·L) | O(n) | n = chromosomes, L = window length (hashing substrings) |
| `CalculateIHS` | O(n·h) | O(n·h) | n = chromosomes, h = markers; EHH recomputed per marker outward |
| `StandardizeIHS` | O(N) | O(N) | N = SNPs |
| `ScanForSelection` | O(N) | O(1) per window | single pass over standardized scores |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [PopulationGeneticsAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Population/PopulationGeneticsAnalyzer.cs)

- `PopulationGeneticsAnalyzer.CalculateEhh(IReadOnlyList<string>)`: EHH for one core-allele subset.
- `PopulationGeneticsAnalyzer.CalculateIHS(IReadOnlyList<string>, IReadOnlyList<int>, int)`: full iHS pipeline → `IhsResult` (unstandardized ln(iHH_A/iHH_D)).
- `PopulationGeneticsAnalyzer.StandardizeIHS(IReadOnlyList<(double, double)>, int)`: frequency-binned z-standardization.
- `PopulationGeneticsAnalyzer.ScanForSelection(IReadOnlyList<double>, int)`: windowed proportion of |iHS|>2.

### 5.2 Current Behavior

EHH is computed from the literal substrings of the supplied haplotype matrix, so any alphabet is
accepted for flanking markers (only the core character is constrained to `'0'`/`'1'`).
Integration recomputes EHH per outward marker rather than incrementally; this is O(n·h) and exact
but not optimized for genome-scale inputs. The implementation does **not** use the repository
suffix tree: this is not a substring-search problem (EHH counts exact whole-window haplotype
multiplicities by hashing, not occurrences of a query pattern in a text), so the suffix tree does
not apply. A separate, pre-existing overload `CalculateIHS(ehh0, ehh1, positions)` (EHH curves
supplied directly) and the region-threshold `ScanForSelection(...)` remain for the MCP layer and
are out of scope for this unit.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- EHH_c = Σ C(n_h,2)/C(n_c,2) [2][3].
- iHH via trapezoidal rule, both directions, truncated at EHH < 0.05 [1][3][4].
- Unstandardized iHS = ln(iHH_A/iHH_D), Voight sign convention [1].
- Frequency-binned standardization (x − E_p)/SD_p [1].
- Genome-wide window summary by proportion of |iHS| > 2 [1].

**Intentionally simplified:**

- Standardization SD uses the sample (N−1) estimator; **consequence:** standardized magnitudes
  differ slightly from a population-SD implementation, with no effect on sign, ordering, or the
  unstandardized score.
- Frequency-bin width defaults to 0.05 (parameterized); **consequence:** different binning
  changes which empirical neighbours define E/SD but not the unstandardized iHS.

**Not implemented:**

- Genetic-map distance scaling and the ad-hoc 20–200 kb gap correction of Voight §M&M;
  **users should rely on:** supplying positions already in the desired distance units (e.g. 4Nr).
- Cross-population XPEHH; **users should rely on:** no current alternative in this class.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | SD estimator = sample (N−1) | Assumption | Magnitude only | accepted | ASM-02; Voight does not specify N vs N−1 |
| 2 | Default bin width 0.05 | Assumption | Choice of empirical neighbours | accepted | matches rehh `freqbin` [4] |
| 3 | Sign = Voight ln(iHH_A/iHH_D) | Deviation (vs selscan) | Sign flip vs selscan | accepted | selscan reciprocal documented [3] |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Single chromosome carries core | EHH = 1 | one chromosome is trivially homozygous [3] |
| Empty EHH sample | EHH = 0 | no pairs; boundary of C(n,2) formula |
| Balanced decay | unstandardized iHS = 0 | INV-03 [1] |
| Monomorphic focal SNP | `ArgumentException` | iHS defined only for polymorphic SNPs [1] |
| Core allele ∉ {'0','1'} | `ArgumentException` | EHH defined on polarized biallelic core [1] |
| Singleton frequency bin | standardized = 0 | SD undefined for one element (ASM-02) |

### 6.2 Limitations

Requires correctly phased, ancestral-state-polarized haplotypes; no phasing or polarization is
performed. Standardization needs a genome-scale set of SNPs to estimate stable per-bin E/SD;
small inputs give noisy standardized scores. Does not model recombination-rate-aware genetic
distances, large-gap masking, or XPEHH. O(n·h) integration is exact but not tuned for very large
panels.

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
var haplotypes = new[] { "AA1GG", "AA1GG", "AA1GG", "TC0TC", "GA0AG", "CT0CA" };
var positions = new[] { 0, 10, 20, 30, 40 };
var ihs = PopulationGeneticsAnalyzer.CalculateIHS(haplotypes, positions, coreIndex: 2);
// ihs.IhhAncestral == 10.0, ihs.IhhDerived == 40.0
// ihs.UnstandardizedIHS == ln(10/40) == -1.386294361
```

**Numerical walk-through:** The 3 derived chromosomes are identical, so EHH_D = 1 at every flank
marker; with positions spaced by 10, each direction contributes ½(1+1)·10 + ½(1+1)·10 = 20, so
iHH_D = 40. The 3 ancestral chromosomes are all distinct, so EHH_A = 0 at the first flank marker;
each direction contributes ½(1+0)·10 = 5 then truncates (EHH < 0.05), so iHH_A = 10. Hence
unstandardized iHS = ln(10/40) = ln(0.25) = −1.386294361, a strong negative score indicating a
long derived haplotype (candidate positive selection on the derived allele) [1].

### 7.2 Applications and Use Cases

- **Genome-wide sweep scans:** ranking SNPs by |iHS| and flagging windows with an excess of
  |iHS| > 2 to localize recent positive selection [1].

### 7.3 Related Tests, Evidence, or Documents

- Tests: [PopulationGeneticsAnalyzer_SelectionSignature_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Population/PopulationGeneticsAnalyzer_SelectionSignature_Tests.cs) — covers `INV-01`…`INV-06`
- Evidence: [POP-SELECT-001-Evidence.md](../../../docs/Evidence/POP-SELECT-001-Evidence.md)

### 7.4 Change History

| Date | Version | Changes |
|------|---------|---------|
| 2026-06-13 | 1.0 | Initial canonical iHS pipeline (POP-SELECT-001). |

## 8. References

1. Voight BF, Kudaravalli S, Wen X, Pritchard JK. 2006. A Map of Recent Positive Selection in the Human Genome. PLoS Biology 4(3):e72. https://doi.org/10.1371/journal.pbio.0040072
2. Sabeti PC, Reich DE, Higgins JM, et al. 2002. Detecting recent positive selection in the human genome from haplotype structure. Nature 419:832–837. https://pubmed.ncbi.nlm.nih.gov/12397357/
3. Szpiech ZA, Hernandez RD. 2014. selscan: an efficient multi-threaded program to perform EHH-based scans for positive selection. Molecular Biology and Evolution 31(10):2824–2827. https://doi.org/10.1093/molbev/msu211
4. Gautier M, Klassmann A, Vitalis R. rehh package vignette. CRAN. https://cran.r-project.org/web/packages/rehh/vignettes/rehh.html
