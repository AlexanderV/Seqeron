# Alpha Diversity

| Field | Value |
|-------|-------|
| Algorithm Group | Metagenomics |
| Test Unit ID | META-ALPHA-001 |
| Related Projects | N/A |
| Implementation Status | Simplified |
| Last Reviewed | 2026-04-30 |

## 1. Overview

Alpha diversity measures diversity within a single sample or site by combining richness and evenness into one set of summary statistics. This document covers observed species richness, Shannon entropy, Simpson concentration, inverse Simpson diversity, Pielou's evenness, and the Chao1 richness estimator. The underlying formulas are deterministic summary statistics over one abundance profile rather than a heuristic search procedure. In this repository, the implementation accepts one taxon-to-abundance map and reports all six metrics in one result record.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

Whittaker introduced the alpha-beta-gamma diversity framework to distinguish within-sample diversity from between-sample turnover and regional diversity (Whittaker, 1960). Alpha diversity therefore answers a local question: how many taxa are present in one sample, and how evenly are their abundances distributed within that sample.

### 2.2 Core Model

Observed species richness counts taxa with positive abundance:

$$
S_{obs} = \left|\{i : p_i > 0\}\right|
$$

Shannon diversity measures uncertainty in the identity of a randomly selected individual (Shannon, 1948):

$$
H = -\sum_{i=1}^{S} p_i \ln(p_i)
$$

Simpson concentration measures the probability that two randomly selected individuals belong to the same species (Simpson, 1949):

$$
\lambda = \sum_{i=1}^{S} p_i^2
$$

Inverse Simpson diversity converts Simpson concentration into an effective-number measure (Hill, 1973):

$$
D = \frac{1}{\lambda}
$$

Pielou's evenness rescales Shannon diversity by the maximum entropy at the observed richness (Pielou, 1966):

$$
J = \frac{H}{\ln(S)}
$$

Chao1 estimates unseen richness from singleton and doubleton counts (Chao, 1984):

$$
\hat{S}_{Chao1} = S_{obs} + \frac{f_1^2}{2f_2}
$$

where $f_1$ is the number of singleton taxa and $f_2$ is the number of doubleton taxa.

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | `ObservedSpecies` equals the number of taxa whose abundance is strictly greater than zero | This is the definition of observed richness and the implementation filters non-positive values before counting |
| INV-02 | `ShannonIndex` is `0` for a single-species sample | Shannon entropy is zero when one category has probability 1 (Shannon, 1948) |
| INV-03 | `SimpsonIndex` is `1` for a single-species sample and decreases as abundance is spread across taxa | The Simpson concentration formula sums squared proportions (Simpson, 1949) |
| INV-04 | `InverseSimpson` is `1 / SimpsonIndex` when `SimpsonIndex > 0` | The inverse Simpson definition is the reciprocal of Simpson concentration (Hill, 1973) |
| INV-05 | `PielouEvenness` is `0` when `ObservedSpecies <= 1` in the current implementation | Pielou's evenness is undefined when `S <= 1`, and the implementation returns `0` in that branch |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `abundances` | `IReadOnlyDictionary<string, double>` | required | Mapping from taxon name to abundance value | `null` or empty input returns an all-zero result; values `<= 0` are ignored for all reported metrics |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `ShannonIndex` | `double` | Shannon diversity computed with the natural logarithm |
| `SimpsonIndex` | `double` | Simpson concentration index |
| `InverseSimpson` | `double` | Reciprocal of Simpson concentration when the Simpson value is positive |
| `Chao1Estimate` | `double` | Chao1 richness estimate for integer count data, or observed richness for non-integer abundance data |
| `ObservedSpecies` | `double` | Count of taxa with strictly positive abundance, stored in the record as `double` |
| `PielouEvenness` | `double` | Shannon diversity divided by `ln(S)` when `S > 1`, otherwise `0` |

### 3.3 Preconditions and Validation

`null` or empty input returns `0` for every field. The implementation filters out abundance values less than or equal to zero before computing any metric. Shannon and Simpson are computed after internal normalization by the total positive abundance, so callers may pass either counts or proportional abundances. Chao1 only uses singleton and doubleton logic when every positive abundance is integer-valued within the method's tolerance; otherwise the implementation returns observed richness.

## 4. Algorithm

### 4.1 High-Level Steps

1. Return an all-zero `AlphaDiversity` record when the input is `null` or empty.
2. Collect only strictly positive abundance values and count them as `ObservedSpecies`.
3. Return an all-zero record if no positive abundances remain.
4. Compute Shannon diversity from normalized positive abundances.
5. Compute Simpson concentration, then derive inverse Simpson when the Simpson value is positive.
6. Compute Pielou's evenness when more than one species is observed.
7. Compute Chao1 from singleton and doubleton counts for integer count data, or fall back to observed richness for non-integer abundance data.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

The implementation distinguishes between count data and proportional data before computing Chao1. If every positive abundance is effectively an integer value, it counts singleton taxa (`f_1`) and doubleton taxa (`f_2`) and uses the Chao1 estimator. If `f_2 = 0`, it uses the bias-corrected branch documented in the source code: `S_obs + f_1 (f_1 - 1) / 2`. If any positive abundance is non-integer, it returns `ObservedSpecies` for `Chao1Estimate` instead of applying the count-data formula.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `CalculateAlphaDiversity` | `O(n)` | `O(n)` | `n` = number of taxa; the implementation materializes the positive abundance values into a list |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [MetagenomicsAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Metagenomics/MetagenomicsAnalyzer.cs)

- `MetagenomicsAnalyzer.CalculateAlphaDiversity(IReadOnlyDictionary<string, double>)`: Computes all alpha-diversity fields and returns an `AlphaDiversity` record.

### 5.2 Current Behavior

The repository implementation uses `Math.Log`, so Shannon diversity is reported with the natural logarithm. Non-positive abundances are ignored. `ObservedSpecies` is populated from the number of remaining positive entries, even though the public record field type is `double`. Chao1 is more precise than the older document text implied: the code uses the standard singleton/doubleton formula for integer count data and a bias-corrected `f_2 = 0` branch, but it falls back to observed richness when the input behaves like proportional abundance data rather than counts.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Observed richness as the count of taxa with positive abundance.
- Shannon diversity using $-\sum p_i \ln(p_i)$.
- Simpson concentration using $\sum p_i^2$.
- Inverse Simpson as the reciprocal of Simpson concentration.
- Pielou's evenness as $H / \ln(S)$ when more than one species is observed.

**Intentionally simplified:**

- Chao1 on proportional or otherwise non-integer abundance data falls back to `ObservedSpecies`; **consequence:** unseen-richness correction is not applied unless the input behaves like count data.

**Not implemented:**

- (none)

### 5.4 Deviations and Assumptions (Optional)

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Chao1 fallback for non-integer abundance data | Deviation | Users supplying relative abundances receive observed richness instead of a singleton/doubleton-based unseen-richness estimate | accepted | The code only applies Chao1 when positive abundances are integer-valued |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Empty input | All reported fields are `0` | The method returns an all-zero record when the input is `null` or empty |
| Single species with positive abundance | `ShannonIndex = 0`, `SimpsonIndex = 1`, `InverseSimpson = 1`, `PielouEvenness = 0`, `ObservedSpecies = 1` | This follows the entropy and Simpson definitions for a one-category distribution |
| Two equally abundant species `(0.5, 0.5)` | `ShannonIndex = ln(2)`, `SimpsonIndex = 0.5`, `InverseSimpson = 2`, `PielouEvenness = 1` | Equal abundances maximize evenness at fixed richness |
| Four equally abundant species `(0.25, 0.25, 0.25, 0.25)` | `ShannonIndex = ln(4)`, `SimpsonIndex = 0.25`, `InverseSimpson = 4`, `PielouEvenness = 1` | Equal abundances maximize evenness at richness `4` |
| Highly uneven abundances such as `(0.99, 0.01)` | Shannon decreases, Simpson approaches `1`, inverse Simpson approaches `1`, and evenness stays below `1` | Dominance by one taxon reduces evenness and effective diversity |

### 6.2 Limitations

Chao1 is only meaningful in this implementation when the input behaves like count data with singleton and doubleton counts. When callers provide proportional abundances, the method preserves the output shape by returning observed richness for `Chao1Estimate` instead of estimating unseen taxa.

## 7. Examples and Related Material

- [META-PROF-001](../../../tests/TestSpecs/META-PROF-001.md) documents taxonomic profile behavior that uses Shannon and Simpson diversity values.
- [Beta_Diversity.md](./Beta_Diversity.md) covers cross-sample diversity instead of within-sample diversity.

## 8. References

1. Shannon, C. E. 1948. A Mathematical Theory of Communication. Bell System Technical Journal.
2. Simpson, E. H. 1949. Measurement of Diversity. Nature 163(4148):688.
3. Hill, M. O. 1973. Diversity and Evenness: A Unifying Notation and Its Consequences. Ecology 54(2):427-432.
4. Chao, A. 1984. Non-parametric estimation of the number of classes in a population. Scandinavian Journal of Statistics 11:265-270.
5. Pielou, E. C. 1966. The measurement of diversity in different types of biological collections. Journal of Theoretical Biology 13:131-144.
6. Whittaker, R. H. 1960. Vegetation of the Siskiyou Mountains, Oregon and California. Ecological Monographs 30(3):279-338.
7. Wikipedia contributors. Diversity index. Wikipedia. https://en.wikipedia.org/wiki/Diversity_index
8. Wikipedia contributors. Alpha diversity. Wikipedia. https://en.wikipedia.org/wiki/Alpha_diversity
