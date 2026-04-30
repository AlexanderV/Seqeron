# Beta Diversity Analysis

| Field | Value |
|-------|-------|
| Algorithm Group | Metagenomics |
| Test Unit ID | META-BETA-001 |
| Related Projects | N/A |
| Implementation Status | Simplified |
| Last Reviewed | 2026-04-30 |

## 1. Overview

Beta diversity measures the dissimilarity between two ecological communities rather than diversity within one sample. This document covers two classical beta-diversity measures implemented together in the repository: Bray-Curtis dissimilarity for abundance-sensitive comparison and Jaccard distance for presence/absence comparison. Both measures are deterministic functions of two sample profiles and quantify different aspects of community turnover. The current implementation returns both metrics in one `BetaDiversity` result, along with shared and sample-specific species counts and a placeholder `UniFracDistance` field.

## 2. Scientific / Formal Basis

> A = Bray-Curtis dissimilarity, B = Jaccard distance

### 2.A Bray-Curtis dissimilarity

#### Domain Context

Bray-Curtis dissimilarity is used when the comparison should retain abundance information instead of collapsing taxa to simple presence/absence. In ecology and metagenomics, it highlights how much abundance mass is shared between two samples and how much is unique to each sample (Bray and Curtis, 1957).

#### Core Model

The Bray-Curtis dissimilarity coefficient is defined as (Bray and Curtis, 1957):

$$
BC_{jk} = 1 - \frac{2 C_{jk}}{S_j + S_k}
$$

where $C_{jk}$ is the sum of the lesser abundance for each taxon shared by samples $j$ and $k$, and $S_j$ and $S_k$ are the total abundances in the two samples.

#### Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-BRAY-01 | `0 <= BrayCurtis <= 1` when abundances are non-negative | The Bray-Curtis formula is a normalized ratio of shared to total abundance (Bray and Curtis, 1957) |
| INV-BRAY-02 | `BrayCurtis = 0` for identical abundance profiles | Shared abundance then equals half of the total abundance sum |
| INV-BRAY-03 | `BrayCurtis = 1` when no positive-abundance taxon is shared and total abundance is positive | Then `C_jk = 0` in the defining formula while the guarded denominator remains non-zero |

### 2.B Jaccard distance

#### Domain Context

Jaccard distance is appropriate when the comparison should ignore abundance magnitude and focus only on whether taxa are present in each sample. It measures turnover in species membership rather than quantitative abundance balance (Jaccard, 1901).

#### Core Model

The Jaccard distance is the complement of Jaccard similarity (Jaccard, 1901):

$$
J_{distance} = 1 - \frac{|A \cap B|}{|A \cup B|}
$$

For presence/absence data this can also be written as:

$$
J_{distance} = \frac{u_A + u_B}{s + u_A + u_B}
$$

where $s$ is the number of shared taxa, $u_A$ is the number unique to sample $A$, and $u_B$ is the number unique to sample $B$.

#### Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-JACCARD-01 | `0 <= JaccardDistance <= 1` | Jaccard distance is the complement of a similarity ratio over set cardinalities |
| INV-JACCARD-02 | `JaccardDistance = 0` for identical species sets | Then `|A ∩ B| = |A ∪ B|` |
| INV-JACCARD-03 | `JaccardDistance = 1` when no taxon is shared and at least one taxon is present | Then `|A ∩ B| = 0` while `|A ∪ B| > 0` |

#### Comparison with Related Methods

| Aspect | Bray-Curtis dissimilarity | Jaccard distance |
|--------|---------------------------|------------------|
| Input interpretation | Uses abundance values | Uses presence/absence only |
| Range | `[0, 1]` | `[0, 1]` |
| Sensitivity to abundance shifts | Yes | No |
| Metric property | Not a true distance metric | True distance metric |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `sample1Name` | `string` | required | Label for the first sample | Stored verbatim in the result |
| `sample1` | `IReadOnlyDictionary<string, double>` | required | Taxon-to-abundance mapping for the first sample | Positive values are treated as present; missing keys default to zero during comparison |
| `sample2Name` | `string` | required | Label for the second sample | Stored verbatim in the result |
| `sample2` | `IReadOnlyDictionary<string, double>` | required | Taxon-to-abundance mapping for the second sample | Positive values are treated as present; missing keys default to zero during comparison |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `Sample1` | `string` | First sample label |
| `Sample2` | `string` | Second sample label |
| `BrayCurtis` | `double` | Bray-Curtis dissimilarity from the two abundance profiles |
| `JaccardDistance` | `double` | Jaccard distance from taxon presence/absence |
| `UniFracDistance` | `double` | Placeholder phylogenetic distance field; the current implementation sets it to `0` |
| `SharedSpecies` | `int` | Number of taxa present with positive abundance in both samples |
| `UniqueToSample1` | `int` | Number of taxa present only in sample 1 |
| `UniqueToSample2` | `int` | Number of taxa present only in sample 2 |

### 3.3 Preconditions and Validation

The method compares the union of taxon keys from both samples. A taxon is treated as present only when its abundance is strictly greater than zero. Missing taxa are treated as abundance `0`. If the union has no positively present taxa, Jaccard distance returns `0`. If the summed abundance over both samples is `0`, Bray-Curtis also returns `0`. The current API does not accept a phylogenetic tree, so `UniFracDistance` is always `0`.

## 4. Algorithm

### 4.A Bray-Curtis dissimilarity

#### High-Level Steps

1. Build the union of species keys from both input dictionaries.
2. Read each species abundance from both samples, defaulting missing taxa to `0`.
3. Sum the lesser abundance across all species to obtain shared abundance.
4. Sum both sample abundances across the same union.
5. Return `1 - 2 * shared / total`, or `0` when the total abundance is `0`.

#### Decision Rules / Reference Tables

The current implementation does not normalize abundance values before applying the Bray-Curtis formula. It applies the formula directly to the supplied values, so the caller controls whether they are raw counts, relative abundances, or another non-negative abundance scale.

#### Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| Bray-Curtis calculation | `O(u)` | `O(u)` | `u` = number of taxa in the union of the two samples |

### 4.B Jaccard distance

#### High-Level Steps

1. Build the union of species keys from both input dictionaries.
2. Mark each species as present in a sample only when its abundance is greater than `0`.
3. Count shared taxa, taxa unique to sample 1, and taxa unique to sample 2.
4. Return `1 - shared / (shared + unique1 + unique2)`, or `0` when the denominator is `0`.

#### Decision Rules / Reference Tables

Presence/absence is derived from abundance values by the rule `abundance > 0`. Species with zero or negative abundance are treated as absent for Jaccard counting.

#### Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| Jaccard calculation | `O(u)` | `O(u)` | `u` = number of taxa in the union of the two samples |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [MetagenomicsAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Metagenomics/MetagenomicsAnalyzer.cs)

- `MetagenomicsAnalyzer.CalculateBetaDiversity(string, IReadOnlyDictionary<string, double>, string, IReadOnlyDictionary<string, double>)`: Computes shared/unique species counts plus Bray-Curtis and Jaccard outputs for one sample pair.

### 5.2 Current Behavior

The repository calculates both metrics in one pass over the union of taxon keys. Bray-Curtis uses the supplied abundance values directly. Jaccard first converts those values to presence/absence using the rule `> 0`. Species with zero abundance are treated as absent. The result record also includes a `UniFracDistance` field, but the current implementation hard-codes it to `0` because no phylogenetic tree is part of the method contract.

### 5.3 Conformance to Theory / Spec

#### 5.3.A Bray-Curtis dissimilarity

**Implemented (verbatim from the cited theory/spec):**

- Bray-Curtis dissimilarity as a normalized function of shared abundance and total abundance.
- The expected `0` value for identical samples and `1` value for samples with no shared positive-abundance taxa when the total-abundance denominator is non-zero.

**Intentionally simplified:**

- (none)

**Not implemented:**

- (none)

#### 5.3.B Jaccard distance

**Implemented (verbatim from the cited theory/spec):**

- Jaccard distance as `1 - |A ∩ B| / |A ∪ B|` over taxon presence/absence sets.
- The expected `0` value for identical sets and `1` value for disjoint non-empty sets.

**Intentionally simplified:**

- (none)

**Not implemented:**

- (none)

### 5.4 Deviations and Assumptions (Optional)

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | `UniFracDistance` placeholder | Deviation | Users receive `0` for the phylogenetic distance field even though UniFrac requires a phylogenetic tree | accepted | The current API computes only Bray-Curtis and Jaccard |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Identical samples | `BrayCurtis = 0`, `JaccardDistance = 0` | Both formulas report zero dissimilarity when the profiles match exactly |
| No shared positive-abundance taxa | `BrayCurtis = 1`, `JaccardDistance = 1` | Shared abundance and set intersection both become zero |
| Empty dictionaries or all-zero abundances | `BrayCurtis = 0`, `JaccardDistance = 0` | Both implementations guard their denominators and return zero when no informative comparison exists |
| Taxon present in one sample with value `0` in the other | Counted as unique to the positive-abundance sample | Presence is defined by `abundance > 0` |

### 6.2 Limitations

These standard beta-diversity measures do not incorporate phylogenetic relationships between taxa. The current implementation also does not model sampling effort explicitly, and presence/absence calculations may be dominated by rare taxa because abundance magnitude is discarded for the Jaccard branch.

## 7. Examples and Related Material

### 7.2 Applications and Use Cases

- Comparing microbial community composition between samples.
- Assessing habitat similarity in ecological studies.
- Monitoring temporal changes in community structure.
- Evaluating treatment effects on species composition.

## 8. References

1. Bray, J. R., and J. T. Curtis. 1957. An ordination of the upland forest communities of southern Wisconsin. Ecological Monographs 27(4):325-349.
2. Jaccard, P. 1901. Étude comparative de la distribution florale dans une portion des Alpes et des Jura. Bulletin de la Société Vaudoise des Sciences Naturelles 37(142):547-579.
3. Whittaker, R. H. 1960. Vegetation of the Siskiyou Mountains, Oregon and California. Ecological Monographs 30(3):279-338.
