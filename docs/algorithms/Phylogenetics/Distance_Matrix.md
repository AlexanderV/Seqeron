# Phylogenetic Distance Matrix Calculation

> **Baseline / reference method.** The Jukes-Cantor and Kimura-2P corrections here are canonical baseline substitution models with simplifying assumptions; they are not always best for real datasets. See [Legacy / Baseline Methods](../CANONICAL_MAP.md).

| Field | Value |
|-------|-------|
| Algorithm Group | Phylogenetics |
| Test Unit ID | PHYLO-DIST-001 |
| Related Projects | N/A |
| Implementation Status | Simplified |
| Last Reviewed | 2026-04-30 |

## 1. Overview

Phylogenetic distance matrices summarize evolutionary divergence between pairs of aligned sequences and provide the numeric substrate for distance-based tree construction. The repository supports four pairwise distance methods in one surface: raw Hamming distance, p-distance, Jukes-Cantor correction, and Kimura 2-parameter correction. All four methods share the same pairwise scan over aligned sequences and differ only in how mismatches are converted into a distance value. The current implementation also applies pairwise deletion for gaps and ambiguous bases.

## 2. Scientific / Formal Basis

> A = Hamming distance, B = p-distance, C = Jukes-Cantor distance, D = Kimura 2-parameter distance

### 2.A Hamming distance

#### Domain Context

Hamming distance is the simplest mismatch count on aligned sequences. It treats every comparable site independently and does not attempt to correct for multiple substitutions.

#### Core Model

$$
d_H(s_1, s_2) = \sum_i \mathbb{1}[s_1[i] \ne s_2[i]]
$$

over the set of comparable positions.

#### Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-HAMMING-01 | Hamming distance is a non-negative integer count | It counts differing comparable positions |
| INV-HAMMING-02 | Identical comparable sequences have distance `0` | No positions contribute to the count |

### 2.B p-distance

#### Domain Context

p-distance converts raw mismatches into a proportion of differing comparable sites. It is commonly used for closely related sequences when no correction for multiple substitutions is desired.

#### Core Model

$$
p = \frac{\text{number of differences}}{\text{number of comparable sites}}
$$

#### Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-PDIST-01 | `0 <= p <= 1` | The numerator cannot exceed the denominator |
| INV-PDIST-02 | p-distance underestimates divergence when multiple substitutions occur at the same site | It treats each observed difference as a single event |

### 2.C Jukes-Cantor distance

#### Domain Context

Jukes-Cantor (JC69) corrects p-distance for unseen multiple substitutions under a one-parameter nucleotide-substitution model (Jukes and Cantor, 1969).

#### Core Model

$$
d_{JC} = -\frac{3}{4} \ln\left(1 - \frac{4p}{3}\right)
$$

where $p$ is the p-distance.

#### Modeling Assumptions

| ID | Assumption | Consequence if Violated |
|----|------------|--------------------------|
| ASM-JC-01 | All nucleotides have equal equilibrium frequencies | JC69 can be a poor correction when base composition is biased |
| ASM-JC-02 | All substitutions occur at the same rate | Transition/transversion bias is not represented |

#### Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-JC-01 | `d_JC >= p` when the correction is defined | JC69 corrects upward for hidden substitutions |
| INV-JC-02 | The correction becomes undefined when `p >= 0.75` | The logarithm argument becomes non-positive |

### 2.D Kimura 2-parameter distance

#### Domain Context

Kimura's 2-parameter model (K80/K2P) refines JC69 by distinguishing transitions from transversions, which often occur at different rates in real data (Kimura, 1980).

#### Core Model

$$
d_{K2P} = -\frac{1}{2} \ln\left((1 - 2S - V) \sqrt{1 - 2V}\right)
$$

where $S$ is the proportion of transitions and $V$ is the proportion of transversions.

#### Modeling Assumptions

| ID | Assumption | Consequence if Violated |
|----|------------|--------------------------|
| ASM-K2P-01 | Transition and transversion rates are each homogeneous within their class | More complex substitution asymmetries are not modeled |

#### Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-K2P-01 | K2P distinguishes `A<->G` and `C<->T` transitions from transversions | The model is parameterized by separate `S` and `V` terms |
| INV-K2P-02 | The correction becomes undefined when either logarithm argument is non-positive | The implementation returns positive infinity in that case |

#### Comparison with Related Methods

| Aspect | Hamming | p-distance | Jukes-Cantor | Kimura 2-parameter |
|--------|---------|------------|--------------|--------------------|
| Output scale | Raw count | Proportion | Corrected evolutionary distance | Corrected evolutionary distance |
| Multiple-substitution correction | No | No | Yes | Yes |
| Transition/transversion distinction | No | No | No | Yes |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `[CalculateDistanceMatrix] alignedSequences` | `IReadOnlyList<string>` | required | Equal-length aligned sequences | Required for meaningful pairwise comparison |
| `[CalculateDistanceMatrix] method` | `DistanceMethod` | `JukesCantor` | Distance method used for every pairwise comparison | One of `Hamming`, `PDistance`, `JukesCantor`, or `Kimura2Parameter` |
| `[CalculatePairwiseDistance] seq1` | `string` | required | First aligned sequence | Must have the same length as `seq2` |
| `[CalculatePairwiseDistance] seq2` | `string` | required | Second aligned sequence | Must have the same length as `seq1` |
| `[CalculatePairwiseDistance] method` | `DistanceMethod` | `JukesCantor` | Distance method for the pairwise comparison | Same enum as the matrix builder |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `[CalculateDistanceMatrix] return value` | `double[,]` | Symmetric `n x n` matrix with zero diagonal |
| `[CalculatePairwiseDistance] return value` | `double` | Pairwise distance under the selected method |

### 3.3 Preconditions and Validation

Pairwise comparison is case-insensitive because the implementation uppercases each site before inspection. Only standard DNA bases `A/C/G/T` are comparable; gaps and ambiguous IUPAC symbols are skipped by pairwise deletion. If no comparable sites remain, the pairwise distance returns `0`. `CalculatePairwiseDistance` throws `ArgumentException` when the two sequences have different lengths, and `CalculateDistanceMatrix` relies on that pairwise routine for its sequence pairs.

## 4. Algorithm

### 4.A Hamming distance

#### High-Level Steps

1. Uppercase both sequences.
2. Skip sites with gaps or non-standard bases.
3. Count the remaining sites where the bases differ.
4. Return the mismatch count.

#### Decision Rules / Reference Tables

Only positions where both bases are one of `A/C/G/T` contribute to the count.

#### Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| Hamming distance | `O(m)` | `O(1)` | `m` = sequence length |

### 4.B p-distance

#### High-Level Steps

1. Perform the shared pairwise scan over comparable sites.
2. Count differences and comparable positions.
3. Return `differences / comparableSites`, or `0` when no comparable sites remain.

#### Decision Rules / Reference Tables

The same comparable-site filter used for Hamming distance defines the denominator of p-distance.

#### Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| p-distance | `O(m)` | `O(1)` | `m` = sequence length |

### 4.C Jukes-Cantor distance

#### High-Level Steps

1. Compute p-distance on comparable sites.
2. Evaluate the JC69 correction formula.
3. Return positive infinity if the logarithm argument is non-positive.

#### Decision Rules / Reference Tables

The correction is applied only after the shared comparable-site scan; ambiguous or gap-containing sites do not contribute to `p`.

#### Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| Jukes-Cantor distance | `O(m)` | `O(1)` | Pairwise scan plus constant-time correction |

### 4.D Kimura 2-parameter distance

#### High-Level Steps

1. Compute counts of transitions, transversions, and comparable sites.
2. Convert those counts into `S` and `V` proportions.
3. Evaluate the K2P correction formula.
4. Return positive infinity if either correction argument is non-positive.

#### Decision Rules / Reference Tables

The implementation treats `A<->G` and `C<->T` as transitions; all other standard-base mismatches are transversions.

#### Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| Kimura 2-parameter distance | `O(m)` | `O(1)` | Pairwise scan plus constant-time correction |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [PhylogeneticAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Phylogenetics/PhylogeneticAnalyzer.cs)

- `PhylogeneticAnalyzer.CalculateDistanceMatrix(IReadOnlyList<string>, DistanceMethod)`: Builds a full symmetric matrix from aligned sequences.
- `PhylogeneticAnalyzer.CalculatePairwiseDistance(string, string, DistanceMethod)`: Computes one pairwise distance.

### 5.2 Current Behavior

The repository performs one shared pairwise scan regardless of the requested method. Gaps and ambiguous characters are skipped, so distance is computed only on positions where both sequences contain standard DNA bases. `CalculateDistanceMatrix` initializes the diagonal to zero and mirrors pairwise distances across the matrix. The JC69 and K2P helpers return `double.PositiveInfinity` when the correction formula becomes undefined at high divergence.

### 5.3 Conformance to Theory / Spec

#### 5.3.A Hamming distance

**Implemented (verbatim from the cited theory/spec):**

- Raw mismatch counting over comparable sites.

**Intentionally simplified:**

- (none)

**Not implemented:**

- (none)

#### 5.3.B p-distance

**Implemented (verbatim from the cited theory/spec):**

- Uncorrected proportion of differing comparable sites.

**Intentionally simplified:**

- (none)

**Not implemented:**

- (none)

#### 5.3.C Jukes-Cantor distance

**Implemented (verbatim from the cited theory/spec):**

- JC69 correction from p-distance.

**Intentionally simplified:**

- (none)

**Not implemented:**

- (none)

#### 5.3.D Kimura 2-parameter distance

**Implemented (verbatim from the cited theory/spec):**

- K2P correction using separate transition and transversion proportions.

**Intentionally simplified:**

- (none)

**Not implemented:**

- (none)

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Identical comparable sequences | Distance `0` for all methods | No mismatches are observed |
| All sites skipped because of gaps or ambiguity | Distance `0` | The implementation returns `0` when `comparableSites = 0` |
| `p >= 0.75` in JC69 | Returns positive infinity | The JC69 logarithm argument becomes non-positive |
| Unequal sequence lengths | `ArgumentException` | Pairwise comparison requires aligned sequences |

### 6.2 Limitations

The current implementation supports only four simple nucleotide-distance models and uses pairwise deletion for gaps and ambiguous symbols. It does not provide richer substitution models, codon-aware distances, or ambiguity-code matching beyond site exclusion.

## 7. Examples and Related Material

- [Tree_Construction.md](./Tree_Construction.md) documents the downstream distance-based tree builders that consume these matrices.
- [../../../tests/TestSpecs/PHYLO-DIST-001.md](../../../tests/TestSpecs/PHYLO-DIST-001.md) captures the repository acceptance scenarios for pairwise distance calculation.

## 8. References

1. Jukes, T. H., and C. R. Cantor. 1969. Evolution of Protein Molecules. In Mammalian Protein Metabolism. Academic Press.
2. Kimura, M. 1980. A simple method for estimating evolutionary rates of base substitutions through comparative studies of nucleotide sequences. Journal of Molecular Evolution 16:111-120.
3. Felsenstein, J. 2004. Inferring Phylogenies. Sinauer Associates.
4. Wikipedia contributors. Models of DNA evolution. Wikipedia. https://en.wikipedia.org/wiki/Models_of_DNA_evolution
5. Wikipedia contributors. Distance matrices in phylogeny. Wikipedia. https://en.wikipedia.org/wiki/Distance_matrices_in_phylogeny

