# K-mer Euclidean Distance

| Field | Value |
|-------|-------|
| Algorithm Group | K-mer / Alignment-free sequence comparison |
| Test Unit ID | KMER-DIST-001 |
| Related Projects | Seqeron.Genomics.Analysis |
| Implementation Status | Production |
| Last Reviewed | 2026-06-13 |

## 1. Overview

K-mer Euclidean distance is an alignment-free dissimilarity measure between two biological
sequences. Each sequence is summarized by the frequencies of its length-*k* substrings
(k-mers); the distance is the Euclidean (L2) distance between the two frequency vectors over
the union of observed k-mers [1][2]. It is exact (deterministic), requires no alignment, and
is used for fast whole-genome phylogeny, clustering, and database screening where alignment
is impractical [1][3]. Identical sequences have distance 0 and more dissimilar sequences have
larger values [1].

## 2. Scientific / Formal Basis

### 2.1 Domain Context

Alignment-free methods map each sequence into a fixed-length vector of word (k-mer)
statistics and compare those vectors, sidestepping the cost and the recombination/shuffling
assumptions of alignment [3]. Word-composition methods map a sequence into a 4^k-dimensional
vector of k-word frequencies (for a 4-letter nucleotide alphabet), and a dissimilarity score
is obtained by a vector metric such as the Euclidean distance, Pearson correlation,
Kullback–Leibler discrepancy, or cosine distance [3].

### 2.2 Core Model

For a sequence *s* of length *L*, the number of overlapping k-mer windows is *L − k + 1*. The
k-mer **count** vector counts how many times each word *w* occurs; the k-mer **frequency**
vector normalizes each count by the total number of k-mers in the sequence (i.e. the sequence
length minus the k-mer length) [2]:

```
f_s(w) = count_s(w) / (L_s − k + 1)
```

Given the union *W* of k-mers occurring in either sequence *x* or *y*, the Euclidean distance is

```
d(x, y) = sqrt( Σ_{w ∈ W} ( f_x(w) − f_y(w) )² )
```

where a word absent from a sequence contributes a 0 component [1]. The difference between the
two word vectors "is very commonly computed by the Euclidean distance" [1], applied to the
relative-frequency vectors [4].

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | d(x, x) = 0 | Equal frequency vectors give a zero sum of squares; "identical sequences yield a distance of 0" [1] |
| INV-02 | d(x, y) = d(y, x) | (f_x − f_y)² = (f_y − f_x)²; Euclidean distance is a metric [3] |
| INV-03 | d(x, y) ≥ 0 | Square root of a sum of squares is non-negative [1][3] |
| INV-04 | Two sequences each consisting of a single distinct k-mer (frequency 1) with disjoint word sets give d = √2 | Frequency vectors are (1,0) and (0,1); √((1−0)²+(0−1)²)=√2 [2] |

### 2.5 Comparison with Related Methods

| Aspect | K-mer Euclidean (frequency) | D2 / cosine |
|--------|-----------------------------|-------------|
| Operates on | normalized frequency vectors [2][4] | raw counts (D2 is an uncentered correlation of counts) [3] |
| Range | ≥ 0, length-normalized | unbounded (D2) / [−1,1]-derived (cosine) |
| Metric | yes (Euclidean) [3] | D2 is a similarity statistic, not a metric [3] |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| seq1 | string | required | First sequence | nucleotide/protein text; null/empty allowed (empty vector) |
| seq2 | string | required | Second sequence | same as seq1 |
| k | int | required | K-mer length | must be > 0 |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| (return) | double | Non-negative Euclidean distance between the two frequency vectors; 0 for identical or both-empty inputs |

### 3.3 Preconditions and Validation

- `k ≤ 0` throws `ArgumentOutOfRangeException`.
- Inputs are upper-cased before counting, so comparison is case-insensitive (ASM-01).
- A null/empty sequence, or one shorter than *k*, produces an empty frequency vector, treated
  as the zero vector (ASM-02); the distance then equals the L2 norm of the other sequence's
  frequency vector, and 0 when both are empty.
- No alphabet restriction is enforced; any character may form a k-mer.

## 4. Algorithm

### 4.1 High-Level Steps

1. Validate `k > 0`.
2. Compute the normalized k-mer frequency vector of each sequence (count ÷ total windows) [2].
3. Form the union of k-mers present in either vector.
4. Sum the squared per-word frequency differences over the union (absent words = 0).
5. Return the square root of that sum [1].

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| KmerDistance | O(n + m) | O(u) | n, m = sequence lengths; u = number of distinct k-mers in the union; one linear pass per sequence to count, one pass over the union |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [KmerAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/KmerAnalyzer.cs)

- `KmerAnalyzer.KmerDistance(string seq1, string seq2, int k)`: returns the Euclidean distance between the two normalized k-mer frequency vectors.
- `KmerAnalyzer.GetKmerFrequencies(string, int)`: builds the per-sequence frequency vector (count ÷ sum of counts).
- `KmerAnalyzer.CountKmers(string, int)`: underlying k-mer counter (upper-cases input; throws for k ≤ 0).

### 5.2 Current Behavior

- Uses the **frequency** variant: `GetKmerFrequencies` divides each count by the sum of counts,
  which equals the number of k-mer windows (L − k + 1) for inputs whose every position forms a
  k-mer, matching the normalization in [2].
- Distance is taken over the union of observed k-mers only (sparse representation), which is
  equivalent to the full 4^k-dimensional vector since unobserved words have frequency 0 on
  both sides and contribute nothing.
- **Suffix tree not used:** this is not a substring-search / occurrence-enumeration task; it is
  a single linear k-mer counting pass per sequence followed by a vector difference. A suffix
  tree would add construction overhead without changing the O(n + m) cost, so the dictionary
  counter is the correct structure.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Word-vector model and union over both sequences with 0 components for absent words [1].
- Frequency normalization: count ÷ (L − k + 1) [2].
- Euclidean distance √(Σ (f_x − f_y)²) over the frequency vectors [1][4].
- Identity property d(x, x) = 0 [1].

**Intentionally simplified:**

- (none)

**Not implemented:**

- Count-based (un-normalized) Euclidean and other metrics (Manhattan, Canberra, Chebyshev, cosine, D2); **users should rely on:** no current alternative in this class — only the frequency Euclidean variant is provided.
- Spaced k-mers / spaced-word frequencies [4]; **users should rely on:** contiguous k-mers only.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Case-insensitive (upper-casing) | Assumption | Mixed-case inputs treated as identical k-mers | accepted | ASM-01; benign for canonical upper-case input |
| 2 | L < k ⇒ empty (zero) vector | Assumption | Defines distance for too-short inputs not covered by sources | accepted | ASM-02 |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Identical sequences | 0 | INV-01 [1] |
| k ≤ 0 | ArgumentOutOfRangeException | Input validation |
| One sequence shorter than k | distance = L2 norm of the other's frequency vector | ASM-02 (empty = zero vector) |
| Both sequences empty | 0 | ASM-02 |
| Disjoint single-k-mer sequences | √2 | INV-04 |
| Lower-case vs upper-case | same as upper-case result | ASM-01 |

### 6.2 Limitations

- Frequency normalization removes overall length information, so two sequences with the same
  relative composition but different lengths can have distance 0.
- Uses contiguous k-mers only (no spaced words); does not account for reverse-complement
  strands or for statistical background correction (D2*, D2S). The raw Euclidean value is a
  dissimilarity, not a calibrated phylogenetic distance [2].

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
double d = KmerAnalyzer.KmerDistance("ATGTGTG", "CATGTG", k: 3);
// d ≈ 0.33166247903553997
```

**Numerical walk-through (Zielezinski et al. 2017, Fig. 1):**

For x = "ATGTGTG" and y = "CATGTG", k = 3, the union of words is {ATG, CAT, GTG, TGT} with
counts c_x = (1, 0, 2, 2) and c_y = (1, 1, 1, 1) [1]. Normalizing by the window counts (5 for
x, 4 for y) gives frequencies f_x = (0.2, 0, 0.4, 0.4) and f_y = (0.25, 0.25, 0.25, 0.25).
The squared differences are (0.0025, 0.0625, 0.0225, 0.0225), summing to 0.11, so
d = √0.11 ≈ 0.33166247903553997.

### 7.3 Related Tests, Evidence, or Documents

- Tests: [KmerAnalyzer_KmerDistance_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Analysis/KmerAnalyzer_KmerDistance_Tests.cs) — covers `INV-01`–`INV-04`
- Evidence: [KMER-DIST-001-Evidence.md](../../../docs/Evidence/KMER-DIST-001-Evidence.md)

## 8. References

1. Zielezinski A, Vinga S, Almeida J, Karlowski WM. 2017. Alignment-free sequence comparison: benefits, applications, and tools. Genome Biology 18:186. https://pmc.ncbi.nlm.nih.gov/articles/PMC5627421/ (DOI: 10.1186/s13059-017-1319-7)
2. Lau AK, et al. 2022. Interpreting alignment-free sequence comparison: what makes a score a good score? NAR Genomics and Bioinformatics. https://pmc.ncbi.nlm.nih.gov/articles/PMC9442500/
3. Vinga S, Almeida J. 2003. Alignment-free sequence comparison—a review. Bioinformatics 19(4):513–523. https://academic.oup.com/bioinformatics/article/19/4/513/218529 (DOI: 10.1093/bioinformatics/btg005)
4. Boden M, et al. 2014. Fast alignment-free sequence comparison using spaced-word frequencies. Bioinformatics 30(14):1991–1999. https://pmc.ncbi.nlm.nih.gov/articles/PMC4080745/
