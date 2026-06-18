# K-mer Statistics

| Field | Value |
|-------|-------|
| Algorithm Group | K-mer |
| Test Unit ID | KMER-STATS-001 |
| Related Projects | Seqeron.Genomics.Analysis |
| Implementation Status | Production |
| Last Reviewed | 2026-06-14 |

## 1. Overview

K-mer statistics summarize the composition of a sequence by the multiset of its overlapping length-k substrings. `AnalyzeKmers` reports the total number of k-mers, the number of distinct k-mers, the maximum/minimum/average multiplicity, and the Shannon entropy of the k-mer frequency distribution. These quantities characterize sequence diversity and repetitiveness: high entropy with a high distinct/total ratio indicates a diverse sequence, while low entropy indicates repetitive composition [1][3]. The computation is exact (not heuristic): every value is determined directly from the k-mer count table.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

A *k-mer* is a substring of length k. For a sequence of length L the overlapping k-mers are the L − k + 1 windows starting at positions 0 … L − k [1][2]. Counting these k-mers yields a multiplicity (frequency) for each distinct k-mer; the resulting distribution is the basis for diversity and complexity measures used in genome analysis, assembly, and alignment-free comparison [3].

### 2.2 Core Model

Let the count table be `mult(α)` for each distinct k-mer α occurring in the sequence.

- **Total k-mers:** `T = L − k + 1` [1][2]; equivalently `T = Σ_α mult(α)`.
- **Distinct k-mers:** `D = |{α : mult(α) > 0}|` (each different k-mer counted once) [1][2].
- **Max / Min multiplicity:** `max_α mult(α)`, `min_α mult(α)`.
- **Average multiplicity:** `T / D`.
- **Shannon entropy:** `E_k = − Σ_α p(α) log₂ p(α)` with `p(α) = mult(α) / T`, the relative frequency of α over the T windows [3]. The same single-sequence form `H_k(s) = − Σ_i p_i log₂ p_i` (p_i = relative frequency of the i-th k-mer) is used as a sequence complexity measure in bits [4]. The convention `0 · log 0 = 0` applies [4].

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | `TotalKmers = L − k + 1` for L ≥ k | number of overlapping length-k windows [1][2] |
| INV-02 | `TotalKmers = Σ_α mult(α)` | each window contributes exactly one k-mer count |
| INV-03 | `UniqueKmers` = number of distinct k-mers | distinct count of the k-mer table [1][2] |
| INV-04 | `MinCount ≤ AverageCount ≤ MaxCount` and `AverageCount = TotalKmers / UniqueKmers` | arithmetic mean of multiplicities |
| INV-05 | `0 ≤ Entropy ≤ log₂(UniqueKmers)`; Entropy = 0 iff one distinct k-mer; Entropy = log₂(D) iff all multiplicities equal | Shannon entropy bounds for a D-symbol distribution [3][4] |
| INV-06 | empty sequence or k > L ⇒ all fields = 0 | `L − k + 1 ≤ 0` ⇒ no k-mers [1] |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `sequence` | `string` | required | Sequence to analyze | null/empty ⇒ all-zero result; upper-cased internally (case-insensitive) |
| `k` | `int` | required | K-mer length | Must be > 0; k > L ⇒ all-zero result |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `TotalKmers` | `int` | Total overlapping k-mers, L − k + 1 |
| `UniqueKmers` | `int` | Number of distinct k-mers |
| `MaxCount` | `int` | Maximum k-mer multiplicity |
| `MinCount` | `int` | Minimum k-mer multiplicity |
| `AverageCount` | `double` | Mean multiplicity (TotalKmers / UniqueKmers), rounded to 2 decimals |
| `Entropy` | `double` | Shannon entropy of the k-mer frequencies, −Σ p log₂ p, in bits |

### 3.3 Preconditions and Validation

Input is upper-cased (case-insensitive); no alphabet restriction (any character may form a k-mer). 0-based windows. `k ≤ 0` throws `ArgumentOutOfRangeException` (via `CountKmers`). Null/empty sequence and `k > L` return `KmerStatistics(0,0,0,0,0,0)` because no k-mers exist (L − k + 1 ≤ 0) [1].

## 4. Algorithm

### 4.1 High-Level Steps

1. Build the k-mer count table with `CountKmers(sequence, k)` (one pass over the L − k + 1 windows).
2. If the table is empty, return the all-zero `KmerStatistics`.
3. Compute `TotalKmers` = sum of counts, `UniqueKmers` = table size, `MaxCount`/`MinCount` = extremes, `AverageCount` = mean (rounded to 2 decimals).
4. Compute `Entropy` = −Σ (count/total) log₂(count/total) over the table.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `AnalyzeKmers` | O(L·k) | O(D·k) | L − k + 1 windows, each k-mer materialized as a length-k string; D distinct k-mers stored in the count map |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [KmerAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/KmerAnalyzer.cs)

- `KmerAnalyzer.AnalyzeKmers(string, int)`: returns the `KmerStatistics` record.
- `KmerAnalyzer.CountKmers(string, int)`: builds the count table (reused).
- `KmerAnalyzer.CalculateKmerEntropy(string, int)`: computes the Shannon entropy term used for `Entropy`.

### 5.2 Current Behavior

`AnalyzeKmers` delegates counting to `CountKmers` and entropy to `CalculateKmerEntropy` (which normalizes via `GetKmerFrequencies`, dividing each count by the total = L − k + 1). `AverageCount` is rounded to two decimals via `Math.Round(x, 2)` for display; the underlying ratio is exact. The unit is not a search/matching operation (it aggregates a precomputed count table), so the repository **suffix tree was not used** — there is no occurrence-enumeration or pattern-location subproblem here; the linear count-table scan in `CountKmers` is optimal.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- `TotalKmers = L − k + 1` [1][2].
- `UniqueKmers` = distinct k-mer count [1][2].
- `Entropy = −Σ p(α) log₂ p(α)`, p(α) = mult(α)/(L − k + 1) [3][4].
- Max/Min/Average multiplicity over the count table.

**Intentionally simplified:**

- `AverageCount` is rounded to 2 decimals for display; **consequence:** the reported mean is a rounded view of the exact ratio TotalKmers/UniqueKmers.

**Not implemented:**

- Canonical (reverse-complement-collapsed) k-mer statistics; **users should rely on:** `KmerAnalyzer.CountKmersBothStrands` for strand-aware counting (KMER-BOTH-001).

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Empty / null sequence | all fields 0 | no k-mers |
| k > L | all fields 0 | L − k + 1 ≤ 0 [1] |
| k ≤ 0 | `ArgumentOutOfRangeException` | k-mer length must be positive |
| Single distinct k-mer (homopolymer) | Entropy = 0; Max = Min = Total | one-component distribution [3] |
| All windows distinct | Entropy = log₂(Total); Max = Min = 1 | uniform distribution [3] |
| Lower-case input | identical to upper-case stats | upper-cased internally |

### 6.2 Limitations

The `UniqueKmers` field name denotes the **distinct** k-mer count, not the count-1 "unique" set computed by `FindUniqueKmers` (KMER-UNIQUE-001) — a documented naming caveat to avoid the distinct/unique confusion [2]. No IUPAC-degenerate handling: ambiguous symbols form ordinary k-mers.

## 7. Examples and Related Material

### 7.1 Worked Example

**Numerical walk-through (GTAGAGCTGT, Wikipedia K-mer example [1]):**

- k=1: counts G=4, T=3, A=2, C=1 ⇒ Total=10, Distinct=4, Max=4, Min=1, Avg=2.5, Entropy = −(0.4·log₂0.4 + 0.3·log₂0.3 + 0.2·log₂0.2 + 0.1·log₂0.1) = 1.846439… bits.
- k=3: all 8 windows distinct ⇒ Total=8, Distinct=8, Max=Min=1, Avg=1.0, Entropy = log₂8 = 3.0 bits.

**API usage example:**

```csharp
var stats = KmerAnalyzer.AnalyzeKmers("GTAGAGCTGT", 1);
// stats.TotalKmers == 10, stats.UniqueKmers == 4, stats.MaxCount == 4,
// stats.MinCount == 1, stats.AverageCount == 2.5, stats.Entropy ≈ 1.84644
```

### 7.3 Related Tests, Evidence, or Documents

- Tests: [KmerAnalyzer_AnalyzeKmers_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/KmerAnalyzer_AnalyzeKmers_Tests.cs) — covers `INV-01`–`INV-06`
- Evidence: [KMER-STATS-001-Evidence.md](../../../docs/Evidence/KMER-STATS-001-Evidence.md)
- Related algorithms: [Unique_And_MinCount_Kmers](../K-mer/Unique_And_MinCount_Kmers.md), [Both_Strand_Kmer_Counting](../K-mer/Both_Strand_Kmer_Counting.md)

## 8. References

1. Wikipedia contributors. 2026. K-mer. Wikipedia. https://en.wikipedia.org/wiki/K-mer
2. Clavijo, B. 2018. k-mer counting, part I: Introduction. BioInfoLogics. https://bioinfologics.github.io/post/2018/09/17/k-mer-counting-part-i-introduction/
3. Manca, V. et al. 2021. Spectral concepts in genome informational analysis. arXiv preprint. arXiv:2106.15351. https://arxiv.org/abs/2106.15351
4. Entropy–Rank Ratio: A Novel Entropy–Based Perspective for DNA Complexity and Classification. 2025. arXiv preprint. arXiv:2511.05300. https://arxiv.org/html/2511.05300
