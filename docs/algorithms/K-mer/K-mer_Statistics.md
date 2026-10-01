# K-mer Statistics

| Field | Value |
|-------|-------|
| Algorithm Group | K-mer |
| Test Unit ID | KMER-STATS-001 |
| Related Projects | Seqeron.Genomics.Analysis |
| Implementation Status | Production |
| Last Reviewed | 2026-10-01 (B06 review) |

## 1. Overview

K-mer statistics summarize the composition of a sequence by the multiset of its overlapping length-k substrings. `AnalyzeKmers` reports the Jellyfish `stats` fields — Total, Distinct, Unique (count-1, here `SingletonKmers`) and Max_count [5] — plus the minimum/mean multiplicity and the Shannon entropy of the k-mer frequency distribution. An overload applies Jellyfish's `-L/--lower-count` and `-U/--upper-count` filters [5]. These quantities characterize sequence diversity and repetitiveness: high entropy with a high distinct/total ratio indicates a diverse sequence, while low entropy indicates repetitive composition [1][3]. The computation is exact (not heuristic): every value is determined directly from the k-mer count table.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

A *k-mer* is a substring of length k. For a sequence of length L the overlapping k-mers are the L − k + 1 windows starting at positions 0 … L − k [1][2]. Counting these k-mers yields a multiplicity (frequency) for each distinct k-mer; the resulting distribution is the basis for diversity and complexity measures used in genome analysis, assembly, and alignment-free comparison [3].

### 2.2 Core Model

Let the count table be `mult(α)` for each distinct k-mer α occurring in the sequence.

- **Total k-mers:** `T = L − k + 1` [1][2]; equivalently `T = Σ_α mult(α)`.
- **Distinct k-mers:** `D = |{α : mult(α) > 0}|` (each different k-mer counted once) [1][2]; Jellyfish "Distinct" [5].
- **Singleton k-mers:** `S = |{α : mult(α) = 1}|`; Jellyfish "Unique" (`uniq += val == 1`) [5]; BioInfoLogics "unique" [2].
- **Max / Min multiplicity:** `max_α mult(α)`, `min_α mult(α)`.
- **Average multiplicity:** `T / D` (exact).
- **Count filter (optional):** only k-mers with `lower ≤ mult(α) ≤ upper` are retained (Jellyfish `compute_stats`: `if(val < low || val > high) continue;`) [5]; every statistic, entropy included (p = mult/T over retained), is computed over the retained k-mers.
- **Shannon entropy:** `E_k = − Σ_α p(α) log₂ p(α)` with `p(α) = mult(α) / T`, the relative frequency of α over the T windows [3]. The same single-sequence form `H_k(s) = − Σ_i p_i log₂ p_i` (p_i = relative frequency of the i-th k-mer) is used as a sequence complexity measure in bits [4]. The convention `0 · log 0 = 0` applies [4].

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | `TotalKmers = L − k + 1` for L ≥ k | number of overlapping length-k windows [1][2] |
| INV-02 | `TotalKmers = Σ_α mult(α)` | each window contributes exactly one k-mer count |
| INV-03 | `UniqueKmers` = number of distinct k-mers | distinct count of the k-mer table [1][2] |
| INV-04 | `MinCount ≤ AverageCount ≤ MaxCount` and `AverageCount = TotalKmers / UniqueKmers` exactly | arithmetic mean of multiplicities |
| INV-07 | `0 ≤ SingletonKmers ≤ DistinctKmers`; `SingletonKmers = |FindUniqueKmers|` | Jellyfish `stats` Unique ≤ Distinct [5] |
| INV-05 | `0 ≤ Entropy ≤ log₂(UniqueKmers)`; Entropy = 0 iff one distinct k-mer; Entropy = log₂(D) iff all multiplicities equal | Shannon entropy bounds for a D-symbol distribution [3][4] |
| INV-06 | empty sequence or k > L ⇒ all fields = 0 | `L − k + 1 ≤ 0` ⇒ no k-mers [1] |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `sequence` | `string` | required | Sequence to analyze | null/empty ⇒ all-zero result; upper-cased internally (case-insensitive) |
| `k` | `int` | required | K-mer length | Must be > 0; k > L ⇒ all-zero result |
| `lowerCount` | `int` | 0 | Jellyfish `-L`: ignore k-mers with count < lowerCount (overload) | ≥ 0 |
| `upperCount` | `int` | `int.MaxValue` | Jellyfish `-U`: ignore k-mers with count > upperCount (overload) | ≥ 0; upper < lower ⇒ all-zero |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `TotalKmers` | `int` | Total k-mers including multiplicity (Jellyfish Total), L − k + 1 when unfiltered |
| `UniqueKmers` | `int` | Number of **distinct** k-mers (legacy name, kept for compatibility; = `DistinctKmers`) |
| `DistinctKmers` | `int` | Number of distinct k-mers (Jellyfish Distinct) |
| `SingletonKmers` | `int` | Number of k-mers with count 1 (Jellyfish Unique) |
| `MaxCount` | `int` | Maximum k-mer multiplicity |
| `MinCount` | `int` | Minimum k-mer multiplicity |
| `AverageCount` | `double` | Exact mean multiplicity, TotalKmers / UniqueKmers |
| `Entropy` | `double` | Shannon entropy of the k-mer frequencies, −Σ p log₂ p, in bits |

### 3.3 Preconditions and Validation

Input is upper-cased (case-insensitive); no alphabet restriction (any character may form a k-mer). 0-based windows. `k ≤ 0` throws `ArgumentOutOfRangeException` (via `CountKmers`). Null/empty sequence and `k > L` return `KmerStatistics(0,0,0,0,0,0)` because no k-mers exist (L − k + 1 ≤ 0) [1].

## 4. Algorithm

### 4.1 High-Level Steps

1. Build the k-mer count table with `CountKmers(sequence, k)` (one pass over the L − k + 1 windows).
2. If the table is empty, return the all-zero `KmerStatistics`.
3. In one pass over the table (skipping counts outside [lower, upper]): `TotalKmers` = sum, `UniqueKmers`/`DistinctKmers` = number retained, `SingletonKmers` = number with count 1, `MaxCount`/`MinCount` = extremes; `AverageCount` = Total/Distinct (exact).
4. `Entropy` = canonical `StatisticsHelper.ShannonIndex(retained counts) / ln 2` — the same computation as `SequenceComplexity.CalculateKmerEntropy` and `KmerAnalyzer.CalculateKmerEntropy` (bit-identical when unfiltered); the sequence is counted only once.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `AnalyzeKmers` | O(L·k) | O(D·k) | L − k + 1 windows, each k-mer materialized as a length-k string; D distinct k-mers stored in the count map |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [KmerAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/KmerAnalyzer.cs)

- `KmerAnalyzer.AnalyzeKmers(string, int)`: returns the `KmerStatistics` record.
- `KmerAnalyzer.CountKmers(string, int)`: builds the count table (reused).
- `KmerAnalyzer.AnalyzeKmers(string, int, int lowerCount, int upperCount = int.MaxValue)`: Jellyfish `-L/-U` filtered statistics. The range predicate is the private `SelectByCountRange` helper shared with `FindKmersWithMinCount(string, int, int, int)` (`jellyfish dump -L/-U`), so the k-mers it summarises are exactly those that method lists (KMER-UNIQUE-001).
- `StatisticsHelper.ShannonIndex`: canonical Shannon entropy (nats; ÷ ln 2 for bits), also behind `KmerAnalyzer.CalculateKmerEntropy` → `SequenceComplexity.CalculateKmerEntropy`.

### 5.2 Current Behavior

`AnalyzeKmers` builds the count table once with `CountKmers` and derives every statistic from it (until B06 the entropy re-counted the sequence through `CalculateKmerEntropy`, and `AverageCount` was rounded to two decimals — both fixed, F6). The unit is not a search/matching operation (it aggregates a precomputed count table), so the repository **suffix tree was not used** — there is no occurrence-enumeration or pattern-location subproblem here; the linear count-table scan in `CountKmers` is optimal.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- `TotalKmers = L − k + 1` [1][2].
- `UniqueKmers` = distinct k-mer count [1][2].
- `Entropy = −Σ p(α) log₂ p(α)`, p(α) = mult(α)/(L − k + 1) [3][4].
- Max/Min/Average multiplicity over the count table.

- Jellyfish `stats` fields Unique/Distinct/Total/Max_count and the `-L/-U` filters [5].

**Intentionally simplified:** none.

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

The `UniqueKmers` field name denotes the **distinct** k-mer count — it is *not* Jellyfish's "Unique" (count 1) [5]. It is kept for backward compatibility; use `DistinctKmers` (same value) and `SingletonKmers` (Jellyfish Unique, = size of `FindUniqueKmers`) [2][5]. Jellyfish counts canonical k-mers with `-C` and skips non-ACGT windows; this method is single-strand with no alphabet filtering (see KMER-COUNT-001 / KMER-BOTH-001). No IUPAC-degenerate handling: ambiguous symbols form ordinary k-mers.

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

### 7.2 Reference cross-check (B06, 2026-10-01)

Python `collections.Counter` replica of Jellyfish `compute_stats` + scipy 1.17.1 `entropy(retained, base=2)`:

| Input | k | L/U | Unique | Distinct | Total | Max | Min | Mean | Entropy (bits) |
|---|---|---|---|---|---|---|---|---|---|
| GTAGAGCTGT | 1 | – | 1 | 4 | 10 | 4 | 1 | 2.5 | 1.8464393446710154 |
| GTAGAGCTGT | 2 | – | 5 | 7 | 9 | 2 | 1 | 1.2857142857142858 | 2.7254805569978684 |
| ATCGATCAC | 3 | – | 5 | 6 | 7 | 2 | 1 | 1.1666666666666667 | 2.521640636343318 |
| AAAA | 2 | – | 0 | 1 | 3 | 3 | 3 | 3.0 | 0 |
| ACGTTGCATGTCGCATGATGCATGAGAGCT | 4 | – | 17 | 21 | 27 | 3 | 1 | 1.2857142857142858 | 4.254525464966174 |
| GTAGAGCTGT | 1 | L=2 | 0 | 3 | 9 | 4 | 2 | 3.0 | 1.5304930567574826 |
| GTAGAGCTGT | 1 | U=3 | 1 | 3 | 6 | 3 | 1 | 2.0 | 1.4591479170272446 |
| GTAGAGCTGT | 2 | L=U=2 | 0 | 2 | 4 | 2 | 2 | 2.0 | 1.0 |
| ACGTTGCATGTCGCATGATGCATGAGAGCT | 4 | L=2 | 0 | 4 | 10 | 3 | 2 | 2.5 | 1.970950594454669 |

### 7.3 Related Tests, Evidence, or Documents

- Tests: [KmerAnalyzer_AnalyzeKmers_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Analysis/KmerAnalyzer_AnalyzeKmers_Tests.cs) — covers `INV-01`–`INV-06`
- Evidence: [KMER-STATS-001-Evidence.md](../../../docs/Evidence/KMER-STATS-001-Evidence.md)
- Related algorithms: [Unique_And_MinCount_Kmers](../K-mer/Unique_And_MinCount_Kmers.md), [Both_Strand_Kmer_Counting](../K-mer/Both_Strand_Kmer_Counting.md)

## 8. References

1. Wikipedia contributors. 2026. K-mer. Wikipedia. https://en.wikipedia.org/wiki/K-mer
2. Clavijo, B. 2018. k-mer counting, part I: Introduction. BioInfoLogics. https://bioinfologics.github.io/post/2018/09/17/k-mer-counting-part-i-introduction/
3. Manca, V. et al. 2021. Spectral concepts in genome informational analysis. arXiv preprint. arXiv:2106.15351. https://arxiv.org/abs/2106.15351
4. Entropy–Rank Ratio: A Novel Entropy–Based Perspective for DNA Complexity and Classification. 2025. arXiv preprint. arXiv:2511.05300. https://arxiv.org/html/2511.05300
5. Marçais G, Kingsford C. 2011. A fast, lock-free approach for efficient parallel counting of occurrences of k-mers. Bioinformatics 27:764–770. Jellyfish source `sub_commands/stats_main.cc` (`compute_stats`) and `stats_main_cmdline.yaggo` (field definitions, `-L/-U`), https://github.com/gmarcais/Jellyfish (opened via raw.githubusercontent.com, 2026-10-01).

