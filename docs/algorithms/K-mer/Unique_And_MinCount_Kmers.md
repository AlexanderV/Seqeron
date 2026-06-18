# Unique K-mers and K-mers with Minimum Count

| Field | Value |
|-------|-------|
| Algorithm Group | K-mer |
| Test Unit ID | KMER-UNIQUE-001 |
| Related Projects | Seqeron.Genomics.Analysis |
| Implementation Status | Production |
| Last Reviewed | 2026-06-14 |

## 1. Overview

These two operations filter the k-mer frequency spectrum of a DNA/RNA sequence. `FindUniqueKmers` returns the **unique** k-mers — those appearing exactly once (occurrence count = 1) [2]. `FindKmersWithMinCount` returns **recurrent** k-mers — those whose overlapping occurrence count is at least a threshold `minCount` (Count(Text, Pattern) ≥ t) [3] — paired with their counts and ordered by count descending. Both are exact, deterministic, combinatorial operations derived directly from the k-mer definition; neither is heuristic or probabilistic.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

A k-mer is a substring of length k contained within a biological sequence [1]. K-mers are extracted with a sliding window of step 1, so adjacent k-mers overlap by k−1 characters and a sequence of length L yields **L − k + 1** total k-mers [1]. The frequency spectrum distinguishes *total* k-mers (with duplicates), *distinct* k-mers (each different string once), and *unique* k-mers (frequency exactly 1) [2].

### 2.2 Core Model

For a sequence `Text` and k-mer `Pattern`, `Count(Text, Pattern)` is the number of overlapping occurrences of `Pattern` in `Text` [3]. Define the multiset of all length-k substrings and group by string to obtain counts c(P) = Count(Text, P).

- **Unique k-mers:** `{ P : c(P) = 1 }` — "Unique k-mers are those that appear only once" [2].
- **K-mers with minimum count:** `{ (P, c(P)) : c(P) ≥ minCount }`, the recurrent k-mers (Count ≥ t) [3].

Worked example (ATCGATCAC, k=3) [2]: 7 total, 6 distinct, 5 unique = {TCG, CGA, GAT, TCA, CAC}; ATC occurs twice (c=2) so it is distinct but not unique.

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | Each k-mer in `FindUniqueKmers` has c(P) = 1 | Filter predicate is `c(P) = 1` [2] |
| INV-02 | unique ⊆ distinct, so &#124;unique&#124; ≤ &#124;distinct&#124; ≤ L − k + 1 | Total count bound [1]; unique is a frequency-1 subset of distinct [2] |
| INV-03 | Each pair in `FindKmersWithMinCount` has c(P) ≥ minCount and the reported count equals the overlapping occurrence count | Filter predicate `c(P) ≥ t` over exact counts [3] |
| INV-04 | `FindKmersWithMinCount` output is ordered by count non-increasing | `OrderByDescending` on count (recurrent-first ranking [3]) |
| INV-05 | With minCount ≤ 1 the returned keys equal the distinct k-mer set | Every observed k-mer has c(P) ≥ 1 [2][3] |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| sequence | string | required | Sequence to analyze | null/empty → empty result; upper-cased internally (case-insensitive) |
| k | int | required | K-mer length | Must be > 0; k > L → empty result |
| minCount | int | required (`FindKmersWithMinCount`) | Inclusive minimum occurrence threshold | ≤ 1 selects all distinct k-mers |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| (FindUniqueKmers) | IEnumerable\<string\> | K-mers with occurrence count = 1; order unspecified |
| (FindKmersWithMinCount) | IEnumerable\<(string Kmer, int Count)\> | Recurrent k-mers with their counts, ordered by Count descending |

### 3.3 Preconditions and Validation

Input is upper-cased (T/U not normalised; standard string k-mers). Indexing is 0-based over the character string. Null or empty sequence returns an empty result; k > sequence length returns empty (L − k + 1 ≤ 0). k ≤ 0 throws `ArgumentOutOfRangeException`.

## 4. Algorithm

### 4.1 High-Level Steps

1. Build the k-mer count map via `CountKmers` (single O(n) pass, overlapping window, step 1) [1].
2. `FindUniqueKmers`: select keys with count == 1 [2].
3. `FindKmersWithMinCount`: select pairs with count ≥ minCount, order by count descending [3].

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| FindUniqueKmers | O(n·k) | O(d·k) | n = length, d = distinct k-mers; substring extraction is O(k) per position |
| FindKmersWithMinCount | O(n·k + d log d) | O(d·k) | dominated by counting; sort over d distinct k-mers |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [KmerAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/KmerAnalyzer.cs)

- `KmerAnalyzer.FindUniqueKmers(string, int)`: returns k-mers with count = 1.
- `KmerAnalyzer.FindKmersWithMinCount(string, int, int)`: returns (k-mer, count) with count ≥ minCount, ordered by count descending.

### 5.2 Current Behavior

Both methods delegate counting to `KmerAnalyzer.CountKmers`, inheriting its null/empty/k-bounds handling and upper-casing. `FindUniqueKmers` returns results in dictionary enumeration order (unspecified); callers needing a stable order must sort. The repository **suffix tree** was evaluated and not used: both methods need the full frequency spectrum (count of every distinct k-mer), which a single linear hash-map pass over L − k + 1 windows yields in O(n·k); a suffix tree adds construction overhead without benefit for this exhaustive-count workload (it shines for many exact-match queries against one text, not for enumerating all k-mer frequencies).

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Unique k-mers = occurrence count exactly 1 [2].
- Min-count filter = overlapping occurrence count ≥ t [3].
- Overlapping, step-1 counting (L − k + 1 total) [1].
- Count-descending ordering of recurrent k-mers (most-frequent-first) [3].

**Intentionally simplified:**

- (none)

**Not implemented:**

- Canonical (reverse-complement-merged) k-mer counting; users should rely on `KmerAnalyzer.CountKmersBothStrands` for strand-aware counts.

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| null / empty sequence | empty result | No k-mers exist |
| k > L | empty result | L − k + 1 ≤ 0 [1] |
| k ≤ 0 | ArgumentOutOfRangeException | k-mer length must be positive [1] |
| Homopolymer (AAAAA, k=3) | no unique k-mers | AAA has count 3 > 1 [2] |
| minCount ≤ 1 | all distinct k-mers | every k-mer has count ≥ 1 [2][3] |
| minCount > max count | empty | no k-mer meets the threshold [3] |

### 6.2 Limitations

Counts forward-strand string k-mers only (no reverse-complement merging); does not normalise RNA/DNA (T vs U treated as different characters). Order of `FindUniqueKmers` is not guaranteed.

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
var unique = KmerAnalyzer.FindUniqueKmers("ATCGATCAC", 3).ToHashSet();
// { "TCG", "CGA", "GAT", "TCA", "CAC" }  — ATC (count 2) excluded

var recurrent = KmerAnalyzer.FindKmersWithMinCount("ACGTACGT", 4, 2).ToList();
// [ ("ACGT", 2) ]
```

### 7.3 Related Tests, Evidence, or Documents

- Tests: [KmerAnalyzer_FindUniqueAndMinCount_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/KmerAnalyzer_FindUniqueAndMinCount_Tests.cs) — covers `INV-01`..`INV-05`
- Evidence: [KMER-UNIQUE-001-Evidence.md](../../../docs/Evidence/KMER-UNIQUE-001-Evidence.md)

## 8. References

1. Wikipedia contributors. 2026. *K-mer*. Wikipedia. https://en.wikipedia.org/wiki/K-mer
2. Clavijo B, et al. 2018. *k-mer counting, part I: Introduction*. BioInfoLogics. https://bioinfologics.github.io/post/2018/09/17/k-mer-counting-part-i-introduction/
3. Compeau P, Pevzner P. 2015. *Bioinformatics Algorithms: An Active Learning Approach* (2nd ed.). Active Learning Publishers. https://www.amazon.com/BIOINFORMATICS-ALGORITHMS-Phillip-Compeau/dp/0990374637
