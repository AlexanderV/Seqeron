# K-mer Search Algorithms

| Field | Value |
|-------|-------|
| Algorithm Group | K-mer Analysis |
| Test Unit ID | KMER-FIND-001 |
| Related Projects | N/A |
| Implementation Status | N/A |
| Last Reviewed | 2026-04-30 |

## 1. Overview

K-mer search algorithms identify k-mers of special interest within a sequence rather than returning the full count map. In this repository, the documented search surface includes the most frequent k-mers, unique k-mers, and `(L, t)` clumps. All three operations are built on exact k-mer counts, with clump finding using a sliding-window update strategy and a deduplicating result set.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

A most frequent k-mer is any k-mer attaining the maximum count in a sequence, and multiple k-mers may tie for that maximum. A unique k-mer occurs exactly once. A pattern forms an `(L, t)` clump if some window of length `L` contains at least `t` occurrences of that pattern. The original document also notes that clumps can indicate biologically interesting regions such as origins of replication or other motif-rich segments. Sources: Rosalind BA1B, Rosalind BA1E, Wikipedia (K-mer).

### 2.2 Core Model

`FindMostFrequentKmers(...)` and `FindUniqueKmers(...)` are filters over exact k-mer counts:

$$
MostFrequent = \{kmer : Count(kmer) = \max_j Count(j)\}
$$

$$
Unique = \{kmer : Count(kmer) = 1\}
$$

For clumps, the algorithm maintains a mutable count map over a sliding window of length `L` and adds any k-mer whose count reaches or exceeds `t` in any visited window.

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | Every k-mer returned by `FindUniqueKmers(...)` has count exactly `1` | The method filters the count dictionary by `kvp.Value == 1` |
| INV-02 | Every k-mer returned by `FindMostFrequentKmers(...)` has the maximum observed count | The method filters by `counts.Values.Max()` |
| INV-03 | `FindClumps(...)` returns each qualifying k-mer at most once | Results are stored in a `HashSet<string>` |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `sequence` | `string` | required | Sequence to search | Empty or null input yields an empty result |
| `k` | `int` | required | K-mer length | `FindMostFrequentKmers(...)` and `FindUniqueKmers(...)` inherit `CountKmers(...)` validation; `FindClumps(...)` returns empty when `k <= 0` |
| `windowSize` | `int` | required for clumps | Window length `L` for clump detection | `FindClumps(...)` returns empty when `windowSize < k` or `windowSize > sequence.Length` |
| `minOccurrences` | `int` | required for clumps | Minimum number of occurrences `t` in a window | `FindClumps(...)` returns empty when `minOccurrences <= 0` |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `mostFrequent` | `IEnumerable<string>` | All k-mers tied for the maximum observed count |
| `unique` | `IEnumerable<string>` | All k-mers observed exactly once |
| `clumps` | `IEnumerable<string>` | Unique set of k-mers that satisfy the `(L, t)` clump condition |

### 3.3 Preconditions and Validation

`FindMostFrequentKmers(...)` and `FindUniqueKmers(...)` use `CountKmers(...)`, so `k <= 0` raises `ArgumentOutOfRangeException` there. `FindClumps(...)` instead treats invalid parameters as empty-result conditions: null or empty sequence, `k <= 0`, `windowSize < k`, `windowSize > sequence.Length`, or `minOccurrences <= 0` all yield no clumps. All string-based searches uppercase the sequence internally.

## 4. Algorithm

### 4.1 High-Level Steps

1. Normalize the input sequence to uppercase.
2. For most-frequent and unique searches, compute the full k-mer count map.
3. Return either the keys tied at the maximum count or the keys with count `1`.
4. For clumps, initialize counts in the first window, record k-mers meeting the threshold, then slide the window by one position while updating counts incrementally.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `FindMostFrequentKmers` | `O(n)` | `O(u)` | Builds exact counts and filters maxima |
| `FindUniqueKmers` | `O(n)` | `O(u)` | Builds exact counts and filters singletons |
| `FindClumps` | `O(n × (L - k + 1))` worst case | `O(u)` | The original document describes this bound and notes that it is typically near-linear with efficient data structures |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [KmerAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/KmerAnalyzer.cs)

- `KmerAnalyzer.FindMostFrequentKmers(string, int)`: Returns all maxima from the count map.
- `KmerAnalyzer.FindUniqueKmers(string, int)`: Returns all singleton k-mers.
- `KmerAnalyzer.FindClumps(string, int, int, int)`: Returns deduplicated clump-forming k-mers.

### 5.2 Current Behavior

All three methods uppercase the input sequence. `FindMostFrequentKmers(...)` and `FindUniqueKmers(...)` reuse the exact counting surface, while `FindClumps(...)` manages a mutable per-window count dictionary and a `HashSet<string>` of discovered clumps. `FindClumps(...)` returns empty rather than throwing on invalid window or threshold parameters.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Identification of all k-mers tied at the maximum count.
- Identification of singleton k-mers.
- Sliding-window `(L, t)` clump detection.

**Intentionally simplified:**

- `FindClumps(...)` returns only the set of qualifying k-mers, not the windows in which they qualified; **consequence:** callers learn which patterns form clumps but not where each supporting window occurred.

**Not implemented:**

- Reporting the full set of clump-supporting windows or multiplicity traces; **users should rely on:** downstream custom analysis if those details are needed.

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Empty sequence | Returns an empty result | No valid windows exist |
| `k <= 0` | Throws for count-based searches; returns empty for `FindClumps(...)` | Different validation paths in source |
| `k > sequence.Length` | Returns an empty result | No valid k-mers exist |
| `windowSize > sequence.Length` | `FindClumps(...)` returns empty | No full window exists |
| `windowSize < k` | `FindClumps(...)` returns empty | A window cannot contain a full k-mer |
| All k-mers equally frequent | `FindMostFrequentKmers(...)` returns all observed k-mers | Every k-mer is tied at the maximum count |

### 6.2 Limitations

The current implementation reports only pattern-level results, not locations for all qualifying clump windows. As with the underlying count-based helpers, memory usage still depends on the number of unique k-mers maintained in dictionaries and sets.

## 7. Examples and Related Material

### 7.1 Worked Example

**Numerical / biological walk-through (optional):**

The original document cites `ACGTTGCATGTCGCATGATGCATGAGAGCT` as a standard frequent-word example in which the most frequent 4-mers are `CATG` and `GCAT`, each appearing 3 times. It also gives `TGCA` as a `(25, 3)` clump in:

```text
gatcagcataagggtcccTGCAATGCATGACAAGCCTGCAgttgttttac
```

### 7.2 Applications and Use Cases (Optional)

- Unique k-mers for marker discovery and genomic fingerprinting.
- Clump finding for motif-rich regions such as origins of replication.

## 8. References

1. Rosalind BA1B - Find the Most Frequent Words in a String. https://rosalind.info/problems/ba1b/
2. Rosalind BA1E - Find Patterns Forming Clumps in a String. https://rosalind.info/problems/ba1e/
3. Wikipedia (K-mer). https://en.wikipedia.org/wiki/K-mer
