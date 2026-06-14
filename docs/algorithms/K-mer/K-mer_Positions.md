# K-mer Positions

| Field | Value |
|-------|-------|
| Algorithm Group | K-mer |
| Test Unit ID | KMER-POSITIONS-001 |
| Related Projects | Seqeron.Genomics.Analysis |
| Implementation Status | Production |
| Last Reviewed | 2026-06-14 |

## 1. Overview

K-mer Positions reports every starting index at which a given k-mer (a fixed pattern) occurs in a sequence. It solves the classical Pattern Matching Problem — "find all occurrences of a pattern in a string" [1] — returning all 0-based start positions, including overlapping occurrences. The result is exact (not heuristic): it enumerates every position `p` where `sequence[p .. p+|kmer|)` equals the k-mer.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

A *k-mer* is a substring of length *k* drawn from a biological sequence [2]. Locating where a specific k-mer occurs underlies motif scanning, repeat localization, primer/probe placement, and read-mapping primitives.

### 2.2 Core Model

Given a text *T* (the sequence) of length *L* and a pattern *P* (the k-mer) of length *k*, the set of occurrences is

`Occ(P, T) = { i ∈ [0, L − k] : T[i .. i+k) = P }` [1]

reported in ascending order using **0-based** indexing [1]. There are at most `L − k + 1` candidate start positions [2]. Occurrences may overlap, and every overlapping start is included [1].

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | Each returned position `p` satisfies `T[p .. p+k) = P` (case-folded). | Direct from the matching predicate [1]. |
| INV-02 | Positions are 0-based and strictly ascending. | Single left-to-right scan; 0-based per spec [1]. |
| INV-03 | The count of returned positions equals the overlapping occurrence count of `P` in `T`. | Every candidate index is tested; overlaps reported [1]. |
| INV-04 | All positions lie in `[0, L − k]`; the result is empty when `k > L`. | Loop bound `i ≤ L − k`; `L − k + 1` candidates [2]. |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| sequence | string | required | Text to search. | Case-insensitive (upper-cased internally); null/empty → empty result. |
| kmer | string | required | Pattern to locate. | Case-insensitive; null/empty → empty result. |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| (return) | IEnumerable&lt;int&gt; | Ascending 0-based start positions of every (possibly overlapping) occurrence. |

### 3.3 Preconditions and Validation

Indexing is **0-based** [1]; the position range is inclusive of `0` and `L − k`. The accepted alphabet is unrestricted text (no DNA/RNA validation); matching is case-insensitive (both arguments are upper-cased with `ToUpperInvariant`). Null/empty `sequence` or `kmer`, and `|kmer| > |sequence|`, all yield an empty sequence — no exception is thrown.

## 4. Algorithm

### 4.1 High-Level Steps

1. If `sequence` or `kmer` is null/empty, return empty.
2. Upper-case both inputs (case-insensitive matching).
3. For each start index `i` from `0` to `L − k`, compare `sequence[i .. i+k)` with the k-mer.
4. Yield `i` for every match, in ascending order.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| FindKmerPositions | O(L·k) | O(L) | Naive forward scan over `L − k + 1` windows, each compared in O(k). The O(L) space is the upper-cased copy of the inputs. A single query against one text does not amortize suffix-tree construction (see §5.2). |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [KmerAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/KmerAnalyzer.cs)

- `KmerAnalyzer.FindKmerPositions(string sequence, string kmer)`: returns ascending 0-based start positions of all overlapping occurrences of `kmer` in `sequence`.

### 5.2 Current Behavior

Implemented as a lazy `IEnumerable<int>` (`yield return`) over a single left-to-right scan, so positions are emitted in ascending order without an explicit sort, and comparison uses `ReadOnlySpan<char>.SequenceEqual` to avoid per-window substring allocation. Both inputs are upper-cased once via `ToUpperInvariant`, making matching case-insensitive (consistent with sibling `KmerAnalyzer.CountKmers`).

**Search-reuse decision (suffix tree evaluated).** The repository `SuffixTree.FindAllOccurrences` ([SuffixTree.Search.cs](../../../src/SuffixTree/Algorithms/SuffixTree/SuffixTree.Search.cs)) was evaluated. It correctly counts overlapping occurrences via leaf collection, but it returns positions in unordered leaf-collection order (the algorithm here requires ascending order) and amortizes only when many patterns are queried against one preprocessed text. For a single k-mer query against one sequence, the O(n) suffix-tree construction is not repaid, and the unordered output would require an extra sort. The naive forward scan is therefore retained; it is simpler, allocation-light, and yields ascending order directly. Correctness (overlapping, 0-based) is unaffected by this choice.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- All 0-based start positions where the pattern occurs as a substring, overlapping occurrences included [1].
- At most `L − k + 1` candidate positions scanned [2].

**Intentionally simplified:**

- (none)

**Not implemented:**

- Approximate / mismatch-tolerant matching; **users should rely on:** an alignment or approximate-search routine (out of scope for exact k-mer location).

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| `kmer` longer than `sequence` | empty | `L − k + 1 ≤ 0`, no candidate positions [2]. |
| `kmer` absent | empty | Only matching starts are reported [1]. |
| `kmer` equals whole `sequence` | `[0]` | One occurrence at index 0 [1]. |
| Self-overlapping (`AA` in `AAAA`) | `[0,1,2]` | Overlapping occurrences all reported [1]. |
| null/empty `sequence` or `kmer` | empty | Repository convention (no spec mandate). |

### 6.2 Limitations

Exact matching only — no mismatches, gaps, or IUPAC-degenerate codes. No DNA/RNA alphabet validation; any characters are matched literally (after case-folding).

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
var positions = KmerAnalyzer.FindKmerPositions("GATATATGCATATACTT", "ATAT").ToList();
// positions == [1, 3, 9]   (Rosalind BA1D sample)
```

**Numerical walk-through:** `GATATATGCATATACTT` (indices 0..16). `ATAT` matches at i=1 (`ATAT`), i=3 (`ATAT`), and i=9 (`ATAT`) → `1 3 9` [1].

### 7.3 Related Tests, Evidence, or Documents

- Tests: [KmerAnalyzer_FindKmerPositions_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/KmerAnalyzer_FindKmerPositions_Tests.cs) — covers `INV-01`–`INV-04`
- Evidence: [KMER-POSITIONS-001-Evidence.md](../../../docs/Evidence/KMER-POSITIONS-001-Evidence.md)

## 8. References

1. Rosalind. 2026 (accessed). Find All Occurrences of a Pattern in a String (Problem BA1D). https://rosalind.info/problems/ba1d/
2. Wikipedia contributors. 2026 (accessed). k-mer. https://en.wikipedia.org/wiki/K-mer
3. Compeau, P., Pevzner, P. 2015. Bioinformatics Algorithms: An Active Learning Approach (Pattern Matching Problem). Active Learning Publishers. https://gerdos.web.elte.hu/edu/bioinformatics_algorithms/week1.pdf
