# Exact Pattern Search (Suffix Tree)

| Field | Value |
|-------|-------|
| Algorithm Group | Pattern Matching |
| Test Unit ID | PAT-EXACT-001 |
| Related Projects | N/A |
| Implementation Status | N/A |
| Last Reviewed | 2026-04-30 |

## 1. Overview

Exact pattern matching finds every occurrence of a pattern `P` within a text `T`. In this repository, the core implementation uses a suffix tree and exposes existence checks, occurrence counting, and occurrence listing. Higher-level genomics wrappers normalize motifs to uppercase and, in one case, sort the returned positions.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

In a suffix tree for text `T`, every suffix of `T` is represented as a path from the root to a leaf. A pattern `P` occurs at position `i` if and only if `P` is a prefix of the suffix beginning at `i`. Exact pattern matching on a suffix tree therefore reduces to following the pattern across tree edges and collecting the leaves below the matched locus. Sources: Gusfield (1997), Ukkonen (1995), Wikipedia (Suffix tree), Rosalind SUBS.

### 2.2 Core Model

The search logic is:

1. Start at the root of the suffix tree.
2. Match the pattern against edge labels.
3. If a mismatch occurs before the pattern is exhausted, the pattern is absent.
4. If the pattern is exhausted, every leaf below the current match point corresponds to an occurrence.

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | Every reported position `i` satisfies `T[i..i+m-1] = P` | Leaves are collected only after a full edge-by-edge match |
| INV-02 | `CountOccurrences(P) = FindAllOccurrences(P).Count` on the core suffix-tree API | `CountOccurrences(...)` returns the matched node's leaf count |
| INV-03 | `Contains(P)` is equivalent to at least one exact occurrence | `Contains(...)` is a matched/not-matched specialization of the same traversal |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `pattern` | `string` or `ReadOnlySpan<char>` | required | Pattern to search in an already built suffix tree | Null string input throws `ArgumentNullException` |
| `sequence` | `DnaSequence` | required for genomics wrappers | DNA sequence whose cached suffix tree is searched | Wrappers return empty for null or empty motif |
| `motif` | `string` | required for genomics wrappers | Exact motif to search | Uppercased before wrapper-level search |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `found` | `bool` | Whether the pattern occurs in the indexed text |
| `count` | `int` | Number of occurrences |
| `positions` | `IReadOnlyList<int>` or `IEnumerable<int>` | Zero-based occurrence positions |

### 3.3 Preconditions and Validation

Core suffix-tree string APIs throw `ArgumentNullException` on null pattern input. Empty-string patterns return all valid start positions `[0..n-1]` from the core tree and a count equal to text length. `MotifFinder.FindExactMotif(...)` and `GenomicAnalyzer.FindMotif(...)` both return empty for null or empty motifs and uppercase the motif before calling the suffix tree.

## 4. Algorithm

### 4.1 High-Level Steps

1. Traverse the suffix tree edge by edge using the pattern characters.
2. Stop immediately on the first mismatch.
3. If the pattern fully matches, collect or count the leaves under the matched node.
4. In genomics wrappers, uppercase the motif before searching.
5. In `MotifFinder.FindExactMotif(...)`, sort the resulting positions before yielding them.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

The core suffix-tree implementation uses:

- Hybrid SIMD comparisons for edge fragments of length at least 8.
- Scalar comparisons for shorter edge fragments.
- Precomputed `LeafCount` values so `CountOccurrences(...)` is `O(m)` after traversal.
- Thread-static buffers to reduce repeated allocations during traversal and leaf collection.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| Build suffix tree | `O(n)` | `O(n)` | Ukkonen-style suffix-tree construction as cited in the original document |
| `Contains` | `O(m)` | `O(1)` auxiliary | Stops after pattern traversal |
| `CountOccurrences` | `O(m)` | `O(1)` auxiliary | Uses precomputed leaf counts |
| `FindAllOccurrences` | `O(m + z)` | `O(z)` | `z` is the number of matches collected from leaves |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [SuffixTree.Search.cs](../../../src/SuffixTree/Algorithms/SuffixTree/SuffixTree.Search.cs), [MotifFinder.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/MotifFinder.cs), [GenomicAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/GenomicAnalyzer.cs)

- `SuffixTree.Contains(string|ReadOnlySpan<char>)`: Exact existence check.
- `SuffixTree.FindAllOccurrences(string|ReadOnlySpan<char>)`: Returns all occurrence positions.
- `SuffixTree.CountOccurrences(string|ReadOnlySpan<char>)`: Returns the occurrence count via `LeafCount`.
- `MotifFinder.FindExactMotif(DnaSequence, string)`: Uppercases the motif and yields sorted positions.
- `GenomicAnalyzer.FindMotif(DnaSequence, string)`: Uppercases the motif and returns the suffix-tree positions directly.

### 5.2 Current Behavior

The core suffix-tree API returns all valid start positions for an empty pattern, returns `true` from `Contains(...)` on an empty pattern, and returns the text length from `CountOccurrences(...)` on an empty pattern. `MotifFinder.FindExactMotif(...)` sorts the positions before yielding them, while `GenomicAnalyzer.FindMotif(...)` returns the underlying suffix-tree list directly. Both wrappers normalize motifs to uppercase for case-insensitive DNA matching.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Suffix-tree traversal for exact pattern matching.
- Occurrence counting from leaf counts.
- Full occurrence collection from leaves under the matched node.

**Intentionally simplified:**

- Wrapper-level exact motif search is restricted to exact uppercase-normalized DNA motifs; **consequence:** approximate or ambiguity-code matching is out of scope for these entry points.

**Not implemented:**

- Approximate or degenerate matching in the exact-pattern search surface; **users should rely on:** `Approximate_Matching_Hamming.md`, `Edit_Distance.md`, or `IUPAC_Degenerate_Matching.md` for those cases.

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Null string pattern | Throws `ArgumentNullException` | Explicit core-string API contract |
| Empty pattern on the core tree | Returns all positions `[0..n-1]` | Shared empty-pattern contract in source |
| Empty motif in genomics wrappers | Returns an empty result | Wrapper-level guard |
| Pattern longer than the text | Returns no matches | No full edge traversal is possible |

### 6.2 Limitations

The documented surface is exact-only. Ambiguity-code matching and approximate alignment are separate algorithms in the repository, and wrapper-level result ordering differs between `MotifFinder` and `GenomicAnalyzer`.

## 7. Examples and Related Material

### 7.1 Worked Example

**Numerical / biological walk-through (optional):**

Examples preserved from the original document:

| Text | Pattern | Occurrences (0-indexed) |
|------|---------|-------------------------|
| `banana` | `ana` | `[1, 3]` |
| `mississippi` | `issi` | `[1, 4]` |
| `GATATATGCATATACTT` | `ATAT` | `[1, 3, 9]` |

Related documentation: [Suffix_Tree.md](Suffix_Tree.md)

## 8. References

1. Gusfield, D. (1997). *Algorithms on Strings, Trees and Sequences: Computer Science and Computational Biology*. Cambridge University Press.
2. Ukkonen, E. (1995). "On-line construction of suffix trees." *Algorithmica*, 14(3), 249-260.
3. Wikipedia contributors. "Suffix tree." *Wikipedia, The Free Encyclopedia*. https://en.wikipedia.org/wiki/Suffix_tree
4. Rosalind Team. "Finding a Motif in DNA." *Rosalind Bioinformatics*. https://rosalind.info/problems/subs/
