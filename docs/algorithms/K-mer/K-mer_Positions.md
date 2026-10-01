# K-mer Positions

| Field | Value |
|-------|-------|
| Algorithm Group | K-mer |
| Test Unit ID | KMER-POSITIONS-001 |
| Related Projects | Seqeron.Genomics.Analysis |
| Implementation Status | Production |
| Last Reviewed | 2026-10-01 |

## 1. Overview

K-mer Positions reports every starting index at which a given k-mer (a fixed pattern) occurs in a sequence. It solves the classical Pattern Matching Problem — "find all occurrences of a pattern in a string" [1] — returning all 0-based start positions, including overlapping occurrences. The result is exact (not heuristic): it enumerates every position `p` where `sequence[p .. p+|kmer|)` equals the k-mer.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

A *k-mer* is a substring of length *k* drawn from a biological sequence [2]. Locating where a specific k-mer occurs underlies motif scanning, repeat localization, primer/probe placement, and read-mapping primitives.

### 2.2 Core Model

Given a text *T* (the sequence) of length *L* and a pattern *P* (the k-mer) of length *k*, the set of occurrences is

`Occ(P, T) = { i ∈ [0, L − k] : T[i .. i+k) = P }` [1]

reported in ascending order using **0-based** indexing [1]. There are at most `L − k + 1` candidate start positions [2]. Occurrences may overlap, and every overlapping start is included [1].

**Indexing conventions.** Rosalind SUBS ("Finding a Motif in DNA") asks for the same set in **1-based** positions (sample `GATATATGCATATACTT` / `ATAT` → `2 4 10`) [4]; BA1D and this method use 0-based (`1 3 9`). Add 1 to convert. Biopython `Bio.SeqUtils.nt_search` (forward strand, 0-based) and Python `re.finditer("(?=P)")` return the 0-based set [5].

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

1. If `sequence` or `kmer` is null/empty, or `k > L`, return empty.
2. Upper-case both inputs (case-insensitive matching).
3. Build the KMP prefix function π of the k-mer (COMPUTE-PREFIX-FUNCTION) [6][7].
4. Read the text once left to right, keeping `q` = number of matched pattern characters; on mismatch fall back `q ← π[q−1]`; when `q = k`, yield `i − k + 1` and set `q ← π[k−1]` so overlapping occurrences are found (KMP-MATCHER) [6][7].

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| FindKmerPositions | O(L + k) | O(L + k) | Knuth–Morris–Pratt: O(k) prefix function + one O(L) pass over the text (amortized; each fallback undoes an earlier advance) [6][7]. Space: upper-cased copies of the inputs + π. Before the 2026-09 review the method compared every window (Θ((L−k+1)·k) worst case: `A`^10⁶ vs `A`^99 999`C` took 5.35 s; KMP 9–26 ms). |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [KmerAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/KmerAnalyzer.cs)

- `KmerAnalyzer.FindKmerPositions(string sequence, string kmer)`: returns ascending 0-based start positions of all overlapping occurrences of `kmer` in `sequence`.

### 5.2 Current Behavior

Implemented as a lazy `IEnumerable<int>` (`yield return`) over a single left-to-right Knuth–Morris–Pratt scan [6][7], so positions are emitted in ascending order without an explicit sort and no per-window substring is allocated. Both inputs are upper-cased once via `ToUpperInvariant`, making matching case-insensitive (consistent with sibling `KmerAnalyzer.CountKmers`). No other single-query linear-time exact matcher exists in the repository (searched for KMP / Z / Boyer–Moore / Aho–Corasick), so the private `ComputeKmpPrefixFunction` is not a duplicate.

**Search-reuse decision (suffix tree evaluated).** The repository `SuffixTree.FindAllOccurrences` ([SuffixTree.Search.cs](../../../src/SuffixTree/Algorithms/SuffixTree/SuffixTree.Search.cs)) was evaluated. It correctly counts overlapping occurrences via leaf collection, but it returns positions in unordered leaf-collection order (the algorithm here requires ascending order) and amortizes only when many patterns are queried against one preprocessed text. For a single k-mer query against one sequence, the O(n) suffix-tree construction is not repaid (measured 2026-10-01, Release: random 10 Mbp DNA, build + one query + sort = 23.3 s and ≈1.7 GB, vs 71 ms for the KMP scan), and the unordered output would require an extra sort. `MotifFinder.FindExactMotif` / `GenomicAnalyzer.FindMotif` take a `DnaSequence` (validated alphabet, cached tree), whereas this method accepts any text, so it cannot delegate to them. Correctness (overlapping, 0-based) is unaffected by this choice.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- All 0-based start positions where the pattern occurs as a substring, overlapping occurrences included [1].
- At most `L − k + 1` candidate positions reported [2].
- Linear-time exact matching (KMP) [6][7].

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

### 7.2 Reference cross-check (2026-10-01)

| Sequence | k-mer | Python `re` lookahead / Biopython 1.88 `nt_search` (0-based) | Rosalind (1-based) | This method |
|---|---|---|---|---|
| GATATATGCATATACTT | ATAT | 1 3 9 | SUBS/BA1D: 2 4 10 | 1 3 9 |
| ACGTACGTACGTACGT | GTA | 2 6 10 | SUBS (legacy sample): 3 7 11 | 2 6 10 |
| AAAAAAAAAA | AAA | 0 1 2 3 4 5 6 7 | — | same |
| ACGTACGTACGT | CGTA | 1 5 | — | same |
| ABABZABABYABABX | ABABX | 10 | — | same |
| AABAACAADAABAABA | AABA | 0 9 12 | — | same |

5000 random cases (alphabets AC / ACGT / ACGTN / mixed case / A; L ≤ 60, k ≤ 8; 35 249 positions) vs Python `re.finditer("(?=P)")` on upper-cased input: 0 mismatches.

### 7.3 Related Tests, Evidence, or Documents

- Tests: [KmerAnalyzer_FindKmerPositions_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Analysis/KmerAnalyzer_FindKmerPositions_Tests.cs) — covers `INV-01`–`INV-04`
- Evidence: [KMER-POSITIONS-001-Evidence.md](../../../docs/Evidence/KMER-POSITIONS-001-Evidence.md)

## 8. References

1. Rosalind. 2026 (accessed). Find All Occurrences of a Pattern in a String (Problem BA1D). https://rosalind.info/problems/ba1d/
2. Wikipedia contributors. 2026 (accessed). k-mer. https://en.wikipedia.org/wiki/K-mer
3. Compeau, P., Pevzner, P. 2015. Bioinformatics Algorithms: An Active Learning Approach (Pattern Matching Problem). Active Learning Publishers. https://gerdos.web.elte.hu/edu/bioinformatics_algorithms/week1.pdf
4. Rosalind. Finding a Motif in DNA (Problem SUBS), 1-based positions. https://rosalind.info/problems/subs/ (statement/samples read from GitHub mirrors: breezedu/rosalind `FindingaMotifinDNA.java`; mtarbit/Rosalind-Problems `e009-subs.py`).
5. Cock, P.J.A. et al. 2009. Biopython. *Bioinformatics* 25:1422–1423. `Bio.SeqUtils.nt_search` (v1.88 source inspected).
6. Knuth, D.E., Morris, J.H., Pratt, V.R. 1977. Fast pattern matching in strings. *SIAM J. Comput.* 6(2):323–350.
7. Cormen, T.H. et al. *Introduction to Algorithms*, §32.4 (KMP-MATCHER, COMPUTE-PREFIX-FUNCTION); failure-array reference code: TheAlgorithms/Python `strings/knuth_morris_pratt.py`.
