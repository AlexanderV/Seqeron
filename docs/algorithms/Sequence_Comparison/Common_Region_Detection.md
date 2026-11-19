# Common Region Detection (Longest Common Substring)

| Field | Value |
|-------|-------|
| Algorithm Group | Analysis / Sequence Comparison |
| Test Unit ID | GENOMIC-COMMON-001 |
| Related Projects | Seqeron.Genomics.Analysis, SuffixTree |
| Implementation Status | Production |
| Last Reviewed | 2026-11-18 (B09 review) |

## 1. Overview

Common region detection finds the longest contiguous stretch of DNA shared by two sequences — the **longest common substring (LCS)** — and, optionally, every shared contiguous substring above a minimum length. It is an exact, deterministic, combinatorial algorithm (no scoring, no gaps): the answer is a *substring*, not a gapped subsequence [3]. It is used to locate conserved blocks, exact homology anchors, and shared repeats between two sequences. The implementation builds a generalized suffix tree and reads the deepest node common to both strings [1][4].

## 2. Scientific / Formal Basis

### 2.1 Domain Context

Given two sequences over the DNA alphabet {A,C,G,T}, a *common substring* is a contiguous string that appears in both. The longest such string is the longest common substring. This differs from the longest common *subsequence*, which permits gaps: "Unlike the longest common subsequence problem … the longest common substring problem seeks a contiguous substring shared by both texts." [3]

### 2.2 Core Model

Formal definition: "Given two strings, S of length m and T of length n, find a longest string which is substring of both S and T." [3]

Suffix-tree solution: build a generalized suffix tree over both strings; the LCS is the path label of the deepest internal node whose subtree contains a leaf from *each* string [1][4]. As GeeksforGeeks states: "The path label from root to the deepest node marked as XY will give the LCS of X and Y." [4] Lengths and starting positions of the longest common substrings can be found "in Θ(n + m) time with the help of a generalized suffix tree." [3] (attributed there to Gusfield 1997 [1]).

When several distinct substrings share the maximal length, multiple answers exist — e.g. "BADANAT" and "CANADAS" "share the maximal-length substrings 'ADA' and 'ANA'." [3] A deterministic representative must therefore be chosen (see INV-03).

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | The returned substring is contiguous and occurs in both sequences at the reported 0-based positions | LCS is a substring of both [3]; positions come from suffix-tree leaves |
| INV-02 | No common contiguous substring strictly longer than the returned one exists | "a longest string which is substring of both" [3] |
| INV-03 | On a length tie, the substring whose first match in `sequence2` ends earliest is returned (deterministic); both positions are the **first (lowest) 0-based occurrence** in their sequence, independent of the suffix-tree layout | documented tie-break of `SuffixTree.FindAllLongestCommonSubstrings` [7]; positions = `PositionsInText[0]` / `PositionsInOther[0]` (ascending lists) |
| INV-04 | Every region from `FindCommonRegions(minLength)` is contiguous, has length ≥ max(1, minLength), occurs in `sequence1` (first occurrence reported) and equals the matching statistic ms(i) = the longest prefix of `sequence2[i..]` occurring in `sequence1` | Chang & Lawler 1994 matching statistics (Gusfield 1997 ch. 7) |
| INV-05 | Empty input or no shared character → `CommonRegion.None` (empty, length 0, positions −1) | only the empty string qualifies → length-0 LCS [3] |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| sequence1 | DnaSequence | required | first sequence (suffix tree built over it) | alphabet A/C/G/T (DnaSequence-validated), 0-based |
| sequence2 | DnaSequence | required | second sequence | alphabet A/C/G/T, 0-based |
| minLength | int | required (`FindCommonRegions` only) | minimum region length | values < 1 treated as 1 |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `CommonRegion.Sequence` | string | the shared substring |
| `CommonRegion.PositionInFirst` | int | 0-based start in `sequence1` (−1 when none) |
| `CommonRegion.PositionInSecond` | int | 0-based start in `sequence2` (−1 when none) |
| `CommonRegion.Length` | int | substring length |
| `CommonRegion.IsEmpty` | bool | true when no shared substring |

### 3.3 Preconditions and Validation

Inputs are `DnaSequence` (case-normalized to upper, alphabet A/C/G/T enforced at construction). Coordinates are 0-based. Empty sequences are accepted and yield `CommonRegion.None` (for `FindLongestCommonRegion`) or an empty enumeration (`FindCommonRegions`). `null` arguments throw `ArgumentNullException` eagerly (also for the lazily enumerated `FindCommonRegions`). Lower-case input is folded to upper case by `DnaSequence` (matching is case-insensitive); IUPAC ambiguity codes (`N`, …), `U` and any other non-ACGT symbol are rejected by the `DnaSequence` constructor and never reach the algorithm.

## 4. Algorithm

### 4.1 High-Level Steps

1. Build the suffix tree of `sequence1` (cached on the `DnaSequence`).
2. **FindLongestCommonRegion:** call `FindAllLongestCommonSubstrings(sequence2)`, which streams `sequence2` through the tree with suffix links (matching statistics) and tracks the deepest common match; return the substring with the first entries of its (ascending) position lists, or `CommonRegion.None` if empty.
3. **FindCommonRegions:** compute the matching statistics ms(i) of `sequence2` against the tree of `sequence1`:
   - *Short matches* (ms(i) < K, K = ⌈log₄(n·m)⌉ + 2): ms(i) ≥ ms(i−1) − 1, so the match is carried over from the previous start and extended one base at a time with `Contains` probes of at most K characters (O(m·K) in total).
   - *Long matches* (ms(i) ≥ max(K, minLength)): from the canonical suffix-link maximal-exact-match stream `FindMaximalExactMatches(sequence2, max(K, minLength))` (MUMmer `-maxmatch`): the longest match at `i` lies on a diagonal whose maximal run is a MEM, hence ms(i) = max over MEMs with query start ≤ i of (start + length − i). Between unrelated random sequences the expected number of MEMs of length ≥ K is O(1).
   - Report each distinct `sequence2[i..i+ms(i)]` (ms(i) ≥ max(1, minLength)) once, in order of first start, with the lowest occurrence in `sequence1` and start `i`.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| FindLongestCommonRegion | O(n + m) | O(n) | suffix-tree construction + suffix-link streaming [3] |
| FindCommonRegions | O(n + m·K + R) + Σ\|region\| | O(n + m + R) | K = ⌈log₄(n·m)⌉ + 2; R = number of maximal matches of length ≥ max(K, minLength); Σ\|region\| is the size of the reported strings. 1 Mbp random pair: ≈ 4 s (Debug); the former per-start binary search over `Contains` copied up to m characters per probe, O(m²), ≈ 2.6 s at 100 kbp growing quadratically |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [GenomicAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/GenomicAnalyzer.cs)

- `GenomicAnalyzer.FindLongestCommonRegion(seq1, seq2)`: longest common substring with positions.
- `GenomicAnalyzer.FindCommonRegions(seq1, seq2, minLength)`: for each start position in `sequence2`, the single longest common substring of length ≥ max(1, minLength) (right-maximal matches keyed by start position, deduplicated) — **not** every common substring (a shorter prefix sharing a start position with a longer match is omitted).
- Underlying LCS: [SuffixTree.Algorithms.cs](../../../src/SuffixTree/Algorithms/SuffixTree/SuffixTree.Algorithms.cs) `LongestCommonSubstringInfo`.

### 5.2 Current Behavior

**Suffix-tree reuse (decision):** the repository's generalized-suffix-tree LCS (`DnaSequence.SuffixTree.FindAllLongestCommonSubstrings`) is used directly for `FindLongestCommonRegion`, giving the Θ(n+m) result of the cited theory [1][3] rather than an O(n·m) DP scan; this is the algorithmically appropriate use of the suffix tree (single text, streamed query). `FindCommonRegions` computes the matching statistics of `sequence2` with the suffix tree of `sequence1` (short matches by `Contains` extension, long matches from the canonical MEM stream — see §4.1), never copying a substring longer than the reported region. Tie-break for `FindLongestCommonRegion` is "first found in `sequence2`" as documented on the suffix-tree method [7]; the reported positions are the lowest occurrences (the suffix tree's own `LongestCommonSubstringInfo` / `FindAllOccurrences` return one occurrence in internal DFS order, which B09 review F15 showed is not the first).

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Longest common *substring* (contiguous) of two strings via generalized suffix tree, Θ(n+m) [1][3][4].
- Deterministic representative on ties (first-in-`sequence2`) [7].
- Right-maximal common substrings of length ≥ minLength, one per start position in `sequence2`, deduplicated (`FindCommonRegions`): the matching statistics of Chang & Lawler (1994) / Gusfield (1997) ch. 7 reported as strings. This is a derived helper, not a notion taken verbatim from [1][3][4]; it is **not** the full set of every common substring (see 6.2).

**Intentionally simplified:**

- `FindCommonRegions` reports only the longest match per start position in `sequence2` (right-maximal), deduplicated — **not** every common substring. For example, for `ACGTACGT` vs `TTACGTGG` with minLength 3 it returns `{TACGT, ACGT, CGT}`; the prefixes `TAC`, `TACG`, `ACG` (which are also common substrings of length ≥ 3) are intentionally omitted. This keeps the output to a useful set of anchor regions instead of the O(n²) set of all substrings; **consequence:** callers needing the exhaustive set must expand each region's prefixes themselves.

**Variants (audit A2-30 / A2-35; the former "not implemented" list):**

- **Modes** (`FindCommonRegions(seq1, seq2, minLength, CommonRegionMode, StrandSelection)`): `RightMaximal` (the default above), `MaximalExactMatches` (every maximal exact match, one entry per pair of occurrences; `mummer -maxmatch`, canonical `SuffixTree.FindMaximalExactMatches`), `MaximalUniqueReference` (string unique in `sequence1`; `-mumreference`) and `MaximalUnique` (unique in both; `-mum`, Delcher et al. 1999) — the B05 suffix-tree engines, not re-implemented.
- **Strands** (`StrandSelection.ReverseComplement` / `Both`; `mummer -r` / `-b`): the second sequence is reverse-complemented and the same engine run; a reverse region has `IsReverseComplement = true` and `PositionInSecond` is the **forward-strand** start of the window that spells the reverse complement of `Sequence` (the `-c` convention). `FindLongestCommonRegion(seq1, seq2, strands)` returns the longer strand (tie: forward).
- **k mismatches** (`FindLongestCommonRegionWithMismatches(seq1, seq2, k, strands)`): the longest pair of equal-length windows with Hamming distance ≤ k (Flouri et al. 2015; Thankachan et al. 2016), by a two-pointer scan of every diagonal, O(n·m). Ties: smallest start in `sequence2`, then in `sequence1`; k = 0 is `FindLongestCommonRegion`. A symbol other than A/C/G/T is a mismatch against everything.
- **N sequences** (`FindLongestCommonRegion(IReadOnlyList<…> sequences, minSupport)`): the longest substring shared by ≥ `minSupport` (2 … N, default N) sequences, every tie, with its first occurrence in each sequence (−1 where absent) — the canonical generalized suffix tree `SuffixTree.FindLongestCommonSubstrings` (Gusfield 1997 §7.6; Rosalind LCSM).
- **String / IUPAC input** (`FindLongestCommonRegion(string, string)`, `FindCommonRegions(string, string, …)`, …): only A/C/G/T match; every maximal run of another symbol (N, IUPAC, U, gap) becomes a unique separator of the suffix tree (query side: a symbol that equals nothing), so no region contains it and an ordinary region next to an N run is unchanged. On an A/C/G/T string the result equals the `DnaSequence` overload.
- Verification: `mummer -maxmatch / -mumreference / -mum -b -l L` 3.23 on 4716 (reference, query, L, mode) rows == an independent Python brute force (all ACGT rows; the `-c` position convention asserted), the k-mismatch LCS on 1680 brute-force rows, the N-sequence LCS on 694 brute-force rows (`scripts/analysis/generate_common_region_strand_oracle.py`, `GenomicAnalyzer_CommonRegionStrands_Oracle_Tests`).

**Not implemented:**

- Gapped / scored common regions (indels); **users should rely on:** the alignment algorithms in `Seqeron.Genomics.Alignment` (`LocalAlign`, `nucmer`-style anchoring in `AnchorBasedAligner`).

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Empty sequence1 or sequence2 | `CommonRegion.None` / empty enumeration | only empty string qualifies → length 0 [3] |
| No shared character | `CommonRegion.None` | length-0 LCS [3] |
| Identical sequences | whole sequence, positions 0/0 | a string is a substring of itself [3] |
| Length tie (e.g. two distinct len-3 substrings) | first-in-`sequence2` representative | documented tie-break [7] |
| minLength < 1 | treated as 1 | a region is non-empty [3] |

### 6.2 Limitations

Exact matching only — no mismatches, gaps, or scoring. Operates on the DNA alphabet enforced by `DnaSequence`. `FindCommonRegions` reports the longest substring per start position in `sequence2`; shorter substrings nested inside a longer reported one are not separately reported unless they begin at a different `sequence2` position.

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
var a = new DnaSequence("ACGTACGT");
var b = new DnaSequence("TTACGTGG");
CommonRegion r = GenomicAnalyzer.FindLongestCommonRegion(a, b);
// r.Sequence == "TACGT", r.Length == 5, r.PositionInFirst == 3, r.PositionInSecond == 1
```

**Numerical walk-through:** `ACGTACGT` vs `TTACGTGG`. The shared contiguous runs are `TACGT` (length 5, at index 3 of the first, index 1 of the second), which is maximal: no length-6 substring of the first appears in the second. Reference length analogue from GeeksforGeeks: "xabxac"/"abcabxabcd" → "abxa", length 4 [4].

### 7.3 Related Tests, Evidence, or Documents

- Tests: [GenomicAnalyzer_FindCommonRegion_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Analysis/GenomicAnalyzer_FindCommonRegion_Tests.cs) — covers INV-01..INV-05
- Evidence: [GENOMIC-COMMON-001-Evidence.md](../../../docs/Evidence/GENOMIC-COMMON-001-Evidence.md)
- Oracle: [generate_common_region_oracle.py](../../../scripts/analysis/generate_common_region_oracle.py) (independent DP + suffix automaton + brute-force matching statistics; 283 pairs)
- Related algorithms: [Repeat_Detection](../Repeat_Analysis/Repeat_Detection.md)

### 7.4 Review notes (B09, 2026-11-18)

- F15: `PositionInFirst` was one occurrence in suffix-tree DFS order, not the first (e.g. `GA` in `AGGAAGAAGGAGGCCGATTACT` reported 15, first occurrence 2).
- F16: `FindCommonRegions` was quadratic (substring copies of up to m characters in a per-start binary search); now matching statistics, linear-ish (1 Mbp pair in seconds).

## 8. References

1. Gusfield, Dan. 1997. Algorithms on Strings, Trees and Sequences: Computer Science and Computational Biology. Cambridge University Press. ISBN 0-521-58519-8. (citation retrieved from https://en.wikipedia.org/wiki/Longest_common_substring on 2026-06-13)
2. Charalampopoulos P., Kociumaka T., Pissis S.P., Radoszewski J. 2021. Faster Algorithms for Longest Common Substring. ESA 2021, LIPIcs 204, Schloss Dagstuhl. https://doi.org/10.4230/LIPIcs.ESA.2021.30
3. Wikipedia contributors. Longest common substring. https://en.wikipedia.org/wiki/Longest_common_substring (accessed 2026-06-13)
4. GeeksforGeeks. Suffix Tree Application 5 — Longest Common Substring. https://www.geeksforgeeks.org/dsa/suffix-tree-application-5-longest-common-substring-2/ (accessed 2026-06-13)
5. Chang W.I., Lawler E.L. 1994. Sublinear approximate string matching and biological applications. Algorithmica 12:327-344 (matching statistics computed in one left-to-right scan with suffix links; confirmed from search snippets/lecture notes, paper not opened).
6. Kurtz S. et al. 2004. Versatile and open software for comparing large genomes. Genome Biol 5:R12 (maximal exact matches; as implemented by `SuffixTree.FindMaximalExactMatches`, MUMmer-identical per B05 F20).
7. Seqeron repository. SuffixTree.Algorithms.cs — `LongestCommonSubstringInfo` (documented tie-break: first found in 'other'). ../../../src/SuffixTree/Algorithms/SuffixTree/SuffixTree.Algorithms.cs
