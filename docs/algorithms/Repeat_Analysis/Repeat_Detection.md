# Repeat Detection (Longest Repeated Substring, Repeated Substrings, Maximal Repeats)

| Field | Value |
|-------|-------|
| Algorithm Group | Repeat Analysis (Genomic Analysis) |
| Test Unit ID | GENOMIC-REPEAT-001 |
| Related Projects | Seqeron.Genomics.Analysis, SuffixTree |
| Implementation Status | Production |
| Last Reviewed | 2026-11-19 (B09 review, audit round 1: line-by-line pass WP26) |

## 1. Overview

Repeat detection finds substrings that recur within a single DNA sequence — the building blocks of tandem repeats, transposable elements, and other repetitive genomic structures. This unit covers four exact operations: `FindLongestRepeat` / `FindAllLongestRepeats`, which return the longest substring(s) occurring at least twice (the Longest Repeated Substring, LRS), `FindRepeats`, which enumerates every distinct substring occurring at least twice with length ≥ a given minimum, and `FindMaximalRepeats`, the subset of those that are *maximal* (left- and right-diverse, Gusfield §7.12 [5]). All are exact (not heuristic): they realise the classical suffix-tree characterisation of repeats [1][2][5]. Occurrences may overlap. `GenomicAnalyzer.FindPalindromes` (reverse-complement palindromes) is a thin facade over `RepeatFinder.FindPalindromes` and is documented in [Palindrome_Detection.md](Palindrome_Detection.md); its contract is summarised in §5.5. Since audit A2-33 / A2-35 the same engine also searches both strands (`FindRepeatsBothStrands`, `FindRepeatPairs`; a string and its reverse complement form one class) and accepts nucleotide / protein strings with N / IUPAC symbols (`FindRepeats(string, …)`, …).

## 2. Scientific / Formal Basis

### 2.1 Domain Context

A substring that appears two or more times in a sequence is a *repeat*. In a suffix tree of the text, every leaf corresponds to a distinct suffix (a position), so a substring that occurs k times is the path label of an internal node having k leaves below it; a substring occurring once ends inside a leaf edge [1]. Repeat finding is therefore a traversal problem on the suffix tree, a family of applications attributed to Gusfield, *Algorithms on Strings, Trees, and Sequences* (ch. 5–7) [4].

### 2.2 Core Model

**Longest Repeated Substring.** "Find the longest string `r` such [that] `r` occurs at least twice in `T`: Find the deepest node that has ≥ 2 leaves under it" [1, §2.1]. "Deepest" is measured by string depth — the number of characters on the root-to-node path — which equals the length of the repeated substring [2][3]. Equivalently, build the suffix tree of `T` (with a terminal sentinel) and return the path label of the internal node of maximum string depth [2].

**All repeats ≥ minLength.** Every substring occurring ≥ 2 times is the longest common prefix (LCP) of two suffixes of `T`; restricting to adjacent suffixes in lexicographic order captures all maximal-length such prefixes [1]. A candidate qualifies iff its LCP length ≥ `minLength` and it occurs ≥ 2 times (i.e. corresponds to an internal node, never a leaf) [1][3].

**Maximal repeats (Gusfield §7.12 [5]).** A substring α is a *maximal repeat* if it occurs at least twice and no occurrence can be extended on both sides in the same way: α is *left-diverse* (the characters immediately before its occurrences are not all equal; the start of the text counts as a distinct symbol) and *right-diverse* (likewise after; the end of the text is a distinct symbol). In the suffix tree of `T$` every internal node is right-diverse (its children start with different characters, `$` being one of them), so the maximal repeats are exactly the internal nodes whose leaf set is left-diverse, and there are at most n of them. Left-diversity is propagated bottom-up as a bit set of preceding characters (leaf = character before its suffix; node = union over its subtree), which is linear. Every longest repeated substring is maximal; maximal ⊂ all repeated substrings; a *supermaximal* repeat (`RepeatFinder.FindSupermaximalRepeats`) is a maximal repeat that is not a substring of another maximal repeat.

**What is a "repeat" here (semantic decision record).** `FindRepeats` deliberately reports *every* distinct repeated substring (not only maximal ones): all lengths on the edge into an internal node, overlapping occurrences counted, no minimum count other than 2, no non-overlap requirement, case-insensitive through `DnaSequence` (upper-cased, ACGT only; the string overloads, §5.3, accept any alphabet). It therefore lists every prefix and suffix of a longer repeat (`ACGTACGT` → `A, AC, ACG, ACGT, C, CG, CGT, G, GT, T`); use `FindMaximalRepeats` for the non-redundant representation. It is *not* the same as `RepeatFinder.FindDirectRepeats` (maximal repeated **pairs**: one result per pair of occurrences, with spacing and length-cap options) nor `SuffixTree.FindMaximalRepeatedPairs` (the same pairs from the tree): a call-through was evaluated and rejected because the result shape differs (strings with occurrence sets vs. pairs) and the sets differ (pairs enumerate prefixes/suffixes only as part of maximal pairs). The relation is exact and tested: the strings of `FindMaximalRepeats` = the distinct substrings of the maximal pairs of `FindDirectRepeats` without spacing/length limits.

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | Every returned repeat occurs ≥ 2 times (`Count ≥ 2`). | Repeats map to internal nodes (≥ 2 leaves); leaves (single occurrence) are excluded [1, §2.1]. |
| INV-02 | `RepeatInfo.Length == Sequence.Length` = root-to-node string depth. | String depth equals substring length [2]. |
| INV-03 | Each `Positions[i]` is a true 0-based start of `Sequence` in the text. | Positions come from `SuffixTree.FindAllOccurrences`, the leaves below the node [1]. |
| INV-04 | `FindLongestRepeat` returns a substring no shorter than any repeated substring. | It is the deepest internal node by definition [1][2]. |
| INV-05 | `FindRepeats` returns only substrings with `Length ≥ minLength`. | Node/length filter in the enumeration loop. |
| INV-06 | `Positions` is sorted ascending. | Implementation sorts via `OrderBy` (output-shape convention; see 5.4). |
| INV-07 | `FindMaximalRepeats` ⊂ `FindRepeats`; every result is left- and right-diverse; at most n results; every longest repeat is maximal. | Gusfield §7.12 [5]; verified against a brute-force definition and a Python oracle (§5.6). |
| INV-08 | `FindLongestRepeat` = first element of `FindAllLongestRepeats` (the tied substring whose first occurrence is leftmost). | Deterministic tie-break (§5.4 #2). |

### 2.5 Comparison with Related Methods

| Aspect | LRS / repeat enumeration (this unit) | Tandem repeat finding (`FindTandemRepeats`) |
|--------|--------------------------------------|----------------------------------------------|
| Repeat type | Any recurring substring (dispersed or overlapping) | Consecutive repeating unit only |
| Structure used | Suffix tree (deepest internal node) | Direct windowed scan |
| Occurrence layout | Positions anywhere, may overlap | Adjacent copies from one start |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `sequence` | `DnaSequence` or `string` | required | Sequence to search. | `DnaSequence`: alphabet A/C/G/T (uppercased & validated); `string` overloads (A2-35): case-insensitive, only the alphabet's symbols repeat (default A/C/G/T; optional `alphabet` for proteins), any other symbol is a unique separator; empty allowed. |
| `maximalOnly` | `bool` | `false` | `FindRepeatsBothStrands`: keep only the classes left- and right-diverse over both strands. | |
| `strands` | `StrandSelection` | `Both` | `FindRepeatPairs`: `Forward` (`repeat-match -f`), `ReverseComplement`, `Both`. | |
| `minLength` | `int` | required | Minimum repeat length for `FindRepeats` / `FindMaximalRepeats`. | Repeats with `Length ≥ minLength` returned; values ≤ 0 still yield only substrings occurring ≥ 2× (no zero-length results). |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `RepeatInfo.Sequence` | `string` | The repeated substring (empty for `None`). |
| `RepeatInfo.Positions` | `IReadOnlyList<int>` | 0-based start positions, ascending. |
| `RepeatInfo.Length` | `int` | `Sequence.Length`. |
| `RepeatInfo.Count` | `int` | `Positions.Count` (≥ 2 for any non-`None` result). |
| `RepeatInfo.IsEmpty` | `bool` | True when no repeat (the `None` value). |

`FindLongestRepeat` returns one `RepeatInfo` (or `RepeatInfo.None`); `FindAllLongestRepeats` an `IReadOnlyList<RepeatInfo>` of all length ties; `FindRepeats` and `FindMaximalRepeats` return an `IEnumerable<RepeatInfo>` in suffix order. The both-strand methods return `StrandedRepeatInfo` (`Sequence` = the ordinal-smaller of the string and its reverse complement, `Positions` = forward starts of the string, `ReversePositions` = forward starts of the windows that spell its reverse complement) and `RepeatPair(FirstPosition, SecondPosition, Length, IsReverseComplement)`.

### 3.3 Preconditions and Validation

Indexing is 0-based. `DnaSequence` normalises to uppercase and rejects non-ACGT characters at construction (`ArgumentException`); these methods themselves throw only `ArgumentNullException` for a `null` sequence (eagerly, at the call; before the 2026-11 review a `null` surfaced as a `NullReferenceException` on enumeration). An empty sequence and a sequence with no repeat both yield `RepeatInfo.None` from `FindLongestRepeat` and an empty enumeration from `FindRepeats`. Occurrences may overlap.

## 4. Algorithm

### 4.1 High-Level Steps

**FindLongestRepeat**
1. Build/obtain the suffix tree of the sequence (`DnaSequence.SuffixTree`, Ukkonen).
2. Collect all deepest internal nodes via the canonical `SuffixTree.FindAllLongestRepeatedSubstrings()` (B05; ordered by first occurrence, positions ascending). `FindAllLongestRepeats` returns them all; `FindLongestRepeat` returns the first (`RepeatInfo.None` when empty).

**FindRepeats**
1. Traverse the suffix tree (`ISuffixTreeDiagnostics.Traverse`, pre-order, children in ascending key order) and collect every internal non-root node with its parent's path-label length `p` and its own path-label length `d`.
2. For each such node with `d ≥ max(1, minLength)`, resolve its occurrence set once (`FindAllOccurrences` of its path label), sort it ascending, and emit the repeats of lengths `max(p + 1, max(1, minLength))` … `d` (all share that occurrence set). (A substring occurs ≥ 2 times iff it ends on an edge into an internal node, so every length on the edge — not just the full node label — is a repeat.)
**FindMaximalRepeats** — the same traversal, but the visitor additionally keeps a stack of bit sets: a leaf pushes the bit of the character before its suffix (a dedicated bit for suffix 0), an exiting internal node ORs its set into its parent and is kept when it has ≥ 2 bits (left-diverse). Only the full node label is emitted (shorter lengths on the edge are right-extendable). O(n) plus output.

3. The emission order equals the order of the former sorted-suffix implementation (by the smallest suffix starting with the repeat, then by length); the former implementation is kept as a reference oracle in `tests/Seqeron/Seqeron.Genomics.Tests/Differential/GenomicAnalyzerFindRepeatsOracleTests.cs`.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `FindLongestRepeat` / `FindAllLongestRepeats` | O(n) suffix-tree build (cached) + one O(n) traversal for the deepest nodes + O(occ) | O(n) | Linear via the cached suffix tree [1][2][3]. |
| `FindMaximalRepeats` | O(n + output) | O(n) | Left-diversity bit sets bottom-up; output ≤ n strings. 1 Mbp random DNA, `minLength` 20: `FindRepeats` + `FindMaximalRepeats` ≈ 4.7 s including the tree build (test `FindRepeats_And_FindMaximalRepeats_OneMegabase_…`). |
| `FindRepeats` | O(n + output) | O(n + node-local occurrence list) | One tree traversal plus, per internal node, one occurrence query O(d + occ); the total is bounded by the size of the result (which can itself be O(n²) for a homopolymer). Replaced the O(n²) `GetAllSuffixes` method (ST-S0 WP 0.4; 20 kbp: 400 MB allocated before, a few MB now). |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [GenomicAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/GenomicAnalyzer.cs)

- `GenomicAnalyzer.FindLongestRepeat(DnaSequence)`: returns the LRS as a `RepeatInfo`.
- `GenomicAnalyzer.FindAllLongestRepeats(DnaSequence)`: every longest repeat (all ties).
- `GenomicAnalyzer.FindRepeats(DnaSequence, int minLength)`: enumerates all repeats ≥ minLength.
- `GenomicAnalyzer.FindMaximalRepeats(DnaSequence, int minLength)`: the maximal (left- and right-diverse) repeats.
- `GenomicAnalyzer.FindPalindromes(DnaSequence, int, int)`: facade over `RepeatFinder.FindPalindromes` (§5.5).
- Both strands (A2-33; [GenomicAnalyzer.Repeats.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/GenomicAnalyzer.Repeats.cs)): `FindRepeatsBothStrands(DnaSequence | string, minLength, maximalOnly)`, `FindAllLongestRepeatsBothStrands(DnaSequence | string)`, `FindLongestRepeatBothStrands(DnaSequence)`, `FindRepeatPairs(DnaSequence | string, minLength = 5, StrandSelection = Both)`.
- String / IUPAC overloads (A2-35): `FindLongestRepeat(string, alphabet?)`, `FindAllLongestRepeats(string, …)`, `FindRepeats(string, minLength, …)`, `FindMaximalRepeats(string, minLength, …)`.
- MCP: `find_repeats` (strand / maximalOnly / IUPAC input), `find_all_longest_repeats` (A2-33 / A2-36).
- Backing structure: [SuffixTree.LongestRepeatedSubstring / FindAllOccurrences / GetAllSuffixes](../../../src/SuffixTree/Algorithms/SuffixTree/SuffixTree.Algorithms.cs).

### 5.2 Current Behavior

Both methods delegate the core string work to the repository suffix tree. `FindLongestRepeat` uses the suffix tree's cached deepest-internal-node bookkeeping (`LongestRepeatedSubstring`) directly. `FindRepeats` traverses the suffix tree and, for every internal node, emits every length on its incoming edge (length ≥ `max(1, minLength)`), resolving positions with `FindAllOccurrences`. Positions are returned ascending. Overlapping occurrences are counted (consistent with the definition [3]).

**Search-reuse decision (mandatory record):** the repository `SuffixTree` (namespace `SuffixTree`) is the algorithmically appropriate structure here — repeat finding is *the* canonical suffix-tree application (deepest internal node) [1][2][3], and many occurrence queries run against one fixed text (O(m) lookups after O(n) build). The existing implementation already builds on it; **no naive scan was introduced**. A naive O(n²)/O(n³) all-pairs substring scan would be strictly worse and was rejected.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- LRS = longest substring occurring ≥ 2 times via the deepest internal node [1, §2.1][2].
- Repeats correspond to internal nodes (≥ 2 occurrences), never leaves [1][3].
- Overlapping occurrences counted (`AAAAAAAAAA`→`AAAAAAAAA`; `ABABABA`→`ABABA`) [3].

**Intentionally simplified:** none. (Until 2026-11 the tie-break of `FindLongestRepeat` was "whichever the tree's child order records"; it is now deterministic — the leftmost first occurrence — see 5.4 #2.)

**Both strands and string input (audit A2-33 / A2-35; formerly "not implemented"):**

- `FindRepeatsBothStrands(sequence, minLength, maximalOnly)` / `FindAllLongestRepeatsBothStrands` / `FindLongestRepeatBothStrands`: the repeats of `S # revcomp(S)` — a string and its reverse complement form one class (reported as the ordinal-smaller of the two, `StrandedRepeatInfo.Sequence`), `Positions` are the forward starts of the string, `ReversePositions` the forward-strand starts of the windows that spell its reverse complement; a palindromic site (GAATTC) is a repeat on its own (`repeat-match` prints the pair `i = k`). `maximalOnly` keeps the classes left- and right-diverse over all their occurrences on both strands, equivalently the strings of at least one maximal repeated pair (Gusfield 1997 §7.12 on the doubled text). One suffix tree of the doubled text is built and traversed with the same engine as the forward search (non-matching symbols and the strand boundary are unique separators).
- `FindRepeatPairs(sequence, minLength, StrandSelection)`: the unfiltered MUMmer `repeat-match` list (`maxLength = int.MaxValue`, `minSpacing = int.MinValue` of the canonical `RepeatFinder.FindDirectRepeats` / `FindReverseComplementRepeats`) as `RepeatPair(FirstPosition, SecondPosition, Length, IsReverseComplement)`; `StrandSelection.Forward` is `repeat-match -f`.
- String / IUPAC overloads of `FindLongestRepeat`, `FindAllLongestRepeats`, `FindRepeats`, `FindMaximalRepeats` (optional `alphabet` for a protein: `"ACDEFGHIKLMNPQRSTVWY"`): only the alphabet's symbols repeat, every other symbol (N, IUPAC, U, gap) is a unique separator that never repeats and makes an adjacent repeat non-extendable; positions are input positions; an A/C/G/T string equals the `DnaSequence` result. More than 55 000 separate non-matching runs are refused (split the sequence).
- Verification: `repeat-match` 3.23 (both strands and `-f`, run on `S + 'N'` because without the sentinel S and revcomp(S) share the final terminator and the tool loses palindromic pairs, B04 F25) on 1900 (sequence, L) rows — the pair lists equal and the class set of its output equals the maximal classes on every ACGT row — plus an independent Python class brute force over `S # revcomp(S)` and an in-test brute force of the string overloads with N runs (`scripts/analysis/generate_strand_repeat_oracle.py`, `GenomicAnalyzer_StrandRepeats_Oracle_Tests`).

**Not implemented:** gapped / approximate repeats (`RepeatFinder.FindApproximateDirectRepeats`, `FindInvertedRepeats`, B04); supermaximal repeats (`RepeatFinder.FindSupermaximalRepeats`).

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | `Positions` sorted ascending | Assumption | Output ordering only; the position *set* is fixed | accepted | INV-06; output-shape convention |
| 2 | Tie-break among equal-length LRS | Decision | Which equal-length winner is returned | fixed (B09 F29) | Spec requires only "a" longest [2]; now the tied substring with the leftmost first occurrence (= first of `FindAllLongestRepeats`, canonical `SuffixTree.FindAllLongestRepeatedSubstrings`), identical for every tree implementation; before: tree child-storage order (unspecified, differs between tree implementations) |
| 4 | Null sequence | Deviation (fix) | `NullReferenceException` raised lazily on enumeration | fixed (B09 F28) | `ArgumentNullException` at the call for all four methods |
| 3 | `minLength ≤ 0` clamped to 1 in `FindRepeats`/`FindMaximalRepeats` | Deviation (fix) | Without the clamp the empty string was emitted as a zero-length "repeat" | fixed | A repeat is a non-empty substring (internal node, not the root) [1]; effective minimum is `max(1, minLength)` |

### 5.5 `GenomicAnalyzer.FindPalindromes` (B09 F27)

Formerly a verbatim copy of `RepeatFinder.FindPalindromesCore` without its validation (DUP_MAP §17; B04 cross-batch request): `minLength = 2` reported every AT/CG/GC/TA dinucleotide, an odd `minLength` scanned only odd lengths and silently returned nothing, `null` threw a `NullReferenceException`, the scan was O(n·L²) (a substring + reverse complement per window and length) and the output was length-major. It now delegates to the canonical `RepeatFinder.FindPalindromes(DnaSequence, int, int)` and maps `PalindromeResult` → `PalindromeInfo(Sequence, Position)`. **Contract change (documented, intended):** `minLength` must be even and ≥ 4, `maxLength` ≥ `minLength` (`ArgumentOutOfRangeException` at the call); `null` → `ArgumentNullException`; results are ordered by position, then length (REVP sample order) instead of length, then position; cost O(n·maxLength) (1 Mbp, 4–12: ≈ 0.4 s). Everything else (default 4–12, 0-based positions, all nested/overlapping windows, odd `maxLength` rounded down) is unchanged. `DnaSequence` is ACGT-only, so IUPAC symbols cannot occur here (the string overload `RepeatFinder.FindPalindromes(string, …)` and MCP `find_palindromes` treat non-ACGT symbols as never pairing). Oracle: Biopython `Seq.reverse_complement` over every window, 910 (sequence, minLength, maxLength) rows, plus the Rosalind REVP sample.

### 5.6 Independent oracles (B09 F28–F30)

`scripts/analysis/generate_repeat_oracle.py` computes, for 628 (sequence, minLength) rows (random ACGT/AC/ACG, homopolymers, Fibonacci, Thue–Morse, planted repeats, the cited worked examples), the repeated substrings, maximal repeats and all longest repeats **twice** — brute force over every substring, and suffix array + Kasai LCP + lcp-intervals (maximal repeats via left characters of the interval) — asserting equality, and writes `tests/.../TestData/Analysis/repeat_oracle.tsv` (count + SHA-256 prefix of the ordered repeat list, explicit maximal-repeat and longest-repeat lists). The C# suffix-tree methods match it on all rows (set, positions, order). An in-test brute force on 480 random strings (4 alphabets × 120 trials × 4 minLengths) additionally pins the definition including emission order, and a 1 Mbp run verifies reported repeats directly against the text.

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Empty sequence | `FindLongestRepeat` → `None`; `FindRepeats` → empty | No substring occurs twice in ε |
| No repeat (e.g. `ACGT`) | `FindLongestRepeat` → `None` | No internal node with ≥ 2 leaves [1][3] |
| Overlapping run (`AAAAAAAAAA`) | LRS `AAAAAAAAA`, positions {0,1} | Overlap allowed [3] |
| `minLength` > all repeats | `FindRepeats` → empty | minLength filter |
| `minLength` ≤ 0 | `FindRepeats` → only substrings occurring ≥ 2×; no zero-length results | INV-01 still gates on ≥ 2 occurrences |

### 6.2 Limitations

The forward-strand methods do not search the reverse complement; the both-strand methods (§5.3) do (the doubled text `S # revcomp(S)` costs about twice the tree memory). `DnaSequence` input is restricted to ACGT; the string overloads treat every other symbol as a unique separator (at most 55 000 separate non-matching runs per call). `FindRepeats` is linear in the input plus the size of its output; the output itself can be quadratic (a homopolymer of length n has n − 1 repeats with up to n positions each), and the in-memory suffix tree needs about 170 bytes per input character, so whole-genome-scale inputs still need the persistent tree or a size limit.

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
var dna = new DnaSequence("ATCGATCGA");
RepeatInfo lrs = GenomicAnalyzer.FindLongestRepeat(dna);
// lrs.Sequence == "ATCGA", lrs.Length == 5, lrs.Count == 2, lrs.Positions == [0, 4]

foreach (var r in GenomicAnalyzer.FindRepeats(new DnaSequence("ACGTACGTTTTTACGT"), 3))
{
    // Full set (8 distinct substrings, each occurring >=2 times, length >=3):
    //   ACG@{0,4,12}, ACGT@{0,4,12}, CGT@{1,5,13},
    //   TAC@{3,11}, TACG@{3,11}, TACGT@{3,11}, TTT@{7,8,9}, TTTT@{7,8}
}
```

**Numerical walk-through:** for `ATCGATCGA`, suffixes starting at 0 and 4 share the prefix `ATCGA` (length 5); no longer prefix is shared by two suffixes, so the deepest internal node spells `ATCGA` — matching the Wikipedia example `ATCGATCGA$` → `ATCGA` [2].

### 7.3 Related Tests, Evidence, or Documents

- Both strands / string input: [GenomicAnalyzer_StrandRepeats_Oracle_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Analysis/GenomicAnalyzer_StrandRepeats_Oracle_Tests.cs) (`repeat-match` 3.23, 1900 rows)
- Tests: [GenomicAnalyzer_FindRepeats_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Analysis/GenomicAnalyzer_FindRepeats_Tests.cs) — `INV-01`..`INV-06`; [GenomicAnalyzer_Repeat_Oracle_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Analysis/GenomicAnalyzer_Repeat_Oracle_Tests.cs) — Python oracle, brute force, `INV-07`/`INV-08`; [GenomicAnalyzer_FindPalindromes_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Analysis/GenomicAnalyzer_FindPalindromes_Tests.cs); [GenomicAnalyzerFindRepeatsOracleTests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Differential/GenomicAnalyzerFindRepeatsOracleTests.cs) — the pre-WP-0.4 implementation as differential oracle
- Evidence: [GENOMIC-REPEAT-001-Evidence.md](../../../docs/Evidence/GENOMIC-REPEAT-001-Evidence.md)

## 8. References

1. Sleator, D. (course staff). 2017. *15-451/651 Algorithm Design and Analysis, Lecture #10: Suffix Trees and Arrays* (§2.1 Longest repeat). Carnegie Mellon University. https://www.cs.cmu.edu/~15451-f17/lectures/lec10-sufftree.pdf
2. Wikipedia contributors. 2026. *Longest repeated substring problem.* https://en.wikipedia.org/wiki/Longest_repeated_substring_problem
3. GeeksforGeeks. *Suffix Tree Application 3 – Longest Repeated Substring.* https://www.geeksforgeeks.org/dsa/suffix-tree-application-3-longest-repeated-substring/
4. Langmead, B. *Suffix Trees (lecture notes)* (cites Gusfield 5.4). Johns Hopkins University. https://www.cs.jhu.edu/~langmea/resources/lecture_notes/08_suffix_trees_v2.pdf
5. Gusfield, D. 1997. *Algorithms on Strings, Trees, and Sequences.* Cambridge University Press. https://doi.org/10.1017/CBO9780511574931
