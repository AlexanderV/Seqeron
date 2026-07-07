# Repeat Detection (Longest Repeated Substring)

| Field | Value |
|-------|-------|
| Algorithm Group | Repeat Analysis (Genomic Analysis) |
| Test Unit ID | GENOMIC-REPEAT-001 |
| Related Projects | Seqeron.Genomics.Analysis, SuffixTree |
| Implementation Status | Production |
| Last Reviewed | 2026-06-15 |

## 1. Overview

Repeat detection finds substrings that recur within a single DNA sequence — the building blocks of tandem repeats, transposable elements, and other repetitive genomic structures. This unit covers two exact operations: `FindLongestRepeat`, which returns the longest substring occurring at least twice (the Longest Repeated Substring, LRS), and `FindRepeats`, which enumerates every distinct substring occurring at least twice with length ≥ a given minimum. Both are exact (not heuristic): they realise the classical suffix-tree characterisation of repeats [1][2]. Occurrences may overlap.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

A substring that appears two or more times in a sequence is a *repeat*. In a suffix tree of the text, every leaf corresponds to a distinct suffix (a position), so a substring that occurs k times is the path label of an internal node having k leaves below it; a substring occurring once ends inside a leaf edge [1]. Repeat finding is therefore a traversal problem on the suffix tree, a family of applications attributed to Gusfield, *Algorithms on Strings, Trees, and Sequences* (ch. 5–7) [4].

### 2.2 Core Model

**Longest Repeated Substring.** "Find the longest string `r` such [that] `r` occurs at least twice in `T`: Find the deepest node that has ≥ 2 leaves under it" [1, §2.1]. "Deepest" is measured by string depth — the number of characters on the root-to-node path — which equals the length of the repeated substring [2][3]. Equivalently, build the suffix tree of `T` (with a terminal sentinel) and return the path label of the internal node of maximum string depth [2].

**All repeats ≥ minLength.** Every substring occurring ≥ 2 times is the longest common prefix (LCP) of two suffixes of `T`; restricting to adjacent suffixes in lexicographic order captures all maximal-length such prefixes [1]. A candidate qualifies iff its LCP length ≥ `minLength` and it occurs ≥ 2 times (i.e. corresponds to an internal node, never a leaf) [1][3].

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | Every returned repeat occurs ≥ 2 times (`Count ≥ 2`). | Repeats map to internal nodes (≥ 2 leaves); leaves (single occurrence) are excluded [1, §2.1]. |
| INV-02 | `RepeatInfo.Length == Sequence.Length` = root-to-node string depth. | String depth equals substring length [2]. |
| INV-03 | Each `Positions[i]` is a true 0-based start of `Sequence` in the text. | Positions come from `SuffixTree.FindAllOccurrences`, the leaves below the node [1]. |
| INV-04 | `FindLongestRepeat` returns a substring no shorter than any repeated substring. | It is the deepest internal node by definition [1][2]. |
| INV-05 | `FindRepeats` returns only substrings with `Length ≥ minLength`. | LCP-length filter in the enumeration loop. |
| INV-06 | `Positions` is sorted ascending. | Implementation sorts via `OrderBy` (output-shape convention; see 5.4). |

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
| `sequence` | `DnaSequence` | required | Sequence to search. | Alphabet A/C/G/T (uppercased & validated by `DnaSequence`); empty allowed. |
| `minLength` | `int` | required | Minimum repeat length for `FindRepeats`. | Repeats with `Length ≥ minLength` returned; values ≤ 0 still yield only substrings occurring ≥ 2× (no zero-length results). |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `RepeatInfo.Sequence` | `string` | The repeated substring (empty for `None`). |
| `RepeatInfo.Positions` | `IReadOnlyList<int>` | 0-based start positions, ascending. |
| `RepeatInfo.Length` | `int` | `Sequence.Length`. |
| `RepeatInfo.Count` | `int` | `Positions.Count` (≥ 2 for any non-`None` result). |
| `RepeatInfo.IsEmpty` | `bool` | True when no repeat (the `None` value). |

`FindLongestRepeat` returns one `RepeatInfo` (or `RepeatInfo.None`); `FindRepeats` returns an `IEnumerable<RepeatInfo>`.

### 3.3 Preconditions and Validation

Indexing is 0-based. `DnaSequence` normalises to uppercase and rejects non-ACGT characters at construction (`ArgumentException`); these methods themselves do not throw for valid `DnaSequence` inputs. An empty sequence and a sequence with no repeat both yield `RepeatInfo.None` from `FindLongestRepeat` and an empty enumeration from `FindRepeats`. Occurrences may overlap.

## 4. Algorithm

### 4.1 High-Level Steps

**FindLongestRepeat**
1. Build/obtain the suffix tree of the sequence (`DnaSequence.SuffixTree`, Ukkonen).
2. Query the deepest internal node via `SuffixTree.LongestRepeatedSubstring()`.
3. If empty, return `RepeatInfo.None`; otherwise enumerate all occurrences (`FindAllOccurrences`), sort ascending, return.

**FindRepeats**
1. Obtain all suffixes (`GetAllSuffixes`) and sort them.
2. For each adjacent pair, compute the LCP length `L`.
3. For every prefix length `ℓ` from `max(1, minLength)` to `L`, take the length-`ℓ` prefix; if it is new, enumerate its occurrences and emit when it occurs ≥ 2 times. (A substring occurs ≥ 2 times iff it is a *prefix* of some adjacent-pair LCP, so all prefixes — not just the full LCP — must be considered for completeness.)

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `FindLongestRepeat` | O(n) deepest-node query after O(n) suffix-tree build; +O(occ) for occurrences | O(n) | Linear via the cached suffix tree [1][2][3]. |
| `FindRepeats` | O(n²) worst case | O(n²) worst case | O(n) suffixes of up to O(n) length sorted + per-pair LCP + occurrence enumeration. |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [GenomicAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/GenomicAnalyzer.cs)

- `GenomicAnalyzer.FindLongestRepeat(DnaSequence)`: returns the LRS as a `RepeatInfo`.
- `GenomicAnalyzer.FindRepeats(DnaSequence, int minLength)`: enumerates all repeats ≥ minLength.
- Backing structure: [SuffixTree.LongestRepeatedSubstring / FindAllOccurrences / GetAllSuffixes](../../../src/SuffixTree/Algorithms/SuffixTree/SuffixTree.Algorithms.cs).

### 5.2 Current Behavior

Both methods delegate the core string work to the repository suffix tree. `FindLongestRepeat` uses the suffix tree's cached deepest-internal-node bookkeeping (`LongestRepeatedSubstring`) directly. `FindRepeats` uses `GetAllSuffixes` + adjacent-pair LCP, enumerating every prefix (length ≥ `max(1, minLength)`) of each LCP for completeness, and resolves positions with `FindAllOccurrences`. Positions are returned ascending. Overlapping occurrences are counted (consistent with the definition [3]).

**Search-reuse decision (mandatory record):** the repository `SuffixTree` (namespace `SuffixTree`) is the algorithmically appropriate structure here — repeat finding is *the* canonical suffix-tree application (deepest internal node) [1][2][3], and many occurrence queries run against one fixed text (O(m) lookups after O(n) build). The existing implementation already builds on it; **no naive scan was introduced**. A naive O(n²)/O(n³) all-pairs substring scan would be strictly worse and was rejected.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- LRS = longest substring occurring ≥ 2 times via the deepest internal node [1, §2.1][2].
- Repeats correspond to internal nodes (≥ 2 occurrences), never leaves [1][3].
- Overlapping occurrences counted (`AAAAAAAAAA`→`AAAAAAAAA`; `ABABABA`→`ABABA`) [3].

**Intentionally simplified:**

- Tie-breaking among equal-length longest repeats: the spec requires only *a* longest substring [2]; the implementation returns whichever the suffix tree's deepest-node bookkeeping records. **Consequence:** for inputs with multiple equal-length maximal repeats, the specific winner is unspecified but always a valid longest repeat.

**Not implemented:**

- Maximal-repeat / supermaximal-repeat classification (left-diverse nodes [4]): not distinguished here; **users should rely on** dedicated repeat-classification tooling — no current in-repo alternative.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | `Positions` sorted ascending | Assumption | Output ordering only; the position *set* is fixed | accepted | INV-06; output-shape convention |
| 2 | Tie-break among equal-length LRS | Assumption | Which equal-length winner is returned | accepted | Spec requires only "a" longest [2]; see 5.3 |
| 3 | `minLength ≤ 0` clamped to 1 in `FindRepeats` | Deviation (fix) | Without the clamp the empty string was emitted as a zero-length "repeat" | fixed | A repeat is a non-empty substring (internal node, not the root) [1]; effective minimum is `max(1, minLength)` |

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

Does not classify maximal vs supermaximal repeats and does not search the reverse complement (dispersed inverted repeats); restricted to the ACGT alphabet enforced by `DnaSequence`. `FindRepeats` is worst-case O(n²) in time and space, so it is unsuitable for whole-genome-scale inputs.

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

- Tests: [GenomicAnalyzer_FindRepeats_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Analysis/GenomicAnalyzer_FindRepeats_Tests.cs) — covers `INV-01`..`INV-06`
- Evidence: [GENOMIC-REPEAT-001-Evidence.md](../../../docs/Evidence/GENOMIC-REPEAT-001-Evidence.md)

## 8. References

1. Sleator, D. (course staff). 2017. *15-451/651 Algorithm Design and Analysis, Lecture #10: Suffix Trees and Arrays* (§2.1 Longest repeat). Carnegie Mellon University. https://www.cs.cmu.edu/~15451-f17/lectures/lec10-sufftree.pdf
2. Wikipedia contributors. 2026. *Longest repeated substring problem.* https://en.wikipedia.org/wiki/Longest_repeated_substring_problem
3. GeeksforGeeks. *Suffix Tree Application 3 – Longest Repeated Substring.* https://www.geeksforgeeks.org/dsa/suffix-tree-application-3-longest-repeated-substring/
4. Langmead, B. *Suffix Trees (lecture notes)* (cites Gusfield 5.4). Johns Hopkins University. https://www.cs.jhu.edu/~langmea/resources/lecture_notes/08_suffix_trees_v2.pdf
5. Gusfield, D. 1997. *Algorithms on Strings, Trees, and Sequences.* Cambridge University Press. https://doi.org/10.1017/CBO9780511574931
