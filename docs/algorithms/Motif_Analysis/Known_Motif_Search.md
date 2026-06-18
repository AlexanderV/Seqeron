# Known Motif Search

| Field | Value |
|-------|-------|
| Algorithm Group | Motif Analysis (Analysis) |
| Test Unit ID | GENOMIC-MOTIFS-001 |
| Related Projects | Seqeron.Genomics.Analysis, SuffixTree |
| Implementation Status | Production |
| Last Reviewed | 2026-06-13 |

## 1. Overview

Known Motif Search locates a *set* of pre-specified motifs (short DNA patterns of biological
interest, e.g. restriction sites or transcription-factor binding consensus) within a subject
sequence, returning for each motif the start positions of **all** of its occurrences. It is the
classical exact set-matching problem — "find all occurrences of each pattern P in text T" [1] —
and is exact (not heuristic): a position is reported iff the motif matches the subject there.
Overlapping occurrences are all reported [1][2]. It is the right tool when many fixed query motifs
are checked against one (re-usable) sequence index.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

A *motif* here is a fixed nucleotide string; "known" means the queries are supplied by the caller
(a curated set), as opposed to *de novo* motif discovery. A canonical biological example is the
EcoRI restriction recognition site `GAATTC` [3]. Exact matching reports every alignment of the
motif against the subject.

### 2.2 Core Model

For text `T` (length `m`) and pattern `P` (length `n`), the exact-matching problem is to "find all
occurrences of a pattern string P ... in a text string T" [1]; the answer is the set
`{ i : T[i..i+n-1] = P }` of 0-based start positions. Occurrences may overlap: for `P = aaa`,
`T = aaaaa` there are "three (overlapping) occurrences" [1] (positions 0, 1, 2). Biopython's
`Seq.count_overlap` confirms the overlap-aware count: `Seq("AAAA").count_overlap("AA") == 3`,
versus the non-overlapping `2` [2]. For a *set* of motifs, each motif is matched independently and
maps to its own position set, mirroring Biopython `Seq.search`, which "Search[es] the substrings
subs in self and yield[s] the index and substring found" for multiple queries [2].

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | Every reported position `p` for motif `m` satisfies `T[p..p+|m|-1] == m` (0-based) | Exact-matching definition [1] |
| INV-02 | All occurrences are reported, including overlapping ones | Exact-matching counts every alignment [1]; Biopython `count_overlap` [2] |
| INV-03 | Each motif's positions are sorted strictly ascending and distinct | Positions form a set [1]; implementation sorts the suffix-tree (DFS-order) output |
| INV-04 | A motif with zero occurrences is omitted from the result | The occurrence set is empty [1]; per-motif hits only when found [2] |
| INV-05 | Result keys are upper-cased motif strings | DNA processed upper-cased [2]; `DnaSequence` upper-cases on construction |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `sequence` | `DnaSequence` | required | Subject sequence (text T) | Upper-cased ACGT (validated by `DnaSequence`); empty allowed |
| `motifs` | `IEnumerable<string>` | required | Set of query motifs | Non-null; individual motifs upper-cased; empty/whitespace motifs skipped |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| (return) | `Dictionary<string, IReadOnlyList<int>>` | Map from upper-cased motif → ascending list of 0-based start positions of all its occurrences. Motifs with no occurrence are absent. |

### 3.3 Preconditions and Validation

`motifs` null → `ArgumentNullException`. Empty/whitespace motifs are skipped (the empty string is
not a motif — a suffix-tree search for `""` would return every position). Motifs are upper-cased
(case-insensitive matching; result keyed by the upper-cased form). Duplicate motifs collapsing to
the same upper-cased key produce a single entry. Indexing is 0-based; positions are starts of the
matched substring. Empty subject sequence → empty result (no motif occurs).

## 4. Algorithm

### 4.1 High-Level Steps

1. Validate `motifs` is non-null; obtain the subject's suffix tree (built once, O(n)).
2. For each motif: skip if empty/whitespace; upper-case it; skip if its key was already processed.
3. Query `SuffixTree.FindAllOccurrences(motif)` to enumerate all start positions (overlaps included).
4. If non-empty, sort positions ascending and store under the upper-cased motif key.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

Core data structure: the repository **suffix tree** (`SuffixTree`), built once for the subject and
reused across all motif queries. `FindAllOccurrences` enumerates leaf start positions in DFS order
(not sorted), so the implementation sorts each motif's positions for a deterministic contract.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| Suffix-tree construction | O(n) | O(n) | Built once for the subject (Ukkonen); cached on `DnaSequence` |
| Per-motif search + sort | O(\|m\| + occ·log occ) | O(occ) | `FindAllOccurrences` is O(\|m\|+occ); sort dominates when occ large |
| Total (k motifs) | O(n + Σ(\|mᵢ\| + occᵢ·log occᵢ)) | O(n + occ) | Index reuse across motifs is the advantage over per-motif naive scan |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [GenomicAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/GenomicAnalyzer.cs)

- `GenomicAnalyzer.FindKnownMotifs(DnaSequence, IEnumerable<string>)`: multi-motif exact search returning per-motif sorted 0-based position lists.

### 5.2 Current Behavior

**Suffix-tree reuse decision (search/matching unit):** this is exact-match occurrence enumeration
of *multiple* fixed queries against *one* reusable text — the textbook ideal use of a suffix tree
(O(n) build amortized across all queries, O(\|m\|+occ) per query). The repository `SuffixTree`
(`DnaSequence.SuffixTree.FindAllOccurrences`) is therefore used directly rather than a naive
per-motif scan; no approximate/scored matching is involved that would disqualify it.

`FindAllOccurrences` returns occurrences in DFS order, which is not sorted; `FindKnownMotifs` sorts
each motif's positions ascending (INV-03). Empty/whitespace motifs are skipped and duplicate
upper-cased keys collapse to one entry. These two behaviors corrected a prior version that returned
DFS-order positions and did not guard empty motifs.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- All start positions where `T[i..i+n-1] == P` reported, for each motif (exact-matching set) [1].
- Overlapping occurrences all reported (e.g. `AAA` in `AAAAA` → {0,1,2}) [1][2].
- Multi-motif set: each motif mapped to its own occurrence set [2].
- Upper-cased processing of DNA / motifs [2].

**Intentionally simplified:**

- (none)

**Not implemented:**

- Degenerate / IUPAC-coded motifs and reverse-complement matching: feature: ambiguous-base or
  both-strand motif search; **users should rely on:** PROSITE/regex motif methods
  (`ProteinMotifFinder`) for patterns, or reverse-complement the motif explicitly before searching.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Empty/whitespace motif skipped | Assumption | A `""` query returns no entry rather than "all positions" | accepted | No authoritative source defines an empty motif's occurrence set as meaningful; mirrors `FindMotif` |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Empty motif set | Empty dictionary | No motif to search [1] |
| Motif absent in subject | Omitted from result | Empty occurrence set [1] (INV-04) |
| Overlapping motif (`AAA` in `AAAAA`) | {0,1,2} | All overlapping occurrences reported [1][2] (INV-02) |
| Lower/mixed-case motif | Matched; key upper-cased | Upper-cased processing [2] (INV-05) |
| Empty/whitespace motif | Skipped | Empty string is not a motif (deviation #1) |
| `motifs` null | `ArgumentNullException` | Input validation |
| Empty subject sequence | Empty dictionary | No motif occurs |

### 6.2 Limitations

Exact matching only — no mismatches, gaps, IUPAC degeneracy, or reverse-complement strand search.
Motif strings are validated only via upper-casing; non-ACGT characters in a motif simply fail to
match (return no positions) rather than raising.

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
var seq = new DnaSequence("GAATTCAAAGAATTC");
var hits = GenomicAnalyzer.FindKnownMotifs(seq, new[] { "GAATTC" });
// hits["GAATTC"] == [0, 9]  (EcoRI sites)
```

**Numerical walk-through:** For `T = "AAAAA"`, motif `"AAA"` (length 3): alignments start at
i = 0 (`AAA`), 1 (`AAA`), 2 (`AAA`); i = 3 would need indices 3..5 but T has length 5 (max start 2).
Result: positions {0, 1, 2}, three overlapping occurrences [1].

### 7.3 Related Tests, Evidence, or Documents

- Tests: [GenomicAnalyzer_FindKnownMotifs_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/GenomicAnalyzer_FindKnownMotifs_Tests.cs) — covers INV-01..INV-05
- Evidence: [GENOMIC-MOTIFS-001-Evidence.md](../../../docs/Evidence/GENOMIC-MOTIFS-001-Evidence.md)
- Related algorithms: [Repeat_Detection](../Repeat_Analysis/Repeat_Detection.md)

## 8. References

1. Gusfield, D. 1997. *Algorithms on Strings, Trees and Sequences: Computer Science and Computational Biology*. Cambridge University Press. ISBN 0-521-58519-8. Exact-matching definition and overlapping-occurrence example (Tufts COMP 150GEN): https://www.cs.tufts.edu/comp/150GEN/classpages/exact.html
2. Cock, P.J.A. et al. Biopython `Bio.Seq` module — `search` and `count_overlap` (master, accessed 2026-06-13): https://raw.githubusercontent.com/biopython/biopython/master/Bio/Seq.py
3. Wikipedia contributors. "Restriction site" (accessed 2026-06-13). EcoRI recognizes GAATTC: https://en.wikipedia.org/wiki/Restriction_site
