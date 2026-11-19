# Known Motif Search

| Field | Value |
|-------|-------|
| Algorithm Group | Motif Analysis (Analysis) |
| Test Unit ID | GENOMIC-MOTIFS-001 |
| Related Projects | Seqeron.Genomics.Analysis, SuffixTree |
| Implementation Status | Production |
| Last Reviewed | 2026-11-18 (B09 review) |

## 1. Overview

Known Motif Search locates a *set* of pre-specified motifs (short DNA patterns of biological
interest, e.g. restriction sites or transcription-factor binding consensus) within a subject
sequence, returning for each motif the start positions of **all** of its occurrences. It is the
classical exact set-matching problem — "find all occurrences of each pattern P in text T" [1] —
and is exact (not heuristic): a position is reported iff the motif matches the subject there.
Overlapping occurrences are all reported [1][2]. Optional extensions, both verified against EMBOSS
`fuzznuc` [4]: IUPAC-degenerate motifs (`TATAWAW`) and the reverse strand. The motif set is supplied
by the caller; there is no built-in catalogue here (the fixed regulatory-element library — TATA box,
Kozak, CAAT, GC box, … — is `MotifFinder.FindRegulatoryElements`, unit MOTIF-REGULATORY-001).

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
| INV-06 | `FindMotif` returns the same ascending list as `MotifFinder.FindExactMotif` (one engine) | B05/B06 cross-batch requests; oracle = `str.find` scan |
| INV-07 | Reverse-strand hit at `p` ⇔ window `T[p..p+n-1]` equals the (IUPAC) reverse complement of the motif; `p` is a forward-strand coordinate; a palindromic site is a hit on both strands at the same `p` | fuzznuc `-complement` [4] (2448-pair oracle) |
| INV-08 | With IUPAC codes a window matches iff every base is in the base set of the motif code at that position; without them IUPAC symbols are literal and match nothing in a DNA sequence | NC-IUB 1984; fuzznuc [4] |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `sequence` | `DnaSequence` | required | Subject sequence (text T) | Upper-cased ACGT (validated by `DnaSequence`); empty allowed |
| `motifs` | `IEnumerable<string>` | required | Set of query motifs | Non-null; individual motifs upper-cased; empty/whitespace motifs skipped |
| `bothStrands` | `bool` | `false` | `FindKnownMotifHits` only: also report reverse-strand hits | |
| `allowIupac` | `bool` | `false` | `FindKnownMotifHits` only: IUPAC codes in the motifs | Non-IUPAC symbols → `ArgumentException` |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| (return) | `Dictionary<string, IReadOnlyList<int>>` | Map from upper-cased motif → ascending list of 0-based start positions of all its occurrences. Motifs with no occurrence are absent. |
| (return, `FindKnownMotifHits`) | `Dictionary<string, IReadOnlyList<KnownMotifHit>>` | Same keys; `KnownMotifHit(Position, IsReverseStrand)` sorted by position (forward before reverse at a tie). `Position` is always the forward-strand window start. |

### 3.3 Preconditions and Validation

`sequence` or `motifs` null → `ArgumentNullException` (a null `sequence` used to raise `NullReferenceException`). Empty/whitespace motifs are skipped (the empty string is
not a motif — a suffix-tree search for `""` would return every position). Motifs are upper-cased
(case-insensitive matching; result keyed by the upper-cased form). Duplicate motifs collapsing to
the same upper-cased key produce a single entry. Indexing is 0-based; positions are starts of the
matched substring. Empty subject sequence → empty result (no motif occurs).

## 4. Algorithm

### 4.1 High-Level Steps

1. Validate `sequence` and `motifs` are non-null.
2. For each motif: skip if empty/whitespace; upper-case it; skip if its key was already processed.
3. Enumerate all start positions (overlaps included) with one of two exact engines, chosen once per call:
   fewer than 128 distinct queries (motifs × strands) → one KMP scan of the text per query
   (`KmerAnalyzer.FindKmerPositions`, ascending by construction); 128 or more → the sequence's cached suffix
   tree (`FindAllOccurrences`, DFS order, then sorted). With `allowIupac`, `MotifFinder.FindDegenerateMotif`.
4. With `bothStrands`, repeat for the IUPAC reverse complement of the motif (`DnaSequence.GetReverseComplementString`)
   and merge by position.
5. If a motif has a hit, store it under the upper-cased key.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

Two interchangeable exact engines. The repository **suffix tree** (`SuffixTree`) is an index: built once for
the subject (cached on the `DnaSequence`) and reused across all queries; `FindAllOccurrences` enumerates leaf
start positions in DFS order (not sorted), so positions are sorted for the deterministic contract. The
**KMP scan** (Knuth–Morris–Pratt 1977; CLRS §32.4) needs no index. Measured on random DNA (Release):

| n | suffix-tree build | one KMP scan (6-mer) | break-even |
|---|---|---|---|
| 100 kb | 62 ms | 0.85 ms | 73 scans |
| 1 Mbp | 924 ms (171 MB) | 6.9 ms | 135 scans |
| 4 Mbp | 4.8 s (684 MB) | 24.5 ms | 197 scans |
| 10 Mbp (B06) | 23.3 s (1.7 GB) | 71 ms | ≈ 330 scans |

The query threshold 128 distinct queries sits at the measured break-even; both engines return identical positions
(tested against `string.IndexOf`).

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| k < 128 queries: KMP scans | O(k·n + Σ\|mᵢ\| + occ) | O(max \|m\| + occ) | No index; one pass per query |
| k ≥ 128 queries: suffix-tree construction | O(n) | O(n) | Built once for the subject (Ukkonen); cached on `DnaSequence` |
| k ≥ 128: per-motif search + sort | O(\|m\| + occ·log occ) | O(occ) | `FindAllOccurrences` is O(\|m\|+occ); sort dominates when occ large |
| `GenomicAnalyzer.FindMotif` (= `MotifFinder.FindExactMotif`) | O(n) build on first use, then O(\|m\| + occ·log occ) | O(n) | Always the suffix-tree index; for a one-off query on a very large sequence call `KmerAnalyzer.FindKmerPositions` (same positions, no index) |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [GenomicAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/GenomicAnalyzer.cs)

- `GenomicAnalyzer.FindKnownMotifs(DnaSequence, IEnumerable<string>)`: multi-motif exact search returning per-motif sorted 0-based position lists (forward strand, exact).
- `GenomicAnalyzer.FindKnownMotifHits(DnaSequence, IEnumerable<string>, bool bothStrands = false, bool allowIupac = false)`: the strand-aware / IUPAC generalisation returning `KnownMotifHit(Position, IsReverseStrand)`.
- `GenomicAnalyzer.FindMotif(DnaSequence, string)`: single exact motif, ascending; delegates to `MotifFinder.FindExactMotif` (B05/B06 dedup, F18).

### 5.2 Current Behavior

**Engine decision.** This is exact-match occurrence enumeration of *multiple* fixed queries against *one*
text. A suffix tree amortises its O(n) build over many queries, but its build is heavy (see §4.2): for the
common case of a handful of motifs the linear scans are 70–330× cheaper in time and need no index memory, so
`FindKnownMotifs` uses the index only from 128 distinct queries on. Behaviour does not depend on the engine.

**Ordering.** `SuffixTree.FindAllOccurrences` returns occurrences in DFS order (Rosalind SUBS sample
`ATAT` in `GATATATGCATATACTT` enumerates as 1, 9, 3); every public entry point returns ascending positions
(INV-03, INV-06). `GenomicAnalyzer.FindMotif` used to return the tree order and duplicated
`MotifFinder.FindExactMotif`; it now delegates (F18).

**Strands and IUPAC (F19).** `FindKnownMotifHits` reproduces EMBOSS fuzznuc 6.6.0 (`-pmismatch 0`,
with/without `-complement`): reverse-strand hits are windows equal to the reverse complement of the motif,
reported in forward-strand coordinates (EcoRI `GAATTC` → both strands at the same position); IUPAC codes
use the NC-IUB base sets and the IUPAC reverse complement (R↔Y, K↔M, B↔V, D↔H).

**Mismatches (A1-20).** `FindKnownMotifHits(…, maxMismatches: k)` is fuzznuc `-pmismatch k`: a window is a hit when
at most k of its m positions disagree with the motif (Hamming distance, substitutions only). The scan is the canonical
`ApproximateMatcher.FindWithMismatches` (one scan per motif and strand, O(n·m), no index; the exact k = 0 route keeps the
KMP / suffix-tree engines). IUPAC codes: with `allowIupac` a position agrees when the base is in the code's base set (the
distance is the Hamming distance to the nearest concrete expansion), implemented as the additive
`ApproximateMatcher.FindWithMismatches(…, iupacPattern: true)` overload (per-position membership =
`IupacHelper.MatchesIupac`); without it every motif symbol is literal (a degenerate symbol costs one mismatch at each
position). Reverse strand: the IUPAC reverse-complement motif with the same bound; each hit carries its own distance
(`KnownMotifHit.Mismatches`, = fuzznuc's `*mismatch n` note; a non-palindromic motif can hit on the two strands with
different distances). A bound ≥ m makes every window a hit (fuzznuc refuses such values); a sequence shorter than the
motif has no window and no hit. fuzznuc 6.6.0 does print "hits" there (`AAA` with `-pmismatch 2` on `AC` gives `1..3`, `*mismatch 2`), but they are
artefacts of its matcher, decided against the binary (audit A2-18, 300 runs on sequences of 1 .. m-1 bases, `fuzznuc_short_sequence_oracle.tsv`): only 49 of the 300 runs print
anything, 41 of the printed hits end beyond the end of the sequence (GFF `end` > sequence length), 10 reverse-strand lines have `end` < `start` (e.g. `2..-1`), and patterns with IUPAC
codes print nothing where the same window would be within the bound; Seqeron reports none (`GenomicAnalyzer_KnownMotifs_ShortSequence_Tests`). Oracle: fuzznuc 6.6.0 `-pmismatch 1..3 -complement Y` (positions and `*mismatch` counts) = an
independent Python Hamming scan on 6576 (sequence, pattern, k) rows (20 patterns, 3932 rows with forward and 3785 with
reverse hits, 4183 with mismatching hits).

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- All start positions where `T[i..i+n-1] == P` reported, for each motif (exact-matching set) [1].
- Overlapping occurrences all reported (e.g. `AAA` in `AAAAA` → {0,1,2}) [1][2].
- Multi-motif set: each motif mapped to its own occurrence set [2].
- Upper-cased processing of DNA / motifs [2].

**The fuzznuc pattern language (A2-32).** `FuzznucPattern.Parse` reads EMBOSS 6.6.0's `-pattern` syntax (`fuzznuc -help`): IUPAC letters, classes `[ACG]`, exclusions `{AG}`, repeats `N(3)`, ranges `A(2,4)` (n may be 0), optional `-` separators, `<` (start of the sequence), `>` (end) and a final `.`;
`[CG](5)TG{A}N(1,5)C` is the help's own example. Classes and exclusions list plain bases (the binary rejects `[RT]`, `{N}`); `{ACGT}` matches nothing; malformed spellings that the binary accepts in unspecified ways (`A<C`, `{}`, `A(0,0)`) are rejected with an `ArgumentException`. Every
non-empty base set is exactly one of the 15 IUPAC codes, so a pattern expands to a finite set of fixed-length IUPAC strings (`Expand()`: `AC(2,3)G` → `ACCG`, `ACCCG`) that the canonical `ApproximateMatcher` Hamming scan (IUPAC pattern codes) runs over; hits are united per window
(`FindPatternHits(sequence, pattern, maxMismatches, bothStrands)` → `PatternHit(Position, Length, IsReverseStrand, Mismatches)`; a range reports every matching length; the reverse strand is the reversed, complemented pattern with swapped anchors; `FuzznucPattern.ParsePatternFile` and the
named-pattern overload cover `-pattern @file`: one pattern per line or FASTA `>name`). No new scanning loop exists. The binary's own search classes are reproduced exactly (EMBOSS `embpat.c`, read on raw.githubusercontent.com, and confirmed on 3300 runs): with mismatches a range pattern is searched by the
recursive `patOUBrute`, which tolerates mismatches only in the mandatory part (the first n repeats of `X(n,m)` and the fixed elements — the extra repeats must match) and fails when pattern elements remain at the end of the scanned strand (a path whose last element is an empty optional element
never ends on the last base); an exact range pattern starting with `<` is searched by a regex that returns only the longest match at the anchor. A sequence symbol other than A/C/G/T matches nothing (the B05 contract, also for the string overloads of `FindMotif` / `FindKnownMotifs` / `FindKnownMotifHits`);
the binary reads a sequence ambiguity code as a base set and lets a pattern `N` match it — not reproduced. Oracle: `scripts/analysis/generate_fuzznuc_pattern_oracle.py` → `fuzznuc_pattern_oracle.tsv`, 3300 random syntax-tree (sequence, pattern, k) rows, fuzznuc 6.6.0 `-pmismatch k -complement Y` == an independent Python reference, 0 disagreements
(`GenomicAnalyzer_PatternHits_Oracle_Tests`).

**Intentionally simplified:**

- (none)

**Not implemented:**

- (Implemented by audit A2-32: the fuzznuc pattern language is parsed by `FuzznucPattern` and searched by `GenomicAnalyzer.FindPatternHits`; see "The fuzznuc pattern language" below.)
  Edit-distance matching is `ApproximateMatcher.FindWithEdits` (PAT-APPROX-002) and position-weight scanning is `MotifFinder` PWM (PAT-PWM-001).
- Motif *discovery* and the fixed regulatory-element library: other units (MOTIF-DISCOVER-001, MOTIF-REGULATORY-001).

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

Exact matching by default; `FindKnownMotifHits(maxMismatches)` adds substitution-tolerant matching (no gaps, O(n·m) per motif and strand, no index). `FindKnownMotifs` / `FindMotif` search the forward strand and treat
IUPAC symbols as literals (they match nothing in a DNA sequence); use `FindKnownMotifHits` for degenerate motifs
and the reverse strand. In exact mode non-ACGT characters in a motif simply fail to match (no positions) rather
than raising; in IUPAC mode a non-IUPAC symbol throws `ArgumentException`. The subject is a `DnaSequence`
(ACGT only), so an ambiguity symbol in the *sequence* cannot occur.

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
var seq = new DnaSequence("GAATTCAAAGAATTC");
var hits = GenomicAnalyzer.FindKnownMotifs(seq, new[] { "GAATTC" });
// hits["GAATTC"] == [0, 9]  (EcoRI sites)
```

```csharp
var hits = GenomicAnalyzer.FindKnownMotifHits(
    new DnaSequence("GAATTCAAAGAATTCTTTATATAAAGTTCAAGAAC"), new[] { "TATAWAW" }, bothStrands: true, allowIupac: true);
// hits["TATAWAW"] == [ (16, reverse), (17, forward) ]   (fuzznuc: 18-24 '+', 17-23 '-', 1-based)
```

**Numerical walk-through:** For `T = "AAAAA"`, motif `"AAA"` (length 3): alignments start at
i = 0 (`AAA`), 1 (`AAA`), 2 (`AAA`); i = 3 would need indices 3..5 but T has length 5 (max start 2).
Result: positions {0, 1, 2}, three overlapping occurrences [1].

### 7.3 Related Tests, Evidence, or Documents

- Oracle: `scripts/analysis/generate_known_motif_oracle.py` → `TestData/Analysis/known_motif_oracle.tsv` (EMBOSS fuzznuc 6.6.0 == Python `re`, 2448 pairs; 816 with forward hits, 776 with reverse hits), tests in [GenomicAnalyzer_KnownMotifs_Oracle_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Analysis/GenomicAnalyzer_KnownMotifs_Oracle_Tests.cs)
- Tests: [GenomicAnalyzer_FindKnownMotifs_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Analysis/GenomicAnalyzer_FindKnownMotifs_Tests.cs) — covers INV-01..INV-05
- Evidence: [GENOMIC-MOTIFS-001-Evidence.md](../../../docs/Evidence/GENOMIC-MOTIFS-001-Evidence.md)
- Related algorithms: [Repeat_Detection](../Repeat_Analysis/Repeat_Detection.md)

## 8. References

1. Gusfield, D. 1997. *Algorithms on Strings, Trees and Sequences: Computer Science and Computational Biology*. Cambridge University Press. ISBN 0-521-58519-8. Exact-matching definition and overlapping-occurrence example (Tufts COMP 150GEN): https://www.cs.tufts.edu/comp/150GEN/classpages/exact.html
2. Cock, P.J.A. et al. Biopython `Bio.Seq` module — `search` and `count_overlap` (master, accessed 2026-06-13): https://raw.githubusercontent.com/biopython/biopython/master/Bio/Seq.py
3. Wikipedia contributors. "Restriction site" (accessed 2026-06-13). EcoRI recognizes GAATTC: https://en.wikipedia.org/wiki/Restriction_site
4. EMBOSS 6.6.0 `fuzznuc` (Rice, Longden & Bleasby 2000, Trends Genet 16:276; Debian `emboss` package, run locally with `-pmismatch 0 -complement Y -rformat gff`): IUPAC patterns, both strands, overlapping hits, forward-strand coordinates for '-' hits.
