# Common Motif Finding

| Field | Value |
|-------|-------|
| Algorithm Group | ProteinMotif |
| Test Unit ID | PROTMOTIF-COMMON-001 |
| Related Projects | Seqeron.Genomics.Analysis |
| Implementation Status | Production |
| Last Reviewed | 2026-06-14 |

## 1. Overview

Common Motif Finding scans a protein sequence against a curated library of well-known
functional motifs expressed as PROSITE-style patterns and reports every occurrence of
every pattern. It answers "which of the common functional signatures (N-glycosylation,
phosphorylation sites, P-loop, RGD, …) occur in this protein, and where?". The result is
specification-driven (exact pattern matching), not probabilistic: a window either satisfies
a PROSITE pattern or it does not [1][2].

## 2. Scientific / Formal Basis

### 2.1 Domain Context

Short functional motifs (active sites, post-translational-modification sites, binding
signatures) are catalogued in PROSITE as regular-expression-like patterns over the 20-letter
amino-acid alphabet [1]. A pattern is a sequence of positional elements; each element fixes
the residue(s) allowed at one sequence position [2].

### 2.2 Core Model

PROSITE pattern element syntax (verbatim from the ScanProsite documentation) [2]:

- A single IUPAC letter (e.g. `G`) — that residue only.
- `[..]` — an *allowed set*: any one of the listed residues (e.g. `[ST]` = Ser or Thr).
- `{..}` — an *excluded set*: any residue except those listed (e.g. `{P}` = any but Pro).
- `x` — any residue.
- `x(n)` — exactly `n` arbitrary residues (`x(2)` = `x-x`).
- `x(n,m)` — between `n` and `m` arbitrary residues.
- `-` separates consecutive elements.

The library patterns used here, each verified against its official ExPASy entry:

| Accession | Name | Pattern | Source |
|-----------|------|---------|--------|
| PS00001 | ASN_GLYCOSYLATION | `N-{P}-[ST]-{P}` | [3] |
| PS00005 | PKC_PHOSPHO_SITE | `[ST]-x-[RK]` | [4] |
| PS00006 | CK2_PHOSPHO_SITE | `[ST]-x(2)-[DE]` | [5] |
| PS00016 | RGD | `R-G-D` | [6] |
| PS00017 | ATP_GTP_A | `[AG]-x(4)-G-K-[ST]` | [7] |

(The full library additionally includes PS00004, PS00007–PS00009, PS00018, PS00028, PS00029
and several literature-derived motifs; all are cited inline in the source file.)

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | Each reported match's substring equals the residues at its reported coordinates | The match is a literal window of the input scanned left to right [2] |
| INV-02 | `0 ≤ Start ≤ End < length` for every match | Matches are in-sequence hits [2] |
| INV-03 | Overlapping occurrences are all reported (unless one is fully contained in another) | PROSITE/ScanProsite default "overlaps, no includes" [2] |
| INV-04 | The scan is deterministic | Regex matching over a fixed pattern set is deterministic |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| proteinSequence | `string` | required | Protein sequence to scan | Upper-cased internally; non-residue characters simply fail to match patterns |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| (each) `MotifMatch.Start` | `int` | 0-based, inclusive start of the matched window |
| `MotifMatch.End` | `int` | 0-based, inclusive end of the matched window |
| `MotifMatch.Sequence` | `string` | The matched residues |
| `MotifMatch.MotifName` | `string` | PROSITE entry name (e.g. `RGD`) |
| `MotifMatch.Pattern` | `string` | PROSITE accession (e.g. `PS00016`) |
| `MotifMatch.Score` / `EValue` | `double` | Information-content score / expected random count (defined under PROTMOTIF-PATTERN-001) |

### 3.3 Preconditions and Validation

Null or empty input yields an empty sequence (no exception). Coordinates are 0-based and
inclusive (repository convention; PROSITE itself reports 1-based [2]). Input is upper-cased
before matching (case-insensitive). Alphabet is protein single-letter codes; the patterns are
PROSITE-defined over the 20 standard amino acids [1].

## 4. Algorithm

### 4.1 High-Level Steps

1. If the sequence is null/empty, return no matches.
2. Upper-case the sequence.
3. For each pattern in the `CommonMotifs` library, find every occurrence in the sequence.
4. Yield each occurrence with the originating pattern's name and accession.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

The reference table is the `CommonMotifs` dictionary (PROSITE accession → pattern). Each entry
stores both the PROSITE pattern string and the equivalent .NET regex. Per-pattern matching
(including overlapping-occurrence discovery via a `(?=(…))` lookahead wrapper) is delegated to
`FindMotifByPattern` (PROTMOTIF-PATTERN-001) [2].

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `FindCommonMotifs` | O(p · n) | O(occurrences) | p = number of library patterns (constant), n = sequence length; each pattern is a bounded-width regex scan |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [ProteinMotifFinder.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/ProteinMotifFinder.cs)

- `ProteinMotifFinder.FindCommonMotifs(string)`: iterates `CommonMotifs.Values` and yields all matches per pattern.
- `ProteinMotifFinder.CommonMotifs`: the curated PROSITE-style pattern library.

### 5.2 Current Behavior

Returns matches grouped by pattern in dictionary-iteration order; within a pattern, matches
are in increasing start order. Matching is case-insensitive (input upper-cased). Overlapping
occurrences are reported via the lookahead engine in `FindMotifByPattern`.

**Search-infrastructure decision (suffix tree):** the repository `SuffixTree`
(`FindAllOccurrences`/`CountOccurrences`) performs *exact substring* matching only and cannot
express PROSITE pattern elements — allowed sets `[ST]`, excluded sets `{P}`, wildcards `x`, or
variable gaps `x(n,m)`. Common-motif finding is therefore class-based pattern matching, not
literal substring search, so the suffix tree is **not used**; matching is delegated to the
regex-based `FindMotifByPattern`. Using the suffix tree would require enumerating every literal
expansion of each pattern, which is exponential and incorrect for variable-length elements.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- PROSITE element semantics `[..]`, `{..}`, `x`, `x(n)`, `x(n,m)`, `-` [2].
- Library patterns PS00001/PS00005/PS00006/PS00016/PS00017 (and the rest of `CommonMotifs`) matching their official ExPASy entries [3][4][5][6][7].
- Overlapping-occurrence reporting (PROSITE default) [2].

**Intentionally simplified:**

- Coordinates are reported 0-based inclusive rather than PROSITE's 1-based; **consequence:** positions are shifted by −1 versus a ScanProsite report, but matched residues and relative positions are identical [2].

**Not implemented:**

- N-/C-terminal anchors (`<`/`>`) inside the curated library entries here (none of the bundled patterns are anchored); **users should rely on:** `FindMotifByProsite` / `ConvertPrositeToRegex` (PROTMOTIF-PROSITE-001) for arbitrary anchored PROSITE patterns.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | 0-based vs PROSITE 1-based coordinates | Assumption | Coordinate origin differs from ScanProsite output | accepted | Matched substring/relative positions unchanged; consistent with sibling units |
| 2 | `FindAllKnownMotifs` listed in Registry but absent in code | Deviation | Registry method-table naming only | accepted | Canonical method is `FindCommonMotifs`; no alias invented |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| null / empty input | empty result | guard clause |
| Pro at a `{P}` position | no N-glycosylation match | PS00001 excluded set [3] |
| two adjacent RGD sites | two matches reported | overlap/multi-occurrence reporting [2][6] |
| lowercase input | matched (case-insensitive) | input upper-cased |
| sequence with no motif | empty result | no window satisfies any pattern |

### 6.2 Limitations

The library is a fixed curated subset of PROSITE plus a few literature motifs; it is not the
full PROSITE database. Patterns are exact (no scoring profile / HMM); biologically these motifs
have high false-positive rates and PROSITE supplies skip-rules and profiles for some entries
that are not modelled here. For arbitrary patterns use `FindMotifByProsite`.

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
foreach (var m in ProteinMotifFinder.FindCommonMotifs("AARGDKK"))
    Console.WriteLine($"{m.MotifName} {m.Pattern} {m.Start}-{m.End} {m.Sequence}");
// RGD PS00016 2-4 RGD
```

**Numerical / biological walk-through:**

For `AAAANFTAAAA` and PS00001 `N-{P}-[ST]-{P}` (regex `N[^P][ST][^P]`): the only window that
satisfies all four elements is positions 4–7 = `N`,`F`(≠P),`T`(∈{S,T}),`A`(≠P) → match
`(Start=4, End=7, "NFTA")` [3].

### 7.3 Related Tests, Evidence, or Documents

- Tests: [ProteinMotifFinder_FindCommonMotifs_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Analysis/ProteinMotifFinder_FindCommonMotifs_Tests.cs) — covers INV-01..INV-04
- Evidence: [PROTMOTIF-COMMON-001-Evidence.md](../../../docs/Evidence/PROTMOTIF-COMMON-001-Evidence.md)
- Related algorithms: [PROSITE_Pattern_Matching](../ProteinMotif/PROSITE_Pattern_Matching.md)

## 8. References

1. Sigrist CJA, de Castro E, Cerutti L, et al. 2013. New and continuing developments at PROSITE. Nucleic Acids Research 41(D1):D344–D347. https://doi.org/10.1093/nar/gks1067
2. ExPASy. ScanProsite documentation — pattern syntax, overlap and coordinate reporting. https://prosite.expasy.org/scanprosite/scanprosite_doc.html (accessed 2026-06-14).
3. ExPASy PROSITE. PS00001 — N-glycosylation site (ASN_GLYCOSYLATION). https://prosite.expasy.org/PS00001 (accessed 2026-06-14).
4. ExPASy PROSITE. PS00005 — PKC phosphorylation site. https://prosite.expasy.org/PS00005 (accessed 2026-06-14).
5. ExPASy PROSITE. PS00006 — CK2 phosphorylation site. https://prosite.expasy.org/PS00006 (accessed 2026-06-14).
6. ExPASy PROSITE. PS00016 — Cell attachment sequence (RGD). https://prosite.expasy.org/PS00016 (accessed 2026-06-14).
7. ExPASy PROSITE. PS00017 — ATP/GTP-binding site motif A (P-loop). https://prosite.expasy.org/PS00017 (accessed 2026-06-14).
