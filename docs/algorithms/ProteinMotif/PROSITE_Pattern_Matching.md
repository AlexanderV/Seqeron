# PROSITE Pattern Matching

| Field | Value |
|-------|-------|
| Algorithm Group | ProteinMotif |
| Test Unit ID | PROTMOTIF-PROSITE-001 |
| Related Projects | Seqeron.Genomics |
| Implementation Status | Simplified |
| Last Reviewed | 2026-04-30 |

## 1. Overview

PROSITE pattern matching in this repository converts PROSITE pattern syntax into .NET regular expressions and then searches protein sequences with the converted regex. The focus is the formal PROSITE pattern language itself, including anchors, exclusions, repetitions, and the rare `[G>]` C-terminal bracket form documented by PROSITE. Matching is deterministic, case-insensitive, and overlap-aware because the converted regex is delegated to the repository's lookahead-based motif search helper. This document is limited to pattern syntax and matching behavior, not PROSITE profile or database services.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

PROSITE defines a formal pattern notation for biologically meaningful protein sequence motifs and describes that notation in the PA-line specification of the PROSITE User Manual. A PROSITE pattern is a declarative description of accepted residues, exclusions, anchors, and repeat counts, and pattern matching asks whether a sequence contains a subsequence satisfying that specification (PROSITE User Manual; Hulo et al., 2007).

### 2.2 Core Model

The core model has two stages. First, parse the PROSITE pattern from left to right and replace each syntax element with an equivalent regular-expression construct. Second, search the protein sequence for every occurrence of the converted regex. The current source and tests cover the standard PROSITE elements `x`, `[ABC]`, `{ABC}`, repetition markers, terminus anchors, and the rare `[G>]` bracket form noted in PROSITE examples such as PS00267 and PS00539 (PROSITE User Manual; ScanProsite documentation).

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | Literal residues, wildcards, character classes, exclusions, repetitions, and terminus anchors preserve their declared PROSITE meaning after conversion. | The converter maps each syntax element to a corresponding regex construct. |
| INV-02 | A period terminates the PROSITE pattern. | The PROSITE PA-line specification treats `.` as the end of the pattern. |
| INV-03 | `[G>]` denotes either a literal `G` or the end of the sequence. | The PROSITE User Manual defines `>` inside the bracketed final element as a C-terminal alternative. |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `prositePattern` | `string` | required | PROSITE-format pattern string | Null or empty input converts to the empty regex string |
| `proteinSequence` | `string` | required for `FindMotifByProsite(...)` | Protein sequence searched with the converted regex | Null or empty input yields no matches |
| `motifName` | `string` | `Custom` | Name stored in each returned `MotifMatch` | Caller-supplied label |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `ConvertPrositeToRegex(...)` | `string` | Converted .NET regex string |
| `Start` | `int` | Inclusive 0-based start index of a PROSITE match |
| `End` | `int` | Inclusive 0-based end index of a PROSITE match |
| `Sequence` | `string` | Uppercased matched substring |
| `MotifName` | `string` | Caller-supplied motif label |
| `Pattern` | `string` | Original PROSITE pattern string |
| `Score` | `double` | Repository motif score from the delegated regex search |
| `EValue` | `double` | Repository E-value from the delegated regex search |

### 3.3 Preconditions and Validation

`ConvertPrositeToRegex(...)` returns the empty string for null or empty input. `FindMotifByProsite(...)` returns no matches for null or empty sequence or pattern input. Matching is case-insensitive because the repository uppercases the sequence before scanning. Coordinates are inclusive 0-based indexes. If conversion produced a malformed regex, the delegated matcher would return no hits rather than throw, because `FindMotifByPattern(...)` catches regex compilation failures.

## 4. Algorithm

### 4.1 High-Level Steps

1. Read the PROSITE pattern from left to right.
2. Translate each syntax element into its regex equivalent, dropping separators and stopping at the first terminating period.
3. Pass the converted regex to `FindMotifByPattern(...)` together with the input sequence and motif name.
4. Return every captured match span from the delegated regex scan.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

| Syntax Element | Meaning | Regex Equivalent |
|---------------|---------|------------------|
| `A` | Specific amino acid | `A` |
| `x` | Any amino acid | `.` |
| `[ABC]` | Any of the listed residues | `[ABC]` |
| `{ABC}` | Any amino acid except the listed residues | `[^ABC]` |
| `-` | Element separator | dropped |
| `(n)` | Repeat the preceding element `n` times | `{n}` |
| `x(n,m)` | Wildcard repeated from `n` to `m` positions | `.{n,m}` |
| `<` | N-terminus anchor | `^` |
| `>` | C-terminus anchor | `$` |
| `.` | Pattern terminator | parsing stops |
| `[G>]` | `G` or end-of-sequence in the final element | `(?:G|$)` |

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `ConvertPrositeToRegex(...)` | `O(p)` | `O(p)` | `p` = PROSITE pattern length; the parser is a single left-to-right pass |
| `FindMotifByProsite(...)` | `O(p + scan)` | `O(k)` | `k` = number of matches; after `O(p)` conversion, matching is delegated to the regex search helper |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [ProteinMotifFinder.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/ProteinMotifFinder.cs)

- `ProteinMotifFinder.ConvertPrositeToRegex(string)`: translates PROSITE notation into a .NET regex string.
- `ProteinMotifFinder.FindMotifByProsite(string, string, string)`: converts the PROSITE pattern and delegates matching to `FindMotifByPattern(...)`.
- `ProteinMotifFinder.FindMotifByPattern(string, string, string, string)`: supplies the overlap-aware regex matcher used by the end-to-end PROSITE helper.

### 5.2 Current Behavior

Repository-specific behavior confirmed by source and tests:

- `ConvertPrositeToRegex(...)` returns `""` for null or empty input.
- The parser drops `-` separators, converts exclusions to negated character classes, emits regex quantifiers for repetitions, and stops processing at the first `.`.
- `ConvertPrositeToRegex(...)` contains explicit handling for `>` inside a bracketed terminal element and emits a non-capturing alternation such as `(?:G|$)` for `[G>]`.
- `FindMotifByProsite(...)` stores the original PROSITE pattern in `MotifMatch.Pattern` and delegates the converted regex to `FindMotifByPattern(...)`.
- End-to-end PROSITE matching is overlap-aware because the delegated matcher uses a lookahead wrapper.
- Returned coordinates are inclusive 0-based indexes, and matching is case-insensitive after uppercasing the input sequence.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Standard PROSITE letters, wildcard `x`, character classes, exclusion classes, separators, repetitions, and terminus anchors.
- Termination of pattern parsing at the first `.`.
- Rare `[G>]` handling for patterns such as PS00267 and PS00539.

**Intentionally simplified:**

- The helper covers PROSITE pattern syntax only; **consequence:** PROSITE profile or matrix entries are outside the scope of this implementation.
- End-to-end matching reuses the repository's generic motif score and E-value calculation; **consequence:** `Score` and `EValue` are repository-defined outputs rather than ScanProsite statistics.
- Conversion targets .NET regex constructs directly; **consequence:** final matching behavior follows the generated regex engine semantics rather than a standalone PROSITE runtime.

**Not implemented:**

- PROSITE profile matching or database-backed annotation workflows; **users should rely on:** external ScanProsite or related PROSITE tooling.
- Alternate output conventions such as 1-based coordinates; **users should rely on:** no current alternative in this repository.

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Empty or null pattern | Conversion returns `""`; matching returns no hits | Both entry points guard the empty-pattern case |
| Trailing period | Content after `.` is ignored | The parser stops at the first terminator |
| `[G>]` at the sequence end | Match may end exactly at the sequence boundary | `[G>]` maps to `G` or end-of-sequence |
| `[G>]` in the middle of a sequence without `G` | No match | The `$` branch cannot match mid-sequence |
| Empty protein sequence | Returns no matches | The delegated motif search short-circuits empty input |

### 6.2 Limitations

This helper is a syntax converter plus regex search wrapper. It does not implement PROSITE profile scanning, does not fetch entries from an external catalog, and does not provide ScanProsite-specific result metadata. It should therefore be used when the task is explicit PROSITE pattern matching rather than full PROSITE annotation.

## 8. References

1. PROSITE User Manual. SIB Swiss Institute of Bioinformatics. https://prosite.expasy.org/prosuser.html
2. ScanProsite Documentation. https://prosite.expasy.org/scanprosite/scanprosite_doc.html
3. Hulo N, Bairoch A, Bulliard V, et al. The 20 years of PROSITE. Nucleic Acids Research. https://doi.org/10.1093/nar/gkm977
4. De Castro E, Sigrist CJA, Gattiker A, et al. ScanProsite. Nucleic Acids Research. https://doi.org/10.1093/nar/gkl124
5. PROSITE PS00267. Tachykinin family signature. https://prosite.expasy.org/PS00267
6. PROSITE PS00539. Pyrokinins signature. https://prosite.expasy.org/PS00539
