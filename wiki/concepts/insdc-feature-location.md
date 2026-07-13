---
type: concept
title: "INSDC feature-table location descriptors (complement/join/order/partial/remote)"
tags: [file-io, algorithm]
sources:
  - docs/Evidence/PARSE-GENBANK-001-Evidence.md
  - docs/Evidence/PARSE-EMBL-001-Evidence.md
  - docs/Validation/reports/PARSE-GENBANK-001.md
  - docs/algorithms/FileIO/EMBL_Parsing.md
  - docs/algorithms/FileIO/GenBank_Parsing.md
source_commit: 0a11c83c145d84b916029acb9bf18c179115ba0e
created: 2026-07-10
updated: 2026-07-13
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: parse-genbank-001-evidence
      evidence: "PARSE-GENBANK-001 and PARSE-EMBL-001 ... Feature Location Syntax (from INSDC/NCBI): n..m / complement(n..m) / join(ranges) / <n..m / n..>m"
      confidence: high
      status: current
---

# INSDC feature-table location descriptors

The **INSDC** (International Nucleotide Sequence Database Collaboration — DDBJ, ENA/EMBL, GenBank)
publishes one shared **Feature Table Definition**: the same feature keys, qualifiers, and
**location-descriptor grammar** underlie both the GenBank and the EMBL flat-file formats. The two
formats differ only in *line syntax* (GenBank's keyword sections vs EMBL's two-character line-type
codes) — a feature's **location string is identical** in both. This page is the shared reference
for that location grammar, so it does not have to be re-explained per dialect. It is traced in two
`PARSE-*` [[test-unit-registry|test units]]: [[parse-genbank-001-evidence]] (GenBank) and
[[parse-embl-001-evidence]] (EMBL), the two flat-file members of the file-parsing family anchored
at [[bed-format-parsing]].

## Why a location is not just a start/end

A GenBank/EMBL feature does not carry a plain `[start, end)` interval like [[bed-format-parsing]].
Its location is a small **expression** that may reverse strand, stitch discontiguous spans, mark
incompleteness, and even point into another record. Coordinates in a descriptor are **1-based,
inclusive** (the first base is position 1), which is the opposite of BED's 0-based half-open
convention — a coordinate-conversion hazard when moving features between formats.

## The descriptor grammar

| Form | Example | Meaning |
|------|---------|---------|
| `n..m` | `100..200` | Contiguous range, positions n through m inclusive |
| `n` | `467` | A single base |
| `n^m` | `123^124` | A site *between* two adjacent bases (EMBL/INSDC; e.g. an insertion point) |
| `<n..m` | `<1..206` | 5′-partial — the feature starts before position n |
| `n..>m` | `4821..>5028` | 3′-partial — the feature extends past position m |
| `<n..>m` | `<1..>500` | Partial at both ends |
| `complement(loc)` | `complement(3300..4037)` | `loc` is on the minus (reverse) strand |
| `join(loc,…)` | `join(1..50,60..100)` | One feature spanning several spans, concatenated end-to-end in listed order |
| `order(loc,…)` | `order(1..50,60..100)` | Same elements, but order/adjacency is **not** asserted |
| `complement(join(…))` | `complement(join(1..50,60..100))` | Operators nest |

**Nesting rule (from the INSDC spec):** `complement` may wrap `join`/`order`, but a `join`/`order`
may **not** be nested inside another `join`/`order` in the same location.

**Remote references** — a span may name another entry: `accession[.version]:span`, e.g.
`J00194.1:100..202`, and may appear inside an operator (`join(1..100,J00194.1:100..202)`). Remote
spans are 1-based inclusive on the *remote* entry.

## Assembling the sequence (operator semantics)

To materialize a feature's nucleotide sequence from its location, the operators compose as (rules
taken verbatim from the INSDC feature-table definition):

- `join(a,b,…)` / `order(a,b,…)` — concatenate the element subsequences in listed order.
- `complement(loc)` — reverse-complement the assembled string for `loc`.
- The two interact so that `complement(join(a,b)) == join(complement(b), complement(a))` — the
  outer complement **reverses element order** as well as strand.
- Worked oracles from the EMBL evidence: `join(1..3,7..9)` over `ACGTACGTAC` → `ACGGTA`;
  `complement(Y.1:1..4)` with remote `Y → AAAC` → `GTTT`.

Because the library is **offline-first and does no network I/O**, remote references are resolved by
a **caller-supplied resolver delegate** (`FeatureLocationHelper.ResolveLocationSequence` /
`EmblParser.ResolveLocationSequence`); a missing/null resolver contributes an empty segment while
local segments are still assembled in place. See [[parse-embl-001-evidence]] for the resolver
mechanics and the remote-reference per-segment prefix-strip fix, and [[parse-genbank-001-evidence]]
for the GenBank `LOCUS`/`ORIGIN`/division specifics.

The GenBank parser's two-stage validation verdict is recorded in [[parse-genbank-001-report]] —
**Stage A PASS (the 1-based-inclusive / `complement` / `join`-vs-`order` / partial-marker grammar is
faithful to INSDC/NCBI), Stage B PASS-WITH-NOTES, State CLEAN**. The one code defect was outside the
location grammar — multi-line qualifier reconstruction inserted a spurious space in `/translation`
and left an unstripped opening quote on wrapped values — fixed to the Biopython reference behaviour
(wrap→single space; `/translation`→no space) and locked with two RED→GREEN tests. Full suite 18213/0.

## The GenBank parser surface (`GenBankParser`)

The GenBank dialect's parser is specified in `docs/algorithms/FileIO/GenBank_Parsing.md`
(test unit **PARSE-GENBANK-001**, Implementation Status *Simplified*). Beyond the shared
location grammar above, it handles the GenBank flat-file **keyword-section record shape**: sections
begin with a keyword in **column 1** and the record terminates with `//` (INV-01). `GenBankParser.Parse(string)`
/ `ParseFile(string)` split input on `\n//`, keep only trimmed blocks that begin with `LOCUS`, then
parse each top-level section by its header keyword (`LOCUS`, `DEFINITION`, `ACCESSION`, `VERSION`,
`KEYWORDS`, `SOURCE`+nested `ORGANISM`, `REFERENCE`, `FEATURES`, `ORIGIN`), merging continuation
lines into one logical string per section. Non-core sections land in `AdditionalFields`.

Unlike EMBL's single delimited `ID` line, the **`LOCUS` line is fixed-column positional** — parsed
into locus name (cols 13–28), length + `bp`/`aa` (30–40), molecule type (45–47), topology
`linear`/`circular` (56–63), GenBank division (three-letter, 65–67), and date `DD-MMM-YYYY` (69–79,
also accepts `DD-MMM-YY`). The division vocabulary is GenBank's 18-code set (PRI, ROD, MAM, VRT, INV,
PLN, BCT, VRL, PHG, SYN, UNA, EST, PAT, STS, GSS, HTG, HTC, ENV) — overlapping but **not** identical
to EMBL's. `KEYWORDS` of just `.` normalizes to an empty list; the organism name comes from `SOURCE`
and the taxonomy lineage from the indented `ORGANISM` subsection. The `ORIGIN` block stores the
sequence as position-numbered 60-base lines; parsing strips digits/whitespace and **uppercases** the
remaining letters.

Key entry points (all `O(n)` in text length; `ParseLocation`/`ExtractSequence`/`TranslateCDS` linear
in the location/subsequence length):

- `ParseLocation(string)` — parses one INSDC location string (empty → zeroed `Location` with no parts).
- `GetFeatures(...)` / `GetCDS(...)` / `GetGenes(...)` — feature selection helpers.
- `ExtractSequence(...)` — 1-based-inclusive local-parts subsequence extraction (`realStart=partStart-1`,
  `realEnd=partEnd`), join concatenation, and IUPAC-aware reverse-complement for `complement(...)`, via
  the shared `FeatureLocationHelper`.
- `GetQualifier(...)` — reads one qualifier from a parsed feature.
- `TranslateCDS(...)` — **GenBank-specific** (no EMBL analog): returns an existing `/translation`
  qualifier verbatim when present, otherwise extracts the CDS nucleotides for the location and
  translates them with a **built-in standard codon table** (unknown codons → `X`). It therefore does
  **not** vary by alternative genetic-code metadata.

The location parser is the same `SequenceFormatHelper.ParseLocationParts(...)` reused by the EMBL
path. **Simplified by design:** feature qualifiers are flattened to key/value strings assembled from
continuation lines, so original quoting/formatting is not round-tripped; there is no validation that
the declared `LOCUS` length matches the parsed `ORIGIN` sequence; sequence case is normalized to
uppercase. GenBank-specific format/testing facts (18 division codes, U49845 canonical record,
defensive contracts) live in the [[parse-genbank-001-evidence]] source page.

## The EMBL parser surface (`EmblParser`)

The EMBL dialect's parser is specified in `docs/algorithms/FileIO/EMBL_Parsing.md`
(test unit **PARSE-EMBL-001**, Implementation Status *Simplified*). Beyond the shared
location grammar above, it handles the EMBL flat-file **line-oriented record shape**:
each line is a two-character code plus content, and records terminate with `//` (INV-01).
`EmblParser.Parse(string)` / `ParseFile(string)` split input on `\n//`, keep only trimmed
blocks that begin with `ID`, then `GroupLinesByPrefix(...)` concatenates repeated same-prefix
lines before field parsing. The `ID` line is decoded as
`accession; SV n; topology; mol_type; data_class; division; length BP`; when `SequenceVersion`
is absent from `ID` it falls back to the `SV` line, and `Accession` falls back to the first
`AC` line. Recognised line codes: `ID`, `AC`, `SV`, `DE`, `KW`, `OS`, `OC`, `RN`, `RA`, `RT`,
`RL`, `FT`, `SQ`. Non-consumed groups (e.g. `DT`, `DR`, `CC`, `OG`) are preserved in
`AdditionalFields` rather than dropped. Sequence letters are extracted from the `SQ` section,
stripped of digits/spaces and uppercased.

Key entry points (all `O(n)` in text length; `ToGenBank` is `O(features + references)`):

- `ParseLocation(string)` — parses one INSDC location string (empty → zeroed `Location`).
- `GetFeatures(...)` / `GetCDS(...)` / `GetGenes(...)` — feature selection helpers.
- `ExtractSequence(...)` — local-parts subsequence extraction, 1-based inclusive, via the
  shared `FeatureLocationHelper` (joins + reverse-complement).
- `ResolveLocationSequence(record, location, resolver)` — the opt-in remote-aware assembly
  described above (caller-supplied `RemoteSequenceResolver`; no network I/O in the library).
- `ToGenBank(...)` — converts an `EmblRecord` into the repository's GenBank record shape,
  so the two INSDC dialects interoperate through one in-memory model.

**Simplified by design:** the parser preserves the main parsed fields but is **not** a
full-fidelity round-trip serializer — feature qualifiers are flattened to a string dictionary
(bare qualifiers → `"true"`), original quoting/line-layout is not reproduced, and malformed
record separators or non-`ID` preambles are skipped rather than repaired. Full EMBL
occurrence-count validation and `SQ` composition cross-checking are not implemented. The
shared location parser is `SequenceFormatHelper.ParseLocationParts(...)`, reused by both the
GenBank and EMBL paths.

## Partial-flag slicing (assumption)

The `<` / `>` partial markers assert biological incompleteness but carry **no extra coordinate** —
the parser slices the stated number verbatim and the partial flag does not move the slice bounds.
This is an implementation assumption recorded in the EMBL evidence, not an INSDC-prescribed
coordinate rule.
