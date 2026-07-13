---
type: concept
title: "FASTA parsing (>defline + sequence lines; default strict-DNA, opt-in alphabets)"
tags: [file-io, algorithm]
sources:
  - docs/algorithms/FileIO/FASTA_Parsing.md
  - docs/Evidence/PARSE-FASTA-001-Evidence.md
source_commit: a84ee65448c2719984cef2d0213ab95bfbdc428a
created: 2026-07-13
updated: 2026-07-13
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: parse-fasta-001-evidence
      evidence: "Test Unit ID: PARSE-FASTA-001 ... Algorithm Group: File I/O (FastaParser)"
      confidence: high
      status: current
---

# FASTA parsing (>defline + sequence lines; default strict-DNA, opt-in alphabets)

**FASTA** is the simplest and most ubiquitous sequence format: a `>` **description line
(defline)** followed by one or more **sequence lines**, repeated for multiple records. Seqeron's
`FastaParser` (test unit **PARSE-FASTA-001**, status *Simplified*) reads such content into typed
`FastaEntry` values (`Id`, `Description`, sequence) and writes them back with configurable line
wrapping. This concept is a member of the file-parsing **FileIO** (`PARSE-*`) family anchored by
[[bed-format-parsing]]; the literature-traced format facts and character-set tables live in the
Evidence source page [[parse-fasta-001-evidence]]. [[test-unit-registry]] tracks the unit,
[[algorithm-validation-evidence]] describes the evidence-artifact pattern, and [[fuzzing]] explains
why parsers are the campaign's highest-priority malformed-input target. Its FASTQ sibling (4-line
records with a per-base quality string) is [[parse-fastq-001-evidence]].

## The format (what the parser consumes)

- **Defline:** begins with `>`; the **first whitespace-delimited token** (up to the first space
  or tab) is the sequence **identifier** (`Id`); the remainder is the optional **description**.
- **Sequence lines:** follow the defline; may be **multi-line (interleaved)** or **single-line**;
  **whitespace inside a sequence line is ignored** (only non-whitespace characters are collected).
- **Multi-FASTA:** several records concatenated in one input, each starting with a fresh `>`.

Canonical single record:

```text
>seq1 First sequence
ATGCATGC
```

## The parser as a state machine

`FastaParser.Parse` is a **line-oriented state machine** with two pieces of state: the *current
header* and a *current sequence buffer* (`StringBuilder`). For each input line:

1. **Line starts with `>`** — emit the previously buffered record **iff** it has both a header
   *and* at least one sequence character, then store the new header and clear the buffer.
2. **Otherwise** — append every non-whitespace character of the line to the sequence buffer.
3. **End of input** — emit the final buffered record under the same both-header-and-sequence rule.

### Invariants (from the spec)

| ID | Invariant |
|----|-----------|
| INV-01 | Every emitted entry comes from a `>`-defline; its `Id` is the first whitespace-delimited token. |
| INV-02 | An entry's sequence is the concatenation of all non-whitespace characters between its defline and the next defline (or EOF). |
| INV-03 | `ToFasta` writes each record as `>` + header followed by sequence chunks no longer than `lineWidth`. |

## Contract & surface (`FastaParser`, `FastaEntry`, `DnaSequence`)

Implementation: `src/Seqeron/Algorithms/Seqeron.Genomics.IO/FastaParser.cs` and
`.../Seqeron.Genomics.Core/DnaSequence.cs`.

| Entry point | Returns | Notes |
|-------------|---------|-------|
| `Parse(string)` | `IEnumerable<FastaEntry>` | Lazy line-oriented iterator over in-memory text. |
| `ParseFile(string)` | `IEnumerable<FastaEntry>` | Opens the file, delegates to the shared parser. |
| `ParseFileAsync(string)` | `IAsyncEnumerable<FastaEntry>` | Async streaming read of a file. |
| `ToFasta(entries, lineWidth=80)` | `string` | Serializes; wraps sequence at `lineWidth` (default **80**). |
| `WriteFile(path, entries, lineWidth=80)` | `void` | `ToFasta` then `File.WriteAllText`. |

`FastaEntry` exposes `Id`, `Description` (`null` when the defline has no remainder),
`Sequence` (a `DnaSequence`), and a reconstructed `Header` (`Id` or `Id + " " + Description`). On the
default path, entry materialization is `new DnaSequence(sequence)`, which **uppercases** the payload
and **rejects** any non-`A/C/G/T` character with `ArgumentException` after whitespace removal — so the
default parser documents a **DNA-oriented subset** of the broader FASTA convention. Header parsing is
**case-preserving** and imposes no identifier-syntax validation beyond the first space/tab split.

### Complexity

| Operation | Time | Space |
|-----------|------|-------|
| `Parse` / `ParseFile` / `ParseFileAsync` | `O(n)` (n = input chars) | `O(m)` (m = longest buffered sequence) |
| `ToFasta` / `WriteFile` | `O(n)` (n = emitted chars) | `O(n)` |

## Opt-in `SequenceAlphabet` overloads (`FastaRecord`)

The default `Parse`/`ParseFile`/`ParseFileAsync` (no alphabet argument) remain **strict-DNA and
byte-for-byte unchanged**. Overloads taking a `SequenceAlphabet` validate the uppercased sequence
against a selectable alphabet and yield **`FastaRecord`** (the **raw sequence string is preserved**,
not wrapped in `DnaSequence`), throwing `ArgumentException` on the first out-of-alphabet character:

- **`StrictDna`** — `A C G T` (same as the default).
- **`IupacNucleotide`** — `A C G T U R Y S W K M B D H V N` + gap `-` (per NC-IUB 1985).
- **`Rna`** — `A C G U`.
- **`Protein`** — 20 standard residues + `B Z J X U O` + stop `*`.

`ToFasta`/`WriteFile` serialization is **not implemented for the non-DNA `FastaRecord`**: users
serialize via the `FastaEntry`/`DnaSequence` formatter or assemble output strings directly from
`FastaRecord.Sequence`.

## Edge cases & intentional simplifications

- **Null / empty / whitespace-only input** → no entries (`Parse` checks `IsNullOrWhiteSpace` and
  `yield break`s).
- **Header without sequence** → the record is **not yielded** (emitted only when buffered length
  `> 0`); malformed header-only entries silently disappear rather than being preserved as empty.
- **Blank lines** → effectively skipped (contribute no characters).
- **Whitespace inside sequence lines** → discarded before materialization; round-trip therefore
  normalizes whitespace-bearing input to a contiguous string.
- **Lowercase DNA** → accepted; output `DnaSequence` is uppercase.
- **Multi-space deflines** → `Description` keeps a single leading space (by-design header-split
  contract).
- **Round-trip** — `Parse → ToFasta → Parse` consistency for the Id/Description split is
  verified by tests (`FastaParserTests`, `FastaParser_Alphabet_Tests`, and the
  `Seqeron.Mcp.Parsers` parse/write/format tests).
