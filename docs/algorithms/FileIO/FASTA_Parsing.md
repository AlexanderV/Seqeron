# FASTA Parsing

| Field | Value |
|-------|-------|
| Algorithm Group | File I/O |
| Test Unit ID | PARSE-FASTA-001 |
| Related Projects | Seqeron.Genomics; Seqeron.Mcp.Parsers |
| Implementation Status | Simplified |
| Last Reviewed | 2026-06-24 |

## 1. Overview

FASTA is a text-based format for representing biological sequences and originated from the FASTA software package (Lipman and Pearson, 1985). The format is specification-driven: each record begins with a defline starting with `>`, followed by one or more sequence lines. In this repository, `FastaParser` implements sequential parsing and formatting for multi-entry FASTA content, with linear-time scans over the input or output length. The implementation materializes entries as `DnaSequence` objects, so the repository documents a DNA-oriented subset of the broader FASTA convention.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

FASTA is a de facto interchange format for nucleotide and protein sequences in bioinformatics workflows. The governing format rules used in this document come from the FASTA format description summarized by Wikipedia and the NCBI BLAST FASTA guidance: a record begins with a single-line defline whose first token is the identifier, and subsequent lines contain sequence characters that may be split across multiple lines.

The canonical single-record layout is:

```text
>ID description
SEQUENCE_LINE_1
SEQUENCE_LINE_2
```

Multiple FASTA entries are formed by concatenating records, each starting with `>`.

### 2.2 Core Model

The parser is a line-oriented state machine with two states: current header and current sequence buffer. For each input line:

- If the line starts with `>`, the previous buffered record is emitted when both a header and at least one sequence character have been collected, then the new header is stored.
- Otherwise, non-whitespace characters from the line are appended to the current sequence buffer.
- At end of input, the final buffered record is emitted when it has both header and sequence content.

Header parsing follows the common FASTA convention already documented in the repository evidence and tests: the first space- or tab-delimited token is the sequence identifier and the remainder of the defline is the optional description.

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | Every emitted entry comes from a defline that starts with `>` and yields an identifier from the first whitespace-delimited token | FASTA defline rule documented in the existing document, evidence, and `CreateEntry` |
| INV-02 | Sequence content for an emitted entry is the concatenation of all non-whitespace characters collected after its defline and before the next defline or end of file | `ParseReader` and `ParseFileAsync` append only non-whitespace characters to a `StringBuilder` |
| INV-03 | `ToFasta` writes each output record as `>` + header followed by sequence chunks no longer than `lineWidth` | `ToFasta` prepends `>` and iterates in steps of `lineWidth` |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `fastaContent` | `string` | required | FASTA text to parse with `Parse` | Null, empty, or whitespace-only content yields no entries |
| `filePath` | `string` | required | Path to a FASTA file for `ParseFile` or `ParseFileAsync` | Must identify a readable file when parsing or writable path when writing |
| `entries` | `IEnumerable<FastaEntry>` | required | Entries to serialize with `ToFasta` or `WriteFile` | Each entry must already contain a valid `DnaSequence` |
| `lineWidth` | `int` | `80` | Maximum sequence characters written per output line | Must be positive for meaningful wrapping; formatter loops in increments of this value |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `Parse` / `ParseFile` | `IEnumerable<FastaEntry>` | Lazily yields parsed FASTA entries |
| `ParseFileAsync` | `IAsyncEnumerable<FastaEntry>` | Asynchronously yields parsed FASTA entries |
| `FastaEntry.Id` | `string` | First token of the defline |
| `FastaEntry.Description` | `string?` | Remaining defline text after the identifier, or `null` when absent |
| `FastaEntry.Sequence` | `DnaSequence` | Parsed DNA sequence, normalized and validated by `DnaSequence` |
| `FastaEntry.Header` | `string` | Reconstructed header as `Id` or `Id + " " + Description` |
| `ToFasta` | `string` | FASTA-formatted text for the supplied entries |
| `WriteFile` | `void` | Writes FASTA-formatted text to disk |

### 3.3 Preconditions and Validation

`Parse` returns an empty sequence for null, empty, or whitespace-only input. Both synchronous and asynchronous parsing require a defline before any emitted sequence data; a header without sequence content is not yielded. Sequence lines are treated as raw character streams with whitespace removed before entry creation. Entry materialization uses `new DnaSequence(sequence)`, so parsed content is normalized to uppercase and rejected with `ArgumentException` if any non-DNA character remains after whitespace removal. Header parsing is case-preserving and does not impose additional validation on identifier syntax beyond splitting on the first space or tab.

## 4. Algorithm

### 4.1 High-Level Steps

1. Read FASTA content line by line.
2. When a line starts with `>`, emit the previous buffered entry if both header and sequence are present, then store the new header and clear the sequence buffer.
3. For non-header lines, append every non-whitespace character to the current sequence buffer.
4. After the scan completes, emit the final buffered entry when it has both a header and at least one sequence character.
5. For formatting, write `>` plus the entry header, then emit the sequence in chunks of `lineWidth` characters.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

The implementation uses a `StringBuilder` as the sequence buffer for both synchronous and asynchronous parsing. Header parsing splits on the first space or tab into at most two parts, preserving the first token as `Id` and the remainder as `Description`. Formatting uses `Substring` over the normalized DNA string in fixed-width chunks.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `Parse` / `ParseFile` / `ParseFileAsync` | `O(n)` | `O(m)` | `n` = total input characters, `m` = longest buffered sequence |
| `ToFasta` | `O(n)` | `O(n)` | `n` = total emitted sequence characters plus headers in the returned string |
| `WriteFile` | `O(n)` | `O(n)` | Delegates to `ToFasta` before `File.WriteAllText` |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [FastaParser.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.IO/FastaParser.cs), [DnaSequence.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Core/DnaSequence.cs)

- `FastaParser.Parse(string)`: Parses FASTA text with a line-oriented iterator.
- `FastaParser.ParseFile(string)`: Opens a file and delegates to the shared parser logic.
- `FastaParser.ParseFileAsync(string)`: Asynchronously reads a FASTA file and yields entries.
- `FastaParser.ToFasta(IEnumerable<FastaEntry>, int)`: Formats entries as FASTA with configurable line wrapping.
- `FastaParser.WriteFile(string, IEnumerable<FastaEntry>, int)`: Writes serialized FASTA text to disk.
- `FastaEntry.Header`: Reconstructs the header line from `Id` and `Description`.
- `DnaSequence.DnaSequence(string)`: Normalizes parsed sequence data to uppercase and rejects characters outside `A/C/G/T`.
- `FastaParser.Parse(string, SequenceAlphabet)` / `ParseFile(string, SequenceAlphabet)` / `ParseFileAsync(string, SequenceAlphabet)`: **Opt-in** overloads that validate the uppercased sequence against a selectable alphabet and yield `FastaRecord` (raw sequence string preserved). Alphabets: `StrictDna` (A/C/G/T — same as default), `IupacNucleotide` (A C G T U R Y S W K M B D H V N + `-`, per NC-IUB 1985), `Rna` (A C G U), `Protein` (20 IUPAC residues + B Z J X U O + `*`).

### 5.2 Current Behavior

Repository-specific behavior confirmed by source and tests:

- Whitespace inside sequence lines is discarded before sequence validation.
- Blank lines are effectively skipped because they contribute no characters to the sequence buffer.
- Header-only entries are dropped because entries are yielded only when `sequenceBuilder.Length > 0`.
- Lowercase nucleotide input is accepted because `DnaSequence` uppercases before validation.
- `ToFasta` uses a default output width of 80 characters per line.
- Parsed entries preserve the identifier and description split used by the defline parser, and round-trip tests verify `Parse -> ToFasta -> Parse` consistency for that split.
- The default `Parse`/`ParseFile`/`ParseFileAsync` (no alphabet argument) remain strict DNA-only and byte-for-byte unchanged. The opt-in `SequenceAlphabet` overloads additionally accept RNA (`U`), protein residues (incl. `*` and ambiguity codes), and IUPAC nucleotide ambiguity/gap codes, validating against the selected alphabet and throwing `ArgumentException` on the first out-of-alphabet character.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Defline parsing based on `>`-prefixed records with identifier as the first token and optional trailing description.
- Multi-entry FASTA parsing by concatenating multiple `>`-prefixed records.
- Multi-line sequence concatenation and configurable output line wrapping.

**Intentionally simplified:**

- The default (no-alphabet) `Parse` path materializes payload as `DnaSequence`; **consequence:** protein FASTA, RNA `U`, ambiguity codes such as `N`, and gap-containing sequence lines are rejected on the default path (preserved intentionally). The opt-in `SequenceAlphabet` overloads remove this restriction for the selected alphabet.
- Header-only records are not yielded; **consequence:** malformed FASTA entries without sequence content disappear from parsed output instead of being preserved as empty records.
- Whitespace is removed from sequence lines before materialization; **consequence:** formatting normalizes whitespace-bearing input to contiguous sequence strings on round trip.

**Not implemented:**

- `ToFasta`/`WriteFile` serialization for the non-DNA `FastaRecord` type; **users should rely on:** the `FastaEntry`/`DnaSequence` formatter for serialization, or assemble output strings directly from `FastaRecord.Sequence`.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | DNA-only on the default `Parse` path | Deviation | The default (no-alphabet) path rejects non-`A/C/G/T` payloads | accepted | Caused by `FastaEntry.Sequence` using `DnaSequence`; the opt-in `SequenceAlphabet` overloads (`Rna`, `Protein`, `IupacNucleotide`) lift this for the selected alphabet |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Null, empty, or whitespace-only input to `Parse` | Returns no entries | `Parse` checks `string.IsNullOrWhiteSpace` and `yield break`s |
| Header without description | `Description` is `null` | `CreateEntry` stores only the first token when no remainder exists |
| Header without sequence | Record is not yielded | Entries are emitted only when buffered sequence length is greater than zero |
| Multi-line sequence | Lines are concatenated | Parser appends sequence characters across successive non-header lines |
| Whitespace inside sequence lines | Whitespace is ignored | Parser appends only characters for which `!char.IsWhiteSpace(c)` |
| Lowercase DNA sequence | Output `DnaSequence` is uppercase | `DnaSequence` normalizes with `ToUpperInvariant()` |
| Long output sequence | Wrapped at `lineWidth` characters per line | `ToFasta` writes fixed-width chunks |

### 6.2 Limitations

The default `Parse`/`ParseFile`/`ParseFileAsync` path documents a DNA-specific subset of FASTA (parsed records instantiated as `DnaSequence`). Opt-in `SequenceAlphabet` overloads accept RNA, protein, and IUPAC-ambiguous nucleotide alphabets (returning `FastaRecord`). In all modes the parser does not preserve empty (header-only) records and does not retain whitespace formatting inside sequence lines. On multi-space deflines the `FastaEntry`/`FastaRecord` `Description` keeps a single leading space (by-design header-split contract; not changed here).

## 7. Examples and Related Material

### 7.1 Worked Example

**Numerical / biological walk-through:**

Input:

```text
>seq1 First sequence
ATGCATGC
>seq2 Second sequence
GGCCTTAA
```

Parsed result:

- Entry 1: `Id = "seq1"`, `Description = "First sequence"`, `Sequence = "ATGCATGC"`
- Entry 2: `Id = "seq2"`, `Description = "Second sequence"`, `Sequence = "GGCCTTAA"`

### 7.3 Related Tests, Evidence, or Documents

- Tests: [FastaParserTests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/IO/FastaParserTests.cs) — covers `INV-01`, `INV-02`, and `INV-03`
- Tests: [FastaParser_Alphabet_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/IO/FastaParser_Alphabet_Tests.cs) — covers the opt-in `SequenceAlphabet` (RNA / protein / IUPAC-nucleotide) overloads and the default strict-DNA regression
- Tests: [FastaParseTests.cs](../../../tests/Seqeron/Seqeron.Mcp.Parsers.Tests/FastaParseTests.cs) — covers parser binding behavior
- Tests: [FastaWriteTests.cs](../../../tests/Seqeron/Seqeron.Mcp.Parsers.Tests/FastaWriteTests.cs) — covers output file generation
- Tests: [FastaFormatTests.cs](../../../tests/Seqeron/Seqeron.Mcp.Parsers.Tests/FastaFormatTests.cs) — covers line wrapping and formatter binding behavior
- Evidence: [PARSE-FASTA-001-Evidence.md](../../../docs/Evidence/PARSE-FASTA-001-Evidence.md)
- Related algorithms: [FASTQ_Parsing.md](FASTQ_Parsing.md)

## 8. References

1. Lipman DJ, Pearson WR. 1985. Rapid and sensitive protein similarity searches. Science 227(4693):1435-1441.
2. Pearson WR, Lipman DJ. 1988. Improved tools for biological sequence comparison. Proceedings of the National Academy of Sciences 85(8):2444-2448.
3. Wikipedia contributors. 2026. FASTA format. Wikipedia. https://en.wikipedia.org/wiki/FASTA_format
4. NCBI. 2026. BLAST topics and FASTA input guidance. NCBI BLAST Help. https://blast.ncbi.nlm.nih.gov/doc/blast-topics/
5. NCBI. 2026. FASTA input format. NCBI BLAST. https://www.ncbi.nlm.nih.gov/blast/fasta.shtml
