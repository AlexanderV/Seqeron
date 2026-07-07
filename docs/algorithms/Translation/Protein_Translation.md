# Protein Translation

| Field | Value |
|-------|-------|
| Algorithm Group | Translation |
| Test Unit ID | TRANS-PROT-001 |
| Related Projects | N/A |
| Implementation Status | Simplified |
| Last Reviewed | 2026-04-30 |

## 1. Overview

Protein translation converts DNA or RNA sequences into amino-acid sequences using a selected genetic code. In this repository, the `Translator` class provides single-sequence translation, six-frame translation, and ORF discovery on one or both strands. The implementation supports frame offsets `0`, `1`, and `2`, uses the `GeneticCode` abstraction for codon mapping, and can optionally stop translation at the first stop codon. ORF finding is heuristic and sequence-based rather than a full annotation engine, but it is sufficient for start/stop scanning within the provided DNA sequence.[1][2][3]

## 2. Scientific / Formal Basis

### 2.1 Domain Context

Translation reads a nucleotide sequence in triplets, producing a protein sequence from the 5' to 3' reading direction.[1] The original document recorded the standard three forward reading frames and the six-frame convention for double-stranded DNA:[2]

| Frame Convention | Meaning |
|------------------|---------|
| `0`, `1`, `2` | Public `Translate` frame offsets in the repository. |
| `+1`, `+2`, `+3` | Forward-strand keys returned by `TranslateSixFrames`. |
| `-1`, `-2`, `-3` | Reverse-complement keys returned by `TranslateSixFrames`. |

An ORF is used here in the practical repository sense of a start-codon-initiated translated region that continues until a stop codon or the end of the scanned sequence, subject to a minimum amino-acid length filter.[3]

### 2.2 Core Model

Translation proceeds by normalizing DNA to RNA, skipping the configured frame offset, reading complete codons in order, and translating each codon with the chosen `GeneticCode`. `TranslateSixFrames` applies this process to the forward strand and to the reverse complement. `FindOrfs` scans each requested strand and frame, starting a new ORF when it encounters a start codon and emitting a result when a stop codon is reached or the sequence ends with an ORF still open.

### 2.3 Modeling Assumptions

| ID | Assumption | Consequence if Violated |
|----|------------|--------------------------|
| ASM-01 | The provided sequence is meaningful in the requested reading frame(s). | The returned protein can be mathematically correct for the frame yet biologically irrelevant. |
| ASM-02 | Start and stop codons from the chosen `GeneticCode` are appropriate for the organism and compartment. | ORF calling can miss valid coding regions or report false ORFs. |

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | `Translate(seq, frame: 0)` matches `TranslateSixFrames(seq)[1]`. | `TranslateSixFrames` calls the same translation helper with frame `0` on the forward strand. |
| INV-02 | Translated protein length is at most `floor((sequenceLength - frame) / 3)`. | Translation consumes only complete codons from the selected frame. |
| INV-03 | `TranslateSixFrames` returns exactly six frame keys: `1`, `2`, `3`, `-1`, `-2`, `-3`. | The method populates each key explicitly. |
| INV-04 | `OrfResult.NucleotideLength == EndPosition - StartPosition + 1`. | `OrfResult` computes the property directly from stored coordinates. |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `[Translate] dna` | `DnaSequence` | required | DNA sequence wrapper. | `null` throws `ArgumentNullException`. |
| `[Translate] rna` | `RnaSequence` | required | RNA sequence wrapper. | `null` throws `ArgumentNullException`. |
| `[Translate] sequence` | `string` | required | DNA or RNA sequence string. | `null` or empty returns an empty protein. |
| `[Translate*] geneticCode` | `GeneticCode` | `GeneticCode.Standard` | Translation table. | Supported built-ins are documented in [Codon_Translation.md](Codon_Translation.md). |
| `[Translate*] frame` | `int` | `0` | Reading-frame offset. | Must be `0`, `1`, or `2`; otherwise `ArgumentOutOfRangeException`. |
| `[Translate*] toFirstStop` | `bool` | `false` | Whether to stop before appending the first translated stop codon. | Applied during codon iteration. |
| `[FindOrfs] minLength` | `int` | `100` | Minimum ORF length in amino acids. | ORFs shorter than this are filtered out. |
| `[FindOrfs] searchBothStrands` | `bool` | `true` | Whether to scan the reverse complement in addition to the forward strand. | Reverse-strand ORFs are reported with negative frame values. |

### 3.2 Output / Return Value

| Name | Type | Description |
|------|------|-------------|
| `Translate` result | `ProteinSequence` | Protein produced from the selected frame and genetic code. |
| `TranslateSixFrames` result | `IReadOnlyDictionary<int, ProteinSequence>` | Six translated proteins keyed by `1`, `2`, `3`, `-1`, `-2`, `-3`. |
| `FindOrfs` result | `IEnumerable<OrfResult>` | ORFs with start position, end position, frame, and translated protein. |

### 3.3 Preconditions and Validation

All translation paths normalize DNA to RNA by replacing `T` with `U`. The string overload uppercases the input before translation; the wrapper overloads delegate their `Sequence` values to the same core helper. Translation ignores incomplete trailing bases because codons are consumed only while three bases remain. `FindOrfs` requires a non-null `DnaSequence` and uses the selected `GeneticCode`'s start-codon set rather than hard-coding `AUG` only.

## 4. Algorithm

### 4.1 High-Level Steps

1. Validate the input object or return an empty protein for an empty string sequence.
2. Normalize the sequence to RNA notation.
3. For `Translate`, begin at the requested frame offset and iterate over complete codons.
4. Translate each codon using the chosen `GeneticCode`.
5. If `toFirstStop` is true and the translated amino acid is `*`, stop without appending the stop symbol.
6. For `TranslateSixFrames`, translate forward frames `0`, `1`, `2` and reverse-complement frames `0`, `1`, `2`, then map them to keys `±1..±3`.
7. For `FindOrfs`, scan each requested strand and frame for start codons, accumulate amino acids until stop or sequence end, and emit ORFs meeting `minLength`.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

The repository uses two frame-numbering conventions simultaneously:

| API surface | Convention |
|-------------|------------|
| `Translate(..., frame: ...)` | Zero-based offsets `0`, `1`, `2` |
| `TranslateSixFrames(...)` and `OrfResult.Frame` | Signed frame labels `+1`, `+2`, `+3`, `-1`, `-2`, `-3` |

ORF detection starts on any `geneticCode.IsStartCodon(codon)` result, not only on `AUG`.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `Translate` | `O(n)` | `O(n)` | Linear in input length, producing one amino acid per complete codon. |
| `TranslateSixFrames` | `O(n)` | `O(n)` | Six linear translations over the forward and reverse-complement strands. |
| `FindOrfs` | `O(n)` per strand | `O(k)` | Linear scan of each requested strand, where `k` is the number of emitted ORFs. |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [Translator.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Core/Translator.cs)

- `Translator.Translate(DnaSequence, GeneticCode?, int, bool)`
- `Translator.Translate(RnaSequence, GeneticCode?, int, bool)`
- `Translator.Translate(string, GeneticCode?, int, bool)`
- `Translator.TranslateSixFrames(DnaSequence, GeneticCode?)`
- `Translator.FindOrfs(DnaSequence, GeneticCode?, int, bool)`

### 5.2 Current Behavior

`Translate(string, ...)` returns an empty `ProteinSequence` for `null` or empty input, while the `DnaSequence` and `RnaSequence` overloads throw on `null`. `toFirstStop` stops translation before appending the stop codon. `TranslateSixFrames` always returns all six frame keys, even for an empty sequence. `FindOrfs` starts an ORF on any codon accepted by `geneticCode.IsStartCodon`, reports reverse-strand ORFs with negative frames, and emits terminal ORFs that reach the end of the scanned sequence without encountering a stop codon if they meet `minLength`.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Translation reads complete codons in the requested frame and maps them through the selected genetic code.[1][2]
- Six-frame translation covers the three forward and three reverse-complement frames.[2]
- ORF detection uses start codons, stop codons, and a minimum length threshold.[3]

**Intentionally simplified:**

- ORF scanning is a direct start-to-stop sequence scan without higher-level gene-model features; **consequence:** results are candidate ORFs rather than full gene annotations.
- Reverse-strand ORF coordinates are reported in the coordinate system of the reverse-complement sequence that is scanned; **consequence:** callers needing original-strand genomic coordinates must remap them externally.
- When an ORF runs off the end of the sequence, it is still emitted if long enough; **consequence:** ORF results are not limited to stop-terminated regions only.

**Not implemented:**

- Genomic-coordinate remapping for reverse-strand ORFs; **users should rely on:** caller-side remapping.
- Nested ORF reporting that opens a new ORF before the current one closes in the same frame; **users should rely on:** no current alternative.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | ORFs that reach the end of the sequence without a stop codon can still be emitted. | Deviation | Results are broader than a strict stop-terminated ORF definition. | accepted | Confirmed in [Translator.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Core/Translator.cs). |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Empty string input | Returns an empty protein. | The string overload special-cases empty input. |
| Sequence shorter than 3 nt | Returns an empty protein. | No complete codon is available. |
| Invalid frame | Throws `ArgumentOutOfRangeException`. | Frame must be `0`, `1`, or `2`. |
| `toFirstStop = true` | Translation stops before appending the stop codon. | The loop breaks on `*`. |
| No start codon in `FindOrfs` | Returns no ORFs. | ORF accumulation begins only at a start codon. |

### 6.2 Limitations

The repository translation utilities are sequence scanners. They do not include gene-model scoring, splice awareness, coordinate remapping for reverse-strand ORFs, or regulatory-context handling. ORF output should therefore be treated as a practical translation aid rather than a full annotation result.

## 7. Examples and Related Material

### 7.3 Related Tests, Evidence, or Documents

- Tests: [TranslatorTests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Core/TranslatorTests.cs) — covers `INV-01`, `INV-02`, `INV-03`, `INV-04`
- Test specification: [TRANS-PROT-001.md](../../../tests/TestSpecs/TRANS-PROT-001.md)
- Related algorithms: [Codon_Translation.md](Codon_Translation.md)

## 8. References

1. Wikipedia contributors. 2026. Translation (biology). Wikipedia. N/A
2. Wikipedia contributors. 2026. Reading frame. Wikipedia. N/A
3. Wikipedia contributors. 2026. Open reading frame. Wikipedia. N/A
4. NCBI. 2026. The Genetic Codes. https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
5. Lodish H et al. 2007. Molecular Cell Biology. N/A
6. Test specification: [TRANS-PROT-001.md](../../../tests/TestSpecs/TRANS-PROT-001.md)
