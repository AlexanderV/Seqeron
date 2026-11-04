# Protein Translation

| Field | Value |
|-------|-------|
| Algorithm Group | Translation |
| Test Unit ID | TRANS-PROT-001 |
| Related Projects | N/A |
| Implementation Status | Production |
| Last Reviewed | 2026-10-14 |

## 1. Overview

Protein translation converts DNA or RNA sequences into amino-acid sequences using a selected genetic code. In this repository, the `Translator` class provides single-sequence translation, six-frame translation, and ORF discovery on one or both strands. The implementation supports frame offsets `0`, `1`, and `2`, all 27 NCBI translation tables through the `GeneticCode` abstraction, IUPAC-ambiguous codons resolved as Biopython 1.88 does, and can optionally stop translation at the first stop codon. ORF finding is a port of EMBOSS getorf (`-find 1` for the legacy overload; all seven `-find` modes and `-circular` in the getorf-modes overload documented in [Six_Frame_Translation.md](Six_Frame_Translation.md)); it is a sequence scanner, not a gene-annotation engine.[1][2][3]

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

Translation proceeds by normalizing the sequence (ASCII upper case, T/U equivalent), skipping the configured frame offset, reading complete codons in order, and translating each codon with the chosen `GeneticCode`. `TranslateSixFrames` applies this process to the forward strand and to the reverse complement. `FindOrfs` (EMBOSS getorf `-find 1` model) scans each requested strand and frame, opening an ORF at a start codon when none is open (nested in-frame starts do not open a new ORF), and emitting a result when a stop codon is reached or the strand ends with an ORF still open (getorf reports such 3'-open ORFs too).

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
| INV-04 | `OrfResult.NucleotideLength == EndPosition - StartPosition + 1`, a multiple of 3; `== 3·(AminoAcidLength + 1)` when `IsTerminated`, else `3·AminoAcidLength`. | `OrfResult` computes the property directly from stored coordinates; the end is stop-inclusive for terminated ORFs and the last complete codon otherwise (B02 F7). |
| INV-05 | `OrfResult.ToGetOrfRange(L)` equals EMBOSS getorf's `[s - e]` header (`-find 1`, `-table` = the code, `-minsize 3·minLength`), and equals `OrfFeature.ToGetOrfRange(L)` of the same ORF from the getorf-modes overload. | Both call one shared mapping (getorf_WriteORF: STOP excluded, reverse strand `s = L − s + 1`, `e = L − e + 1`). |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `[Translate] dna` | `DnaSequence` | required | DNA sequence wrapper. | `null` throws `ArgumentNullException`. |
| `[Translate] rna` | `RnaSequence` | required | RNA sequence wrapper. | `null` throws `ArgumentNullException`. |
| `[Translate] sequence` | `string` | required | DNA or RNA sequence string (IUPAC ambiguity codes allowed, case-insensitive). | `null` or empty returns an empty protein (after `frame` / `toFirstStop` validation). |
| `[Translate*] geneticCode` | `GeneticCode` | `GeneticCode.Standard` | Translation table. | Supported built-ins are documented in [Codon_Translation.md](Codon_Translation.md). |
| `[Translate*] frame` | `int` | `0` | Reading-frame offset. | Must be `0`, `1`, or `2`; otherwise `ArgumentOutOfRangeException`. |
| `[Translate*] toFirstStop` | `bool` | `false` | Whether to stop before appending the first translated stop codon. | Applied during codon iteration. Throws `ArgumentException` for tables with dual-coding stop codons (27, 28, 31), as Biopython raises `ValueError`.[7] |
| `[Translate string] allowNucleotideX` | `bool` | `false` (5-argument overload `Translate(string, GeneticCode?, int, bool, bool)`) | Opt-in Biopython-compatible nucleotide `X` (= any base; possible-stop X-codons throw). Only the string overload: `DnaSequence`/`RnaSequence` never contain `X`, so `FindOrfs` / `TranslateSixFrames` are unaffected. | `false` = the 4-argument behaviour. See [Codon_Translation.md](Codon_Translation.md) §3.3. |
| `[FindOrfs] minLength` | `int` | `100` | Minimum ORF length in amino acids. | ORFs shorter than this are filtered out. |
| `[FindOrfs] searchBothStrands` | `bool` | `true` | Whether to scan the reverse complement in addition to the forward strand. | Reverse-strand ORFs are reported with negative frame values. |
| `[FindOrfs] geneticCode` | `GeneticCode` | `GeneticCode.Standard` | Start/stop codon sets (NCBI). | Tables 27, 28, 31 (dual-coding stops) → `ArgumentException` (B02 F8). |
| `[FindOrfs string] sequence` | `string` | required | IUPAC DNA/RNA string overload (e.g. genomes with N runs). | Non-IUPAC symbol → `ArgumentException`; an ambiguous codon is never a START (getorf); STOP iff all expansions are stops (B02 F45). |
| `[OrfResult.ToGetOrfRange] sequenceLength` | `int` | required | Length L of the searched sequence. | ≤ 0 → `ArgumentOutOfRangeException`. |

### 3.2 Output / Return Value

| Name | Type | Description |
|------|------|-------------|
| `Translate` result | `ProteinSequence` | Protein produced from the selected frame and genetic code. |
| `TranslateSixFrames` result | `IReadOnlyDictionary<int, ProteinSequence>` | Six translated proteins keyed by `1`, `2`, `3`, `-1`, `-2`, `-3`. |
| `FindOrfs` result | `IEnumerable<OrfResult>` | ORFs with 0-based inclusive start/end (stop-inclusive end for terminated ORFs, INSDC CDS convention; reverse-strand ORFs in reverse-complement coordinates), signed frame, and translated protein (initiator always `M`, getorf `-methionine Y`). `IsTerminated` tells stop-terminated from 3'-open ORFs. |
| `OrfResult.ToGetOrfRange(L)` | `(int Start, int End)` | getorf's 1-based forward-strand `[s - e]` header (STOP excluded; reverse strand `s > e`). |

### 3.3 Preconditions and Validation

All translation paths treat `T` and `U` alike. The string overload folds ASCII lower case to upper case before translation (non-ASCII letters are rejected, not folded; B02 F31); the wrapper overloads delegate their `Sequence` values to the same core helper. Translation ignores incomplete trailing bases because codons are consumed only while three bases remain. `FindOrfs` requires a non-null `DnaSequence` (or an IUPAC string) and uses the selected `GeneticCode`'s start-codon set rather than hard-coding `AUG` only.

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

ORF detection starts on any unambiguous `geneticCode.IsStartCodon(codon)` result, not only on `AUG`, and reports the initiator as `M` whatever the start codon (getorf `-methionine Y`; Biopython `translate(cds=True)`; B02 F6).

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
- `Translator.Translate(string, GeneticCode?, int, bool, bool allowNucleotideX)` (opt-in Biopython nucleotide `X`)
- `Translator.TranslateSixFrames(DnaSequence | string, GeneticCode?[, SixFrameConvention])`
- `Translator.FindOrfs(DnaSequence | string, GeneticCode?, int, bool)` → `OrfResult`; `OrfResult.ToGetOrfRange(int)`
- getorf-modes overload `Translator.FindOrfs(DnaSequence | string, OrfFindMode, …)` → `OrfFeature` (see [Six_Frame_Translation.md](Six_Frame_Translation.md))

### 5.2 Current Behavior

`Translate(string, ...)` returns an empty `ProteinSequence` for `null` or empty input, while the `DnaSequence` and `RnaSequence` overloads throw on `null`. `toFirstStop` stops translation before appending the stop codon. `TranslateSixFrames` always returns all six frame keys, even for an empty sequence. `FindOrfs` starts an ORF on any unambiguous codon accepted by `geneticCode.IsStartCodon`, reports reverse-strand ORFs with negative frames (in reverse-complement coordinates; `OrfResult.ToGetOrfRange(L)` gives getorf's forward-strand header), and emits 3'-open ORFs that reach the end of the scanned strand without a stop codon if they meet `minLength`, as getorf does. Results equal getorf `-find 1` (B02 F6/F44: 0 mismatches over 491 658 reference ORFs).

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Translation reads complete codons in the requested frame and maps them through the selected genetic code.[1][2] Output is identical to Biopython 1.88 `translate(seq[frame:], table, to_stop)` for all 27 NCBI tables, DNA/RNA, lower case and IUPAC-ambiguous input (review 2026-09: 5,200+ randomized cases).[7]
- Six-frame translation covers the three forward and three reverse-complement frames.[2]
- ORF detection is the EMBOSS getorf `-find 1` model: start codons of the selected table, stop codons, initiator read as Met, 3'-open ORFs reported, a minimum length threshold; reverse-strand coordinates map to getorf's forward-strand header with `OrfResult.ToGetOrfRange` (B02 F6–F8, F44, F45, this review F56).[3]
- All 27 NCBI tables, IUPAC ambiguity (B/Z/J), dual-coding-stop tables rejected for `toFirstStop`/ORFs, ASCII-only case folding, opt-in nucleotide `X` (B02 F1–F5, F31, F33); `ProteinSequence` accepts U/O and reports EMBOSS pepstats residue classes (F32, F47); EMBOSS transeq six-frame conventions (F53).

**Intentionally simplified:**

- ORF scanning has no gene-model features (scoring, splicing, regulatory context); **consequence:** results are candidate ORFs rather than full gene annotations.

**Not implemented:**

- Nested ORF reporting (a new ORF opened at an in-frame START before the current one closes) — getorf `-find 1` does not report them either; **users should rely on:** the getorf-modes overload with `OrfFindMode.TranslationStopToStop` (STOP→STOP regions) and caller-side START search.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | `minLength` in amino acids (initiator included), not nucleotides | Convention | getorf `-minsize` is in nucleotides (÷3 for `-find 1`): `minLength = k` ≡ `-minsize 3k` | accepted | [Six_Frame_Translation.md](Six_Frame_Translation.md) §5.4 #2 |
| 2 | End coordinate includes the STOP for terminated ORFs (INSDC CDS convention); reverse-strand ORFs in reverse-complement coordinates | Convention | getorf prints the range without the STOP in forward-strand coordinates; `OrfResult.ToGetOrfRange(L)` converts | accepted | B02 F8, F56; Six_Frame §5.4 #3–#4 |
| 3 | A START in the last complete codon of a frame is reported as a 1-aa open ORF when `minLength ≤ 1` | Divergence | getorf never writes it (branch-order artefact) | accepted | unobservable for `minLength ≥ 2`; Six_Frame §5.4 #5 |
| 4 | Ambiguous codons (string overload) resolved by the all-expansions rule | Divergence | EMBOSS's `trnComplete` table differs for a few codons (e.g. table 1 `TRR`, `TRA`) | accepted | Six_Frame §5.4 #7 |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Empty string input | Returns an empty protein. | The string overload special-cases empty input. |
| Sequence shorter than 3 nt | Returns an empty protein. | No complete codon is available. |
| Invalid frame | Throws `ArgumentOutOfRangeException`. | Frame must be `0`, `1`, or `2`. |
| `toFirstStop = true` | Translation stops before appending the stop codon. | The loop breaks on `*`. |
| `toFirstStop = true` with table 27, 28 or 31 | `ArgumentException`. | These tables have stop codons that are also translated as an amino acid (e.g. table 27 UGA = W or STOP), so "first stop" is undefined; Biopython 1.88 `_translate_str` raises `ValueError`.[7] |
| Ambiguous codon resolving to Asx/Glx/Xle (RAY, SAR, MTH in table 1) | Protein contains `B`, `Z`, `J`. | Biopython 1.88 translates these as B/Z/J; `ProteinSequence` accepts the IUPAC ambiguity codes B, Z, J.[7] |
| Selenocysteine `U` / pyrrolysine `O` in a protein | Accepted by `ProteinSequence` (three-letter Sec/Pyl). Never produced by translation (UGA/UAG → `*`; recoding is not modelled). | Biopython 1.88 `IUPACData.extended_protein_letters` ("ACDEFGHIKLMNPQRSTVWYBXZJUO") / `protein_letters_1to3_extended`; `ProteinSequence.ValidCharacters` = these 26 letters + `*`.[7] |
| `ProteinSequence.TypePercentage` / `TypeCount` residue classes | Nine EMBOSS `pepstats` classes (Tiny, Small, Aliphatic, Aromatic, Non-polar, Polar, Charged, Basic, Acidic) as `AminoAcidType` members; membership from `Eamino.dat` (U as C, O as K, J as I/L; X, `*` in no class); B/Z split to D/N, E/Q (`embIepCompC`); denominator = pepstats "Residues" (one terminal `*` dropped); unrounded Mole%. Insulin A: Tiny 38.095, Non-polar 57.143, Polar 42.857; LACI_ECOLI: Small 54.722, Charged 19.167, Basic 10.278, Acidic 8.889. Supersedes the earlier unsourced 4-class polarity grouping rounded to 2 dp (B02 F47). | EMBOSS 6.6.0 `pepstats`.[8] |
| Non-IUPAC nucleotide symbol (e.g. `-`, `X`) | `ArgumentException`. | Biopython raises `TranslationError` for `-`. |
| Nucleotide `X` with `Translate(string, code, frame, toFirstStop, allowNucleotideX: true)` | `X` = any base: `"ATGGCX"`→`"MA"`; a possible-stop X-codon (`TAX`, `XXX`) throws `ArgumentException` when reached (not after an earlier stop with `toFirstStop`). | Biopython 1.88 `Seq.translate` (`ambiguous_dna_values['X'] == 'GATC'`; `TranslationError` for possible-stop X-codons); opt-in, default rejects `X` (review 2026-09 B02 F33). |
| No start codon in `FindOrfs` | Returns no ORFs. | ORF accumulation begins only at a start codon. |

### 6.2 Limitations

The repository translation utilities are sequence scanners. They do not include gene-model scoring, splice awareness, or regulatory-context handling (reverse-strand coordinates are remapped to getorf's forward-strand header with `OrfResult.ToGetOrfRange` / `OrfFeature.ToGetOrfRange`). ORF output should therefore be treated as a practical translation aid rather than a full annotation result.

## 7. Examples and Related Material

### 7.3 Related Tests, Evidence, or Documents

- Tests: [TranslatorTests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Core/TranslatorTests.cs) — covers `INV-01`, `INV-02`, `INV-03`, `INV-04`; [Translator_OrfResultGetOrfRange_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Core/Translator_OrfResultGetOrfRange_Tests.cs) — covers `INV-05` (getorf 6.6.0 binary output)
- Test specification: [TRANS-PROT-001.md](../../../tests/TestSpecs/TRANS-PROT-001.md)
- Related algorithms: [Codon_Translation.md](Codon_Translation.md), [Six_Frame_Translation.md](Six_Frame_Translation.md)

## 8. References

1. Wikipedia contributors. 2026. Translation (biology). Wikipedia. N/A
2. Wikipedia contributors. 2026. Reading frame. Wikipedia. N/A
3. Wikipedia contributors. 2026. Open reading frame. Wikipedia. N/A
4. NCBI. 2026. The Genetic Codes. https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
5. Lodish H et al. 2007. Molecular Cell Biology. N/A
6. Test specification: [TRANS-PROT-001.md](../../../tests/TestSpecs/TRANS-PROT-001.md)
7. Biopython 1.88 `Bio/Seq.py` `_translate_str` and `Bio/Data/IUPACData.py` (`extended_protein_letters`, `protein_letters_1to3_extended`). https://raw.githubusercontent.com/biopython/biopython/master/Bio/Seq.py
8. EMBOSS 6.6.0 `emboss/pepstats.c` (`pepstats_issize`), `nucleus/embiep.c` (`embIepCompC`, B/Z split) and `data/Eamino.dat`. https://raw.githubusercontent.com/kimrutherford/EMBOSS/master/emboss/pepstats.c
9. EMBOSS 6.6.0 `emboss/getorf.c` (`getorf_FindORFs`, `getorf_WriteORF`) and the `getorf` binary. https://raw.githubusercontent.com/kimrutherford/EMBOSS/master/emboss/getorf.c
