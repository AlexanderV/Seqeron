# Codon Translation

| Field | Value |
|-------|-------|
| Algorithm Group | Translation |
| Test Unit ID | TRANS-CODON-001 |
| Related Projects | N/A |
| Implementation Status | Complete |
| Last Reviewed | 2026-09-28 |

## 1. Overview

Codon translation maps a nucleotide triplet to a single-letter amino-acid code according to a selected genetic code table.[1][4] In this repository, the `GeneticCode` class provides codon translation, start/stop-codon classification, reverse lookup from amino acid to codons, and factory access to the supported NCBI tables. The implementation accepts DNA or RNA codons, is case-insensitive, resolves IUPAC-ambiguous codons exactly as Biopython's ambiguous codon tables do, and throws on non-IUPAC symbols. All 27 NCBI translation tables of `gc.prt` Version 4.6 (1–6, 9–16, 21–33) are supported, built verbatim from the NCBI `ncbieaa`/`sncbieaa` strings; `Standard` (1), `VertebrateMitochondrial` (2), `YeastMitochondrial` (3) and `BacterialPlastid` (11) are also exposed as properties.[4][7][8]

## 2. Scientific / Formal Basis

### 2.1 Domain Context

The genetic code is the rule set that translates codons into amino acids. The original document records these core properties:[1][2][3]

| Property | Meaning |
|----------|---------|
| Triplet code | Each codon contains exactly 3 nucleotides. |
| Non-overlapping | Codons are read sequentially without overlap. |
| Degenerate | 64 codons encode 20 amino acids plus stop signals. |
| Nearly universal | Most organisms use the standard code with limited alternative tables. |

Representative tables (all 27 NCBI tables are supported; the full data is `gc.prt` v4.6):[4][7]

| Table | Name | Key AA Differences | Start Codons |
|-------|------|--------------------|--------------|
| 1 | Standard | Universal default | `AUG`, `UUG`, `CUG` |
| 2 | Vertebrate Mitochondrial | `AGA/AGG = *`, `AUA = M`, `UGA = W` | `AUG`, `AUA`, `AUU`, `AUC`, `GUG` |
| 3 | Yeast Mitochondrial | `CUU/CUC/CUA/CUG = T`, `AUA = M`, `UGA = W` | `AUG`, `AUA`, `GUG` |
| 11 | Bacterial/Plastid | Same amino-acid mapping as standard | `AUG`, `GUG`, `UUG`, `CUG`, `AUU`, `AUC`, `AUA` |

### 2.2 Core Model

Each table is built from the NCBI `gc.prt` strings: codon *i* (bases ordered T, C, A, G at each position) maps to `ncbieaa[i]`; it is a stop codon if `ncbieaa[i]` or `sncbieaa[i]` is `*`, and a start codon if `sncbieaa[i]` is `M`.[7] In tables 27, 28 and 31 some codons are *dual-coding* (amino acid in `ncbieaa`, `*` in `sncbieaa`): `Translate` returns the amino acid and `IsStopCodon` reports a stop, as in Biopython.[8]

Codon translation is a constant-time dictionary lookup after normalization of the input codon to uppercase RNA notation. An IUPAC-ambiguous codon (R, Y, S, W, K, M, B, D, H, V, N) is expanded into all concrete codons it stands for (Biopython `AmbiguousForwardTable` / `Seq._translate_str`):[8]

1. all expansions are stops → `*`;
2. stops mixed with amino acids ("possible stop", e.g. `TAN`) → `X`;
3. one amino acid (e.g. `GCN`) → that amino acid;
4. several amino acids → the most specific IUPAC ambiguous residue covering them: `B` (D/N, e.g. `RAY`), `Z` (E/Q, e.g. `SAR`), `J` (I/L, e.g. `MTH`), else `X`.

`IsStartCodon` / `IsStopCodon` accept an ambiguous codon only when every expansion is a start / stop (Biopython `list_ambiguous_codons`), e.g. `YTG` is a start and `TAR`, `TRA` are stops in table 1.

### 2.3 Modeling Assumptions

| ID | Assumption | Consequence if Violated |
|----|------------|--------------------------|
| ASM-01 | The caller selects a genetic code table appropriate for the organism or compartment being modeled. | Translation results can be correct for the chosen table but wrong biologically. |
| ASM-02 | The input to `Translate` is exactly one codon. | Longer or shorter inputs are rejected because codon translation is defined on triplets only. |

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | Every codon present in a table maps to exactly one amino acid character. | Each `GeneticCode` stores a dictionary from codon to amino acid. |
| INV-02 | Every supported table defines all 64 codons. | Each table is built from a 64-character NCBI `ncbieaa` string. |
| INV-03 | Stop codons translate to `'*'`. | The built-in codon tables store `'*'` for stop codons. |
| INV-04 | `ATG` and `AUG` translate identically. | `Translate` normalizes `T` to `U` before lookup. |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `[Translate] codon` | `string` | required | Single DNA or RNA codon. | Must be exactly 3 characters or `Translate` throws `ArgumentException`. |
| `[IsStartCodon] codon` | `string` | required | Codon to classify as a start codon. | Invalid length or `null` returns `false`. |
| `[IsStopCodon] codon` | `string` | required | Codon to classify as a stop codon. | Invalid length or `null` returns `false`. |
| `[GetCodonsForAminoAcid] aminoAcid` | `char` | required | Amino-acid code for reverse lookup. | Lookup is case-insensitive. |
| `[GetByTableNumber] tableNumber` | `int` | required | NCBI table number. | One of `SupportedTableNumbers` (1–6, 9–16, 21–33); otherwise `ArgumentException`. |

### 3.2 Output / Return Value

| Name | Type | Description |
|------|------|-------------|
| `Translate` result | `char` | Single-letter amino acid, `'*'` for stop, or for an ambiguous IUPAC codon the resolved residue / `B` / `Z` / `J` / `X` / `*` (Section 2.2). |
| `IsStartCodon` result | `bool` | Whether the normalized codon is in the table's start-codon set. |
| `IsStopCodon` result | `bool` | Whether the normalized codon is in the table's stop-codon set. |
| `GetCodonsForAminoAcid` result | `IEnumerable<string>` | All codons in the table that encode the supplied amino acid. |
| `GetByTableNumber` result | `GeneticCode` | The supported genetic code instance for the requested table number. |

### 3.3 Preconditions and Validation

`Translate` throws `ArgumentException` when the input is `null`, empty, or not exactly three characters long. After normalization, a codon found in the table is translated directly. If the codon is not present in the table but consists only of valid IUPAC nucleotide symbols, it is resolved by the ambiguity rule of Section 2.2. Non-IUPAC codons (including `X`) throw `ArgumentException`.[5]

## 4. Algorithm

### 4.1 High-Level Steps

1. Validate that the input codon is non-empty and exactly three characters.
2. Normalize the codon to uppercase RNA notation by replacing `T` with `U`.
3. If the normalized codon exists in the active table, return its amino acid.
4. Otherwise, if all characters are valid IUPAC nucleotide symbols, expand and resolve (Section 2.2; cached per table).
5. Otherwise, throw `ArgumentException`.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

Degeneracy of the standard code (table 1):[4]

| Number of codons | Amino acids |
|------------------|-------------|
| 1 | Met (`M`), Trp (`W`) |
| 2 | Phe, Tyr, His, Gln, Asn, Lys, Asp, Glu, Cys |
| 3 | Ile |
| 4 | Val, Pro, Thr, Ala, Gly |
| 6 | Leu, Ser, Arg |

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `Translate` | `O(1)` | `O(1)` | Constant-time normalization and dictionary lookup per codon. |
| `IsStartCodon` / `IsStopCodon` | `O(1)` | `O(1)` | Set membership after normalization. |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [GeneticCode.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Core/GeneticCode.cs)

- `GeneticCode.Translate(string)`
- `GeneticCode.IsStartCodon(string)`
- `GeneticCode.IsStopCodon(string)`
- `GeneticCode.GetCodonsForAminoAcid(char)`
- `GeneticCode.GetByTableNumber(int)`

### 5.2 Current Behavior

The repository supports all NCBI `gc.prt` v4.6 tables. `Translate` normalizes `T` to `U`, looks up the codon in the active table, and resolves ambiguous IUPAC codons as Biopython does (`GCN`→`A`, `TAR`→`*`, `RAY`→`B`, `ANN`/`NNN`/`TAN`→`X`). Inputs containing non-IUPAC symbols such as `XYZ` or `12G` throw `ArgumentException`.[5] `IsStartCodon` and `IsStopCodon` return `false` rather than throwing when the input is `null`, empty, or not exactly three characters long.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- All 27 NCBI translation tables, codon assignments and start/stop sets, from `gc.prt` Version 4.6.[7]
- IUPAC ambiguity resolution identical to Biopython's ambiguous codon tables.[8]
- Stop codons are represented as `'*'` and DNA codons are treated equivalently to RNA codons through `T -> U` normalization.[3][5]

Cross-check: every one of the 15³ IUPAC codons × 27 tables (translation, start and stop classification) equals Biopython 1.88 (`tests/.../TestData/GeneticCode/biopython_ambiguous_codons.tsv`); gc.prt v4.6 and Biopython's tables were verified identical.

**Documented divergences:**

- `X` as a *nucleotide* symbol is rejected (it is not an IUPAC nucleotide code); Biopython inconsistently treats it as `N` inside some codons (`GCX`→`A`) but rejects `XXX`.
- Dual-coding stop codons (tables 27, 28, 31) translate to their amino acid; context-dependent termination is not modelled (same as Biopython).

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| `AUG` | Returns `M`. | Standard start codon mapping. |
| `UAA`, `UAG`, `UGA` | Return `*` in the standard table. | Standard stop codons. |
| `ATG` | Returns `M`. | DNA is normalized to RNA. |
| Lowercase or mixed case input | Same result as uppercase. | Input is uppercased before lookup. |
| `ANN`, `NNN`, `TAN` | Return `X`. | Expansions mix stops and/or unrelated amino acids (Biopython). |
| `GCN`, `TAR`, `RAY`, `SAR`, `MTH` | Return `A`, `*`, `B`, `Z`, `J`. | Biopython ambiguous codon tables. |
| `XYZ`, `12G` | Throw `ArgumentException`. | They contain invalid non-IUPAC symbols. |

### 6.2 Limitations

`Translate` returns a single character per codon; ambiguity beyond B/Z/J is collapsed to `X`. Codon-context effects (selenocysteine/pyrrolysine recoding, context-dependent stops of tables 27/28/31) are not modelled.

## 7. Examples and Related Material

### 7.3 Related Tests, Evidence, or Documents

- Tests: [GeneticCodeTests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Core/GeneticCodeTests.cs) — covers `INV-01`, `INV-02`, `INV-03`, `INV-04`
- Test specification: [TRANS-CODON-001.md](../../../tests/TestSpecs/TRANS-CODON-001.md)
- Related algorithms: [Protein_Translation.md](Protein_Translation.md)

## 8. References

1. Wikipedia contributors. 2026. Genetic code. Wikipedia. https://en.wikipedia.org/wiki/Genetic_code
2. Wikipedia contributors. 2026. Start codon. Wikipedia. https://en.wikipedia.org/wiki/Start_codon
3. Wikipedia contributors. 2026. Stop codon. Wikipedia. https://en.wikipedia.org/wiki/Stop_codon
4. NCBI. 2026. The Genetic Codes. https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
5. Test specification: [TRANS-CODON-001.md](../../../tests/TestSpecs/TRANS-CODON-001.md)
6. Crick FH. 1968. The origin of the genetic code. Journal of Molecular Biology. N/A
7. NCBI. Genetic code table `gc.prt`, Version 4.6 (Elzanowski A, Ostell J). NCBI C++ Toolkit, `src/objects/seqfeat/gc.prt`. https://raw.githubusercontent.com/ncbi/ncbi-cxx-toolkit-public/master/src/objects/seqfeat/gc.prt
8. Biopython 1.88. `Bio/Data/CodonTable.py` (AmbiguousForwardTable, list_possible_proteins, list_ambiguous_codons) and `Bio/Seq.py` (`_translate_str`). https://raw.githubusercontent.com/biopython/biopython/master/Bio/Data/CodonTable.py
