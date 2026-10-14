# translate_rna

Translate RNA sequence to protein.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Sequence |
| **Tool Name** | `translate_rna` |
| **Method ID** | `Translator.Translate(RnaSequence)` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Translates an RNA sequence to a protein sequence using the standard genetic code. RNA uses U (uracil) instead of T (thymine). The translation reads codons and converts them to amino acids. The NCBI translation table is selected with `geneticCode` (default 1 = Standard; any NCBI gc.prt table 1-6, 9-16, 21-33), as Biopython `Seq.translate(table=n, to_stop=...)`. Input letters are case-insensitive A, C, G, U plus the IUPAC ambiguity codes R, Y, S, W, K, M, B, D, H, V, N (and X when `allowNucleotideX` is true); T (use translate_dna) is rejected, and any other character is an error naming the character and its 0-based position. An ambiguous codon translates to the single amino acid all its expansions encode (e.g. GCN → A), to B/Z/J, or to X; trailing bases that do not form a full codon are ignored. With `allowNucleotideX`, X stands for any base (Biopython `ambiguous_dna_values['X']`), but an X-containing codon that may be a stop is an error (Biopython `TranslationError`). `toFirstStop` is rejected for tables 27, 28, 31 (dual-coding stop codons; Biopython `ValueError`).

## Core Documentation Reference

- Source: [Translator.cs#L97](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Core/Translator.cs#L97)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | The RNA sequence to translate (uses U, not T) |
| `frame` | integer | No | Reading frame (0, 1, or 2, default: 0) |
| `toFirstStop` | boolean | No | Stop at first stop codon (default: false) |
| `geneticCode` | integer | No | NCBI translation table number (default: 1 = Standard; 1-6, 9-16, 21-33) |
| `allowNucleotideX` | boolean | No | Accept nucleotide X as any base, Biopython semantics (default: false) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `protein` | string | Translated protein sequence |
| `frame` | integer | Reading frame used |
| `rnaLength` | integer | Original RNA length |

## Examples

### Example 1: Basic RNA translation

```json
{
  "tool": "translate_rna",
  "arguments": {
    "sequence": "AUGCGAUGA",
    "frame": 0
  }
}
```

**Response:**
```json
{
  "protein": "MR*",
  "frame": 0,
  "rnaLength": 9
}
```

### Example 2: Stop at first stop codon

```json
{
  "tool": "translate_rna",
  "arguments": {
    "sequence": "AUGCGAUGAAAA",
    "toFirstStop": true
  }
}
```

**Response:**
```json
{
  "protein": "MR",
  "frame": 0,
  "rnaLength": 12
}
```

### Example 3: Vertebrate mitochondrial code (table 2)

```json
{
  "tool": "translate_rna",
  "arguments": {
    "sequence": "AUGAGAUGA",
    "geneticCode": 2
  }
}
```

**Response** (Biopython 1.88 `Seq("AUGAGAUGA").translate(table=2)` = `M*W`: AGA is a stop and TGA is Trp under table 2):
```json
{
  "protein": "M*W",
  "frame": 0,
  "rnaLength": 9
}
```

### Example 4: IUPAC ambiguity codes

`AUGGCNUAR` → `MA*` (GCN → Ala; TAR → stop), and `MA` with `toFirstStop: true` (Biopython 1.88).

## See Also

- [translate_dna](translate_dna.md) - Translate DNA to protein
- [rna_from_dna](rna_from_dna.md) - Transcribe DNA to RNA
