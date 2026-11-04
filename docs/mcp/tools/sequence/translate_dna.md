# translate_dna

Translate DNA sequence to protein.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Sequence |
| **Tool Name** | `translate_dna` |
| **Method ID** | `Translator.Translate(DnaSequence)` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Translates a DNA sequence to a protein sequence using the standard genetic code. The translation reads codons (triplets of nucleotides) and converts them to amino acids. Stop codons (*) are included unless `toFirstStop` is true. The NCBI translation table is selected with `geneticCode` (default 1 = Standard; any NCBI gc.prt table 1-6, 9-16, 21-33), as Biopython `Seq.translate(table=n, to_stop=...)`. Input letters are case-insensitive A, C, G, T plus the IUPAC ambiguity codes R, Y, S, W, K, M, B, D, H, V, N (and X when `allowNucleotideX` is true); U (use translate_rna) is rejected, and any other character is an error naming the character and its 0-based position. An ambiguous codon translates to the single amino acid all its expansions encode (e.g. GCN → A), to B/Z/J, or to X; trailing bases that do not form a full codon are ignored. With `allowNucleotideX`, X stands for any base (Biopython `ambiguous_dna_values['X']`), but an X-containing codon that may be a stop is an error (Biopython `TranslationError`). `toFirstStop` is rejected for tables 27, 28, 31 (dual-coding stop codons; Biopython `ValueError`).

## Core Documentation Reference

- Source: [Translator.cs#L97](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Core/Translator.cs#L97)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | The DNA sequence to translate |
| `frame` | integer | No | Reading frame (0, 1, or 2, default: 0) |
| `toFirstStop` | boolean | No | Stop at first stop codon (default: false) |
| `geneticCode` | integer | No | NCBI translation table number (default: 1 = Standard; 1-6, 9-16, 21-33) |
| `allowNucleotideX` | boolean | No | Accept nucleotide X as any base, Biopython semantics (default: false) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `protein` | string | Translated protein sequence |
| `frame` | integer | Reading frame used |
| `dnaLength` | integer | Original DNA length |

## Examples

### Example 1: Basic translation

```json
{
  "tool": "translate_dna",
  "arguments": {
    "sequence": "ATGCGATGA",
    "frame": 0
  }
}
```

**Response:**
```json
{
  "protein": "MR*",
  "frame": 0,
  "dnaLength": 9
}
```

### Example 2: Stop at first stop codon

```json
{
  "tool": "translate_dna",
  "arguments": {
    "sequence": "ATGCGATGAAAA",
    "toFirstStop": true
  }
}
```

**Response:**
```json
{
  "protein": "MR",
  "frame": 0,
  "dnaLength": 12
}
```

### Example 3: Vertebrate mitochondrial code (table 2)

```json
{
  "tool": "translate_dna",
  "arguments": {
    "sequence": "ATGAGATGA",
    "geneticCode": 2
  }
}
```

**Response** (Biopython 1.88 `Seq("ATGAGATGA").translate(table=2)` = `M*W`: AGA is a stop and TGA is Trp under table 2):
```json
{
  "protein": "M*W",
  "frame": 0,
  "dnaLength": 9
}
```

### Example 4: IUPAC ambiguity codes

`ATGGCNTAR` → `MA*` (GCN → Ala; TAR → stop), and `MA` with `toFirstStop: true` (Biopython 1.88).

## See Also

- [translate_rna](translate_rna.md) - Translate RNA to protein
- [rna_from_dna](rna_from_dna.md) - Transcribe DNA to RNA
