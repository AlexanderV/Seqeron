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

Translates a DNA sequence to a protein sequence using the standard genetic code. The translation reads codons (triplets of nucleotides) and converts them to amino acids. Stop codons (*) are included unless `toFirstStop` is true.

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | The DNA sequence to translate |
| `frame` | integer | No | Reading frame (0, 1, or 2, default: 0) |
| `toFirstStop` | boolean | No | Stop at first stop codon (default: false) |

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

## See Also

- [translate_rna](translate_rna.md) - Translate RNA to protein
- [rna_from_dna](rna_from_dna.md) - Transcribe DNA to RNA
