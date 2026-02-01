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

Translates an RNA sequence to a protein sequence using the standard genetic code. RNA uses U (uracil) instead of T (thymine). The translation reads codons and converts them to amino acids.

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | The RNA sequence to translate (uses U, not T) |
| `frame` | integer | No | Reading frame (0, 1, or 2, default: 0) |
| `toFirstStop` | boolean | No | Stop at first stop codon (default: false) |

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

## See Also

- [translate_dna](translate_dna.md) - Translate DNA to protein
- [rna_from_dna](rna_from_dna.md) - Transcribe DNA to RNA
