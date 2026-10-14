# molecular_weight_nucleotide

Calculate molecular weight of a DNA or RNA sequence.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Sequence |
| **Tool Name** | `molecular_weight_nucleotide` |
| **Method ID** | `SequenceStatistics.CalculateNucleotideMolecularWeight` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Calculates the molecular weight of a DNA or RNA sequence in Daltons (Da). Uses average molecular weights for nucleotides as monophosphates. Accounts for different weights between DNA (deoxyribonucleotides) and RNA (ribonucleotides).

## Core Documentation Reference

- Source: [SequenceStatistics.cs#L491](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/SequenceStatistics.cs#L491)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | The DNA or RNA sequence (min length: 1) |
| `isDna` | boolean | No | True for DNA, false for RNA (default: true) |
| `doubleStranded` | boolean | No | Add the complementary strand (Biopython `double_stranded=True`; default false). AGC DNA → 1890.2057 |
| `circular` | boolean | No | Circular molecule, one more water lost per strand (default false). AGC DNA → 931.5942; ds+circular 1854.1751 |
| `monoisotopic` | boolean | No | Monoisotopic masses, water 18.010565 (default false). ATGCGTAC ds+circular → 4940.810128 |
| `strict` | boolean | No | Reject symbols without a mass (N; U in DNA; T in RNA) instead of skipping them (default false) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `molecularWeight` | number | Molecular weight in Daltons |
| `unit` | string | Unit of measurement ("Da") |
| `sequenceType` | string | "DNA" or "RNA" |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1001 | strict: '<c>' is not a valid unambiguous letter for DNA/RNA (strict mode). |

## Examples

### Example 1: DNA molecular weight

**User Prompt:**
> What's the molecular weight of DNA "ATGC"?

**Expected Tool Call:**
```json
{
  "tool": "molecular_weight_nucleotide",
  "arguments": {
    "sequence": "ATGC",
    "isDna": true
  }
}
```

**Response:**
```json
{
  "molecularWeight": 1307.8,
  "unit": "Da",
  "sequenceType": "DNA"
}
```

### Example 2: RNA molecular weight

**User Prompt:**
> Calculate molecular weight of RNA "AUGC"

**Expected Tool Call:**
```json
{
  "tool": "molecular_weight_nucleotide",
  "arguments": {
    "sequence": "AUGC",
    "isDna": false
  }
}
```

**Response:**
```json
{
  "molecularWeight": 1357.2,
  "unit": "Da",
  "sequenceType": "RNA"
}
```

## Performance

- **Time Complexity:** O(n) where n is sequence length
- **Space Complexity:** O(1)

## See Also

- [molecular_weight_protein](molecular_weight_protein.md) - Protein molecular weight
- [nucleotide_composition](nucleotide_composition.md) - DNA/RNA composition
