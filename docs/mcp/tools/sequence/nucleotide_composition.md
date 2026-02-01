# nucleotide_composition

Calculate nucleotide composition of a DNA/RNA sequence.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Sequence |
| **Tool Name** | `nucleotide_composition` |
| **Method ID** | `SequenceStatistics.CalculateNucleotideComposition` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Calculates the nucleotide composition of a DNA or RNA sequence, including counts of each nucleotide (A, T, G, C, U) and the GC content ratio. Works with both DNA (containing T) and RNA (containing U) sequences.

## Core Documentation Reference

- Source: [SequenceStatistics.cs#L48](../../../../src/Seqeron/Seqeron.Genomics/SequenceStatistics.cs#L48)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | The DNA or RNA sequence to analyze (min length: 1) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `length` | integer | Total sequence length |
| `A` | integer | Count of adenine |
| `T` | integer | Count of thymine |
| `G` | integer | Count of guanine |
| `C` | integer | Count of cytosine |
| `U` | integer | Count of uracil (RNA) |
| `other` | integer | Count of other characters |
| `gcContent` | number | GC content ratio (0-1) |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |

## Examples

### Example 1: DNA composition

**User Prompt:**
> What's the nucleotide composition of "ATGCATGC"?

**Expected Tool Call:**
```json
{
  "tool": "nucleotide_composition",
  "arguments": {
    "sequence": "ATGCATGC"
  }
}
```

**Response:**
```json
{
  "length": 8,
  "A": 2,
  "T": 2,
  "G": 2,
  "C": 2,
  "U": 0,
  "other": 0,
  "gcContent": 0.5
}
```

### Example 2: RNA composition

**User Prompt:**
> Analyze the RNA sequence "AUGCAUGC"

**Expected Tool Call:**
```json
{
  "tool": "nucleotide_composition",
  "arguments": {
    "sequence": "AUGCAUGC"
  }
}
```

**Response:**
```json
{
  "length": 8,
  "A": 2,
  "T": 0,
  "G": 2,
  "C": 2,
  "U": 2,
  "other": 0,
  "gcContent": 0.5
}
```

## Performance

- **Time Complexity:** O(n) where n is sequence length
- **Space Complexity:** O(1)

## See Also

- [amino_acid_composition](amino_acid_composition.md) - Get amino acid frequencies
- [dna_validate](dna_validate.md) - Validate DNA sequences
