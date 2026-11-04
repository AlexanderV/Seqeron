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

Calculates the nucleotide composition of a DNA or RNA sequence, including counts of each nucleotide (A, T, G, C, U), of N (any base) and of all other characters, and the GC content ratio. Works with both DNA (containing T) and RNA (containing U) sequences.

## Core Documentation Reference

- Source: [SequenceStatistics.cs#L62](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/SequenceStatistics.cs#L62)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | The DNA or RNA sequence to analyze (min length: 1) |

## Output Schema

Wire names are the MCP serializer's camelCase forms of the result record (`A` → `a`, `GcContent` → `gcContent`, …), as locked by `B03McpDocExamplesTests.NucleotideComposition_WireNames_AreCamelCase`.

| Field | Type | Description |
|-------|------|-------------|
| `length` | integer | Total sequence length |
| `a` | integer | Count of adenine |
| `t` | integer | Count of thymine |
| `g` | integer | Count of guanine |
| `c` | integer | Count of cytosine |
| `u` | integer | Count of uracil (RNA) |
| `other` | integer | Count of other characters (not A/C/G/T/U/N) |
| `n` | integer | Count of N (any base; added FIN-B03 A22) |
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
  "a": 2,
  "t": 2,
  "g": 2,
  "c": 2,
  "u": 0,
  "other": 0,
  "gcContent": 0.5,
  "n": 0
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
  "a": 2,
  "t": 0,
  "g": 2,
  "c": 2,
  "u": 2,
  "other": 0,
  "gcContent": 0.5,
  "n": 0
}
```

## Performance

- **Time Complexity:** O(n) where n is sequence length
- **Space Complexity:** O(1)

## See Also

- [amino_acid_composition](amino_acid_composition.md) - Get amino acid frequencies
- [dna_validate](dna_validate.md) - Validate DNA sequences
