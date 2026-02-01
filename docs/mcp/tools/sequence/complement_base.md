# complement_base

Get the Watson-Crick complement of a single nucleotide base.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Sequence |
| **Tool Name** | `complement_base` |
| **Method ID** | `SequenceExtensions.GetComplementBase` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Returns the Watson-Crick complementary base for a given nucleotide. For DNA: A↔T, C↔G. For RNA: A↔U, C↔G. Input is case-insensitive, output is always uppercase.

## Core Documentation Reference

- Source: [SequenceExtensions.cs#L83](../../../../src/Seqeron/Seqeron.Genomics/SequenceExtensions.cs#L83)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `nucleotide` | string | Yes | The nucleotide base (A, T, G, C, or U), exactly 1 character |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `complement` | string | The complementary base |
| `original` | string | The input base |

## Errors

| Code | Message |
|------|---------|
| 1001 | Must provide exactly one nucleotide character |

## Examples

### Example 1: DNA complement

**User Prompt:**
> What's the complement of "A"?

**Expected Tool Call:**
```json
{
  "tool": "complement_base",
  "arguments": {
    "nucleotide": "A"
  }
}
```

**Response:**
```json
{
  "complement": "T",
  "original": "A"
}
```

### Example 2: RNA complement

**User Prompt:**
> Get complement of "U"

**Expected Tool Call:**
```json
{
  "tool": "complement_base",
  "arguments": {
    "nucleotide": "U"
  }
}
```

**Response:**
```json
{
  "complement": "A",
  "original": "U"
}
```

## Performance

- **Time Complexity:** O(1)
- **Space Complexity:** O(1)

## See Also

- [dna_reverse_complement](dna_reverse_complement.md) - Full sequence reverse complement
- [nucleotide_composition](nucleotide_composition.md) - Sequence composition analysis
