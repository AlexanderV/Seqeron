# gc_content

Calculate GC content of a DNA/RNA sequence.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Sequence |
| **Tool Name** | `gc_content` |
| **Method ID** | `SequenceExtensions.CalculateGcContentFast` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Calculates the GC content (percentage of G and C nucleotides) of a DNA/RNA sequence. GC content is an important metric for primer design, genome analysis, and sequence characterization.

## Core Documentation Reference

- Source: [SequenceExtensions.cs#L41](../../../../src/Seqeron/Seqeron.Genomics/SequenceExtensions.cs#L41)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | The DNA or RNA sequence (min length: 1) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `gcContent` | number | GC content percentage (0-100) |
| `gcCount` | integer | Number of G and C nucleotides |
| `totalCount` | integer | Total sequence length |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |

## Examples

### Example 1: 50% GC content

**User Prompt:**
> What's the GC content of "ATGC"?

**Expected Tool Call:**
```json
{
  "tool": "gc_content",
  "arguments": {
    "sequence": "ATGC"
  }
}
```

**Response:**
```json
{
  "gcContent": 50,
  "gcCount": 2,
  "totalCount": 4
}
```

### Example 2: GC-rich sequence

**User Prompt:**
> Calculate GC content of "GCGCGC"

**Expected Tool Call:**
```json
{
  "tool": "gc_content",
  "arguments": {
    "sequence": "GCGCGC"
  }
}
```

**Response:**
```json
{
  "gcContent": 100,
  "gcCount": 6,
  "totalCount": 6
}
```

## Performance

- **Time Complexity:** O(n) where n is sequence length
- **Space Complexity:** O(1)

## See Also

- [nucleotide_composition](nucleotide_composition.md) - Full composition analysis
- [summarize_sequence](summarize_sequence.md) - Comprehensive statistics
