# linguistic_complexity

Calculate linguistic complexity of a sequence.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Sequence |
| **Tool Name** | `linguistic_complexity` |
| **Method ID** | `SequenceStatistics.CalculateLinguisticComplexity` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Calculates linguistic complexity of a sequence based on k-mer diversity. The measure compares observed k-mers to maximum possible k-mers for each k from 1 to maxK. Values range from 0 (minimal complexity, like homopolymers) to 1 (maximum complexity with diverse k-mers).

## Core Documentation Reference

- Source: [SequenceStatistics.cs#L615](../../../../src/Seqeron/Seqeron.Genomics/SequenceStatistics.cs#L615)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | The sequence to analyze (min length: 1) |
| `maxK` | integer | No | Maximum k-mer length to consider (default: 6) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `complexity` | number | Linguistic complexity (0-1) |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |

## Examples

### Example 1: High complexity sequence

**User Prompt:**
> What's the linguistic complexity of "ATGCATGCATGC"?

**Expected Tool Call:**
```json
{
  "tool": "linguistic_complexity",
  "arguments": {
    "sequence": "ATGCATGCATGC"
  }
}
```

**Response:**
```json
{
  "complexity": 0.95
}
```

### Example 2: Low complexity sequence

**User Prompt:**
> Calculate complexity of "AAAAAAAAAA"

**Expected Tool Call:**
```json
{
  "tool": "linguistic_complexity",
  "arguments": {
    "sequence": "AAAAAAAAAA"
  }
}
```

**Response:**
```json
{
  "complexity": 0.17
}
```

## Performance

- **Time Complexity:** O(n × maxK) where n is sequence length
- **Space Complexity:** O(4^maxK) for k-mer storage

## See Also

- [shannon_entropy](shannon_entropy.md) - Shannon entropy measure
- [nucleotide_composition](nucleotide_composition.md) - Base composition analysis
