# shannon_entropy

Calculate Shannon entropy of a sequence.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Sequence |
| **Tool Name** | `shannon_entropy` |
| **Method ID** | `SequenceStatistics.CalculateShannonEntropy` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Calculates Shannon entropy of a sequence, measuring the information content or randomness. Higher entropy values indicate more complexity. For DNA with equal base frequencies, maximum entropy is 2 bits per nucleotide. A homopolymer (e.g., AAAAAA) has entropy of 0.

## Core Documentation Reference

- Source: [SequenceStatistics.cs#L580](../../../../src/Seqeron/Algorithms/Seqeron.Genomics/SequenceStatistics.cs#L580)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | The sequence to analyze (min length: 1) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `entropy` | number | Shannon entropy value (bits per symbol) |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |

## Examples

### Example 1: DNA sequence entropy

**User Prompt:**
> What's the Shannon entropy of "ATGCATGCAT"?

**Expected Tool Call:**
```json
{
  "tool": "shannon_entropy",
  "arguments": {
    "sequence": "ATGCATGCAT"
  }
}
```

**Response:**
```json
{
  "entropy": 2.0
}
```

### Example 2: Low complexity sequence

**User Prompt:**
> Calculate entropy of "AAAAAAAAAA"

**Expected Tool Call:**
```json
{
  "tool": "shannon_entropy",
  "arguments": {
    "sequence": "AAAAAAAAAA"
  }
}
```

**Response:**
```json
{
  "entropy": 0.0
}
```

## Performance

- **Time Complexity:** O(n) where n is sequence length
- **Space Complexity:** O(k) where k is alphabet size

## See Also

- [linguistic_complexity](linguistic_complexity.md) - Linguistic complexity measure
- [nucleotide_composition](nucleotide_composition.md) - Base composition analysis
