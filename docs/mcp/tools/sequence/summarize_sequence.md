# summarize_sequence

Generate comprehensive summary statistics for a DNA/RNA sequence.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Sequence |
| **Tool Name** | `summarize_sequence` |
| **Method ID** | `SequenceStatistics.SummarizeNucleotideSequence` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Generates comprehensive summary statistics for a DNA/RNA sequence in a single call. Includes length, nucleotide composition, GC content, Shannon entropy, linguistic complexity, and melting temperature.

## Core Documentation Reference

- Source: [SequenceStatistics.cs#L775](../../../../SuffixTree.Genomics/SequenceStatistics.cs#L775)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | The DNA or RNA sequence to analyze (min length: 1) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `length` | integer | Sequence length |
| `gcContent` | number | GC content percentage |
| `entropy` | number | Shannon entropy |
| `complexity` | number | Linguistic complexity |
| `meltingTemperature` | number | Melting temperature (°C) |
| `composition` | object | Nucleotide counts (A, T, G, C, U, N) |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |

## Examples

### Example 1: DNA sequence summary

**User Prompt:**
> Summarize the sequence "ATGCGATCGATCG"

**Expected Tool Call:**
```json
{
  "tool": "summarize_sequence",
  "arguments": {
    "sequence": "ATGCGATCGATCG"
  }
}
```

**Response:**
```json
{
  "length": 13,
  "gcContent": 46.15,
  "entropy": 1.98,
  "complexity": 0.85,
  "meltingTemperature": 38,
  "composition": { "A": 3, "T": 3, "G": 4, "C": 3, "U": 0, "N": 0 }
}
```

## Performance

- **Time Complexity:** O(n × k) where n is sequence length, k is max k-mer size
- **Space Complexity:** O(4^k) for k-mer storage

## See Also

- [nucleotide_composition](nucleotide_composition.md) - Basic composition only
- [shannon_entropy](shannon_entropy.md) - Entropy calculation
- [linguistic_complexity](linguistic_complexity.md) - Complexity calculation
- [melting_temperature](melting_temperature.md) - Tm calculation

