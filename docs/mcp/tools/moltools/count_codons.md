# count_codons

Count ACGT codons in a coding DNA sequence (frame 0).

## Overview

| Property | Value |
|----------|-------|
| **Server** | MolTools |
| **Tool Name** | `count_codons` |
| **Method ID** | `CodonUsageAnalyzer.CountCodons` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Splits the input into non-overlapping frame-0 triplets and tallies how many times each codon occurs. A trailing partial codon (length < 3) and any codon containing a non-ACGT character are skipped silently. Matching is case-insensitive.

## Core Documentation Reference

- Source: [CodonUsageAnalyzer.cs#L29](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/CodonUsageAnalyzer.cs#L29)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | Coding DNA sequence (frame 0), non-empty. |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `counts` | object | Codon → occurrence count. |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |

## Examples

### Example 1: `ATGATGTTT` → `{ ATG: 2, TTT: 1 }`

### Example 2: Skipping invalid/partial codons

`ATGANGTTTC` → `ANG` (contains N) is skipped and trailing `C` is ignored → `{ ATG: 1, TTT: 1 }`.

## Performance

- **Time Complexity:** O(n) in sequence length.
- **Space Complexity:** O(k) in distinct codons.

## See Also

- [rscu](rscu.md), [codon_usage_statistics](codon_usage_statistics.md)
