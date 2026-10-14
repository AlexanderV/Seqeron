# count_codons

Count codons in a coding DNA or RNA sequence (frame 0).

## Overview

| Property | Value |
|----------|-------|
| **Server** | MolTools |
| **Tool Name** | `count_codons` |
| **Method ID** | `CodonUsageAnalyzer.CountCodons` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Splits the input into non-overlapping frame-0 triplets and tallies how many times each codon occurs. Input may be DNA or RNA: U is read as T and the keys are uppercase DNA codons; matching is case-insensitive. A trailing partial codon (length < 3) is ignored, and a triplet containing any other character (e.g. N) is skipped without shifting the reading frame. Counting is independent of the genetic code.

## Core Documentation Reference

- Source: [CodonUsageAnalyzer.cs#L34](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/CodonUsageAnalyzer.cs#L34)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | Coding sequence (frame 0), DNA or RNA, case-insensitive (U read as T; keys are uppercase DNA codons), non-empty. |

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

### Example 3: RNA input

`augaugUUU` → `{ ATG: 2, TTT: 1 }` (U read as T, case-insensitive). The count does not depend on the genetic code.

## Performance

- **Time Complexity:** O(n) in sequence length.
- **Space Complexity:** O(k) in distinct codons.

## See Also

- [rscu](rscu.md), [codon_usage_statistics](codon_usage_statistics.md)
