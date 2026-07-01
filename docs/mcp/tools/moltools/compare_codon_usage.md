# compare_codon_usage

Codon-frequency similarity between two coding sequences.

## Overview

| Property | Value |
|----------|-------|
| **Server** | MolTools |
| **Tool Name** | `compare_codon_usage` |
| **Method ID** | `CodonOptimizer.CompareCodonUsage` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Computes the codon-usage similarity between two coding sequences as `1 − ½·Σ|f1(codon) − f2(codon)|`, where `f` is the relative frequency of a codon within each sequence. Result is in `[0,1]`: 1 = identical codon distribution, 0 = fully disjoint. An input that is empty or has no complete codons contributes 0.

## Core Documentation Reference

- Source: [CodonOptimizer.cs#L929](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/CodonOptimizer.cs#L929)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence1` | string | Yes | First coding sequence (non-null). |
| `sequence2` | string | Yes | Second coding sequence (non-null). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `similarity` | number | Codon-usage similarity in 0..1. |

## Errors

| Code | Message |
|------|---------|
| 1001 | First sequence cannot be null |
| 1002 | Second sequence cannot be null |

## Examples

### Example 1: Identical sequences → 1.0

`ATGATG` vs `ATGATG` → 1.0.

### Example 2: Half overlap → 0.5

`ATGTTT` (AUG:0.5, UUU:0.5) vs `ATGATG` (AUG:1.0) → `½(|0.5−1| + |0.5−0|) = 0.5`, similarity = 0.5.

## Performance

- **Time Complexity:** O(n) in sequence length.
- **Space Complexity:** O(k) in distinct codons.

## See Also

- [codon_usage_statistics](codon_usage_statistics.md), [build_codon_table](build_codon_table.md)
