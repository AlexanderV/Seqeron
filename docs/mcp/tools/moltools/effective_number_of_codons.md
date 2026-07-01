# effective_number_of_codons

Wright's Effective Number of Codons (ENC / Nc).

## Overview

| Property | Value |
|----------|-------|
| **Server** | MolTools |
| **Tool Name** | `effective_number_of_codons` |
| **Method ID** | `CodonUsageAnalyzer.CalculateEnc` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Computes Wright's Effective Number of Codons (Nc), a single measure of how far a gene departs from uniform synonymous-codon usage. It uses per-amino-acid homozygosity `F̂` grouped by degeneracy class (Wright Eq. 1), class averaging (Eq. 4), the isoleucine Eq. 5a fallback, and `Nc = 2 + 9/F̂₂ + 1/F̂₃ + 5/F̂₄ + 3/F̂₆` (Eq. 3). The result is clamped to `[20, 61]` (20 = extreme bias, 61 = no bias).

## Core Documentation Reference

- Source: [CodonUsageAnalyzer.cs#L283](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/CodonUsageAnalyzer.cs#L283)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | Coding DNA sequence (frame 0), non-empty. |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `enc` | number | Effective Number of Codons, 20..61. |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |

## Examples

### Example 1: Fully-populated biased gene

With `F̂₂=0.6, F̂₃=0.2667, F̂₄=0.4333, F̂₆=0.3333`, `Nc = 2 + 9/0.6 + 1/0.2667 + 5/0.4333 + 3/0.3333 = 41.288461538461526` (value verified against an independent Wright/codonW reference).

## See Also

- [rscu](rscu.md), [codon_usage_statistics](codon_usage_statistics.md)
