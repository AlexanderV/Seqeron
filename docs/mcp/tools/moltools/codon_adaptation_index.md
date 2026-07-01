# codon_adaptation_index

Codon Adaptation Index from a reference RSCU table.

## Overview

| Property | Value |
|----------|-------|
| **Server** | MolTools |
| **Tool Name** | `codon_adaptation_index` |
| **Method ID** | `CodonUsageAnalyzer.CalculateCai` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Computes CAI (Sharp & Li 1987) for a coding DNA sequence using a caller-supplied **reference RSCU** table (DNA-alphabet codons). Relative adaptiveness `w = RSCU(codon) / max RSCU(synonymous)`; CAI is the geometric mean of `w` over scored codons. Single-codon amino acids (Met, Trp), stop codons, and codons with `w = 0` are excluded from the mean.

This differs from [cai_from_organism_table](cai_from_organism_table.md), which takes a codon-usage **frequency** table.

## Core Documentation Reference

- Source: [CodonUsageAnalyzer.cs#L134](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/CodonUsageAnalyzer.cs#L134)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | Coding DNA sequence (frame 0), non-empty. |
| `reference_rscu` | object | Yes | Reference RSCU table: DNA codon → RSCU value. |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `cai` | number | Codon Adaptation Index in 0..1. |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1002 | Reference RSCU table cannot be null |

## Examples

### Example 1: Non-optimal codon

RSCU `{TTT:1.0, TTC:0.5}` → w[TTC] = 0.5. Sequence `TTC` → CAI = 0.5.

### Example 2: Two codons

Sequence `TTTTTC` → w = {1.0, 0.5}; CAI = √0.5 ≈ 0.7071.

## Performance

- **Time Complexity:** O(n) in sequence length.
- **Space Complexity:** O(1) beyond the reference table.

## See Also

- [cai_from_organism_table](cai_from_organism_table.md) - CAI from a frequency table.
- [rscu](rscu.md) - Compute RSCU for a sequence.
