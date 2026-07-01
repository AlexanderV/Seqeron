# cai_from_organism_table

Codon Adaptation Index against an organism codon-usage frequency table.

## Overview

| Property | Value |
|----------|-------|
| **Server** | MolTools |
| **Tool Name** | `cai_from_organism_table` |
| **Method ID** | `CodonOptimizer.CalculateCAI` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Computes CAI (Sharp & Li 1987) for a coding sequence relative to an organism codon-usage **frequency** table. For each sense codon the relative adaptiveness is `w = f(codon) / max f(synonymous codons)`; CAI is the geometric mean of `w` over all scored codons. Stop codons are skipped, codons whose amino-acid group has no frequency data are skipped, and `w` is clamped at `1e-6` to avoid `ln(0)` on incomplete tables. Empty sequence → 0.

This is distinct from [codon_adaptation_index](codon_adaptation_index.md), which takes a reference **RSCU** dictionary.

## Core Documentation Reference

- Source: [CodonOptimizer.cs#L473](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/CodonOptimizer.cs#L473)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `coding_sequence` | string | Yes | Coding sequence (DNA or RNA), non-empty. |
| `target_organism` | object | Yes | Preset id (`EColiK12` \| `Yeast` \| `Human`) or an inline custom table with `codonFrequencies` (RNA alphabet). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `cai` | number | Codon Adaptation Index in 0..1. |

## Errors

| Code | Message |
|------|---------|
| 1001 | Coding sequence cannot be null or empty |
| 1002 | Unknown codon-usage preset (expected EColiK12 \| Yeast \| Human) |

## Examples

### Example 1: Non-optimal codon

Table `{UUU:0.8, UUC:0.2}`, sequence `TTC` → `UUC`, `w = 0.2/0.8 = 0.25`, CAI = 0.25.

```json
{ "tool": "cai_from_organism_table", "arguments": { "coding_sequence": "TTC", "target_organism": { "codonFrequencies": { "UUU": 0.8, "UUC": 0.2 } } } }
```
→ `{ "cai": 0.25 }`

### Example 2: Geometric mean of two codons

Sequence `TTTTTC` → `UUU` (w=1), `UUC` (w=0.25); CAI = √(1·0.25) = 0.5.

## Performance

- **Time Complexity:** O(n) in sequence length.
- **Space Complexity:** O(1) beyond the table.

## See Also

- [codon_adaptation_index](codon_adaptation_index.md) - CAI from a reference RSCU dictionary.
- [build_codon_table](build_codon_table.md) - Build a custom frequency table.
