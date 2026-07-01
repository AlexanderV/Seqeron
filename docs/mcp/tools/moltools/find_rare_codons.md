# find_rare_codons

Report codons rarer than a threshold in the target organism.

## Overview

| Property | Value |
|----------|-------|
| **Server** | MolTools |
| **Tool Name** | `find_rare_codons` |
| **Method ID** | `CodonOptimizer.FindRareCodons` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Splits the coding sequence into frame-0 codons (T→U), and reports each codon whose frequency in the target organism's codon-usage table is below `threshold` (default 0.15). Each result carries the codon's 0-based nucleotide position, RNA codon, amino acid, and frequency. Rare codons slow translation and are candidates for optimization.

## Core Documentation Reference

- Source: [CodonOptimizer.cs#L663](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/CodonOptimizer.cs#L663)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `coding_sequence` | string | Yes | Coding sequence (DNA or RNA), non-empty. |
| `target_organism` | object | Yes | Preset id (`EColiK12` \| `Yeast` \| `Human`) or inline custom table. |
| `threshold` | number | No | Frequency threshold (default 0.15). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `rareCodons` | array | Rare codons, each `{position, codon, aminoAcid, frequency}`. |

## Errors

| Code | Message |
|------|---------|
| 1001 | Coding sequence cannot be null or empty |

## Examples

### Example 1: Table `{UUU:0.8, UUC:0.1}`, threshold 0.15

`TTTTTC` → `UUC` at position 3 (frequency 0.1 < 0.15) is the only rare codon.

## See Also

- [optimize_codons](optimize_codons.md), [cai_from_organism_table](cai_from_organism_table.md)
