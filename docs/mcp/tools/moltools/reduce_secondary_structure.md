# reduce_secondary_structure

Reduce mRNA secondary structure by synonymous codon swaps.

## Overview

| Property | Value |
|----------|-------|
| **Server** | MolTools |
| **Tool Name** | `reduce_secondary_structure` |
| **Method ID** | `CodonOptimizer.ReduceSecondaryStructure` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Greedily swaps codons for synonymous alternatives that lower a heuristic local self-complementarity score within a sliding window (default 40 nt), reducing mRNA secondary structure while preserving the protein. Pairing uses the canonical `RnaSecondaryStructure.CanPair` rule (Watson-Crick plus the G·U wobble, ViennaRNA pair set). Output is upper-case RNA trimmed to whole codons; a sequence shorter than `window_size` keeps all its codons. This is a self-complementarity heuristic, not a folding model — use `rna_minimum_free_energy` for thermodynamics.

## Core Documentation Reference

- Source: [CodonOptimizer.cs#L584](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/CodonOptimizer.cs#L584)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `coding_sequence` | string | Yes | Coding sequence (DNA or RNA), non-empty. |
| `target_organism` | object | Yes | Preset id or inline custom table. |
| `window_size` | integer | No | Sliding-window size in nt (default 40, positive). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `optimizedSequence` | string | Sequence with reduced local self-complementarity. |

## Errors

| Code | Message |
|------|---------|
| 1001 | Coding sequence cannot be null or empty |
| 1002 | Window size must be positive |

## Examples

### Example 1: Below window → `ATGATG` → `AUGAUG` (codons unchanged, RNA-normalised).

### Example 2: A 45-nt sequence (> 40-nt window) is processed with synonymous swaps that preserve its length and protein.

## See Also

- [optimize_codons](optimize_codons.md)
