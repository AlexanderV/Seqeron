# three_prime_stability

Primer 3′-end nearest-neighbor stability (ΔG°37).

## Overview

| Property | Value |
|----------|-------|
| **Server** | MolTools |
| **Tool Name** | `three_prime_stability` |
| **Method ID** | `PrimerDesigner.Calculate3PrimeStability` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Computes the SantaLucia (1998) unified nearest-neighbor ΔG°37 (kcal/mol, 1 M NaCl) of a primer's last 5 bases, including initiation terms (terminal G·C = +0.98, terminal A·T = +1.03). This matches Primer3's `PRIMER_MAX_END_STABILITY`. A more negative value means a more stable — and more mispriming-prone — 3′ end. Sequences shorter than 5 bases return 0.

## Core Documentation Reference

- Source: [PrimerDesigner.cs#L427](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/PrimerDesigner.cs#L427)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | Primer sequence (non-empty). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `deltaG` | number | 3′-end ΔG°37 in kcal/mol. |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |

## Examples

### Example 1: `GCGCG` → `−6.86` (GC+CG+GC+CG + 0.98 + 0.98).

### Example 2: `TATAT` → `−0.86` (TA+AT+TA+AT + 1.03 + 1.03).

Only the last 5 bases matter, so `AAAAAGCGCG` also gives `−6.86`.

## See Also

- [primer_dimer](primer_dimer.md), [evaluate_primer](evaluate_primer.md)
