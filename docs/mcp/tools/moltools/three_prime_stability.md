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

Computes Primer3's 3′-end stability (`oligotm.c` `end_oligodg(seq, 5)`): the SantaLucia (1998) nearest-neighbor ΔG°37 (kcal/mol, 1 M NaCl) of the primer's last 5 bases — the whole primer if shorter — with initiation +1.96, +0.05 per terminal A·T and +0.43 for a self-complementary sequence (for a 5-mer identical to terminal G·C +0.98 / A·T +1.03). Primer3 reports the same magnitude with the opposite sign as `PRIMER_*_END_STABILITY`. A more negative value means a more stable — and more mispriming-prone — 3′ end. N is accepted with Primer3's N parameters; any other non-ACGT character in the 3′ window is rejected.

## Core Documentation Reference

- Source: [PrimerDesigner.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/PrimerDesigner.cs) (`Calculate3PrimeStability`)

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
| 1002 | The 3′-terminal 5 bases may contain only A, C, G, T or N |

## Examples

### Example 1: `GCGCG` → `−6.86` (GC+CG+GC+CG + 0.98 + 0.98).

### Example 2: `TATAT` → `−0.86` (TA+AT+TA+AT + 1.03 + 1.03).

Only the last 5 bases matter, so `AAAAAGCGCG` also gives `−6.86`. A 4-mer is scored whole: `ACGT` → `−2.56` (Primer3 end_oligodg 2.56).

## See Also

- [primer_dimer](primer_dimer.md), [evaluate_primer](evaluate_primer.md)
