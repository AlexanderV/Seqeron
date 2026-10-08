# primer_dimer

Primer3 alignment-mode 3′-end primer-dimer check between two primers (PRIMER_PAIR_COMPL_END).

## Overview

| Property | Value |
|----------|-------|
| **Server** | MolTools |
| **Tool Name** | `primer_dimer` |
| **Method ID** | `PrimerDesigner.HasPrimerDimer` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Screens a primer pair for 3′-end dimer formation exactly as Primer3's alignment mode
(`PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT=0`, `PRIMER_PAIR_COMPL_END`): the score is the maximum of the
`dpal` end-anchored (`DPAL_GLOBAL_END`) alignment of `primer1` against the reverse complement of
`primer2` and of `primer2` against the reverse complement of `primer1` (+1 per complementary pair,
−1 per mismatch, −0.25 against N, −2 per single-base gap, floored at 0). Two primers whose 3′-terminal
k bases are reverse complements score k. A dimer is flagged when the score is at least
`min_complementarity` (default 4 = Primer3's default `PRIMER_PAIR_MAX_COMPL_END` 3.00 exceeded).
Identical poly-A primers are not a dimer (A·A cannot pair). Primer3's default thermodynamic check
(ntthal `PRIMER_PAIR_COMPL_ANY_TH` / `_COMPL_END_TH`) is the C# method
`PrimerDesigner.CalculatePrimer3PairComplementarity`.

## Core Documentation Reference

- Source: [PrimerDesigner.cs#L1898](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/PrimerDesigner.cs#L1898) (`HasPrimerDimer`, `CalculatePrimerDimerEndComplementarity`)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `primer1` | string | Yes | First primer sequence (non-empty). |
| `primer2` | string | Yes | Second primer sequence (non-empty). |
| `min_complementarity` | integer | No | Minimum 3′-end complementarity score to flag (default 4). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `hasDimer` | boolean | True if a 3′-dimer is flagged. |
| `complementaryBases` | integer | The compl_end score truncated to an integer (exact for ACGT primers). |
| `complEndScore` | number | Primer3 PRIMER_PAIR_COMPL_END score. |
| `complAnyScore` | number | Primer3 alignment-mode PRIMER_PAIR_COMPL_ANY: dpal local alignment of `primer1` with the reverse complement of `primer2` (same scoring; Primer3 default limit 8.00; `PrimerDesigner.CalculatePrimerDimerAnyComplementarity`). |

## Errors

| Code | Message |
|------|---------|
| 1001 | First primer cannot be null or empty |
| 1002 | Second primer cannot be null or empty |

## Examples

### Example 1: 3′-complementary pair

`AACCGGTTAACCATCGATCG` + `AACCGGTTAACGATCGAT` (3′ ends ATCGATCG / CGATCGAT) → `complementaryBases = 8`, `hasDimer = true` (primer3-py check_primers PRIMER_PAIR_0_COMPL_END = 8.0).

### Example 2: GGCC 3′ ends, threshold 5

`TTCAGTCAGTCAGTGGCC` + `ACTGACTGACTGAGGCC` → score 4; with `min_complementarity = 5` → `hasDimer = false` (true at the default 4).

### Example 3: identical poly-A

`AAAAAAAA` + `AAAAAAAA` → score 0, `hasDimer = false`.

## See Also

- [three_prime_stability](three_prime_stability.md), [evaluate_primer](evaluate_primer.md)
