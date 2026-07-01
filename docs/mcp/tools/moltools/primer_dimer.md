# primer_dimer

Heuristic 3′-end primer-dimer check between two primers.

## Overview

| Property | Value |
|----------|-------|
| **Server** | MolTools |
| **Tool Name** | `primer_dimer` |
| **Method ID** | `PrimerDesigner.HasPrimerDimer` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Screens a primer pair for 3′-end dimer formation. It reverse-complements `primer2`, aligns the up-to-8-bp 3′ ends, and counts complementary positions. A dimer is flagged when the count reaches `min_complementarity` (default 4). The tool returns both the boolean flag and the complementary-base count.

## Core Documentation Reference

- Source: [PrimerDesigner.cs#L398](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/PrimerDesigner.cs#L398)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `primer1` | string | Yes | First primer sequence (non-empty). |
| `primer2` | string | Yes | Second primer sequence (non-empty). |
| `min_complementarity` | integer | No | Minimum complementary 3′-end bases to flag (default 4). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `hasDimer` | boolean | True if a 3′-dimer is flagged. |
| `complementaryBases` | integer | Complementary positions in the 3′-end window. |

## Errors

| Code | Message |
|------|---------|
| 1001 | First primer cannot be null or empty |
| 1002 | Second primer cannot be null or empty |

## Examples

### Example 1: Strong dimer

`AAAAAAAA` + `AAAAAAAA` → revcomp `TTTTTTTT`; all 8 positions complementary → `hasDimer = true`, `complementaryBases = 8`.

### Example 2: Below threshold

`AAAAAAAA` + `AAAACCCC` (revcomp `GGGGTTTT`) → 4 complementary; with `min_complementarity = 5` → `hasDimer = false`.

## See Also

- [three_prime_stability](three_prime_stability.md), [evaluate_primer](evaluate_primer.md)
