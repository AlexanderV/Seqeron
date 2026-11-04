# oligo_extinction_coefficient

260 nm molar extinction coefficient: per-base sum (default) or nearest-neighbour model.

## Overview

| Property | Value |
|----------|-------|
| **Server** | MolTools |
| **Tool Name** | `oligo_extinction_coefficient` |
| **Method ID** | `ProbeDesigner.CalculateExtinctionCoefficient` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Estimates an oligonucleotide's molar extinction coefficient at 260 nm (case-insensitive).

- Default (`nearest_neighbor = false`, `ProbeDesigner.CalculateExtinctionCoefficient`): sum of mononucleotide contributions A = 15400, C = 7400, G = 11500, T = 8700, U = 9900 M⁻¹·cm⁻¹; any other character contributes the fallback constant 10000. Ignores base-stacking hypochromicity.
- `nearest_neighbor = true` (`ProbeDesigner.CalculateExtinctionCoefficientNearestNeighbor`): ε = Σ ε(dinucleotide N<sub>i</sub>N<sub>i+1</sub>) − Σ ε(internal mononucleotides), with the Cantor, Warshaw & Shapiro (1970) DNA table (`is_dna = true`) or the Warshaw & Tinoco (1966) RNA table (`is_dna = false`); a base outside A/C/G/T (DNA) or A/C/G/U (RNA) is an error.

## Core Documentation Reference

- Source: [ProbeDesigner.cs#L3615](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/ProbeDesigner.cs#L3615)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | Oligonucleotide sequence (non-empty). |
| `nearest_neighbor` | boolean | No | Use the nearest-neighbour model (default false). |
| `is_dna` | boolean | No | Nearest-neighbour table: true = DNA (default), false = RNA. |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `extinctionCoefficient` | number | ε₂₆₀ in M⁻¹·cm⁻¹. |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1001 | Nearest-neighbour ε260 needs only A/C/G/T (A/C/G/U) bases |

## Examples

### Example 1: `ACGT` → `15400 + 7400 + 11500 + 8700 = 43000`.

### Example 2: `N` → `10000` (unknown-base fallback).

### Example 3: `ACGT`, `nearest_neighbor = true` → `21200 + 18000 + 20000 − 7400 − 11500 = 40300`.

## See Also

- [oligo_concentration_from_absorbance](oligo_concentration_from_absorbance.md), [analyze_oligo](analyze_oligo.md)
