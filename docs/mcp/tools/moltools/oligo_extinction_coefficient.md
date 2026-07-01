# oligo_extinction_coefficient

Per-base sum of 260 nm molar extinction contributions.

## Overview

| Property | Value |
|----------|-------|
| **Server** | MolTools |
| **Tool Name** | `oligo_extinction_coefficient` |
| **Method ID** | `ProbeDesigner.CalculateExtinctionCoefficient` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Estimates an oligonucleotide's molar extinction coefficient at 260 nm by summing per-base contributions (case-insensitive): A = 15400, C = 7400, G = 11500, T = 8700, U = 9900 M⁻¹·cm⁻¹; any other character contributes the fallback constant 10000.

## Core Documentation Reference

- Source: [ProbeDesigner.cs#L1341](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/ProbeDesigner.cs#L1341)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | Oligonucleotide sequence (non-empty). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `extinctionCoefficient` | number | ε₂₆₀ in M⁻¹·cm⁻¹. |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |

## Examples

### Example 1: `ACGT` → `15400 + 7400 + 11500 + 8700 = 43000`.

### Example 2: `N` → `10000` (unknown-base fallback).

## See Also

- [oligo_concentration_from_absorbance](oligo_concentration_from_absorbance.md), [analyze_oligo](analyze_oligo.md)
