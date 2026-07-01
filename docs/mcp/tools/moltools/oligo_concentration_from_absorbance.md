# oligo_concentration_from_absorbance

Beer–Lambert oligonucleotide concentration from A260.

## Overview

| Property | Value |
|----------|-------|
| **Server** | MolTools |
| **Tool Name** | `oligo_concentration_from_absorbance` |
| **Method ID** | `ProbeDesigner.CalculateConcentration` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Applies the Beer–Lambert law to convert a 260 nm absorbance reading into an oligonucleotide concentration in micromolar: `c (µM) = A₂₆₀ / (ε · path) · 1e6`, where ε is the molar extinction coefficient (M⁻¹·cm⁻¹) and `path` is the cuvette path length in cm.

## Core Documentation Reference

- Source: [ProbeDesigner.cs#L1367](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/ProbeDesigner.cs#L1367)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `absorbance260` | number | Yes | Absorbance at 260 nm. |
| `extinction_coefficient` | number | Yes | ε in M⁻¹·cm⁻¹ (positive). |
| `path_length` | number | No | Path length in cm (default 1.0, positive). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `concentrationMicromolar` | number | Concentration in µM. |

## Errors

| Code | Message |
|------|---------|
| 1001 | Extinction coefficient must be positive |
| 1002 | Path length must be positive |

## Examples

### Example 1: `A260 = 1`, `ε = 10000`, `path = 1` → `1/10000·1e6 = 100` µM.

### Example 2: `path = 2` cm halves the result → `50` µM.

## See Also

- [oligo_extinction_coefficient](oligo_extinction_coefficient.md), [analyze_oligo](analyze_oligo.md)
