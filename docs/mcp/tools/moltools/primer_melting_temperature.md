# primer_melting_temperature

Compute a primer's melting temperature.

## Overview

| Property | Value |
|----------|-------|
| **Server** | MolTools |
| **Tool Name** | `primer_melting_temperature` |
| **Method ID** | `PrimerDesigner.CalculateMeltingTemperature` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Estimates the melting temperature (Tm, °C) of a primer. For fewer than 14 valid (ACGT) bases it uses the **Wallace rule** `Tm = 2·(A+T) + 4·(G+C)`; for 14 or more bases it uses the **Marmur–Doty** formula `Tm = 64.9 + 41·(GC − 16.4)/N`. Non-ACGT characters are ignored when counting.

## Core Documentation Reference

- Source: [PrimerDesigner.cs#L197](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/PrimerDesigner.cs#L197)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `primer` | string | Yes | Primer sequence (non-empty). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `tm` | number | Melting temperature in °C. |

## Errors

| Code | Message |
|------|---------|
| 1001 | Primer cannot be null or empty |

## Examples

### Example 1: Wallace rule

`ACGT` (2 AT, 2 GC) → `2·2 + 4·2 = 12` °C.

### Example 2: Marmur–Doty

`GCGCGCGCGCGCGCGCGCGC` (20 nt, GC=20) → `64.9 + 41·(20−16.4)/20 = 72.28` °C.

## See Also

- [primer_melting_temperature_salt](primer_melting_temperature_salt.md), [analyze_oligo](analyze_oligo.md)
