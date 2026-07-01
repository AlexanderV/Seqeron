# primer_melting_temperature_salt

Salt-corrected primer melting temperature.

## Overview

| Property | Value |
|----------|-------|
| **Server** | MolTools |
| **Tool Name** | `primer_melting_temperature_salt` |
| **Method ID** | `PrimerDesigner.CalculateMeltingTemperatureWithSalt` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Takes the Wallace/Marmur–Doty base Tm (see [primer_melting_temperature](primer_melting_temperature.md)) and adds a Schildkraut–Lifson monovalent-cation correction `16.6·log10([Na+]/1000)` (Na+ in mM), then rounds to one decimal place.

## Core Documentation Reference

- Source: [PrimerDesigner.cs#L227](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/PrimerDesigner.cs#L227)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `primer` | string | Yes | Primer sequence (non-empty). |
| `na_concentration` | number | No | Na+ concentration in mM (default 50, must be positive). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `tm` | number | Salt-corrected Tm (°C, 1 decimal). |

## Errors

| Code | Message |
|------|---------|
| 1001 | Primer cannot be null or empty |
| 1002 | Na+ concentration must be positive |

## Examples

### Example 1: 50 mM Na+

`ACGT` base Tm = 12; correction = `16.6·log10(0.05) ≈ −21.60`; result = `−9.6` °C.

### Example 2: 1 M Na+

At `[Na+] = 1000` mM the correction is `16.6·log10(1) = 0`, so the result is the rounded base Tm = `12.0` °C.

## See Also

- [primer_melting_temperature](primer_melting_temperature.md), [analyze_oligo](analyze_oligo.md)
