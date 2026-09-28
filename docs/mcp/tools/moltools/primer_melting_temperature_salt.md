# primer_melting_temperature_salt

Salt-adjusted primer melting temperature (OligoCalc).

## Overview

| Property | Value |
|----------|-------|
| **Server** | MolTools |
| **Tool Name** | `primer_melting_temperature_salt` |
| **Method ID** | `PrimerDesigner.CalculateMeltingTemperatureWithSalt` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

OligoCalc "Salt Adjusted" Tm (Kibbe 2007, NAR 35:W43), the [Na+]-aware counterpart of the basic Tm in [primer_melting_temperature](primer_melting_temperature.md) (whose formulas already assume 50 mM Na+):

- fewer than 14 valid bases: `Tm = 2·(A+T) + 4·(G+C) − 16.6·log10(0.050) + 16.6·log10([Na+])`
- 14 or more: `Tm = 100.5 + 41·(G+C)/N − 820/N + 16.6·log10([Na+])`

[Na+] in mol/L (the argument is in mM). Rounded to one decimal place. OligoCalc cross-check: the 39-mer `GAGCAGGATCCCTATAGAGTGACAAAAGGATCTTGGTCC` at 50 mM gives 78 °C (basic 67.6 °C).

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
| `tm` | number | Salt-adjusted Tm (°C, 1 decimal). |

## Errors

| Code | Message |
|------|---------|
| 1001 | Primer cannot be null or empty |
| 1002 | Na+ concentration must be positive |

## Examples

### Example 1: 50 mM Na+

`ACGT` (< 14 nt): the Wallace Tm is defined at 50 mM, so the result is `12.0` °C. `ACGTACGTACGTACGTACGT`: `100.5 + 20.5 − 41 − 21.6 = 58.4` °C.

### Example 2: 1 M Na+

`ACGT` at `[Na+] = 1000` mM: `12 + 16.6·log10(1/0.050) = 33.6` °C.

## See Also

- [primer_melting_temperature](primer_melting_temperature.md), [analyze_oligo](analyze_oligo.md)
