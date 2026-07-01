# estimate_cell_divisions_from_telomere_length

Estimate the number of cell divisions from current telomere length.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Chromosome |
| **Tool Name** | `estimate_cell_divisions_from_telomere_length` |
| **Method ID** | `ChromosomeAnalyzer.EstimateCellDivisionsFromTelomereLength` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Assumes telomeres shorten by a fixed amount per division:
`divisions = max(0, (birthLength - currentLength) / lossPerDivision)`.

## Core Documentation Reference

- Source: [ChromosomeAnalyzer.cs#L1675](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/ChromosomeAnalyzer.cs#L1675)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `currentLength` | integer | Yes | Current telomere length in bp (≥ 0). |
| `birthLength` | integer | No | Telomere length at birth (default 15000, > 0). |
| `lossPerDivision` | integer | No | Loss per division in bp (default 50, > 0). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `cellDivisions` | number | Estimated divisions (≥ 0). |

## Errors

| Code | Message |
|------|---------|
| 1001 | Current length cannot be negative |
| 1002 | Birth length must be positive |
| 1003 | Loss per division must be positive |

## Example

`currentLength = 10000` (defaults) → `(15000 - 10000) / 50 = 100`.

## References

- [ChromosomeAnalyzer.EstimateCellDivisionsFromTelomereLength](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/ChromosomeAnalyzer.cs#L1675)
