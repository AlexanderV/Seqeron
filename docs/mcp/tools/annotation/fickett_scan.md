# fickett_scan

Sliding-window Fickett (1982) TESTCODE, as EMBOSS `tcode`.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Annotation |
| **Tool Name** | `fickett_scan` |
| **Method ID** | `GenomeAnnotator.ScanFickettScore` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Windows start at 0, `stepSize`, … while `start + windowSize < length` (EMBOSS loop bound); each window gets the
Fickett 1982 TESTCODE and its call (`< 0.74` NonCoding, `>= 0.95` Coding). `windowSize >= 200` (EMBOSS minimum,
default 200), `stepSize >= 1` (default 3). Coordinates are 0-based half-open (EMBOSS prints `start+1 .. end`).

## Core Documentation Reference

- Source: [GenomeAnnotator.cs#L1263](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/GenomeAnnotator.cs#L1263) — `ScanFickettScore`
- Algorithm: [Coding_Potential_Calculation.md §9](../../../algorithms/Extended_Annotation/Coding_Potential_Calculation.md)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | DNA sequence |
| `windowSize` | integer | No | Window length, >= 200 (default 200) |
| `stepSize` | integer | No | Step, >= 1 (default 3) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `windows` | array | `{ start, end, score, call }` per window |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1001 | windowSize < 200 or stepSize < 1 |

## Example

A 240-nt sequence with the defaults gives 14 windows, `[0,200)` … `[39,239)`.

## See Also

- [fickett_score](fickett_score.md)
