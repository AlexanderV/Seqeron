# find_tandem_repeats

Find approximate tandem repeats with Tandem Repeats Finder (Benson 1999).

## Overview

| Property | Value |
|----------|-------|
| **Server** | Chromosome |
| **Tool Name** | `find_tandem_repeats` |
| **Method ID** | `GenomeAssemblyAnalyzer.FindTandemRepeats` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Runs the Tandem Repeats Finder model (`RepeatFinder.FindApproximateTandemRepeats`; match 2, mismatch −7,
indel −7, PM 80, PI 10; wraparound alignment, consensus, TRF statistics). Reported per repeat: `unit` =
TRF consensus pattern, `copies` = TRF copy number (aligned consensus columns / consensus size, so
fractional), `purity` = TRF "percent matches between adjacent copies" as a fraction. A repeat is kept when
its alignment score ≥ `minScore` (TRF Minscore, default 50), its period lies in
`minUnitLength`–`maxUnitLength` (TRF MaxPeriod ≤ 2000) and `copies` ≥ `minCopies`. Non-ACGT symbols never
match. `start`/`end` are 0-based inclusive.

## Core Documentation Reference

- Source: [GenomeAssemblyAnalyzer.cs#L1367](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/GenomeAssemblyAnalyzer.cs#L1367)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequences` | array | Yes | Sequences `{ id, sequence }`. |
| `minUnitLength` | integer | No | Minimum unit length (default 2, > 0). |
| `maxUnitLength` | integer | No | Maximum unit length (default 50, ≥ min). |
| `minCopies` | integer | No | Minimum copy number (default 3, > 0). |
| `minScore` | integer | No | TRF minimum alignment score (default 50, > 0). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `items[]` | array | `{ sequenceId, start, end, unit, copies, purity }`. |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequences cannot be null |
| 1002 | Minimum unit length must be positive |
| 1003 | Maximum unit length must be >= minimum unit length |
| 1004 | Minimum copies must be positive |
| 1005 | Minimum score must be positive |

## Example

`"CAG"` × 30 (90 bp) → `{ unit: "CAG", copies: 30, start: 0, end: 89, purity: 1.0 }`; `"ATG"`×7 · `ATC` · `"ATG"`×7 → `copies 15`, `purity 40/42` (score 81).

## References

- Benson G. (1999) Tandem repeats finder, Nucleic Acids Res 27:573–580; TRF 4.10.0 README.
- [GenomeAssemblyAnalyzer.FindTandemRepeats](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/GenomeAssemblyAnalyzer.cs#L1367)
