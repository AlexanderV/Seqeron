# identify_whole_chromosome_aneuploidy

Identify whole-chromosome aneuploidies from per-bin copy-number states.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Chromosome |
| **Tool Name** | `identify_whole_chromosome_aneuploidy` |
| **Method ID** | `ChromosomeAnalyzer.IdentifyWholeChromosomeAneuploidy` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

For each chromosome, finds the dominant per-bin copy number. If it covers at least `minFraction` of the
bins and is not 2, the chromosome is reported with an ISCN term (0 Nullisomy, 1 Monosomy, 3 Trisomy,
4 Tetrasomy, 5 Pentasomy, else "Copy number = N"). Diploid (copy number 2) is normal and not reported.

## Core Documentation Reference

- Source: [ChromosomeAnalyzer.cs#L1586](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/ChromosomeAnalyzer.cs#L1586)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `copyNumberStates` | array | Yes | Per-bin copy-number states (from `detect_aneuploidy`). |
| `minFraction` | number | No | Minimum dominant fraction in `(0, 1]` (default 0.8). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `items[]` | array | `{ chromosome, copyNumber, type }`. |

## Errors

| Code | Message |
|------|---------|
| 1001 | Copy-number states cannot be null |
| 1002 | Minimum fraction must be in (0, 1] |

## Example

Three bins all at copy number 3 → `{ chromosome: "chr1", copyNumber: 3, type: "Trisomy" }`.

## References

- [ChromosomeAnalyzer.IdentifyWholeChromosomeAneuploidy](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/ChromosomeAnalyzer.cs#L1586)
- [detect_aneuploidy](detect_aneuploidy.md)
