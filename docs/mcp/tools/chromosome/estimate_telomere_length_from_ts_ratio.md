# estimate_telomere_length_from_ts_ratio

Estimate telomere length in bp from a qPCR T/S ratio.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Chromosome |
| **Tool Name** | `estimate_telomere_length_from_ts_ratio` |
| **Method ID** | `ChromosomeAnalyzer.EstimateTelomereLengthFromTSRatio` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

The Telomere/Single-copy gene (T/S) ratio from qPCR is proportional to telomere length. This tool
scales a known reference: `telomereLength = referenceLength * tsRatio / referenceRatio`.

## Core Documentation Reference

- Source: [ChromosomeAnalyzer.cs#L488](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/ChromosomeAnalyzer.cs#L488)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `tsRatio` | number | Yes | Sample T/S ratio. |
| `referenceRatio` | number | No | Reference T/S ratio (default 1.0, > 0). |
| `referenceLength` | number | No | Telomere length (bp) at the reference ratio (default 7000, > 0). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `telomereLength` | number | Estimated telomere length in bp. |

## Errors

| Code | Message |
|------|---------|
| 1001 | Reference ratio must be positive |
| 1002 | Reference length must be positive |

## Example

`tsRatio = 2.0` (defaults) → `7000 * 2 / 1 = 14000`.

## References

- [ChromosomeAnalyzer.EstimateTelomereLengthFromTSRatio](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/ChromosomeAnalyzer.cs#L488)
