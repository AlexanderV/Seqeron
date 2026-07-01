# arm_ratio

Compute chromosome arm ratio (p/q).

## Overview

| Property | Value |
|----------|-------|
| **Server** | Chromosome |
| **Tool Name** | `arm_ratio` |
| **Method ID** | `ChromosomeAnalyzer.CalculateArmRatio` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Given the centromere position and total length, the short arm is `p = centromerePosition` and the
long arm is `q = chromosomeLength - centromerePosition`; the ratio returned is `p / q`.

## Core Documentation Reference

- Source: [ChromosomeAnalyzer.cs#L1631](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/ChromosomeAnalyzer.cs#L1631)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `centromerePosition` | integer | Yes | Centromere position in bp (`0 < pos < length`). |
| `chromosomeLength` | integer | Yes | Total chromosome length in bp (> 0). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `armRatio` | number | `p/q`. |

## Errors

| Code | Message |
|------|---------|
| 1001 | Centromere position must be positive |
| 1002 | Chromosome length must be positive |
| 1003 | Centromere position must be less than chromosome length |

## Example

`centromerePosition = 40`, `chromosomeLength = 100` → `40 / 60 ≈ 0.667`.

## References

- [ChromosomeAnalyzer.CalculateArmRatio](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/ChromosomeAnalyzer.cs#L1631)
- [classify_chromosome_by_arm_ratio](classify_chromosome_by_arm_ratio.md)
