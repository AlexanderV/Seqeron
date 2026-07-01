# classify_chromosome_by_arm_ratio

Classify a chromosome from its arm ratio per Levan et al. (1964).

## Overview

| Property | Value |
|----------|-------|
| **Server** | Chromosome |
| **Tool Name** | `classify_chromosome_by_arm_ratio` |
| **Method ID** | `ChromosomeAnalyzer.ClassifyChromosomeByArmRatio` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Normalizes the input to `r = long arm / short arm ≥ 1` (accepting either p/q or q/p) and classifies:

| r (long/short) | Type |
|----------------|------|
| ≤ 0 (single arm) | Telocentric |
| ≤ 1.7 | Metacentric |
| ≤ 3.0 | Submetacentric |
| < 7.0 | Subtelocentric |
| ≥ 7.0 | Acrocentric |

## Core Documentation Reference

- Source: [ChromosomeAnalyzer.cs#L1654](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/ChromosomeAnalyzer.cs#L1654)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `armRatio` | number | Yes | Arm ratio (p/q or q/p). Must be finite. |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `classification` | string | Levan (1964) category. |

## Errors

| Code | Message |
|------|---------|
| 1001 | Arm ratio must be a finite number |

## Example

`armRatio = 10.0` → `Acrocentric`. `armRatio = 0.4` normalizes to `r = 2.5` → `Submetacentric`.

## References

- Levan A, Fredga K, Sandberg AA (1964). *Hereditas* 52(2):201–220.
- [ChromosomeAnalyzer.ClassifyChromosomeByArmRatio](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/ChromosomeAnalyzer.cs#L1654)
