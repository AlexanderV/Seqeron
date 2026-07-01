# detect_ploidy

Detect ploidy level from normalized read-depth values.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Chromosome |
| **Tool Name** | `detect_ploidy` |
| **Method ID** | `ChromosomeAnalyzer.DetectPloidy` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Takes the median of the normalized depths, forms `ratio = median / expectedDiploidDepth`, and
reports `ploidy = round(ratio * 2)` clamped to `[1, 8]`. Confidence is `1 - 2 * |fractional part|`
(clamped at 0), so an exactly integer ploidy scores 1.0. An empty depth list yields ploidy 2 with
confidence 0.

## Core Documentation Reference

- Source: [ChromosomeAnalyzer.cs#L360](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/ChromosomeAnalyzer.cs#L360)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `normalizedDepths` | number[] | Yes | Normalized read depths. |
| `expectedDiploidDepth` | number | No | Depth expected for a diploid locus (default 1.0, > 0). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `ploidyLevel` | integer | Detected ploidy, `[1, 8]`. |
| `confidence` | number | `1 - 2*frac`, `[0, 1]`. |

## Errors

| Code | Message |
|------|---------|
| 1001 | Normalized depths cannot be null |
| 1002 | Expected diploid depth must be positive |

## Example

Depths `[1.0, 1.0, 1.0]` with `expectedDiploidDepth = 1.0` → `{ "ploidyLevel": 2, "confidence": 1.0 }`.
Depths `[2.0, 2.0, 2.0]` → `{ "ploidyLevel": 4, "confidence": 1.0 }`.

## References

- [ChromosomeAnalyzer.DetectPloidy](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/ChromosomeAnalyzer.cs#L360)
