# au_n

Compute auN (area under the Nx curve).

## Overview

| Property | Value |
|----------|-------|
| **Server** | Chromosome |
| **Tool Name** | `au_n` |
| **Method ID** | `GenomeAssemblyAnalyzer.CalculateAuN` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

`auN` is a length-weighted contiguity metric that, unlike N50, does not depend on an arbitrary
threshold: `auN = Σ(lᵢ²) / Σ(lᵢ)`. Returns 0 for empty input.

## Core Documentation Reference

- Source: [GenomeAssemblyAnalyzer.cs#L321](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/GenomeAssemblyAnalyzer.cs#L321)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `lengths` | integer[] | Yes | Sequence lengths. |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `auN` | number | `Σ(l²)/Σ(l)`. |

## Errors

| Code | Message |
|------|---------|
| 1001 | Lengths cannot be null |

## Example

`[100,90,80,70,60,50,40,30,20,10]` → `38500 / 550 = 70.0`.

## References

- [GenomeAssemblyAnalyzer.CalculateAuN](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/GenomeAssemblyAnalyzer.cs#L321)
