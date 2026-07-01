# repeat_content

Compute repeat content from repeat annotations.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Chromosome |
| **Tool Name** | `repeat_content` |
| **Method ID** | `GenomeAssemblyAnalyzer.CalculateRepeatContent` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Sums inclusive `(end - start + 1)` lengths across the repeat annotations, computes
`total × 100 / genomeLength`, and groups lengths by `repeatClass`.

## Core Documentation Reference

- Source: [GenomeAssemblyAnalyzer.cs#L816](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/GenomeAssemblyAnalyzer.cs#L816)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `repeats` | array | Yes | Repeat annotations `{ sequenceId, start, end, repeatClass, repeatFamily, divergencePercent, strand }`. |
| `genomeLength` | integer | Yes | Genome length in bp (> 0). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `totalRepeatLength` | integer | Total repeat length. |
| `repeatPercentage` | number | `total × 100 / genomeLength`. |
| `repeatClassLengths` | object | Map of repeat class → summed length. |

## Errors

| Code | Message |
|------|---------|
| 1001 | Repeats cannot be null |
| 1002 | Genome length must be positive |

## Example

`LINE` (0–99, 100 bp) + `SINE` (100–149, 50 bp), genome 1000 →
`totalRepeatLength 150, repeatPercentage 15.0, { LINE: 100, SINE: 50 }`.

## References

- [GenomeAssemblyAnalyzer.CalculateRepeatContent](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/GenomeAssemblyAnalyzer.cs#L816)
