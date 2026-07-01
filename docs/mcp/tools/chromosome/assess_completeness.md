# assess_completeness

Assess assembly completeness with marker genes (BUSCO-like).

## Overview

| Property | Value |
|----------|-------|
| **Server** | Chromosome |
| **Tool Name** | `assess_completeness` |
| **Method ID** | `GenomeAssemblyAnalyzer.AssessCompleteness` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Aligns each marker gene against the assembly (k-mer based) and classifies it as complete
(single-copy or duplicated), fragmented, or missing using the `identityThreshold` and
`coverageThreshold`. Reports the counts plus completeness and duplication percentages.

## Core Documentation Reference

- Source: [GenomeAssemblyAnalyzer.cs#L520](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/GenomeAssemblyAnalyzer.cs#L520)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `assembly` | array | Yes | Assembled sequences `{ id, sequence }`. |
| `markerGenes` | array | Yes | Marker genes `{ id, sequence }`. |
| `identityThreshold` | number | No | Minimum identity in `[0, 1]` (default 0.9). |
| `coverageThreshold` | number | No | Minimum coverage in `[0, 1]` (default 0.9). |

## Output Schema

Fields: `totalGenes`, `complete`, `completeSingleCopy`, `completeDuplicated`, `fragmented`, `missing`,
`completenessPercent`, `duplicationPercent`.

## Errors

| Code | Message |
|------|---------|
| 1001 | Assembly cannot be null |
| 1002 | Marker genes cannot be null |
| 1003 | Identity threshold must be in [0, 1] |
| 1004 | Coverage threshold must be in [0, 1] |

## Example

One marker gene fully contained in the assembly → `totalGenes 1, complete 1, completeSingleCopy 1,
missing 0, completenessPercent 100`.

## References

- [GenomeAssemblyAnalyzer.AssessCompleteness](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/GenomeAssemblyAnalyzer.cs#L520)
