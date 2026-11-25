# repeat_content

Compute RepeatMasker-style repeat content from repeat annotations.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Chromosome |
| **Tool Name** | `repeat_content` |
| **Method ID** | `GenomeAssemblyAnalyzer.CalculateRepeatContent` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Counts every masked base once ("bases masked", RepeatMasker `ProcessRepeats`): per sequence the annotations
are swept in start order and each contributes only the bases beyond the furthest end seen so far, credited to
its `repeatClass`; overlapping or contained annotations are therefore not double counted and the class lengths
sum to the total. `repeatPercentage = total × 100 / genomeLength`. Coordinates are 0-based inclusive; an
annotation with `end < start` or `start < 0`, or more masked bases than `genomeLength`, is rejected.

## Core Documentation Reference

- Source: [GenomeAssemblyAnalyzer.cs#L1418](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/GenomeAssemblyAnalyzer.cs#L1418)

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
| `repeatClassLengths` | object | Map of repeat class → bases credited to the class (sums to the total). |

## Errors

| Code | Message |
|------|---------|
| 1001 | Repeats cannot be null |
| 1002 | Genome length must be positive |

## Example

`LINE` (0–99, 100 bp) + `SINE` (100–149, 50 bp), genome 1000 →
`totalRepeatLength 150, repeatPercentage 15.0, { LINE: 100, SINE: 50 }`; `LINE` 0–99 + `SINE` 50–149 → also 150 (overlap counted once, LINE 100 / SINE 50).

## References

- RepeatMasker `ProcessRepeats` (`generateTableOutput`, "Define length of coverage based on possibly overlapping fragments").
- [GenomeAssemblyAnalyzer.CalculateRepeatContent](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/GenomeAssemblyAnalyzer.cs#L1418)
