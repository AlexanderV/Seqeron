# nx_statistics

Compute Nx and Lx for a single threshold.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Chromosome |
| **Tool Name** | `nx_statistics` |
| **Method ID** | `GenomeAssemblyAnalyzer.CalculateNx` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Walks descending-sorted lengths accumulating coverage until it first reaches `threshold%` of the
total. `Nx` is the length of the sequence that crossed the threshold; `Lx` is how many sequences
that took (inclusive "at least x%" test, matching QUAST). Miller, Koren & Sutton (2010).

## Core Documentation Reference

- Source: [GenomeAssemblyAnalyzer.cs#L237](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/GenomeAssemblyAnalyzer.cs#L237)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sortedLengths` | integer[] | Yes | Lengths sorted descending. |
| `totalLength` | integer | Yes | Sum of all lengths (≥ 0). |
| `threshold` | integer | Yes | Percent in `[0, 100]` (50 for N50). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `threshold` | integer | Echoed threshold. |
| `nx` | integer | Nx length. |
| `lx` | integer | Lx count. |
| `cumulativeLength` | integer | Cumulative length when the threshold was reached. |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sorted lengths cannot be null |
| 1002 | Total length cannot be negative |
| 1003 | Threshold must be in [0, 100] |

## Example

`[100,90,80,70,60,50,40,30,20,10]`, total 550, threshold 50 → cumulative 270 at 80, then 340 at 70
crosses 275 (50% of 550) → `nx = 70, lx = 4, cumulativeLength = 340`.

## References

- [GenomeAssemblyAnalyzer.CalculateNx](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/GenomeAssemblyAnalyzer.cs#L237)
- [nx_curve](nx_curve.md), [au_n](au_n.md)
