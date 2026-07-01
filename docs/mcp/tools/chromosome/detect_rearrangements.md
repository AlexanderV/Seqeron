# detect_rearrangements

Detect chromosomal rearrangements from synteny blocks.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Chromosome |
| **Tool Name** | `detect_rearrangements` |
| **Method ID** | `ChromosomeAnalyzer.DetectRearrangements` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Examines synteny blocks (typically from `find_synteny_blocks`) for signatures of inversions,
translocations, deletions and duplications, and reports each detected event with its type, breakpoint
positions and optional size/description. A single collinear forward block yields no rearrangements.

## Core Documentation Reference

- Source: [ChromosomeAnalyzer.cs#L1426](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/ChromosomeAnalyzer.cs#L1426)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `syntenyBlocks` | array | Yes | Synteny blocks (from `find_synteny_blocks`). Empty → no rearrangements. |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `items[]` | array | `{ type, chromosome1, position1, chromosome2?, position2?, size?, description? }`. |

## Errors

| Code | Message |
|------|---------|
| 1001 | Synteny blocks cannot be null |

## Example

A single collinear forward block → `{ "items": [] }`.

## References

- [ChromosomeAnalyzer.DetectRearrangements](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/ChromosomeAnalyzer.cs#L1426)
- [find_synteny_blocks](find_synteny_blocks.md)
