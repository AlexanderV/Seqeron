# longest_dinucleotide_repeat

Longest dinucleotide tandem-repeat unit count.

## Overview

| Property | Value |
|----------|-------|
| **Server** | MolTools |
| **Tool Name** | `longest_dinucleotide_repeat` |
| **Method ID** | `PrimerDesigner.FindLongestDinucleotideRepeat` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Finds the longest tandem repeat of a 2-nt unit (e.g. `ATATAT` = 3 units of `AT`) and returns the number of repeat units, case-insensitive. Sequences shorter than 4 nt return 0.

## Core Documentation Reference

- Source: [PrimerDesigner.cs#L273](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/PrimerDesigner.cs#L273)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | Nucleotide sequence (non-empty). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `repeats` | integer | Longest dinucleotide repeat unit count. |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |

## Examples

### Example 1: `ATATAT` → `3` (AT × 3).

### Example 2: `AT` → `0` (shorter than 4 nt).

## See Also

- [longest_homopolymer](longest_homopolymer.md)
