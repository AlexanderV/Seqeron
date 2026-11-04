# longest_homopolymer

Length of the longest run of identical consecutive nucleotides.

## Overview

| Property | Value |
|----------|-------|
| **Server** | MolTools |
| **Tool Name** | `longest_homopolymer` |
| **Method ID** | `PrimerDesigner.FindLongestHomopolymer` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Scans the sequence (case-insensitive) and returns the length of its longest homopolymer — the maximal run of the same nucleotide (e.g. `AAAA` → 4). A sequence with no adjacent repeats returns 1. This is the quantity Primer3 limits with `PRIMER_MAX_POLY_X`; as in Primer3 (`_pr_violates_poly_x`), N is a worst-case wildcard (`ANA` → 3, `GNGNG` → 5, `ANGNG` → 4).

## Core Documentation Reference

- Source: [PrimerDesigner.cs#L1673](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/PrimerDesigner.cs#L1673)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | Nucleotide sequence (non-empty). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `length` | integer | Longest homopolymer run length. |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |

## Examples

### Example 1: `AAATTGC` → `3` (the `AAA` run).

### Example 2: `ACGT` → `1` (no adjacent repeats).

## See Also

- [longest_dinucleotide_repeat](longest_dinucleotide_repeat.md)
