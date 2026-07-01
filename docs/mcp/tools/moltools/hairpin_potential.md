# hairpin_potential

Detect whether a sequence can fold into a hairpin.

## Overview

| Property | Value |
|----------|-------|
| **Server** | MolTools |
| **Tool Name** | `hairpin_potential` |
| **Method ID** | `PrimerDesigner.HasHairpinPotential` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Returns true if the sequence contains a self-complementary stem of at least `min_stem_length` separated by a loop of at least `min_loop_length`. A short-sequence O(n²) scan is used below 100 bp; a suffix-tree scan is used at or above 100 bp. Sequences shorter than `2·min_stem_length + min_loop_length` cannot form a hairpin and return false.

## Core Documentation Reference

- Source: [PrimerDesigner.cs#L307](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/PrimerDesigner.cs#L307)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | Nucleotide sequence (non-empty). |
| `min_stem_length` | integer | No | Minimum stem length (default 4, positive). |
| `min_loop_length` | integer | No | Minimum loop length (default 3, non-negative). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `hasHairpin` | boolean | True if a hairpin can form. |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1002 | Minimum stem length must be positive |
| 1003 | Minimum loop length cannot be negative |

## Examples

### Example 1: `GGGGAAACCCC` → `true` (stem `GGGG` / loop `AAA` / stem `CCCC`).

### Example 2: `AAAAAAAAAAA` → `false` (no complementary stem).

## See Also

- [three_prime_stability](three_prime_stability.md), [evaluate_primer](evaluate_primer.md)
