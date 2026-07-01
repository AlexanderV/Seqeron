# can_pair

Whether two RNA bases can pair.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `can_pair` |
| **Method ID** | `RnaSecondaryStructure.CanPair` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Returns `true` when two RNA bases form a canonical Watson-Crick pair (A-U, G-C) or a
G-U wobble pair (Crick 1966), and `false` otherwise. Comparison is case-insensitive.

## Core Documentation Reference

- Source: [RnaSecondaryStructure.cs#L423](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/RnaSecondaryStructure.cs#L423)
- Evidence: `docs/Evidence/RNA-PAIR-001-Evidence.md`

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `base1` | string | Yes | First RNA base (single character) |
| `base2` | string | Yes | Second RNA base (single character) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `result` | boolean | True if the two bases can form a Watson-Crick or wobble pair |

## Errors

| Code | Message |
|------|---------|
| 1001 | Expected a single character (length-1 string) |

## Examples

### Example 1: Watson-Crick pair

**Input:** `{ "base1": "G", "base2": "C" }` → **Response:** `{ "result": true }`

### Example 2: Wobble pair

**Input:** `{ "base1": "G", "base2": "U" }` → **Response:** `{ "result": true }`

### Example 3: Non-pair

**Input:** `{ "base1": "A", "base2": "G" }` → **Response:** `{ "result": false }`

## Performance

- **Time Complexity:** O(1). **Space Complexity:** O(1).

## See Also

- [base_pair_type](base_pair_type.md) — classify the pair (WatsonCrick/Wobble)
- [rna_complement_base](rna_complement_base.md)
