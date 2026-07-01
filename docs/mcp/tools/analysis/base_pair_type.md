# base_pair_type

Classify an RNA base-pair candidate.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `base_pair_type` |
| **Method ID** | `RnaSecondaryStructure.GetBasePairType` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Classifies a pair of RNA bases as `WatsonCrick` (A-U, G-C), `Wobble` (G-U), or
returns `null` when the two bases cannot form a canonical/wobble pair. Comparison
is case-insensitive. This underpins RNA secondary-structure prediction.

## Core Documentation Reference

- Source: [RnaSecondaryStructure.cs#L433](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/RnaSecondaryStructure.cs#L433)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `base1` | string | Yes | First RNA base (single character) |
| `base2` | string | Yes | Second RNA base (single character) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `type` | string \| null | `WatsonCrick`, `Wobble`, or `null` if the bases cannot pair |

## Errors

| Code | Message |
|------|---------|
| 1001 | Expected a single character (length-1 string) |

## Examples

### Example 1: Watson-Crick pair

**User Prompt:**
> What kind of RNA pair is A and U?

**Expected Tool Call:**
```json
{
  "tool": "base_pair_type",
  "arguments": { "base1": "A", "base2": "U" }
}
```

**Response:**
```json
{ "type": "WatsonCrick" }
```

### Example 2: Wobble pair

**User Prompt:**
> Classify the RNA pair G / U.

**Expected Tool Call:**
```json
{
  "tool": "base_pair_type",
  "arguments": { "base1": "G", "base2": "U" }
}
```

**Response:**
```json
{ "type": "Wobble" }
```

### Example 3: Non-pairing bases

**Input:** `{ "base1": "A", "base2": "G" }` → **Response:** `{ "type": null }`

## Performance

- **Time Complexity:** O(1) (lookup table).
- **Space Complexity:** O(1).

## See Also

- [can_pair](can_pair.md) — boolean can-pair test
- [rna_complement_base](rna_complement_base.md) — RNA complement of a base
