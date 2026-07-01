# detect_pseudoknots

Detect pseudoknots as crossing RNA base pairs.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `detect_pseudoknots` |
| **Method ID** | `RnaSecondaryStructure.DetectPseudoknots` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Identifies pseudoknots as **crossing** base pairs. For two pairs normalized to
`(open < close)` and ordered so the earlier-opening pair is `(i, j)`, they cross — and form a
pseudoknot — exactly when `i < k < j < l`. Nested pairs (`i < k < l < j`) and disjoint pairs
(`j < k`) are not pseudoknots. Each crossing pair-of-pairs is reported once. Each input base
pair carries a `type` of `WatsonCrick`, `Wobble`, or `NonCanonical`.

## Core Documentation Reference

- Source: [RnaSecondaryStructure.cs#L1934](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/RnaSecondaryStructure.cs#L1934)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `basePairs` | object[] | Yes | Base pairs: `{ position1, position2, base1, base2, type }` |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `items` | array | One pseudoknot per crossing pair-of-pairs: `start1,end1,start2,end2,crossingPairs` |

## Errors

| Code | Message |
|------|---------|
| 1002 | Unknown base pair type: … |

## Examples

### Example 1: Crossing pairs

**Input:** `{ "basePairs": [ {"position1":0,"position2":5,...}, {"position1":3,"position2":8,...} ] }`
→ 0 < 3 < 5 < 8 crosses → **Response:** `{ "items": [ { "start1":0, "end1":5, "start2":3, "end2":8 } ] }`

### Example 2: Nested pairs (no pseudoknot)

**Input:** pairs `(0,8)` and `(3,5)` → nested, not crossing → **Response:** `{ "items": [] }`

## Performance

- **Time Complexity:** O(p²) for p base pairs. **Space Complexity:** O(1) per report.

## See Also

- [predict_rna_structure](predict_rna_structure.md)
- [parse_dot_bracket](parse_dot_bracket.md)
