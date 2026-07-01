# predict_low_complexity_seg

SEG low-complexity regions in a protein.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `predict_low_complexity_seg` |
| **Method ID** | `DisorderPredictor.PredictLowComplexityRegions` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Finds **low-complexity protein regions** with the SEG algorithm (Wootton & Federhen
1993/1996): a sliding-window Shannon entropy with a K1 trigger window and K2 extension.
Each region is classified by its dominant residue type (e.g. `Q-rich`). This is the
DisorderPredictor SEG variant (see also `find_protein_low_complexity_regions`).

## Core Documentation Reference

- Source: [DisorderPredictor.cs#L682](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/DisorderPredictor.cs#L682)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | Protein sequence (min length 1) |
| `triggerWindow` | integer | No | Trigger window length (default 12) |
| `triggerThreshold` | number | No | K1 trigger entropy threshold in bits (default 2.2) |
| `extensionThreshold` | number | No | K2 extension entropy threshold in bits (default 2.5) |
| `minLength` | integer | No | Minimum reported region length (default 1) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `items` | array | `{ start, end, type }` |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |

## Examples

### Example 1: Poly-Q region

**User Prompt:**
> Find low-complexity regions in a 26-residue poly-glutamine.

**Expected Tool Call:**
```json
{
  "tool": "predict_low_complexity_seg",
  "arguments": { "sequence": "QQQQQQQQQQQQQQQQQQQQQQQQQQ" }
}
```

**Response:**
```json
{ "items": [ { "start": 0, "end": 25, "type": "Q-rich" } ] }
```

### Example 2: High complexity

**User Prompt:**
> Low-complexity regions in a sequence with all 20 residues repeated?

**Expected Tool Call:**
```json
{
  "tool": "predict_low_complexity_seg",
  "arguments": { "sequence": "ACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWY" }
}
```

**Response:**
```json
{ "items": [] }
```

## Performance

- **Time Complexity:** O(n · window).
- **Space Complexity:** O(number of regions).

## See Also

- [find_protein_low_complexity_regions](find_protein_low_complexity_regions.md)
- [predict_disorder](predict_disorder.md)
