# predict_disorder

TOP-IDP intrinsic-disorder prediction for a protein.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `predict_disorder` |
| **Method ID** | `DisorderPredictor.PredictDisorder` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Predicts **intrinsic disorder** using the TOP-IDP scale (Campen et al. 2008): a
sliding-window normalized disorder score per residue, plus contiguous disordered
regions (IDRs) with confidence and subtype classification, and the overall disorder
content. A residue is disordered when its normalized score exceeds `disorderThreshold`;
a region must be at least `minRegionLength` residues.

> Note: the per-residue/per-region confidence is gated in the library's Strict mode; the
> MCP server runs the Genomics Permissive best-effort branch to return real predictions.

## Core Documentation Reference

- Source: [DisorderPredictor.cs#L238](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/DisorderPredictor.cs#L238)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | Protein sequence (min length 1) |
| `windowSize` | integer | No | Sliding window size (default 21) |
| `disorderThreshold` | number | No | Threshold on TOP-IDP normalized score (default 0.542) |
| `minRegionLength` | integer | No | Minimum reported region length (default 5) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `sequence` | string | Echoed sequence |
| `residuePredictions` | array | Per-residue `{ position, residue, disorderScore, isDisordered }` |
| `disorderedRegions` | array | IDRs `{ start, end, meanScore, confidence, regionType }` |
| `overallDisorderContent` | number | Fraction of disordered residues |
| `meanDisorderScore` | number | Mean per-residue disorder score |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |

## Examples

### Example 1: Fully disordered (poly-P)

**User Prompt:**
> Predict disorder for a 30-residue poly-proline.

**Expected Tool Call:**
```json
{
  "tool": "predict_disorder",
  "arguments": { "sequence": "PPPPPPPPPPPPPPPPPPPPPPPPPPPPPP", "minRegionLength": 5 }
}
```

**Response:**
```json
{ "disorderedRegions": [ { "start": 0, "end": 29 } ], "overallDisorderContent": 1.0 }
```
Proline is the most disorder-promoting residue (normalized score 1.0).

### Example 2: Fully ordered (poly-W)

**User Prompt:**
> Predict disorder for a 30-residue poly-tryptophan.

**Expected Tool Call:**
```json
{
  "tool": "predict_disorder",
  "arguments": { "sequence": "WWWWWWWWWWWWWWWWWWWWWWWWWWWWWW", "minRegionLength": 5 }
}
```

**Response:**
```json
{ "disorderedRegions": [], "overallDisorderContent": 0.0 }
```

## Performance

- **Time Complexity:** O(n · windowSize).
- **Space Complexity:** O(n).

## See Also

- [predict_morfs](predict_morfs.md)
- [predict_low_complexity_seg](predict_low_complexity_seg.md)
- [disorder_propensity](disorder_propensity.md)
