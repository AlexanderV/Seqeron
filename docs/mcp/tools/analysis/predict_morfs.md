# predict_morfs

Predict Molecular Recognition Features (MoRFs) within IDRs.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `predict_morfs` |
| **Method ID** | `DisorderPredictor.PredictMoRFs` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Predicts **Molecular Recognition Features** (MoRFs) — short interaction-prone segments
within intrinsically disordered regions — by hydropathy enrichment (a heuristic
inspired by Mohan et al. 2006). Each candidate reports its span and a score. MoRFs are
constrained to `[minLength, maxLength]` residues.

> Note: this unit is gated in the library's Strict mode; the MCP server runs the
> Genomics Permissive best-effort branch to return real predictions.

## Core Documentation Reference

- Source: [DisorderPredictor.cs#L800](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/DisorderPredictor.cs#L800)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | Protein sequence (min length 1) |
| `minLength` | integer | No | Minimum MoRF length (default 10) |
| `maxLength` | integer | No | Maximum MoRF length (default 25) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `items` | array | MoRFs `{ start, end, score }` |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |

## Examples

### Example 1: Hydrophobic dip in disorder

**User Prompt:**
> Predict MoRFs in a poly-P sequence with a central hydrophobic poly-L stretch.

**Expected Tool Call:**
```json
{
  "tool": "predict_morfs",
  "arguments": { "sequence": "PPPP…(25) + LLLL…(30) + PPPP…(25)", "minLength": 10, "maxLength": 25 }
}
```

**Response:**
```json
{ "items": [ { "start": 29, "end": 50, "score": 0.275934 } ] }
```
The hydrophobic L stretch embedded in disorder is a MoRF candidate.

### Example 2: Ordered poly-L

**User Prompt:**
> Predict MoRFs in a 40-residue poly-leucine.

**Expected Tool Call:**
```json
{
  "tool": "predict_morfs",
  "arguments": { "sequence": "LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL" }
}
```

**Response:**
```json
{ "items": [] }
```
A fully ordered sequence has no IDR and therefore no MoRF.

## Performance

- **Time Complexity:** O(n · maxLength).
- **Space Complexity:** O(number of MoRFs).

## See Also

- [predict_disorder](predict_disorder.md)
- [predict_low_complexity_seg](predict_low_complexity_seg.md)
