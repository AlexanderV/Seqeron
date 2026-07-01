# find_protein_low_complexity_regions

Low-complexity regions in a protein via the SEG algorithm.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `find_protein_low_complexity_regions` |
| **Method ID** | `ProteinMotifFinder.FindLowComplexityRegions` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Finds **low-complexity regions** in a protein using the SEG algorithm (Wootton &
Federhen 1993): a sliding-window Shannon entropy (bits/residue) with a two-pass
trigger (`K1`) / extension (`K2`) scheme. Windows below `K1` seed a region that is
then extended over neighbouring windows below `K2`. Each region reports start, end and
its minimum complexity. Compositionally biased tracts (e.g. poly-Q) are the typical
hits.

## Core Documentation Reference

- Source: [ProteinMotifFinder.cs#L1132](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/ProteinMotifFinder.cs#L1132)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `proteinSequence` | string | Yes | Protein sequence (min length 1) |
| `windowSize` | integer | No | Sliding window size W (default 12, ≥ 1) |
| `triggerComplexity` | number | No | Trigger complexity K1 in bits/residue (default 2.2) |
| `extensionComplexity` | number | No | Extension complexity K2 in bits/residue (default 2.5) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `items` | array | Low-complexity regions: `{ start, end, complexity }` |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1002 | Window size must be at least 1 |

## Examples

### Example 1: Poly-Q low-complexity region

**User Prompt:**
> Find low-complexity regions in a protein with a flanked poly-Q tract.

**Expected Tool Call:**
```json
{
  "tool": "find_protein_low_complexity_regions",
  "arguments": { "proteinSequence": "ACDEFGHIKLMN + Q×20 + NMLKIHGFEDCA" }
}
```

**Response:**
```json
{ "items": [ { "start": 6, "end": 37, "complexity": 0.0 } ] }
```
The SEG two-pass extension grows the zero-entropy poly-Q core into residues 6–37.

### Example 2: All distinct residues (high complexity)

**User Prompt:**
> Low-complexity regions in "ACDEFGHIKLMNPQRSTVWY"?

**Expected Tool Call:**
```json
{
  "tool": "find_protein_low_complexity_regions",
  "arguments": { "proteinSequence": "ACDEFGHIKLMNPQRSTVWY" }
}
```

**Response:**
```json
{ "items": [] }
```
All 20 residues distinct ⇒ maximal complexity, no region.

## Performance

- **Time Complexity:** O(n · windowSize).
- **Space Complexity:** O(number of regions).

## See Also

- [predict_low_complexity_seg](predict_low_complexity_seg.md) — DisorderPredictor SEG variant
- [find_low_complexity_regions](find_low_complexity_regions.md) — DNA low-complexity regions
