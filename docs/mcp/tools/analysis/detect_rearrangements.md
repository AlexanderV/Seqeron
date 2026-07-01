# detect_rearrangements

Detect genome rearrangements as breakpoints of the signed gene-order permutation.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `detect_rearrangements` |
| **Method ID** | `ComparativeGenomics.DetectRearrangements` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Detects rearrangements between two gene orders as **breakpoints** of the signed gene-order
permutation (Bafna & Pevzner 1998). Orthologous markers (via `orthologMap`) are read in
genome-1 order, relabelled to genome-2 rank with a sign for relative strand, and extended with
sentinels. Every consecutive pair that is not an identity adjacency (`y ≠ x + 1`) is a
breakpoint, reported as one event and classified `Inversion` (sign flip across the boundary) or
`Transposition` (orientation-preserving discontinuity). Identical gene order yields no events.

## Core Documentation Reference

- Source: [ComparativeGenomics.cs#L582](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/ComparativeGenomics.cs#L582)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `genome1Genes` | object[] | Yes | Genome-1 genes `{ id, genomeId, start, end, strand, sequence? }` |
| `genome2Genes` | object[] | Yes | Genome-2 genes |
| `orthologMap` | object | Yes | Map genome-1 gene id → genome-2 gene id |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `items` | array | Rearrangement events: `type, genomeId, position, length, targetPosition` |

## Errors

| Code | Message |
|------|---------|
| 1001 | Genome 1 must contain at least one gene |
| 1001 | Genome 2 must contain at least one gene |
| 1001 | Ortholog map cannot be null |

## Examples

### Example 1: Identical gene order

Three genes A,B,C in the same order and strand in both genomes → **Response:** `{ "items": [] }`.

### Example 2: Local inversion

Gene B is on the opposite strand in genome 2 → signed permutation `[1, -2, 3]` → two
breakpoints at boundaries `1->-2` and `-2->3`, both classified `Inversion`:

```json
{ "items": [
  { "type": "Inversion", "genomeId": "G1", "position": 20 },
  { "type": "Inversion", "genomeId": "G1", "position": 40 }
] }
```

## Performance

- **Time Complexity:** O(n log n) over n orthologous markers. **Space Complexity:** O(n).

## See Also

- [find_syntenic_blocks](find_syntenic_blocks.md)
- [reversal_distance](reversal_distance.md)
- [compare_genomes](compare_genomes.md)
