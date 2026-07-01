# longest_orfs_per_frame

Return the longest ORF in each of the (up to) six reading frames.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Annotation |
| **Tool Name** | `longest_orfs_per_frame` |
| **Method ID** | `GenomeAnnotator.FindLongestOrfsPerFrame` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Runs [`find_orfs`](find_orfs.md) with `minLength=1` and `requireStartCodon=true`, then keeps the
single longest ORF (by translated protein length) per reading frame. Forward frames are keyed
`1..3`; when `searchBothStrands` is true, reverse-strand frames `-1..-3` are added, so the result
lists three or six frame slots. Frames that contain no start/stop ORF surface as an **empty ORF
sentinel** (all-zero `start`/`end`/`frame` and empty `sequence`/`proteinSequence`), not a positive
hit — this is the natural default of the underlying value-type lookup.

## Core Documentation Reference

- Source: [GenomeAnnotator.cs#L236](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/GenomeAnnotator.cs#L236)
- Algorithm: [ORF_Detection.md](../../../algorithms/Annotation/ORF_Detection.md)

## Input Schema

| Parameter | Type | Required | Default | Description |
|-----------|------|----------|---------|-------------|
| `dnaSequence` | string | Yes | — | DNA sequence to search |
| `searchBothStrands` | boolean | No | true | Whether to also search the reverse complement (frames -1..-3) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `frames` | array | `{ frame, orf }` per reading frame; `orf` is the longest ORF or an empty sentinel ORF |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |

## Examples

### Example 1: Single forward ORF (forward strand only)

`ATGAAAAAAAAATAA`, `searchBothStrands=false` → frame 1 carries the `MKKK*` ORF (positions 0..15);
frames 2 and 3 carry the empty sentinel ORF.

**Response:**
```json
{ "frames": [ { "frame": 1, "orf": { "start": 0, "end": 15, "frame": 1, "isReverseComplement": false, "sequence": "ATGAAAAAAAAATAA", "proteinSequence": "MKKK*" } }, { "frame": 2, "orf": { "start": 0, "end": 0, "sequence": null, "proteinSequence": null } }, { "frame": 3, "orf": { "start": 0, "end": 0, "sequence": null, "proteinSequence": null } } ] }
```

### Example 2: Both strands return six frame slots

The same sequence with `searchBothStrands=true` returns six slots in order `1, 2, 3, -1, -2, -3`.

## Performance

- **Time Complexity:** O(n) across the six frames
- **Space Complexity:** O(k)

## See Also

- [find_orfs](find_orfs.md) - All ORFs across all frames
- [predict_genes](predict_genes.md) - ORF-based gene prediction
