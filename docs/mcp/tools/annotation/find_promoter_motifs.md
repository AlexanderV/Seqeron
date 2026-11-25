# find_promoter_motifs

Find -10 (Pribnow/TATA) and -35 box bacterial promoter motifs.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Annotation |
| **Tool Name** | `find_promoter_motifs` |
| **Method ID** | `GenomeAnnotator.FindPromoterMotifs` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Reports the **-35** and **-10** (Pribnow) boxes of bacterial σ70 promoters. Delegates to
`GenomeAnnotator.FindPromoterMotifs`, which takes every forward-strand pair from the canonical
`MotifFinder.FindSigma70Promoters` (Harley & Reynolds 1987): a -35 hexamer within 2 mismatches of `TTGACA`
followed, after a 15–21 bp spacer, by a -10 hexamer within 2 mismatches of `TATAAT`. Each box of a pair is
reported once (0-based `position`, `type` `-35 box` / `-10 box`, the hexamer `sequence`), ordered by position.
`score` = 1 − mismatches/6 (consensus identity, not promoter strength). An isolated box with no partner is not
reported. For pairs with strand and spacer use `find_sigma70_promoters`; for a free-energy promoter score use
`predict_sigma70_promoters` (Promoter Calculator v1.0).

## Core Documentation Reference

- Source: [GenomeAnnotator.cs#L794](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/GenomeAnnotator.cs#L794)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `dnaSequence` | string | Yes | DNA sequence to scan |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `motifs` | array | `{ position, type, sequence, score }` per promoter-motif hit |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |

## Examples

### Example 1: consensus σ70 promoter

`GGTTGACACCCCCCCCCCCCCCCCCTATAATGG` (TTGACA, 17-bp spacer, TATAAT)

**Response:**
```json
{ "motifs": [
  { "position": 2, "type": "-35 box", "sequence": "TTGACA", "score": 1.0 },
  { "position": 25, "type": "-10 box", "sequence": "TATAAT", "score": 1.0 } ] }
```

### Example 2: unpaired box

`GGGGGTTGACAGGGGG` → `{ "motifs": [] }` (no -10 box at a 15–21 bp spacer).

## Performance

- **Time Complexity:** O(n · 7) (spacer window 15–21)
- **Space Complexity:** O(k)

## See Also

- [find_sigma70_promoters](../analysis/find_sigma70_promoters.md) - σ70 −35/−10 pairs
- [predict_sigma70_promoters](../analysis/predict_sigma70_promoters.md) - Promoter Calculator v1.0
- [find_ribosome_binding_sites](find_ribosome_binding_sites.md) - Shine-Dalgarno RBS motifs
- [predict_genes](predict_genes.md) - ORF-based gene prediction
