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

Scans a DNA sequence for the two canonical bacterial (E. coli σ70) promoter elements: the **-35 box**
consensus `TTGACA` and the **-10 box** (Pribnow box) consensus `TATAAT`, along with their common
sub-variants (5-bp prefix/suffix and 4-bp prefix). Each hit reports a 0-based `position`, the box
`type` (`-35 box` / `-10 box`), the matched `sequence`, and a probability-weighted `score` derived
from E. coli position-specific nucleotide occurrence frequencies (Harley & Reynolds 1987); the full
6-bp consensus scores `1.0`. A full `TTGACA` therefore yields four `-35 box` hits (`TTGACA`, `TTGAC`,
`TGACA`, `TTGA`).

## Core Documentation Reference

- Source: [GenomeAnnotator.cs#L685](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/GenomeAnnotator.cs#L685)

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

### Example 1: Full -35 box consensus

`GGGGGTTGACAGGGGG` → `TTGACA` at position 5, `-35 box`, score 1.0 (plus three sub-variant hits).

**Response:**
```json
{ "motifs": [ { "position": 5, "type": "-35 box", "sequence": "TTGACA", "score": 1.0 } ] }
```

### Example 2: Full -10 box consensus

`CCCCCTATAATCCCCC` → `TATAAT` at position 5, `-10 box`, score 1.0.

## Performance

- **Time Complexity:** O(n · motifs)
- **Space Complexity:** O(k)

## See Also

- [find_ribosome_binding_sites](find_ribosome_binding_sites.md) - Shine-Dalgarno RBS motifs
- [predict_genes](predict_genes.md) - ORF-based gene prediction
