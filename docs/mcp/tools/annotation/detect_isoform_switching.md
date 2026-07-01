# detect_isoform_switching

Detect isoform-usage switching between two conditions.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Annotation |
| **Tool Name** | `detect_isoform_switching` |
| **Method ID** | `TranscriptomeAnalyzer.DetectIsoformSwitching` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

For each gene with at least two isoforms, computes per-condition usage proportions (isoform expression ÷
gene total) and reports a switch when one isoform's usage rises by more than `switchThreshold` while another
falls by more than `switchThreshold`. The result names the decreasing isoform (`transcriptId1`), the
increasing isoform (`transcriptId2`), and a switch score `|deltaUp| + |deltaDown|`.

## Core Documentation Reference

- Source: [TranscriptomeAnalyzer.cs#L1028](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/TranscriptomeAnalyzer.cs#L1028)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `isoformData` | array | Yes | Per-isoform `{ isoform, expression1, expression2 }` |
| `switchThreshold` | number | No | Minimum \|delta-usage\| (default 0.3) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `switches[].geneId` | string | Gene id |
| `switches[].transcriptId1` | string | Isoform whose usage decreased |
| `switches[].transcriptId2` | string | Isoform whose usage increased |
| `switches[].switchScore` | number | \|deltaUp\| + \|deltaDown\| |

## Errors

| Code | Message |
|------|---------|
| 1001 | isoformData cannot be null |

## Examples

### Example 1: Dominant isoform switch

TX1 (0.9→0.1) and TX2 (0.1→0.9) → switch with score 1.6:

**Response:**
```json
{ "switches": [ { "geneId": "GENE1", "transcriptId1": "TX1", "transcriptId2": "TX2", "switchScore": 1.6 } ] }
```

### Example 2: Stable usage

No isoform crosses the threshold:

**Response:**
```json
{ "switches": [] }
```

## Performance

- **Time Complexity:** O(n log n) per gene (usage sort)
- **Space Complexity:** O(n)

## See Also

- [find_dominant_isoforms](find_dominant_isoforms.md) — per-gene dominant isoform
