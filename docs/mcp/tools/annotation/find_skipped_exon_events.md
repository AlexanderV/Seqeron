# find_skipped_exon_events

Compute Percent Spliced In (PSI) for skipped-exon candidates.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Annotation |
| **Tool Name** | `find_skipped_exon_events` |
| **Method ID** | `TranscriptomeAnalyzer.FindSkippedExonEvents` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

For each candidate exon, computes the Percent Spliced In value `PSI = inclusion / (inclusion + skipping)`
(PMC3330053; SUPPA2) and emits a `SkippedExon` event with that inclusion level. Exons with zero total reads
are skipped (undefined PSI). Delta-PSI is 0 because this is a single-sample measurement.

## Core Documentation Reference

- Source: [TranscriptomeAnalyzer.cs#L950](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/TranscriptomeAnalyzer.cs#L950)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `exonData` | array | Yes | Per-exon `{ geneId, exonStart, exonEnd, inclusionReads, skippingReads }` |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `events[].geneId` | string | Gene id |
| `events[].eventType` | string | `SkippedExon` |
| `events[].start` / `.end` | integer | Exon coordinates |
| `events[].inclusionLevel` | number | PSI |
| `events[].deltaPsi` | number | 0 (single-sample) |

## Errors

| Code | Message |
|------|---------|
| 1001 | exonData cannot be null |

## Examples

### Example 1: PSI 0.8

Inclusion 80 / skipping 20 → PSI 0.8:

**Response:**
```json
{ "events": [ { "geneId": "G1", "eventType": "SkippedExon", "inclusionLevel": 0.8, "deltaPsi": 0 } ] }
```

### Example 2: No reads

An exon with zero total reads is skipped:

**Response:**
```json
{ "events": [] }
```

## Performance

- **Time Complexity:** O(n) for n exons
- **Space Complexity:** O(n)

## See Also

- [detect_differential_splicing](detect_differential_splicing.md) — two-condition delta-PSI
