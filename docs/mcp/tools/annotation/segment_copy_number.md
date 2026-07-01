# segment_copy_number

Segment copy-number probe data into log-ratio plateaus.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Annotation |
| **Tool Name** | `segment_copy_number` |
| **Method ID** | `StructuralVariantAnalyzer.SegmentCopyNumber` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Segments position-sorted copy-number probes (CBS-like): a new segment starts at a chromosome change or when
a probe's log-ratio deviates from the running segment mean by more than `changeThreshold`, provided the
current segment already has at least `minProbes` probes. Each emitted segment reports its span, mean
log-ratio, estimated copy number `round(2 · 2^meanLogR)` (clamped to [0, 10]), mean B-allele frequency, and
probe count. Segments with fewer than `minProbes` probes are dropped.

## Core Documentation Reference

- Source: [StructuralVariantAnalyzer.cs#L716](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/StructuralVariantAnalyzer.cs#L716)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `probes` | array | Yes | Probes `{ chromosome, position, logRatio, baf }` |
| `changeThreshold` | number | No | Log-ratio change to start a new segment (default 0.3) |
| `minProbes` | integer | No | Minimum probes per segment (default 5) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `segments[]` | object | Segments (span, logRatio, copyNumber, bAlleleFrequency, probeCount) |

## Errors

| Code | Message |
|------|---------|
| 1001 | probes cannot be null |

## Examples

### Example 1: Diploid plateau

Five probes at log-ratio 0.0 form one segment with copy number 2:

**Response:**
```json
{ "segments": [ { "start": 100, "end": 500, "logRatio": 0.0, "copyNumber": 2, "probeCount": 5 } ] }
```

### Example 2: Too few probes

Fewer than `minProbes` (5) probes yield no segment:

**Response:**
```json
{ "segments": [] }
```

## Performance

- **Time Complexity:** O(n log n) for n probes
- **Space Complexity:** O(n)

## See Also

- [identify_cnvs](identify_cnvs.md) — turn non-baseline segments into CNVs
