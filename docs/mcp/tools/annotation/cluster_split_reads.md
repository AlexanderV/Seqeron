# cluster_split_reads

Cluster split reads into breakpoint candidates.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Annotation |
| **Tool Name** | `cluster_split_reads` |
| **Method ID** | `StructuralVariantAnalyzer.ClusterSplitReads` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Groups split reads on the same chromosome whose primary positions are within `clusterDistance`, and emits a
breakpoint for each cluster meeting `minSupport`. The breakpoint position is the average of the cluster's
primary and supplementary positions, and its quality is `min(support · 15, 100)`.

## Core Documentation Reference

- Source: [StructuralVariantAnalyzer.cs#L495](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/StructuralVariantAnalyzer.cs#L495)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `splitReads` | array | Yes | Split reads to cluster |
| `clusterDistance` | integer | No | Max primary-position distance in a cluster (default 10) |
| `minSupport` | integer | No | Minimum supporting reads per cluster (default 2) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `breakpoints[]` | object | Breakpoints (positions, strands, supportingReads, quality) |

## Errors

| Code | Message |
|------|---------|
| 1001 | splitReads cannot be null |

## Examples

### Example 1: Two nearby split reads

Primaries at 1000 and 1004 (within 10) form one breakpoint at position 1002 with quality 30:

**Response:**
```json
{ "breakpoints": [ { "position1": 1002, "position2": 2000, "supportingReads": 2, "quality": 30 } ] }
```

### Example 2: Below minimum support

A single split read is below the default `minSupport` of 2:

**Response:**
```json
{ "breakpoints": [] }
```

## Performance

- **Time Complexity:** O(n log n) for n split reads
- **Space Complexity:** O(n)

## See Also

- [find_split_reads](find_split_reads.md) — produce the split reads
- [cluster_discordant_pairs](cluster_discordant_pairs.md) — read-pair clustering
