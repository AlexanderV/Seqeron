# cluster_discordant_pairs

Cluster discordant read pairs into structural-variant candidates.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Annotation |
| **Tool Name** | `cluster_discordant_pairs` |
| **Method ID** | `StructuralVariantAnalyzer.ClusterDiscordantPairs` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Groups discordant read pairs that share the same chromosome pair and whose anchor positions are within
`clusterDistance` of one another, then emits one structural variant per cluster with at least `minSupport`
supporting pairs. The SV type is inferred from the cluster's paired-end signature (e.g. interchromosomal →
Translocation).

## Core Documentation Reference

- Source: [StructuralVariantAnalyzer.cs#L300](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/StructuralVariantAnalyzer.cs#L300)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `discordantPairs` | array | Yes | Discordant read-pair signatures |
| `clusterDistance` | integer | No | Max anchor distance within a cluster (default 500) |
| `minSupport` | integer | No | Minimum supporting pairs per cluster (default 3) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `variants[]` | object | SV candidates (id, chromosome, start/end, type, length, quality, supportingReads) |

## Errors

| Code | Message |
|------|---------|
| 1001 | discordantPairs cannot be null |

## Examples

### Example 1: Three interchromosomal pairs

Three chr1→chr2 pairs at nearby anchors form one Translocation with 3 supporting reads.

### Example 2: Below minimum support

Two supporting pairs are below the default `minSupport` of 3:

**Response:**
```json
{ "variants": [] }
```

## Performance

- **Time Complexity:** O(n log n) (sort + single clustering pass) for n pairs
- **Space Complexity:** O(n)

## See Also

- [find_discordant_pairs](find_discordant_pairs.md) — produce the discordant pairs
- [merge_overlapping_svs](merge_overlapping_svs.md) — merge SV calls
