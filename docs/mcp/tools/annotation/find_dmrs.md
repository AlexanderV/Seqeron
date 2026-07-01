# find_dmrs

Identify differentially methylated regions (DMRs) between two methylation samples.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Annotation |
| **Tool Name** | `find_dmrs` |
| **Method ID** | `EpigeneticsAnalyzer.FindDMRs` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Implements the methylKit tiling-window DMR model (Akalin et al. 2012). CpG positions from both samples are
grouped into fixed-size genomic windows (`windowSize`, default 1000 bp). A window with at least
`minCpGCount` covered cytosines is reported when its mean per-site methylation difference
(`sample2 − sample1`, in fraction units) exceeds `minDifference` in absolute value (strict `>`). Each region
carries a two-sided Fisher's-exact p-value on the pooled methylated/unmethylated read counts, and is
labelled **Hypermethylated** when the treatment is higher than the control (positive difference) or
**Hypomethylated** when lower.

## Core Documentation Reference

- Source: [EpigeneticsAnalyzer.cs#L607](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/EpigeneticsAnalyzer.cs#L607)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sample1` | array | Yes | Control sample methylation sites |
| `sample2` | array | Yes | Treatment sample methylation sites |
| `windowSize` | integer | No | Tiling window width in bp (default 1000) |
| `minDifference` | number | No | Minimum absolute mean difference, strict (default 0.25) |
| `minCpGCount` | integer | No | Minimum CpG count per window (default 3) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `regions[].start` | integer | Region start (first window position) |
| `regions[].end` | integer | Region end (last window position) |
| `regions[].meanDifference` | number | mean(sample2 − sample1) over the window |
| `regions[].pValue` | number | Two-sided Fisher's exact test on pooled counts |
| `regions[].cpGCount` | integer | Number of CpG sites in the region |
| `regions[].annotation` | string | `Hypermethylated` or `Hypomethylated` |

## Errors

| Code | Message |
|------|---------|
| 1001 | sample1 cannot be null |
| 1001 | sample2 cannot be null |
| 1001 | Both samples cannot be empty |

## Examples

### Example 1: Hypermethylated region

Control unmethylated, treatment fully methylated at positions 10/20/30 (one window):

**Response:**
```json
{ "regions": [ { "start": 10, "end": 30, "meanDifference": 1.0, "pValue": 0.0, "cpGCount": 3, "annotation": "Hypermethylated" } ] }
```

### Example 2: No difference

Identical samples produce no DMRs:

**Response:**
```json
{ "regions": [] }
```

## Performance

- **Time Complexity:** O(n log n) for n covered positions (sort + single tiling pass)
- **Space Complexity:** O(n)

## See Also

- [methylation_profile](methylation_profile.md) — single-sample methylation summary
- [methylation_from_bisulfite](methylation_from_bisulfite.md) — build per-CpG sites from reads
