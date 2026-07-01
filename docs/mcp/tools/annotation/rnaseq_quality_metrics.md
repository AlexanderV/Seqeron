# rnaseq_quality_metrics

Compute basic RNA-seq quality-control metrics.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Annotation |
| **Tool Name** | `rnaseq_quality_metrics` |
| **Method ID** | `TranscriptomeAnalyzer.CalculateQualityMetrics` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Computes standard RNA-seq QC rates: `mappingRate = mappedReads / totalReads`, `exonicRate = exonicReads /
mappedReads`, `rRnaRate = rRnaReads / mappedReads`, and `detectedGenes` = the number of genes with a
non-zero read count. Rates are 0 when their denominator is 0.

## Core Documentation Reference

- Source: [TranscriptomeAnalyzer.cs#L1225](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/TranscriptomeAnalyzer.cs#L1225)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `totalReads` | number | Yes | Total sequenced reads (≥ 0) |
| `mappedReads` | number | Yes | Reads mapped to reference (≥ 0) |
| `exonicReads` | number | Yes | Reads mapped to exons (≥ 0) |
| `rRnaReads` | number | Yes | Reads mapped to rRNA (≥ 0) |
| `geneCounts` | number[] | Yes | Per-gene read counts |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `mappingRate` | number | mappedReads / totalReads |
| `exonicRate` | number | exonicReads / mappedReads |
| `rRnaRate` | number | rRnaReads / mappedReads |
| `detectedGenes` | integer | genes with count > 0 |

## Errors

| Code | Message |
|------|---------|
| 1001 | geneCounts cannot be null |
| 1001 | read counts cannot be negative. |

## Examples

### Example 1: Typical library

**Response:**
```json
{ "mappingRate": 0.9, "exonicRate": 0.8, "rRnaRate": 0.05, "detectedGenes": 3 }
```

### Example 2: Empty library

```json
{ "mappingRate": 0, "exonicRate": 0, "rRnaRate": 0, "detectedGenes": 0 }
```

## Performance

- **Time Complexity:** O(g) for g genes
- **Space Complexity:** O(1)

## See Also

- [calculate_tpm](calculate_tpm.md) — TPM abundance
- [perform_pca](perform_pca.md) — sample PCA
