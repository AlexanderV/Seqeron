# calculate_tpm

Compute TPM and FPKM normalized expression from raw counts and gene lengths.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Annotation |
| **Tool Name** | `calculate_tpm` |
| **Method ID** | `TranscriptomeAnalyzer.CalculateTPM` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Normalizes per-gene raw counts to length- and depth-adjusted expression:

- **TPM** — for each gene compute a rate `rate_i = count_i / length_i`, then
  `TPM_i = rate_i / Σ rate * 1_000_000`. Across the gene set the TPM values sum to 1,000,000.
- **FPKM** — `FPKM_i = count_i × 10^9 / (length_i × N)` where `N` is the total mapped reads.

When all counts are zero the denominator is undefined, so every TPM and FPKM is emitted as `0`.

## Core Documentation Reference

- Source: [TranscriptomeAnalyzer.cs#L112](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/TranscriptomeAnalyzer.cs#L112)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `geneCounts` | array | Yes | Per-gene `{ geneId, rawCount, length }` |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `expressions` | array | Per-gene `{ geneId, rawCount, tpm, fpkm, length }` |

## Errors

| Code | Message |
|------|---------|
| 1001 | geneCounts cannot be null |
| 1001 | geneCounts cannot be empty |

## Examples

### Example 1: Two equal-length genes

Counts 10 and 30 over length 1000 each → rates 0.01 and 0.03 (sum 0.04). TPM = 250000 / 750000
(sum 1,000,000); FPKM (total reads 40) = 250000 / 750000.

**Response:**
```json
{ "expressions": [
  { "geneId": "G1", "rawCount": 10, "tpm": 250000, "fpkm": 250000, "length": 1000 },
  { "geneId": "G2", "rawCount": 30, "tpm": 750000, "fpkm": 750000, "length": 1000 }
] }
```

### Example 2: All-zero counts

**Response:**
```json
{ "expressions": [
  { "geneId": "G1", "rawCount": 0, "tpm": 0, "fpkm": 0, "length": 1000 },
  { "geneId": "G2", "rawCount": 0, "tpm": 0, "fpkm": 0, "length": 1000 }
] }
```

## Performance

- **Time Complexity:** O(g) over the genes
- **Space Complexity:** O(g)

## See Also

- [quantile_normalize](quantile_normalize.md) - Cross-sample quantile normalization
- [log2_transform](log2_transform.md) - Variance-stabilizing log transform
- [differential_expression](differential_expression.md) - Compare expression between groups
