# predict_imprinted_genes

Predict imprinted genes from allele-specific methylation differences.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Annotation |
| **Tool Name** | `predict_imprinted_genes` |
| **Method ID** | `EpigeneticsAnalyzer.PredictImprintedGenes` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Flags candidate imprinted genes from parent-of-origin methylation. A gene is reported when the absolute
difference between its maternal and paternal methylation is at least `minDifference` (default 0.4). Each
hit records:

- `parentalOrigin` — **Maternal** when maternal methylation is higher, otherwise **Paternal**.
- `imprintingScore` — `|diff| / (maternal + paternal + 0.01)`, capped at 1.
- `hasDMR` — `true` when `|diff| > 0.5` (a differentially methylated region).

## Core Documentation Reference

- Source: [EpigeneticsAnalyzer.cs#L1094](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/EpigeneticsAnalyzer.cs#L1094)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `genes` | array | Yes | Per-gene `{ geneId, start, end, maternalMethylation, paternalMethylation }` |
| `minDifference` | number | No | Minimum allele methylation difference (default 0.4) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `imprinted[].geneId` | string | Gene identifier |
| `imprinted[].start` | integer | Gene start |
| `imprinted[].end` | integer | Gene end |
| `imprinted[].imprintingScore` | number | diff / (maternal + paternal + 0.01), capped at 1 |
| `imprinted[].parentalOrigin` | string | `Maternal` or `Paternal` |
| `imprinted[].hasDMR` | boolean | true when \|diff\| > 0.5 |

## Errors

| Code | Message |
|------|---------|
| 1001 | genes cannot be null |

## Examples

### Example 1: Maternally methylated locus

IGF2 with maternal 0.9 / paternal 0.1 (diff 0.8):

**Response:**
```json
{ "imprinted": [ { "geneId": "IGF2", "start": 1000, "end": 2000, "imprintingScore": 0.7920792079207921, "parentalOrigin": "Maternal", "hasDMR": true } ] }
```

### Example 2: Below difference threshold

A balanced locus (0.5 / 0.5) is not imprinted:

**Response:**
```json
{ "imprinted": [] }
```

## Performance

- **Time Complexity:** O(n) for n genes
- **Space Complexity:** O(k) for k imprinted hits

## See Also

- [find_dmrs](find_dmrs.md) — differentially methylated regions between samples
- [methylation_profile](methylation_profile.md) — methylation summaries
