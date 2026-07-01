# differential_expression

Two-group differential-expression analysis.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Annotation |
| **Tool Name** | `differential_expression` |
| **Method ID** | `TranscriptomeAnalyzer.AnalyzeDifferentialExpression` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

For each gene, computes the log2 fold change between the two groups' mean expression (with a pseudocount),
a t-test p-value, and a Benjamini-Hochberg FDR-adjusted p-value across all genes. A gene is flagged
significant when both `|log2 fold change| ≥ foldChangeThreshold` and `adjustedPValue < pValueThreshold`.
Regulation is **Upregulated** (positive log2FC), **Downregulated** (negative), or **Unchanged**.

## Core Documentation Reference

- Source: [TranscriptomeAnalyzer.cs#L254](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/TranscriptomeAnalyzer.cs#L254)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `expressionData` | array | Yes | Per-gene `{ geneId, group1[], group2[] }` |
| `foldChangeThreshold` | number | No | Minimum \|log2FC\| for significance (default 1.0) |
| `pValueThreshold` | number | No | Adjusted-p-value threshold (default 0.05) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `results[].geneId` | string | Gene identifier |
| `results[].log2FoldChange` | number | log2 fold change (group2 vs group1) |
| `results[].pValue` | number | Raw t-test p-value |
| `results[].adjustedPValue` | number | BH-adjusted p-value |
| `results[].isSignificant` | boolean | Both criteria satisfied |
| `results[].regulation` | string | Upregulated / Downregulated / Unchanged |

## Errors

| Code | Message |
|------|---------|
| 1001 | expressionData cannot be null |
| 1001 | expressionData cannot be empty. |

## Examples

### Example 1: Upregulated gene

GENE1 ~10× higher in group2 → significant, Upregulated.

### Example 2: Downregulated gene

Reversing the groups gives a negative log2 fold change and Downregulated.

## Performance

- **Time Complexity:** O(g · r) for g genes, r replicates (plus O(g log g) BH sort)
- **Space Complexity:** O(g)

## See Also

- [over_representation_analysis](over_representation_analysis.md) — pathway ORA on DE genes
- [enrichment_score](enrichment_score.md) — GSEA-like enrichment
