# predict_gene_structure

Predict exon/intron gene structure.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Annotation |
| **Tool Name** | `predict_gene_structure` |
| **Method ID** | `SpliceSitePredictor.PredictGeneStructure` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Predicts a gene's exon/intron layout by greedily selecting non-overlapping high-scoring introns, then
deriving the intervening exons. Exons are typed **Single**, **Initial**, **Internal**, or **Terminal**
(Gilbert 1978). The tool returns the exons, introns, the spliced (exon-joined) sequence, and an overall
score. Exon and intron lengths together cover the entire input.

## Core Documentation Reference

- Source: [SpliceSitePredictor.cs#L625](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/SpliceSitePredictor.cs#L625)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | RNA/DNA sequence to analyze |
| `minExonLength` | integer | No | Minimum exon length (default 30) |
| `minIntronLength` | integer | No | Minimum intron length (default 60) |
| `minScore` | number | No | Minimum combined intron score (default 0.5) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `exons` | array | Predicted exon segments (with `type`) |
| `introns` | array | Predicted introns |
| `splicedSequence` | string | Exon-joined (intron-removed) sequence |
| `overallScore` | number | Overall structure score |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |

## Examples

### Example 1: Two-exon gene

A two-exon sequence with one GT-AG intron yields an Initial and a Terminal exon plus one U2 intron.

### Example 2: Single-exon gene

A sequence with no GU dinucleotide is a single exon with no introns:

**Response:**
```json
{ "exons": [ { "type": "Single" } ], "introns": [] }
```

## Performance

- **Time Complexity:** O(D·A) for donor/acceptor pairing plus greedy selection
- **Space Complexity:** O(n)

## See Also

- [predict_introns](predict_introns.md) — intron candidates only
- [detect_alternative_splicing](detect_alternative_splicing.md) — alternative splicing events
