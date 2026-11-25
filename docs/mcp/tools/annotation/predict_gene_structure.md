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

Predicts a gene's exon/intron layout. `model = greedy` (default) greedily selects non-overlapping
high-scoring introns (splice-motif model, not coding-aware), then derives the intervening exons.
`model = splice_dp` runs the exact splice-site + open-reading-frame dynamic program
(`SpliceSitePredictor.PredictGeneStructureDp`, GENSCAN / GeneSplicer style): it maximises the sum of
MaxEntScan 5' + 3' scores (bits) minus `intronPenalty` per intron subject to an open reading frame from the
start codon (`startPosition`, default first AUG) to the first in-frame stop, with codons possibly split by
introns; UTRs are intron-free, U12 introns are not modelled. It errors when no open reading frame is
reachable. Measured on Ensembl coding introns: HBB 2/2, TP53 6/9, CFTR exons 2-5 0/3 (greedy 0/14); low precision
(see the algorithm doc). Exons are typed **Single**, **Initial**, **Internal**, or **Terminal**
(Gilbert 1978). The tool returns the exons, introns, the spliced (exon-joined) sequence, and an overall
score. Exon and intron lengths together cover the entire input.

## Core Documentation Reference

- Source: [SpliceSitePredictor.cs#L1337](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/SpliceSitePredictor.cs#L1337) (greedy), [SpliceSitePredictor.GeneStructureDp.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/SpliceSitePredictor.GeneStructureDp.cs) (splice_dp)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | RNA/DNA sequence to analyze |
| `minExonLength` | integer | No | Minimum exon length (default 30 greedy, 3 splice_dp; effective >= 2 for splice_dp) |
| `minIntronLength` | integer | No | Minimum intron length (default 60) |
| `minScore` | number | No | greedy only: minimum combined intron score (default 0.5) |
| `model` | string | No | `greedy` (default) or `splice_dp` |
| `startPosition` | integer | No | splice_dp: 0-based first coding nucleotide (default first AUG) |
| `maxIntronLength` | integer | No | splice_dp: maximum intron length (default 100000) |
| `minDonorScore` / `minAcceptorScore` | number | No | splice_dp: minimum MaxEntScan score5 / score3 in bits (default 0) |
| `intronPenalty` | number | No | splice_dp: bits subtracted per intron (default 12; unpublished, tunable) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `exons` | array | Predicted exon segments (with `type`) |
| `introns` | array | Predicted introns |
| `splicedSequence` | string | Exon-joined (intron-removed) sequence |
| `overallScore` | number | Overall structure score (greedy: mean intron score; splice_dp: objective value) |
| `cdsStart` | integer / null | splice_dp: genomic index of the start codon |
| `codingSequence` | string / null | splice_dp: spliced CDS from the start codon to the stop codon (RNA) |
| `intronScore` | number / null | splice_dp: sum of (score5 + score3 - intronPenalty) |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1001 | splice_dp: no open reading frame reachable / unknown model |

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

- **Time Complexity:** greedy: O(D·A) for donor/acceptor pairing plus greedy selection; splice_dp: O(acceptors x donors in the reading-frame window) (TP53 19 kb ~0.1 s)
- **Space Complexity:** O(n)

## See Also

- [predict_introns](predict_introns.md) — intron candidates only
- [detect_alternative_splicing](detect_alternative_splicing.md) — alternative splicing events
