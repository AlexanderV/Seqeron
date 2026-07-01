# predict_introns

Predict introns by pairing donor and acceptor splice sites.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Annotation |
| **Tool Name** | `predict_introns` |
| **Method ID** | `SpliceSitePredictor.PredictIntrons` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Finds candidate 5' donor and 3' acceptor sites, pairs them into introns whose length falls within
`[minIntronLength, maxIntronLength]` and whose combined score is at least `minScore`, and locates a branch
point where possible. GT-AG introns are classified as **U2** (major spliceosome); AT-AC as **U12**
(Burge et al. 1999). Each intron reports its span, length, flanking sites, optional branch point,
sequence, type, and score.

## Core Documentation Reference

- Source: [SpliceSitePredictor.cs#L625](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/SpliceSitePredictor.cs#L625)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | RNA/DNA sequence to analyze |
| `minIntronLength` | integer | No | Minimum intron length (default 60) |
| `maxIntronLength` | integer | No | Maximum intron length (default 100000) |
| `minScore` | number | No | Minimum combined score (default 0.5) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `introns[].start` / `.end` | integer | Intron span |
| `introns[].length` | integer | Intron length |
| `introns[].donorSite` / `.acceptorSite` | object | Flanking splice sites |
| `introns[].branchPoint` | object \| null | Branch point when found |
| `introns[].sequence` | string | Intron sequence |
| `introns[].type` | string | `U2` (GT-AG) or `U12` |
| `introns[].score` | number | Combined score |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |

## Examples

### Example 1: GT-AG intron classified U2

A two-exon sequence with a `GUAAGU`…`CAG` intron ≥ 60 nt yields a `U2` intron whose sequence starts `GU`
and ends `AG`.

### Example 2: All candidates below minimum length

Raising `minIntronLength` above every candidate length yields no introns:

**Response:**
```json
{ "introns": [] }
```

## Performance

- **Time Complexity:** O(D·A) for D donors and A acceptors
- **Space Complexity:** O(k) for k introns

## See Also

- [find_donor_sites](find_donor_sites.md) / [find_acceptor_sites](find_acceptor_sites.md) — the paired sites
- [predict_gene_structure](predict_gene_structure.md) — full exon/intron structure
