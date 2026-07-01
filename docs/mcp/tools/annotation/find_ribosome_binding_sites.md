# find_ribosome_binding_sites

Locate Shine-Dalgarno (ribosome binding site) motifs upstream of forward-strand ORFs.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Annotation |
| **Tool Name** | `find_ribosome_binding_sites` |
| **Method ID** | `GenomeAnnotator.FindRibosomeBindingSites` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Finds Shine-Dalgarno ribosome-binding motifs (the `AGGAGG` consensus and its sub-motifs `GGAGG`,
`AGGAG`, `GAGG`, `AGGA`) on the **forward strand**, upstream of every forward-strand ORF of at least
30 aa. For each ORF the tool scans the `upstreamWindow` nucleotides preceding the start codon and
keeps a motif only when its aligned spacing to the start codon falls in `[minDistance, maxDistance]`
(default 4–15 nt, the functional range per Shine & Dalgarno 1975 and Chen et al. 1994). Each hit
reports the motif's 0-based `position`, the matched `sequence`, and a `score` equal to
`motif.Length / 6` (so the full `AGGAGG` consensus scores 1.0).

## Core Documentation Reference

- Source: [GenomeAnnotator.cs#L280](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/GenomeAnnotator.cs#L280)

## Input Schema

| Parameter | Type | Required | Default | Description |
|-----------|------|----------|---------|-------------|
| `dnaSequence` | string | Yes | — | DNA sequence to scan |
| `upstreamWindow` | integer | No | 20 | Upstream window size (nt) to search before each ORF start |
| `minDistance` | integer | No | 4 | Minimum aligned distance from motif to start codon (nt) |
| `maxDistance` | integer | No | 15 | Maximum aligned distance from motif to start codon (nt) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `sites` | array | `{ position, sequence, score }` per Shine-Dalgarno hit |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |

## Examples

### Example 1: Canonical AGGAGG at 8 nt spacing

A sequence built as `10×C + AGGAGG + 8×C + (ATG…TAA, 100 aa ORF)` yields the full `AGGAGG` consensus
at position 10 with score 1.0.

**Response:**
```json
{ "sites": [ { "position": 10, "sequence": "AGGAGG", "score": 1.0 } ] }
```

### Example 2: Too close to the start codon

With `AGGAGG` placed at aligned spacing 2 (below `minDistance=4`), the full consensus is filtered
out of the results.

## Performance

- **Time Complexity:** O(ORFs · window · motifs)
- **Space Complexity:** O(k) for the returned hits

## See Also

- [find_orfs](find_orfs.md) - ORF detection (the RBS anchors)
- [find_promoter_motifs](find_promoter_motifs.md) - -10/-35 bacterial promoter boxes
- [predict_genes](predict_genes.md) - ORF-based gene prediction
