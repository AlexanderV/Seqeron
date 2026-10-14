# find_rare_codon_clusters

Detect rare-codon clusters with the Sherlocc rule (Chartier et al. 2012).

## Overview

| Property | Value |
|----------|-------|
| **Server** | MolTools |
| **Tool Name** | `find_rare_codon_clusters` |
| **Method ID** | `CodonOptimizer.FindRareCodonClusters` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Splits the coding sequence into frame-0 codons (T→U) and marks a codon as rare ("pause") when its frequency in the target organism's codon-usage table is strictly below `threshold` (default 0.15) — the same per-codon criterion as [find_rare_codons](find_rare_codons.md); an ambiguous triplet is never rare. A window of `window_size` codons (default 7) holding at least `min_rare_codons` rare codons (default 4) is a cluster — the published Sherlocc rule "a seven position-wide window … containing at least four pause positions out of seven" (Chartier, Gaudreault & Najmanovich 2012, Bioinformatics 28:1438, doi:10.1093/bioinformatics/bts149). Overlapping cluster windows are merged into maximal regions, so a long rare run is reported once.

## Core Documentation Reference

- Source: [CodonOptimizer.cs#L1020](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/CodonOptimizer.cs#L1020)
- Algorithm: [CODON-RARE-001](../../../../tests/TestSpecs/CODON-RARE-001.md), Evidence: [CODON-RARE-001-Evidence](../../../Evidence/CODON-RARE-001-Evidence.md)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `coding_sequence` | string | Yes | Coding sequence (DNA or RNA, frame 0), non-empty. |
| `target_organism` | object | Yes | Preset id (`EColiK12` \| `Yeast` \| `Human`) or inline custom table. |
| `threshold` | number | No | Per-codon rare-frequency cutoff, strict < (default 0.15). |
| `window_size` | integer | No | Cluster window width in codons (default 7); must be ≥ 1. |
| `min_rare_codons` | integer | No | Minimum rare codons per window (default 4); must be ≥ 1. |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `clusters` | array | Maximal, non-overlapping clusters in codon order, each `{startCodon, endCodon, rareCount}` (0-based inclusive codon indices; rare codons inside the region). |

## Errors

| Code | Message |
|------|---------|
| 1001 | Coding sequence cannot be null or empty |
| 1002 | Window size must be at least 1 codon |
| 1003 | Minimum rare codons must be at least 1 |

## Examples

### Example 1: Two separated rare runs (E. coli K12 preset, Sherlocc defaults)

7× `AGA` + 8× `CGC` + 7× `AGA` (22 codons; AGA frequency 0.04 < 0.15) → two clusters: `{startCodon: 0, endCodon: 9, rareCount: 7}` and `{startCodon: 12, endCodon: 21, rareCount: 7}`.

### Example 2: Custom rule

`AGAAGAAGA`, `window_size` 3, `min_rare_codons` 3 → `{startCodon: 0, endCodon: 2, rareCount: 3}`.

## See Also

- [codon_minmax_profile](codon_minmax_profile.md), [find_rare_codons](find_rare_codons.md), [optimize_codons](optimize_codons.md)
