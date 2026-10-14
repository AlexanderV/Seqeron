# find_conserved_rare_codon_clusters

Find rare-codon clusters conserved across a codon alignment of homologous genes (Sherlocc 1.0, Chartier et al. 2012).

## Overview

| Property | Value |
|----------|-------|
| **Server** | MolTools |
| **Tool Name** | `find_conserved_rare_codon_clusters` |
| **Method ID** | `CodonOptimizer.FindConservedRareCodonClusters` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Sherlocc 1.0 stage-3 analysis (Chartier, Gaudreault & Najmanovich 2012, Bioinformatics 28:1438, doi:10.1093/bioinformatics/bts149; `sherlocc.pl`, github.com/mtthchrtr/sherlocc commit b3b4d15) of a protein-family **codon alignment**. Delegates to `CodonOptimizer.FindConservedRareCodonClusters` (single-table or per-sequence-table overload).

1. Every cell (sequence × codon column) gets its codon's usage from the organism table (stage 2, l. 494–522). By default this is the **per-thousand** usage (Sherlocc's CUTG `/1000` column, `files/cut`): the `EColiK12`, `Yeast` and `Human` presets carry it (Kazusa 316407 / 4932 / 9606; E. coli since 2026-10-14); an inline table's `codonFrequencies` values are taken **as** per-thousand usage; `use_codon_frequencies` = true (`SherloccOptions.UseCodonFrequencies`) scores the per-amino-acid fractions / inline values instead (the threshold is then on that scale). A table without per-thousand usage is rejected by the C# API unless `UseCodonFrequencies` is set; no MCP input lacks it. A codon of only `-`/`.` (gap), a partially gapped or ambiguous triplet, or a codon absent from the table has no value (Sherlocc `x`).
2. For every column *i* the valued cells of columns *i*−3 … *i*+3 (`averaging_window`, default 7) of **all** sequences are averaged (l. 780–798).
3. The column is representative only when it has ≥ `min_coverage` × sequences × `averaging_window` valued cells (default 0.75, l. 772/799 — the gap/conservation filter); otherwise its average is `null` (Sherlocc "NA").
4. A representative column with average ≤ `threshold` is a pause (l. 803, inclusive).
5. Every window of `cluster_window` columns (default 7) holding ≥ `min_pauses` pauses (default 4) yields a cluster from its first to its last pause (l. 838–854); a cluster starting inside the previous one is merged with it, abutting clusters are not (l. 894–961).
6. Per cluster: pause count, Sherlocc "Usage Frequency Average" (mean of the column averages over the span; a non-representative column contributes 0, as Perl's undefined value) and the 1-based "Middle Point" `stop − ⌊length/2⌋`; per alignment: total cluster length, fraction of the alignment in clusters, mean cluster usage (l. 875–879, 1131–1155).

Give either `target_organism` (one table for all rows) or `target_organisms_per_sequence` (one table per row, in alignment order — Sherlocc scores each family member with its own organism's usage, l. 497). The analysis is independent of the genetic code (each cell is a per-codon table lookup), so the tool has no `genetic_code` parameter. Not ported: FoldIndex/Uversky disorder prediction and stage 1–2 retrieval (UniProt/ENA queries); no cluster P-values (not computed by `sherlocc.pl`).

## Core Documentation Reference

- Source: [CodonOptimizer.cs#L2035](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/CodonOptimizer.cs#L2035)
- Algorithm: [CODON-RARE-001](../../../../tests/TestSpecs/CODON-RARE-001.md), Evidence: [CODON-RARE-001-Evidence](../../../Evidence/CODON-RARE-001-Evidence.md)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `aligned_coding_sequences` | string[] | Yes | Aligned in-frame coding sequences (DNA or RNA, any case), equal length, a multiple of 3; gaps `-` or `.`. |
| `threshold` | number | Yes | Pause cutoff: column average ≤ threshold (per thousand by default; Sherlocc README examples 14–16). |
| `target_organism` | object | One of | One table for every sequence: preset id (`EColiK12` \| `Yeast` \| `Human`, all with Kazusa per-thousand usage) or inline custom table whose `codonFrequencies` are per-thousand usage. |
| `target_organisms_per_sequence` | object[] | One of | One table per aligned sequence, same order (instead of `target_organism`). |
| `use_codon_frequencies` | boolean | No | `false` (default): per-thousand usage (Sherlocc); `true`: the tables' `codonFrequencies`. |
| `averaging_window` | integer | No | Centred averaging window in columns (default 7); positive and odd. |
| `cluster_window` | integer | No | Cluster window width in columns (default 7); ≥ 1. |
| `min_pauses` | integer | No | Minimum pauses in a cluster window (default 4); ≥ 1. |
| `min_coverage` | number | No | Representativeness filter fraction in [0, 1] (default 0.75). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `sequenceCount` | integer | Number of aligned sequences. |
| `alignmentLength` | integer | Number of codon columns. |
| `columnAverages` | (number \| null)[] | Per column: window-and-family average usage, `null` when not representative. |
| `columnUsageCounts` | integer[] | Per column: number of valued window cells. |
| `pauseColumns` | integer[] | 0-based pause columns. |
| `clusters` | array | Each `{startColumn, endColumn, pauseCount, meanUsage, middlePoint}` — 0-based inclusive columns (Sherlocc prints them 1-based), 1-based middle point. |
| `totalClusterLength` | integer | Σ (endColumn − startColumn + 1). |
| `clusterFraction` | number | totalClusterLength / alignmentLength. |
| `meanClusterUsage` | number \| null | Mean of the clusters' `meanUsage`; `null` without clusters. |

## Errors

| Code | Message |
|------|---------|
| 1001 | At least one aligned coding sequence is required |
| 1002 | Supply exactly one of target_organism or target_organisms_per_sequence |
| 1003 | Sequence length is not a multiple of 3 / all aligned sequences must have the same length |
| 1004 | One codon-usage table per sequence is required |
| 1005 | Table has no per-thousand usage (set use_codon_frequencies; C# API only since 2026-10-14 - every preset and inline table carries per-thousand usage) |
| 1006 | Threshold must be a number; AveragingWindow must be a positive odd number; ClusterWindow/MinPauses must be at least 1; MinCoverage must be in [0, 1] |
| 1007 | Unknown codon-usage preset |

## Examples

### Example 1: Yeast family alignment, one merged cluster (sherlocc.pl-locked)

Three yeast sequences × 20 codon columns (two gaps) with a rare Arg/Pro/Leu/Gly/Ala stretch, `threshold` 15, `target_organism` `{ "preset": "Yeast" }` → pauses at columns 8–11, one cluster `{startColumn: 8, endColumn: 11, pauseCount: 4, meanUsage: 9.839285714285714, middlePoint: 10}` (Perl `clusters=9-12|9.83928571428571`), `clusterFraction` 0.2.

### Example 2: Per-sequence tables

Three rows scored with Yeast / Human / Yeast (`target_organisms_per_sequence`), `threshold` 16 → two clusters, columns 3–9 and 20–26 (Perl `4-10|9.20962052101093,21-27|9.30926559372343`), `meanClusterUsage` 9.25944305736718.

### Example 3: E. coli per thousand and fractions

`["AGAAGAAGAAGAAGAAGAAGA", "AGGAGGAGGAGGAGGAGGAGG"]`, `threshold` 0.15, `{ "preset": "EColiK12" }`, `use_codon_frequencies` true → column 3 average 0.03 (AGA 0.04, AGG 0.02), pauses 2–4, no cluster; without `use_codon_frequencies` (Kazusa 316407 per thousand: AGA 2.0, AGG 1.1), `threshold` 15 → column 3 average 1.55, pauses 2–4, no cluster (before 2026-10-14 this call was rejected: the preset had no per-thousand usage).

## See Also

- [find_rare_codon_clusters](find_rare_codon_clusters.md) (single sequence, per-codon rule), [codon_minmax_profile](codon_minmax_profile.md), [find_rare_codons](find_rare_codons.md)
