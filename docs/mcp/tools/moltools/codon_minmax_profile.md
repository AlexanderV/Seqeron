# codon_minmax_profile

%MinMax codon-usage profile (Clarke & Clark 2008) over a sliding codon window.

## Overview

| Property | Value |
|----------|-------|
| **Server** | MolTools |
| **Tool Name** | `codon_minmax_profile` |
| **Method ID** | `CodonOptimizer.CalculateMinMaxProfile` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Splits the coding sequence into frame-0 codons (T→U) and, for every window of `window_size` consecutive codons (default 18, Clarke & Clark 2008), computes the signed %MinMax value: with `Xij` the usage frequency of the codon used, and `Xavg`/`Xmax`/`Xmin` the mean/max/min frequency of its synonymous family, a window with Σ`Xij` > Σ`Xavg` yields %Max = Σ(`Xij`−`Xavg`)/Σ(`Xmax`−`Xavg`)×100 (positive); Σ`Xij` < Σ`Xavg` yields %Min = Σ(`Xavg`−`Xij`)/Σ(`Xavg`−`Xmin`)×100, returned negative; equality gives 0. +100 = window of the most common synonymous codons, −100 = rarest. Reproduced term-for-term from the Clark lab reference code (`calculateMinMax`, CHARMING.py). Synonymous families (including the stop family UAA/UAG/UGA) are those of the NCBI Standard code; Met/Trp contribute 0 to numerator and denominator; an ambiguous triplet occupies its window position but contributes nothing.

**Frequency scale.** The reference implementation takes per-thousand usage (Kazusa "/1000" column); an inline table with per-thousand values reproduces it exactly. The presets (`EColiK12`, `Yeast`, `Human`) hold per-amino-acid relative fractions, which weight every residue equally (equal to the reference only for windows of a single amino acid).

## Core Documentation Reference

- Source: [CodonOptimizer.cs#L919](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/CodonOptimizer.cs#L919)
- Algorithm: [CODON-RARE-001](../../../../tests/TestSpecs/CODON-RARE-001.md), Evidence: [CODON-RARE-001-Evidence](../../../Evidence/CODON-RARE-001-Evidence.md)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `coding_sequence` | string | Yes | Coding sequence (DNA or RNA, frame 0), non-empty. |
| `target_organism` | object | Yes | Preset id (`EColiK12` \| `Yeast` \| `Human`) or inline custom table. |
| `window_size` | integer | No | Sliding-window width in codons (default 18); must be ≥ 1. |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `windows` | array | One entry per window start (codon indices 0 .. codonCount − window_size), each `{windowStartCodon, percentMinMax}`; empty when the sequence has fewer than `window_size` complete codons. |

## Errors

| Code | Message |
|------|---------|
| 1001 | Coding sequence cannot be null or empty |
| 1002 | Window size must be at least 1 codon |

## Examples

### Example 1: E. coli K12 preset, window 3

`CGCCGCAGAAGAAGA` (CGC CGC AGA AGA AGA) → 3 windows at codons 0, 1, 2; the last (3× AGA, frequency 0.04 in an Arg family with mean 1/6 and min 0.02) is −(0.16667−0.04)/(0.16667−0.02)×100 = **−86.3636…**.

`CGCCGCCGC` (3× most-common Arg codon) → **+100**; `CUGAGA`, window 2 → **36.4706…** (%Max).

## See Also

- [find_rare_codon_clusters](find_rare_codon_clusters.md), [find_rare_codons](find_rare_codons.md), [optimize_codons](optimize_codons.md)
