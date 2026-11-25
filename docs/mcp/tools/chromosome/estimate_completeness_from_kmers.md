# estimate_completeness_from_kmers

Read k-mer spectrum evaluation: Merqury solid threshold, GenomeScope error rate and genome size; with an
assembly, Merqury k-mer completeness and QV.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Chromosome |
| **Tool Name** | `estimate_completeness_from_kmers` |
| **Method ID** | `GenomeAssemblyAnalyzer.EstimateCompletenessFromKmers` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

- **Histogram** h(c): number of distinct k-mers with count c (as `meryl histogram`).
- **Solid threshold** `solidThreshold` (Merqury `filt.sh` → `kmerHistToPloidyDepth`): walking the histogram
  lines, the first line whose slope turns from negative to positive gives `filt = c − 1` (0 when there is
  no such line; a zero slope never qualifies). Solid k-mers have count > filt.
- **Read error rate**: E = Σ_{c ≤ filt} c·h(c), T = Σ c·h(c); `errorRate = 1 − (1 − E/T)^(1/k)`
  (closed form of GenomeScope's `uniroot`).
- **Genome size**: `round((T − E) / peak)`, peak = histogram mode above filt (smallest count on ties) or
  `expectedCoverage`. No GenomeScope model fit: for a heterozygous genome pass the homozygous peak.
- **With `assembly`** (Merqury `spectra-cn.sh`): assembly k-mers counted per sequence (ACGT windows,
  canonical when `canonical`); `completeness = |asm ∩ solid| / |solid|` (fraction; Merqury prints ×100);
  `qualityValue = −10·log10(1 − (1 − ASM_ONLY/TOTAL)^(1/k))` with ASM_ONLY = assembly k-mer instances
  absent from the read spectrum.

Without an assembly, `completeness` and the assembly fields are `null` (k-mer completeness needs an
assembly).

## Core Documentation Reference

- Source: [GenomeAssemblyAnalyzer.cs#L1033](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/GenomeAssemblyAnalyzer.cs#L1033)
- Evidence: [ASSEMBLY-COMPLETE-001](../../../Evidence/ASSEMBLY-COMPLETE-001-Evidence.md)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `kmerSpectrum` | array | Yes | `{ kmer, count }`: each k-mer once, one length k, count ≥ 1; canonical when an assembly is given with `canonical = true`. |
| `expectedCoverage` | integer | No | Peak coverage; 0 detects the peak (default 0, ≥ 0). |
| `assembly` | array | No | Assembly `{ id, sequence }` for completeness and QV. |
| `canonical` | boolean | No | Count assembly k-mers canonically, as meryl (default true). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `completeness` | number / null | Merqury k-mer completeness (fraction); null without assembly. |
| `errorRate` | number | Read error rate 1 − (1 − E/T)^(1/k). |
| `estimatedGenomeSize` | integer | round((T − E) / peak). |
| `k`, `solidThreshold`, `peakCoverage` | integer | k, Merqury filt, peak used. |
| `totalKmers`, `errorKmers`, `solidKmers` | integer | T, E, distinct solid k-mers. |
| `assemblySolidKmers`, `assemblyKmers`, `assemblyOnlyKmers` | integer / null | Merqury intersect / TOTAL / ASM_ONLY. |
| `assemblyErrorRate`, `qualityValue` | number / null | Merqury E and QV (QV null when unbounded, i.e. ASM_ONLY = 0). |

## Errors

| Code | Message |
|------|---------|
| 1001 | K-mer spectrum cannot be null |
| 1002 | Expected coverage cannot be negative |
| 1003 | Invalid spectrum (empty, duplicate k-mer, mixed lengths, count < 1, or non-canonical / non-ACGT k-mer with an assembly) |

## Example

20 canonical 3-mers with histogram 1:5, 2:2, 3:1, 4:3, 5:6, 6:3 and assembly `ACGTTGCAAC`, `GGNNCATTAC` →
`solidThreshold 3, errorRate 0.0589640, estimatedGenomeSize 12, completeness 2/12, qualityValue 8.98186`.

## References

- [GenomeAssemblyAnalyzer.EstimateCompletenessFromKmers](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/GenomeAssemblyAnalyzer.cs#L1033)
- Rhie A et al. (2020) Merqury. Genome Biol 21:245.
- Vurture GW et al. (2017) GenomeScope. Bioinformatics 33:2202.
