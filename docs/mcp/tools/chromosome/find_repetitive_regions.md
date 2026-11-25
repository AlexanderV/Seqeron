# find_repetitive_regions

Identify repetitive regions with the NCBI WindowMasker k-mer-count model.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Chromosome |
| **Tool Name** | `find_repetitive_regions` |
| **Method ID** | `GenomeAssemblyAnalyzer.FindRepetitiveRegions` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Port of NCBI WindowMasker's first masking pass (Morgulis et al. 2006; `CSeqMasker::DoMask`). k-mers
(units) are counted over all sequences pooled, on both strands by default (`count(u) = F(u) + F(rc(u))`,
so a self-complementary unit counts twice per occurrence, as in `win_mask_gen_counts.cpp`); a k-mer
containing any non-ACGT symbol is not counted and a window containing one is skipped. A window of
`windowSize` bases (step 1) scores the integer mean of its unit counts; score ≥ `minCopies`
(T_threshold) starts or extends a region, score < `extendCopies` (T_extend, default `minCopies`) closes
it. Regions are 0-based inclusive; `copies` is the largest unit count inside the region.
T_threshold/T_extend are counts here (WindowMasker derives them from percentiles); the optional merge
pass and discontiguous units are not implemented.

## Core Documentation Reference

- Source: [GenomeAssemblyAnalyzer.cs#L1214](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/GenomeAssemblyAnalyzer.cs#L1214)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequences` | array | Yes | Sequences `{ id, sequence }`. |
| `kmerSize` | integer | No | K-mer size (default 15, > 0). |
| `minCopies` | integer | No | Minimum k-mer copies (default 3, > 0). |
| `windowSize` | integer | No | Window length in bases (default 100, ≥ `kmerSize`). |
| `bothStrands` | boolean | No | Count units on both strands (default true, WindowMasker). |
| `extendCopies` | integer | No | T_extend, 1..`minCopies` (default = `minCopies`). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `items[]` | array | `{ sequenceId, start, end, copies }`. |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequences cannot be null |
| 1002 | K-mer size must be positive |
| 1003 | Minimum copies must be positive |
| 1004 | Window size must be positive |

## Example

`"ACGTACGT"` × 100 (800 bp), `kmerSize = 8`, `windowSize = 50` → one region `[0, 799]`; ACGTACGT is its own reverse complement (199 occurrences × 2 strands), so `copies = 398`.

## References

- Morgulis A. et al. (2006) WindowMasker, Bioinformatics 22:134–141; NCBI C++ Toolkit `src/algo/winmask`.
- [GenomeAssemblyAnalyzer.FindRepetitiveRegions](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/GenomeAssemblyAnalyzer.cs#L1214)
