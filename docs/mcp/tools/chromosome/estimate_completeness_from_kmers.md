# estimate_completeness_from_kmers

Estimate completeness, error rate and genome size from a k-mer spectrum.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Chromosome |
| **Tool Name** | `estimate_completeness_from_kmers` |
| **Method ID** | `GenomeAssemblyAnalyzer.EstimateCompletenessFromKmers` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

From a k-mer count spectrum, separates low-count error k-mers from the coverage peak (auto-detected
when `expectedCoverage = 0`), then estimates completeness, error rate and genome size. A clean spectrum
whose k-mers all sit at the expected coverage has completeness 1 and error rate 0.

## Core Documentation Reference

- Source: [GenomeAssemblyAnalyzer.cs#L628](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/GenomeAssemblyAnalyzer.cs#L628)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `kmerSpectrum` | array | Yes | `{ kmer, count }` entries. |
| `expectedCoverage` | integer | No | Expected coverage; 0 auto-detects the peak (default 0, ≥ 0). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `completeness` | number | Estimated completeness. |
| `errorRate` | number | Estimated error rate. |
| `estimatedGenomeSize` | integer | Estimated genome size. |

## Errors

| Code | Message |
|------|---------|
| 1001 | K-mer spectrum cannot be null |
| 1002 | Expected coverage cannot be negative |

## Example

100 distinct k-mers all at coverage 30, `expectedCoverage = 30` →
`{ completeness: 1.0, errorRate: 0.0, estimatedGenomeSize: 100 }`.

## References

- [GenomeAssemblyAnalyzer.EstimateCompletenessFromKmers](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/GenomeAssemblyAnalyzer.cs#L628)
