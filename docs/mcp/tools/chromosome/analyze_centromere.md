# analyze_centromere

Locate the centromere region of a chromosome sequence and classify its type.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Chromosome |
| **Tool Name** | `analyze_centromere` |
| **Method ID** | `ChromosomeAnalyzer.AnalyzeCentromere` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Scans a chromosome sequence in overlapping windows for alpha-satellite-like content (high k-mer
repeat content combined with low GC variability), anchors the best-scoring window, extends its
boundaries while repeat content stays above 70% of the threshold, and classifies the centromere
position by arm ratio (q/p) per Levan et al. (1964):

| Arm ratio (q/p) | Type |
|-----------------|------|
| p = 0 | Telocentric |
| ≤ 1.7 | Metacentric |
| 1.7 – 3.0 | Submetacentric |
| 3.0 – 7.0 | Subtelocentric |
| ≥ 7.0 | Acrocentric |

If no window exceeds `minAlphaSatelliteContent`, the type is `Unknown` and `start`/`end` are null.

## Core Documentation Reference

- Source: [ChromosomeAnalyzer.cs#L504](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/ChromosomeAnalyzer.cs#L504)
- Reference: Levan A, Fredga K, Sandberg AA (1964). "Nomenclature for centromeric position on chromosomes." *Hereditas* 52(2):201–220.

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `chromosomeName` | string | Yes | Chromosome name (echoed in result). |
| `sequence` | string | Yes | Chromosome nucleotide sequence. |
| `windowSize` | integer | No | Scan window size in bp (default 100000, must be > 0). |
| `minAlphaSatelliteContent` | number | No | Minimum alpha-satellite-like content in [0, 1] (default 0.3). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `chromosome` | string | Chromosome name echoed from input. |
| `start` | integer \| null | Centromere start (bp); null if not detected. |
| `end` | integer \| null | Centromere end (bp); null if not detected. |
| `length` | integer | Centromere length (bp) = end − start; 0 if not detected. |
| `centromereType` | string | Metacentric / Submetacentric / Subtelocentric / Acrocentric / Telocentric / Unknown. |
| `alphaSatelliteContent` | number | Repeat-content score of the best window (≥ 0). |
| `isAcrocentric` | boolean | True when `centromereType` is Acrocentric. |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1002 | Chromosome name cannot be null or empty |
| 1003 | Window size must be positive |
| 1004 | Minimum alpha-satellite content must be in [0, 1] |

## Examples

### Example 1: Uniform highly-repetitive sequence

**User Prompt:**
> Find the centromere in this 300 kb poly-A chromosome (10 kb windows).

**Expected Tool Call:**
```json
{
  "tool": "analyze_centromere",
  "arguments": { "chromosomeName": "chr1", "sequence": "AAAA…(300000×A)", "windowSize": 10000, "minAlphaSatelliteContent": 0.3 }
}
```

**Response:**
```json
{
  "chromosome": "chr1",
  "start": 0,
  "end": 295000,
  "length": 295000,
  "centromereType": "Metacentric",
  "alphaSatelliteContent": 1.0,
  "isAcrocentric": false
}
```

### Example 2: Random non-repetitive sequence

**User Prompt:**
> Is there a centromere in this random 500 kb sequence?

**Response:**
```json
{
  "chromosome": "chr1",
  "start": null,
  "end": null,
  "length": 0,
  "centromereType": "Unknown",
  "alphaSatelliteContent": 0.0,
  "isAcrocentric": false
}
```

## References

- [ChromosomeAnalyzer.AnalyzeCentromere](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/ChromosomeAnalyzer.cs#L504)
- Levan et al. (1964), *Hereditas* 52(2):201–220.
