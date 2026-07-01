# fst

Compute Wright's variance-based Fst between two populations.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Population |
| **Tool Name** | `fst` |
| **Method ID** | `PopulationGeneticsAnalyzer.CalculateFst` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Computes Wright's fixation index `Fst = σ²_S / (p̄(1−p̄))`, summed over loci, where `σ²_S` is the sample-size-weighted variance of allele frequencies between the two populations and `p̄(1−p̄)` is the expected heterozygosity of the pooled allele frequency. Fst is 0 for identical populations (panmixia) and 1 for fixed differences (complete differentiation). Empty inputs return 0. The two populations must have the same per-locus count; a mismatch throws an error.

## Core Documentation Reference

- Source: [PopulationGeneticsAnalyzer.cs#L609](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Population/PopulationGeneticsAnalyzer.cs#L609)
- Algorithm doc: [F_Statistics.md](../../../algorithms/Population_Genetics/F_Statistics.md)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `population1` | object[] | Yes | Per-variant `{alleleFreq, sampleSize}` for population 1 |
| `population2` | object[] | Yes | Per-variant `{alleleFreq, sampleSize}` for population 2 (same locus count) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `fst` | number | Wright's Fst, in [0, 1] |

## Errors

| Code | Message |
|------|---------|
| 1001 | The two populations' per-locus allele frequency counts must match. |

## Examples

### Example 1: Single locus p1 = 0.8, p2 = 0.2

**User Prompt:**
> Fst between two populations with alt-allele frequency 0.8 and 0.2 (100 samples each)?

**Expected Tool Call:**
```json
{
  "tool": "fst",
  "arguments": {
    "population1": [{ "alleleFreq": 0.8, "sampleSize": 100 }],
    "population2": [{ "alleleFreq": 0.2, "sampleSize": 100 }]
  }
}
```

**Response:**
```json
{ "fst": 0.36 }
```

p̄ = 0.5, variance = 0.09, het = 0.25 → Fst = 0.09/0.25 = 0.36.

### Example 2: Fixed differences

**User Prompt:**
> Fst when one population is fixed for the allele and the other lacks it entirely?

**Expected Tool Call:**
```json
{
  "tool": "fst",
  "arguments": {
    "population1": [{ "alleleFreq": 1.0, "sampleSize": 100 }],
    "population2": [{ "alleleFreq": 0.0, "sampleSize": 100 }]
  }
}
```

**Response:**
```json
{ "fst": 1.0 }
```

Complete differentiation.

## Performance

- **Time Complexity:** O(loci)
- **Space Complexity:** O(1)

## References

- Wright, S. (1965). The interpretation of population structure by F-statistics. *Evolution* 19:395–420.
- Wikipedia contributors. [Fixation index](https://en.wikipedia.org/wiki/Fixation_index).

## See Also

- [pairwise_fst](pairwise_fst.md), [f_statistics](f_statistics.md)
