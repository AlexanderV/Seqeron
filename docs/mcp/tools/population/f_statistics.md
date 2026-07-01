# f_statistics

Compute Wright's F-statistics (Fis, Fit, Fst) between two populations.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Population |
| **Tool Name** | `f_statistics` |
| **Method ID** | `PopulationGeneticsAnalyzer.CalculateFStatistics` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Computes the heterozygosity-based F-statistics from per-variant observed heterozygote counts and allele frequencies. With H_I the mean observed heterozygosity, H_S the mean within-population expected heterozygosity, and H_T the total expected heterozygosity (pooled):

- Fis = 1 − H_I/H_S
- Fit = 1 − H_I/H_T
- Fst = 1 − H_S/H_T

The partition identity (1 − Fit) = (1 − Fis)(1 − Fst) holds exactly. Fst lies in [0, 1]; Fis and Fit can be negative under excess heterozygosity. Empty `variantData` returns all zeros.

## Core Documentation Reference

- Source: [PopulationGeneticsAnalyzer.cs#L674](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Population/PopulationGeneticsAnalyzer.cs#L674)
- Algorithm doc: [F_Statistics.md](../../../algorithms/Population_Genetics/F_Statistics.md)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `pop1Name` | string | Yes | Population 1 name |
| `pop2Name` | string | Yes | Population 2 name |
| `variantData` | object[] | Yes | Per-variant `{hetObs1, n1, hetObs2, n2, alleleFreq1, alleleFreq2}` |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `fst` | number | Fst = 1 − H_S/H_T, in [0, 1] |
| `fis` | number | Fis = 1 − H_I/H_S (may be negative) |
| `fit` | number | Fit = 1 − H_I/H_T (may be negative) |
| `population1` | string | Echoed population 1 name |
| `population2` | string | Echoed population 2 name |

## Errors

None (empty `variantData` returns zeros).

## Examples

### Example 1: Two-locus hand-calculated F-statistics

**User Prompt:**
> Compute F-statistics for Pop1 vs Pop2 with two loci.

**Expected Tool Call:**
```json
{
  "tool": "f_statistics",
  "arguments": {
    "pop1Name": "Pop1",
    "pop2Name": "Pop2",
    "variantData": [
      { "hetObs1": 20, "n1": 50, "hetObs2": 25, "n2": 50, "alleleFreq1": 0.4, "alleleFreq2": 0.5 },
      { "hetObs1": 30, "n1": 50, "hetObs2": 15, "n2": 50, "alleleFreq1": 0.5, "alleleFreq2": 0.3 }
    ]
  }
}
```

**Response:**
```json
{ "fis": 0.05263, "fit": 0.07692, "fst": 0.02564, "population1": "Pop1", "population2": "Pop2" }
```

H_I = 0.45, H_S = 0.475, H_T = 0.4875 → Fis = 1/19, Fit = 1/13, Fst = 1/39.

### Example 2: Excess heterozygosity → negative Fis

**User Prompt:**
> Single locus, HetObs 60/100 and 80/100 with p = 0.3 and 0.7.

**Expected Tool Call:**
```json
{
  "tool": "f_statistics",
  "arguments": {
    "pop1Name": "Pop1",
    "pop2Name": "Pop2",
    "variantData": [
      { "hetObs1": 60, "n1": 100, "hetObs2": 80, "n2": 100, "alleleFreq1": 0.3, "alleleFreq2": 0.7 }
    ]
  }
}
```

**Response:**
```json
{ "fis": -0.66667, "fit": -0.4, "fst": 0.16, "population1": "Pop1", "population2": "Pop2" }
```

Fis = −2/3, Fit = −2/5, Fst = 4/25.

## Performance

- **Time Complexity:** O(loci)
- **Space Complexity:** O(1)

## References

- Wright, S. (1965). The interpretation of population structure by F-statistics. *Evolution* 19:395–420.
- Wikipedia contributors. [F-statistics](https://en.wikipedia.org/wiki/F-statistics).

## See Also

- [fst](fst.md), [pairwise_fst](pairwise_fst.md)
