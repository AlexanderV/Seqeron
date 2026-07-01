# pairwise_fst

Compute the symmetric pairwise Wright's Fst matrix for a set of populations.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Population |
| **Tool Name** | `pairwise_fst` |
| **Method ID** | `PopulationGeneticsAnalyzer.CalculatePairwiseFst` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Runs `fst` (Wright's variance-based Fst) on every population pair and returns the population id list together with the symmetric N×N Fst matrix. The diagonal is 0 (self-comparison) and `matrix[i][j] == matrix[j][i]`.

## Core Documentation Reference

- Source: [PopulationGeneticsAnalyzer.cs#L651](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Population/PopulationGeneticsAnalyzer.cs#L651)
- Algorithm doc: [F_Statistics.md](../../../algorithms/Population_Genetics/F_Statistics.md)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `populations` | object[] | Yes | Populations with `populationId` and per-variant `{alleleFreq, sampleSize}` |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `populationIds` | string[] | Population ids in matrix order |
| `matrix` | number[][] | Symmetric N×N Fst matrix (zero diagonal) |

## Errors

None.

## Examples

### Example 1: Three single-locus populations

**User Prompt:**
> Pairwise Fst matrix for populations with allele frequency 0.5, 0.6, 0.9 (100 samples each).

**Expected Tool Call:**
```json
{
  "tool": "pairwise_fst",
  "arguments": {
    "populations": [
      { "populationId": "Pop1", "variants": [{ "alleleFreq": 0.5, "sampleSize": 100 }] },
      { "populationId": "Pop2", "variants": [{ "alleleFreq": 0.6, "sampleSize": 100 }] },
      { "populationId": "Pop3", "variants": [{ "alleleFreq": 0.9, "sampleSize": 100 }] }
    ]
  }
}
```

**Response (key cells):** `matrix[0][1] = 1/99 ≈ 0.0101`, `matrix[0][2] = 4/21 ≈ 0.1905`, `matrix[1][2] = 3/25 = 0.12`; diagonal all 0.

### Example 2: Two identical populations → zero matrix

**User Prompt:**
> Pairwise Fst for two identical populations.

**Expected Tool Call:**
```json
{
  "tool": "pairwise_fst",
  "arguments": {
    "populations": [
      { "populationId": "A", "variants": [{ "alleleFreq": 0.4, "sampleSize": 100 }] },
      { "populationId": "B", "variants": [{ "alleleFreq": 0.4, "sampleSize": 100 }] }
    ]
  }
}
```

**Response:** `matrix = [[0, 0], [0, 0]]`.

## Performance

- **Time Complexity:** O(N²·loci)
- **Space Complexity:** O(N²)

## References

- Wright, S. (1965). The interpretation of population structure by F-statistics. *Evolution* 19:395–420.

## See Also

- [fst](fst.md), [f_statistics](f_statistics.md)
