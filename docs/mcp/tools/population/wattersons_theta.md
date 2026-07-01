# wattersons_theta

Compute Watterson's θ estimator from segregating sites, sample size, and sequence length.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Population |
| **Tool Name** | `wattersons_theta` |
| **Method ID** | `PopulationGeneticsAnalyzer.CalculateWattersonTheta` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Computes Watterson's per-site estimator `θ_W = S / (a₁ · L)`, where `a₁ = Σ_{i=1}^{n-1} 1/i` is the (n−1)th harmonic number, S the number of segregating sites, and L the sequence length. For `n = 2`, `a₁ = 1`, so `θ_W = S / L`. Returns 0 when `n < 2` or `L ≤ 0`.

## Core Documentation Reference

- Source: [PopulationGeneticsAnalyzer.cs#L239](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Population/PopulationGeneticsAnalyzer.cs#L239)
- Algorithm doc: [Diversity_Statistics.md](../../../algorithms/Population_Genetics/Diversity_Statistics.md)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `segregatingSites` | integer | Yes | Number of segregating sites S |
| `sampleSize` | integer | Yes | Sample size n |
| `sequenceLength` | integer | Yes | Sequence length L in nucleotides |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `theta` | number | Watterson's θ per site (0 when n < 2 or L ≤ 0) |

## Errors

None (degenerate cases return 0).

## Examples

### Example 1: S = 10, n = 10, L = 1000

**User Prompt:**
> Watterson's θ for 10 segregating sites in 10 sequences of length 1000.

**Expected Tool Call:**
```json
{
  "tool": "wattersons_theta",
  "arguments": { "segregatingSites": 10, "sampleSize": 10, "sequenceLength": 1000 }
}
```

**Response:**
```json
{ "theta": 0.00353 }
```

a₁ = Σ_{i=1}^{9} 1/i ≈ 2.8290, θ = 10/(2.8290·1000) ≈ 0.00353.

### Example 2: Minimum sample size n = 2 → θ = S/L

**User Prompt:**
> Watterson's θ for S = 5, n = 2, L = 100.

**Expected Tool Call:**
```json
{
  "tool": "wattersons_theta",
  "arguments": { "segregatingSites": 5, "sampleSize": 2, "sequenceLength": 100 }
}
```

**Response:**
```json
{ "theta": 0.05 }
```

For n = 2, a₁ = 1, so θ = 5/100 = 0.05.

## Performance

- **Time Complexity:** O(n) (harmonic sum)
- **Space Complexity:** O(1)

## References

- Watterson, G. A. (1975). On the number of segregating sites in genetical models without recombination. *Theor. Popul. Biol.* 7:256–276.

## See Also

- [diversity_statistics](diversity_statistics.md), [nucleotide_diversity](nucleotide_diversity.md), [tajimas_d](tajimas_d.md)
