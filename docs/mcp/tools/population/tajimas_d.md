# tajimas_d

Compute Tajima's D from k̂, segregating sites, and sample size.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Population |
| **Tool Name** | `tajimas_d` |
| **Method ID** | `PopulationGeneticsAnalyzer.CalculateTajimasD` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Computes `D = (k̂ − S/a₁) / √(e₁·S + e₂·S·(S−1))` (Tajima 1989), where k̂ is the average number of pairwise differences (NOT per-site), S is the number of segregating sites, and a₁ = Σ_{i=1}^{n-1} 1/i. Returns 0 when S = 0, n < 3, or the estimated variance is non-positive. Negative D indicates an excess of rare variants (population expansion / purifying selection); positive D indicates a deficit (balancing selection / contraction).

## Core Documentation Reference

- Source: [PopulationGeneticsAnalyzer.cs#L266](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Population/PopulationGeneticsAnalyzer.cs#L266)
- Algorithm doc: [Diversity_Statistics.md](../../../algorithms/Population_Genetics/Diversity_Statistics.md)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `averagePairwiseDifferences` | number | Yes | k̂ — average pairwise differences (NOT per-site) |
| `segregatingSites` | integer | Yes | Number of segregating sites S |
| `sampleSize` | integer | Yes | Sample size n (requires n ≥ 3 for a non-zero result) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `tajimasD` | number | Tajima's D (0 for S = 0, n < 3, or variance ≤ 0) |

## Errors

None (degenerate cases return 0).

## Examples

### Example 1: Wikipedia example

**User Prompt:**
> Tajima's D for k̂ = 2.0, S = 4, n = 5.

**Expected Tool Call:**
```json
{
  "tool": "tajimas_d",
  "arguments": { "averagePairwiseDifferences": 2.0, "segregatingSites": 4, "sampleSize": 5 }
}
```

**Response:**
```json
{ "tajimasD": 0.273 }
```

a₁ = 25/12, S/a₁ ≈ 1.92, d = 0.08, Var ≈ 0.0856 → D ≈ 0.273 (Wikipedia hand-calculation).

### Example 2: Neutral (k̂ = S/a₁) → D = 0

**User Prompt:**
> Tajima's D when k̂ equals the Watterson expectation.

**Expected Tool Call:**
```json
{
  "tool": "tajimas_d",
  "arguments": { "averagePairwiseDifferences": 22.2, "segregatingSites": 100, "sampleSize": 50 }
}
```

**Response:**
```json
{ "tajimasD": 0.0 }
```

When k̂ = S/a₁ the numerator is 0, so D = 0.

## Performance

- **Time Complexity:** O(n) (harmonic sums)
- **Space Complexity:** O(1)

## References

- Tajima, F. (1989). Statistical method for testing the neutral mutation hypothesis by DNA polymorphism. *Genetics* 123:585–595.

## See Also

- [diversity_statistics](diversity_statistics.md), [nucleotide_diversity](nucleotide_diversity.md), [wattersons_theta](wattersons_theta.md)
