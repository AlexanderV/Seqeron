# diversity_statistics

Compute combined nucleotide-diversity statistics from a set of aligned sequences.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Population |
| **Tool Name** | `diversity_statistics` |
| **Method ID** | `PopulationGeneticsAnalyzer.CalculateDiversityStatistics` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Given an aligned set of equal-length sequences, computes in one pass:

- **π** (nucleotide diversity) — average per-site pairwise differences,
- **Watterson's θ** = S / (a₁·L), a₁ = Σ_{i=1}^{n-1} 1/i,
- **Tajima's D** from k̂ = π·L, S and n (Tajima 1989),
- **S** — number of segregating (polymorphic) sites,
- **n** — sample size,
- **observed heterozygosity** — Nei (1978) unbiased gene diversity per site,
- **expected heterozygosity** — basic gene diversity (1 − Σp²) per site.

Fewer than 2 sequences yields all zeros (with `sampleSize` equal to the input count).

## Core Documentation Reference

- Source: [PopulationGeneticsAnalyzer.cs#L309](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Population/PopulationGeneticsAnalyzer.cs#L309)
- Algorithm doc: [Diversity_Statistics.md](../../../algorithms/Population_Genetics/Diversity_Statistics.md)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequences` | string[] | Yes | Aligned sequences of equal length |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `nucleotideDiversity` | number | π, average per-site pairwise differences |
| `wattersonTheta` | number | Watterson's θ per site |
| `tajimasD` | number | Tajima's D |
| `segregatingSites` | integer | Number of polymorphic sites S |
| `sampleSize` | integer | Number of sequences n |
| `heterozygosityObserved` | number | Nei (1978) unbiased gene diversity per site |
| `heterozygosityExpected` | number | Basic gene diversity (1 − Σp²) per site |

## Errors

None (empty or single-sequence input returns zeros rather than throwing).

## Examples

### Example 1: Wikipedia Tajima's D example

**User Prompt:**
> Compute diversity statistics for these 5 aligned sequences (length 20).

**Expected Tool Call:**
```json
{
  "tool": "diversity_statistics",
  "arguments": {
    "sequences": [
      "00000000000000000000",
      "00100000000010000010",
      "00000000000010000010",
      "00000010000000000010",
      "00000010000010000010"
    ]
  }
}
```

**Response (key fields):**
```json
{
  "nucleotideDiversity": 0.1,
  "wattersonTheta": 0.096,
  "tajimasD": 0.273,
  "segregatingSites": 4,
  "sampleSize": 5
}
```

π = k̂/L = 2.0/20 = 0.1; θ_W = 4/(a₁·20) with a₁ = 25/12 ≈ 0.096; D ≈ 0.273 (Wikipedia hand-calculation).

### Example 2: Three sequences, length 8

**User Prompt:**
> Diversity statistics for ACGTACGT / ACGTATGT / ACGTACGA.

**Expected Tool Call:**
```json
{
  "tool": "diversity_statistics",
  "arguments": { "sequences": ["ACGTACGT", "ACGTATGT", "ACGTACGA"] }
}
```

**Response:**
```json
{
  "nucleotideDiversity": 0.16667,
  "wattersonTheta": 0.16667,
  "tajimasD": 0.0,
  "segregatingSites": 2,
  "sampleSize": 3,
  "heterozygosityObserved": 0.16667,
  "heterozygosityExpected": 0.11111
}
```

S = 2 (positions 5, 7); k̂ = 4/3, π = 4/(3·8) = 1/6; a₁(3) = 3/2, θ_W = 2/(1.5·8) = 1/6; k̂ = S/a₁ so D = 0; H_exp = 1/9, H_obs = (3/2)·(1/9) = 1/6.

## Performance

- **Time Complexity:** O(n²·L) for the pairwise π scan
- **Space Complexity:** O(n·L)

## References

- Tajima, F. (1989). Statistical method for testing the neutral mutation hypothesis by DNA polymorphism. *Genetics* 123:585–595.
- Watterson, G. A. (1975). On the number of segregating sites. *Theor. Popul. Biol.* 7:256–276.
- Nei, M. (1978). Estimation of average heterozygosity and genetic distance. *Genetics* 89:583–590.

## See Also

- [nucleotide_diversity](nucleotide_diversity.md), [wattersons_theta](wattersons_theta.md), [tajimas_d](tajimas_d.md)
