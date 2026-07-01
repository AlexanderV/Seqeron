# alpha_diversity

Compute alpha-diversity metrics for a single sample.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Metagenomics |
| **Tool Name** | `alpha_diversity` |
| **Method ID** | `MetagenomicsAnalyzer.CalculateAlphaDiversity` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Given a per-species abundance vector, computes:

- **Shannon index** `H = -Σ pᵢ ln pᵢ` (Shannon 1948)
- **Simpson index** `λ = Σ pᵢ²` (Simpson 1949)
- **Inverse Simpson** `1/λ` (Hill 1973 number of order 2)
- **Chao1** richness estimate (Chao 1984); requires integer count data, otherwise falls
  back to the observed species count
- **Observed species** — species with positive abundance
- **Pielou's evenness** `J = H / ln(S)` (0 when `S ≤ 1`)

Abundances are normalised internally, so raw counts or fractions are both accepted.
Zero-abundance entries are dropped.

## Core Documentation Reference

- Source: [MetagenomicsAnalyzer.cs#L480](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Metagenomics/MetagenomicsAnalyzer.cs#L480)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `abundances` | AbundanceItem[] | Yes | Per-species `{ name, fraction }` entries. |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `shannonIndex` | number | Shannon entropy H. |
| `simpsonIndex` | number | Simpson concentration λ. |
| `inverseSimpson` | number | 1 / λ. |
| `chao1Estimate` | number | Chao1 richness estimate or S_obs. |
| `observedSpecies` | number | Count of species with positive abundance. |
| `pielouEvenness` | number | Evenness J = H / ln(S). |

## Errors

None. An empty vector yields all-zero metrics.

## Example

**User Prompt:**
> Compute alpha diversity for two equally abundant species.

Input `[{ "name": "s1", "fraction": 0.5 }, { "name": "s2", "fraction": 0.5 }]`:

**Response:**
```json
{
  "shannonIndex": 0.6931471805599453,
  "simpsonIndex": 0.5,
  "inverseSimpson": 2.0,
  "chao1Estimate": 2.0,
  "observedSpecies": 2,
  "pielouEvenness": 1.0
}
```

Shannon = ln 2, Simpson = 2·0.5² = 0.5, inverse Simpson = 2, evenness = 1 (perfectly even).

## References

- [MetagenomicsAnalyzer.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Metagenomics/MetagenomicsAnalyzer.cs) — `CalculateAlphaDiversity`
- Shannon C.E. (1948); Simpson E.H. (1949); Hill M.O. (1973); Chao A. (1984) Scand. J. Stat. 11:265.

## See Also

- [beta_diversity](beta_diversity.md) — pairwise between-sample distances
- [taxonomic_profile](taxonomic_profile.md) — profile with species-level Shannon/Simpson
