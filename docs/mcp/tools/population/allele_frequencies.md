# allele_frequencies

Calculate major and minor allele frequencies from diploid genotype counts.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Population |
| **Tool Name** | `allele_frequencies` |
| **Method ID** | `PopulationGeneticsAnalyzer.CalculateAlleleFrequencies` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Computes the two allele frequencies at a biallelic diploid locus from observed genotype counts. For genotype counts n(AA), n(Aa), n(aa) the frequencies are

- p (major) = (2·n(AA) + n(Aa)) / (2·N)
- q (minor) = (2·n(aa) + n(Aa)) / (2·N)

where N = n(AA) + n(Aa) + n(aa). By construction p + q = 1 for any non-empty sample. When every genotype count is zero the tool returns (0, 0).

## Core Documentation Reference

- Source: [PopulationGeneticsAnalyzer.cs#L138](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Population/PopulationGeneticsAnalyzer.cs#L138)
- Algorithm doc: [Allele_Frequency.md](../../../algorithms/Population_Genetics/Allele_Frequency.md)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `homozygousMajor` | integer | Yes | Count of homozygous-major genotypes (AA); must be ≥ 0 |
| `heterozygous` | integer | Yes | Count of heterozygous genotypes (Aa); must be ≥ 0 |
| `homozygousMinor` | integer | Yes | Count of homozygous-minor genotypes (aa); must be ≥ 0 |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `majorFreq` | number | Major allele frequency p, in [0, 1] |
| `minorFreq` | number | Minor allele frequency q, in [0, 1] |

## Errors

| Code | Message |
|------|---------|
| 1001 | Genotype count cannot be negative. |

## Examples

### Example 1: Wikipedia four-o'clock flower example

**User Prompt:**
> A sample of 100 four-o'clock plants has 49 red (AA), 42 pink (Aa), and 9 white (aa). What are the allele frequencies?

**Expected Tool Call:**
```json
{
  "tool": "allele_frequencies",
  "arguments": { "homozygousMajor": 49, "heterozygous": 42, "homozygousMinor": 9 }
}
```

**Response:**
```json
{ "majorFreq": 0.70, "minorFreq": 0.30 }
```

`p = (2×49 + 42)/200 = 140/200 = 0.70`, `q = (2×9 + 42)/200 = 60/200 = 0.30`.

### Example 2: Wikipedia diploid example

**User Prompt:**
> 10 individuals: 6 AA, 3 AB, 1 BB. Allele frequencies?

**Expected Tool Call:**
```json
{
  "tool": "allele_frequencies",
  "arguments": { "homozygousMajor": 6, "heterozygous": 3, "homozygousMinor": 1 }
}
```

**Response:**
```json
{ "majorFreq": 0.75, "minorFreq": 0.25 }
```

`p = (2×6 + 3)/20 = 15/20 = 0.75`, `q = (2×1 + 3)/20 = 5/20 = 0.25`.

## Performance

- **Time Complexity:** O(1)
- **Space Complexity:** O(1)

## References

- Gillespie, J. H. (2004). *Population Genetics: A Concise Guide*, 2nd ed. Johns Hopkins University Press.
- Wikipedia contributors. [Allele frequency](https://en.wikipedia.org/wiki/Allele_frequency), [Genotype frequency](https://en.wikipedia.org/wiki/Genotype_frequency).

## See Also

- [minor_allele_frequency](minor_allele_frequency.md) — MAF from a genotype vector
- [filter_variants_by_maf](filter_variants_by_maf.md) — filter variants by MAF interval
