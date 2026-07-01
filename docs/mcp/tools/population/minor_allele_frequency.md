# minor_allele_frequency

Compute minor allele frequency (MAF) from a diploid genotype vector.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Population |
| **Tool Name** | `minor_allele_frequency` |
| **Method ID** | `PopulationGeneticsAnalyzer.CalculateMAF` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Sums alternate-allele copies from the `0/1/2` genotype encoding (0 = hom-ref, 1 = het, 2 = hom-alt), forms the alternate-allele frequency `altFreq = Σgenotypes / (2n)`, and returns the folded minor allele frequency `MAF = min(altFreq, 1 − altFreq)`. MAF is always in `[0, 0.5]`. Empty input and monomorphic loci return 0.

## Core Documentation Reference

- Source: [PopulationGeneticsAnalyzer.cs#L164](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Population/PopulationGeneticsAnalyzer.cs#L164)
- Algorithm doc: [Allele_Frequency.md](../../../algorithms/Population_Genetics/Allele_Frequency.md)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `genotypes` | integer[] | Yes | Genotype vector with values in {0, 1, 2} |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `maf` | number | Minor allele frequency in [0, 0.5] |

## Errors

None (empty input returns 0).

## Examples

### Example 1: Alt frequency below 0.5

**User Prompt:**
> MAF for genotypes [0,0,0,0,0,1,1,1,1,2]?

**Expected Tool Call:**
```json
{
  "tool": "minor_allele_frequency",
  "arguments": { "genotypes": [0, 0, 0, 0, 0, 1, 1, 1, 1, 2] }
}
```

**Response:**
```json
{ "maf": 0.3 }
```

Alt alleles = 4·1 + 1·2 = 6, total = 20, altFreq = 0.3, MAF = min(0.3, 0.7) = 0.3.

### Example 2: Balanced polymorphism (MAF = 0.5)

**User Prompt:**
> MAF for [0,1,1,2]?

**Expected Tool Call:**
```json
{
  "tool": "minor_allele_frequency",
  "arguments": { "genotypes": [0, 1, 1, 2] }
}
```

**Response:**
```json
{ "maf": 0.5 }
```

Alt alleles = 4, total = 8, altFreq = 0.5, MAF = 0.5.

## Performance

- **Time Complexity:** O(n)
- **Space Complexity:** O(n)

## References

- Wikipedia contributors. [Minor allele frequency](https://en.wikipedia.org/wiki/Minor_allele_frequency).

## See Also

- [allele_frequencies](allele_frequencies.md), [filter_variants_by_maf](filter_variants_by_maf.md)
