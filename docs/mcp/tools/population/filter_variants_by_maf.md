# filter_variants_by_maf

Filter variants by an inclusive minor-allele-frequency interval.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Population |
| **Tool Name** | `filter_variants_by_maf` |
| **Method ID** | `PopulationGeneticsAnalyzer.FilterByMAF` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Keeps only the variants whose folded minor allele frequency `MAF = min(alleleFrequency, 1 − alleleFrequency)` falls inside the inclusive interval `[minMAF, maxMAF]`. Both bounds are inclusive (compared with `>=` and `<=`). The filter operates on the stored `alleleFrequency` field (it does not recompute frequencies from genotypes) and preserves input order.

## Core Documentation Reference

- Source: [PopulationGeneticsAnalyzer.cs#L182](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Population/PopulationGeneticsAnalyzer.cs#L182)
- Algorithm doc: [Allele_Frequency.md](../../../algorithms/Population_Genetics/Allele_Frequency.md)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `variants` | object[] | Yes | Variants (`id`, `chromosome`, `position`, `referenceAllele`, `alternateAllele`, `alleleFrequency`, `sampleCount`) |
| `minMAF` | number | No | Inclusive lower bound on MAF (default 0.01) |
| `maxMAF` | number | No | Inclusive upper bound on MAF (default 0.5) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `items` | object[] | Variants passing the filter, in input order |

## Errors

None.

## Examples

### Example 1: Drop a rare variant below minMAF

**User Prompt:**
> Filter out variants with MAF below 0.01.

**Expected Tool Call:**
```json
{
  "tool": "filter_variants_by_maf",
  "arguments": {
    "variants": [
      { "id": "V1", "chromosome": "chr1", "position": 100, "referenceAllele": "A", "alternateAllele": "G", "alleleFrequency": 0.005, "sampleCount": 100 },
      { "id": "V2", "chromosome": "chr1", "position": 200, "referenceAllele": "A", "alternateAllele": "G", "alleleFrequency": 0.05, "sampleCount": 100 }
    ],
    "minMAF": 0.01
  }
}
```

**Response:** only `V2` (MAF 0.005 < 0.01 is excluded).

### Example 2: High allele frequency folds to MAF = 0.05

**User Prompt:**
> Does a variant with alleleFrequency 0.95 pass a [0.01, 0.1] MAF filter?

**Expected Tool Call:**
```json
{
  "tool": "filter_variants_by_maf",
  "arguments": {
    "variants": [
      { "id": "V1", "chromosome": "chr1", "position": 100, "referenceAllele": "A", "alternateAllele": "G", "alleleFrequency": 0.95, "sampleCount": 100 }
    ],
    "minMAF": 0.01,
    "maxMAF": 0.1
  }
}
```

**Response:** `V1` passes — MAF = min(0.95, 0.05) = 0.05 ∈ [0.01, 0.1].

## Performance

- **Time Complexity:** O(n)
- **Space Complexity:** O(1) incremental

## References

- Wikipedia contributors. [Minor allele frequency](https://en.wikipedia.org/wiki/Minor_allele_frequency).

## See Also

- [minor_allele_frequency](minor_allele_frequency.md), [allele_frequencies](allele_frequencies.md)
