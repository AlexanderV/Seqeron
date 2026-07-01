# differential_abundance

Test for differential taxon abundance between two condition groups.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Metagenomics |
| **Tool Name** | `differential_abundance` |
| **Method ID** | `MetagenomicsAnalyzer.DifferentialAbundance` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

For every taxon present in either condition, computes the mean abundance in each group, the
**log2 fold-change** `log2(mean2 / mean1)` (reported in `foldChange`), and a Welch's t-test
p-value (normal approximation). A taxon is flagged **significant** when
`p < pValueThreshold` **and** `|log2FC| > 1`.

## Core Documentation Reference

- Source: [MetagenomicsAnalyzer.cs#L1594](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Metagenomics/MetagenomicsAnalyzer.cs#L1594)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `condition1Samples` | AbundanceSample[] | Yes | Per-sample vectors for condition 1. |
| `condition2Samples` | AbundanceSample[] | Yes | Per-sample vectors for condition 2. |
| `pValueThreshold` | number | No | Significance threshold (default 0.05). |

Each `AbundanceSample` wraps an `items` array of `{ name, fraction }`.

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `items[].taxon` | string | Taxon name. |
| `items[].foldChange` | number | log2 fold-change `log2(mean2/mean1)`. |
| `items[].pValue` | number | Welch's t-test p-value. |
| `items[].significant` | boolean | `p < threshold` AND `|log2FC| > 1`. |

## Errors

None. An empty condition yields no items.

## Example

**User Prompt:**
> Is taxon T differentially abundant between my two conditions (1,1,1 vs 4,4,4)?

**Response:**
```json
{ "items": [ { "taxon": "T", "foldChange": 2.0, "pValue": 0.0, "significant": true } ] }
```

log2(4/1) = 2; with zero within-group variance and differing means the Welch p is 0, so the
taxon is significant.

## References

- [MetagenomicsAnalyzer.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Metagenomics/MetagenomicsAnalyzer.cs) — `DifferentialAbundance`
- Welch B.L. (1947) Biometrika 34:28.

## See Also

- [beta_diversity](beta_diversity.md)
- [taxonomic_profile](taxonomic_profile.md)
