# functional_diversity

Compute functional richness and Shannon functional diversity.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Metagenomics |
| **Tool Name** | `functional_diversity` |
| **Method ID** | `MetagenomicsAnalyzer.CalculateFunctionalDiversity` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

From a set of functional annotations, reports:

- **functionalRichness** — the number of distinct non-empty function labels
- **functionalDiversity** — the Shannon entropy `-Σ pᵢ ln pᵢ` of the function-count distribution
- **pathwayCounts** — the number of annotations per non-empty pathway

## Core Documentation Reference

- Source: [MetagenomicsAnalyzer.cs#L1352](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Metagenomics/MetagenomicsAnalyzer.cs#L1352)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `annotations` | FunctionalAnnotation[] | Yes | Annotations (e.g. from `predict_functions`). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `functionalRichness` | number | Distinct function count. |
| `functionalDiversity` | number | Shannon entropy of function counts. |
| `pathwayCounts[]` | array | `{ pathway, count }` per pathway. |

## Errors

None. An empty annotation set yields richness 0, diversity 0, and no pathway counts.

## Example

**User Prompt:**
> Compute functional diversity for annotations F1, F1, F2.

**Response:**
```json
{ "functionalRichness": 2, "functionalDiversity": 0.6365141682948128, "pathwayCounts": [ { "pathway": "P1", "count": 2 }, { "pathway": "P2", "count": 1 } ] }
```

Richness = 2 distinct functions; Shannon = −(2/3·ln 2/3 + 1/3·ln 1/3).

## References

- [MetagenomicsAnalyzer.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Metagenomics/MetagenomicsAnalyzer.cs) — `CalculateFunctionalDiversity`
- Shannon C.E. (1948).

## See Also

- [predict_functions](predict_functions.md)
