# taxonomic_profile

Aggregate per-read classifications into a sample-level taxonomic profile.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Metagenomics |
| **Tool Name** | `taxonomic_profile` |
| **Method ID** | `MetagenomicsAnalyzer.GenerateTaxonomicProfile` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Aggregates the per-read classifications from [classify_reads](classify_reads.md) into relative
abundances at kingdom, phylum, genus, and species level, plus **Shannon** and **Simpson**
diversity computed on the species-level abundances. Reads whose kingdom is empty or
`"Unclassified"` are excluded from `classifiedReads` and from the abundance denominator, while
`totalReads` counts every input read.

## Core Documentation Reference

- Source: [MetagenomicsAnalyzer.cs#L420](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Metagenomics/MetagenomicsAnalyzer.cs#L420)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `classifications` | TaxonomicClassification[] | Yes | Per-read classifications (e.g. from `classify_reads`). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `kingdomAbundance` / `phylumAbundance` / `genusAbundance` / `speciesAbundance` | array | `{ name, fraction }` per rank. |
| `shannonDiversity` | number | Shannon index at species level. |
| `simpsonDiversity` | number | Simpson index at species level. |
| `totalReads` | integer | All input reads. |
| `classifiedReads` | integer | Reads with a real kingdom. |

## Errors

None. Empty input yields an empty profile.

## Example

**User Prompt:**
> Summarise 4 classified reads (2 E. coli, 2 S. aureus) and 1 unclassified read.

**Response:**
```json
{
  "totalReads": 5,
  "classifiedReads": 4,
  "speciesAbundance": [ { "name": "coli", "fraction": 0.5 }, { "name": "aureus", "fraction": 0.5 } ],
  "shannonDiversity": 0.6931471805599453,
  "simpsonDiversity": 0.5
}
```

Two even species → Shannon = ln 2, Simpson = 0.5; the unclassified read is excluded from the
denominator.

## References

- [MetagenomicsAnalyzer.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Metagenomics/MetagenomicsAnalyzer.cs) — `GenerateTaxonomicProfile`
- Shannon C.E. (1948); Simpson E.H. (1949); Segata N. et al. (2012) Nat. Methods 9:811 (MetaPhlAn).

## See Also

- [classify_reads](classify_reads.md)
- [alpha_diversity](alpha_diversity.md)
