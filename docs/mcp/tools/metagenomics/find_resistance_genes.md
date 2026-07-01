# find_resistance_genes

Search genes for antibiotic-resistance markers.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Metagenomics |
| **Tool Name** | `find_resistance_genes` |
| **Method ID** | `MetagenomicsAnalyzer.FindResistanceGenes` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Scans each gene for antibiotic-resistance markers by simple substring (motif) containment
against a resistance database. When a gene's sequence contains a database motif, a hit is
reported with the entry's name, antibiotic class, and an identity of
`motifLength / geneLength`. A gene may produce multiple hits (one per matching motif).

## Core Documentation Reference

- Source: [MetagenomicsAnalyzer.cs#L1562](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Metagenomics/MetagenomicsAnalyzer.cs#L1562)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `genes` | GeneInput[] | Yes | Genes to scan. |
| `resistanceDatabase` | ResistanceDatabaseEntry[] | Yes | `{ motif, name, antibioticClass }` entries. |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `items[].geneId` | string | Gene id of the hit. |
| `items[].resistanceGene` | string | Matched database entry name. |
| `items[].antibioticClass` | string | Antibiotic class. |
| `items[].identity` | number | `motifLength / geneLength`. |

## Errors

None. Genes with no motif match are omitted.

## Example

**User Prompt:**
> Does gene `AAACGTACGT` carry a beta-lactam resistance motif (`CGTACGT`)?

**Response:**
```json
{ "items": [ { "geneId": "g1", "resistanceGene": "blaX-like", "antibioticClass": "beta-lactam", "identity": 0.7 } ] }
```

Identity = 7 / 10 = 0.7.

## References

- [MetagenomicsAnalyzer.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Metagenomics/MetagenomicsAnalyzer.cs) — `FindResistanceGenes`
- Zankari E. et al. (2012) J. Antimicrob. Chemother. 67:2640 (ResFinder).

## See Also

- [predict_functions](predict_functions.md)
