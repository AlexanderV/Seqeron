# classify_reads

Classify metagenomic reads with the Kraken k-mer / LCA algorithm.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Metagenomics |
| **Tool Name** | `classify_reads` |
| **Method ID** | `MetagenomicsAnalyzer.ClassifyReads` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

For each read the tool collects its canonical k-mer hits against the supplied database, builds
the Kraken classification tree over the hit taxa (weighted by k-mer count), scores every
root-to-leaf (RTL) path, and assigns the leaf of the maximum-scoring path. Ties are broken by
the LCA of the tied leaves. Reads with no hits are assigned the root (unclassified).
Confidence reuses Kraken 2's C/Q score = clade k-mers / non-ambiguous k-mers queried
(Wood & Salzberg 2014).

## Core Documentation Reference

- Source: [MetagenomicsAnalyzer.cs#L191](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Metagenomics/MetagenomicsAnalyzer.cs#L191)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `reads` | ReadInput[] | Yes | `{ id, sequence }` reads. |
| `kmerDatabase` | KmerDatabaseEntry[] | Yes | Canonical `{ kmer, taxonId }` entries (e.g. from `build_kmer_database`). |
| `taxonomy` | TaxonNodeInput[] | Yes | `{ id, name, rank, parentId }` nodes. |
| `k` | integer | No | k-mer length (default 31). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `items[].readId` | string | Echoed read id. |
| `items[].taxonId` | integer | Assigned taxon id. |
| `items[].taxonName` / `.rank` | string | Assigned taxon name and rank. |
| `items[].rtlScore` | integer | Maximum RTL path score. |
| `items[].confidence` | number | Kraken2 C/Q confidence. |
| `items[].matchedKmers` / `.totalKmers` | integer | C and Q. |
| `items[].kingdom … species` | string | Rank lineage of the assigned taxon. |

## Errors

| Code | Message |
|------|---------|
| 1002 | `k` must be positive (`ArgumentOutOfRangeException`). |
| 1003 | A required argument is null (`ArgumentNullException`). |

## Example

**User Prompt:**
> Classify a read whose 4 k-mers all hit E. coli.

Database `{AAAA,AAAC,AACA,ACAA} → 100`, read `AAAACAA`, k = 4:

**Response:**
```json
{
  "items": [
    { "readId": "r", "taxonId": 100, "taxonName": "Escherichia coli", "rank": "species", "rtlScore": 4, "matchedKmers": 4, "totalKmers": 4, "confidence": 1.0 }
  ]
}
```

If the k-mers split 2/2 between two sibling species, the read is assigned their genus LCA
(taxon 20, rank genus, RTL score 2).

## References

- [MetagenomicsAnalyzer.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Metagenomics/MetagenomicsAnalyzer.cs) — `ClassifyReads`
- Wood D.E. & Salzberg S.L. (2014) Genome Biology 15:R46 (Kraken).

## See Also

- [build_kmer_database](build_kmer_database.md) — build the database
- [taxonomic_profile](taxonomic_profile.md) — aggregate classifications into a profile
