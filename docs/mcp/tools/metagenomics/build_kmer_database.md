# build_kmer_database

Build a Kraken canonical-k-mer → taxon-id database.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Metagenomics |
| **Tool Name** | `build_kmer_database` |
| **Method ID** | `MetagenomicsAnalyzer.BuildKmerDatabase` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Builds the k-mer → taxon database consumed by [classify_reads](classify_reads.md)
(Wood & Salzberg 2014, Kraken). For every reference sequence, each length-`k` window is
canonicalised (the lexicographic minimum of the window and its reverse complement) and mapped
to the reference's taxon. When a canonical k-mer is contributed by several taxa, its stored
value is collapsed to the **lowest common ancestor** of those taxa in the taxonomy tree.

Ambiguous (non-ACGT) windows are skipped and input is upper-cased. A reference whose taxon id
is not present in the taxonomy tree raises an error.

## Core Documentation Reference

- Source: [MetagenomicsAnalyzer.cs#L363](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Metagenomics/MetagenomicsAnalyzer.cs#L363)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `referenceGenomes` | ReferenceGenomeInput[] | Yes | `{ taxonId, sequence }` references. |
| `taxonomy` | TaxonNodeInput[] | Yes | `{ id, name, rank, parentId }` nodes; root is self-parented. |
| `k` | integer | No | k-mer length (default 31). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `entries[].kmer` | string | Canonical k-mer. |
| `entries[].taxonId` | integer | LCA taxon id of the references sharing that k-mer. |
| `count` | integer | Number of distinct canonical k-mers. |

## Errors

| Code | Message |
|------|---------|
| 4001 | Reference taxon id is not in the taxonomy tree (`KeyNotFoundException`). |

## Example

**User Prompt:**
> Build a k=4 database for two Escherichia species that share the palindrome AGCT.

References `(100, "AGCTAAAA")` and `(101, "AGCTCCCC")` under genus `Escherichia(20)`:

**Response (selected entries):**
```json
{
  "entries": [
    { "kmer": "AGCT", "taxonId": 20 },
    { "kmer": "GCTA", "taxonId": 100 },
    { "kmer": "GAGC", "taxonId": 101 }
  ]
}
```

The shared `AGCT` collapses to `LCA(100, 101) = 20`; species-specific k-mers stay at the species.

## References

- [MetagenomicsAnalyzer.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Metagenomics/MetagenomicsAnalyzer.cs) — `BuildKmerDatabase`
- Wood D.E. & Salzberg S.L. (2014) Genome Biology 15:R46 (Kraken).

## See Also

- [classify_reads](classify_reads.md) — classify reads against this database
