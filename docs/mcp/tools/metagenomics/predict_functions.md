# predict_functions

Predict functional annotations for proteins.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Metagenomics |
| **Tool Name** | `predict_functions` |
| **Method ID** | `MetagenomicsAnalyzer.PredictFunctions` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

For each protein, scans the function database for signatures (motifs) contained in the protein
sequence. Each match is scored with an ungapped BLOSUM62 self-alignment, converted to a BLAST
bit score `(λ·S − ln K)/ln 2` and E-value `K·m·n·e^(−λ·S)` (λ = 0.3176, K = 0.134). The best
hit — the lowest E-value — is transferred, carrying the entry's function, pathway, and KO plus
an inferred COG category.

## Core Documentation Reference

- Source: [MetagenomicsAnalyzer.cs#L1114](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Metagenomics/MetagenomicsAnalyzer.cs#L1114)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `proteins` | ProteinInput[] | Yes | `{ geneId, proteinSequence }` entries. |
| `functionDatabase` | FunctionDatabaseEntry[] | Yes | `{ motif, function, pathway, ko }` entries. |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `items[].geneId` | string | Gene id. |
| `items[].function` / `.pathway` / `.koNumber` | string | Transferred annotation. |
| `items[].cogCategory` | string | Inferred COG category. |
| `items[].eValue` | number | BLAST E-value of the best hit. |
| `items[].bitScore` | number | BLAST bit score of the best hit. |

## Errors

| Code | Message |
|------|---------|
| 1003 | A required argument is null (`ArgumentNullException`). |

## Example

**User Prompt:**
> Annotate protein `WWW` against a database containing motif `WWW` = tryptophanase.

**Response:**
```json
{ "items": [ { "geneId": "g1", "function": "tryptophanase", "pathway": "Amino acid metabolism", "koNumber": "K01667", "bitScore": 18.0202932787533, "eValue": 3.3852730346546e-05 } ] }
```

## References

- [MetagenomicsAnalyzer.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Metagenomics/MetagenomicsAnalyzer.cs) — `PredictFunctions`
- Altschul S.F. et al. (1990) J. Mol. Biol. 215:403 (BLAST); NCBI BLOSUM62.

## See Also

- [functional_diversity](functional_diversity.md)
- [find_resistance_genes](find_resistance_genes.md)
