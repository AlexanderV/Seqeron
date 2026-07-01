# predict_genes

Predict gene annotations from a DNA sequence using ORF-based heuristics.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Annotation |
| **Tool Name** | `predict_genes` |
| **Method ID** | `GenomeAnnotator.PredictGenes` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Finds ORFs of at least `minOrfLength` amino acids on both strands (start codon required), orders them
by start position, and emits one gene annotation per ORF. Every predicted gene is typed `CDS` with
product `hypothetical protein`, a strand of `+` (forward) or `-` (reverse), and a sequential id
`{prefix}_{n:D4}` (e.g. `gene_0001`). The `attributes` map carries `frame`, `protein_length` (the
translated length excluding the stop), and `translation` (the protein, ending with `*`).

## Core Documentation Reference

- Source: [GenomeAnnotator.cs#L402](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/GenomeAnnotator.cs#L402)

## Input Schema

| Parameter | Type | Required | Default | Description |
|-----------|------|----------|---------|-------------|
| `dnaSequence` | string | Yes | — | DNA sequence to annotate |
| `minOrfLength` | integer | No | 100 | Minimum ORF length in amino acids |
| `prefix` | string | No | gene | Gene identifier prefix |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `genes` | array | `{ geneId, start, end, strand, type, product, attributes }` per predicted gene |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |

## Examples

### Example 1: Single forward gene

`ATG` + 297 `A` + `TAA` (a 100-aa ORF), `minOrfLength=50` → one `CDS` gene `gene_0001` on the `+`
strand spanning 0..303.

**Response:**
```json
{ "genes": [ { "geneId": "gene_0001", "start": 0, "end": 303, "strand": "+", "type": "CDS", "product": "hypothetical protein" } ] }
```

### Example 2: Two genes get sequential ids

Two ORFs separated by a spacer, `prefix="test"` → `test_0001`, `test_0002`.

## Performance

- **Time Complexity:** O(n) ORF scan + O(k log k) ordering
- **Space Complexity:** O(k)

## See Also

- [find_orfs](find_orfs.md) - Underlying ORF detection
- [to_gff3](to_gff3.md) - Serialize predicted genes to GFF3
- [predict_gene_structure](predict_gene_structure.md) - Splice-aware gene structure
