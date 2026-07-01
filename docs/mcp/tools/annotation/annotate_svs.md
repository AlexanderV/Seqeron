# annotate_svs

Annotate structural variants with overlapping genes/exons and a coarse functional impact.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Annotation |
| **Tool Name** | `annotate_svs` |
| **Method ID** | `StructuralVariantAnalyzer.AnnotateSVs` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

For each structural variant, finds the same-chromosome genes it overlaps (`sv.Start <= gene.End &&
sv.End >= gene.Start`) and the exons within them it hits (`GENE:exonN`, 1-based). The functional
impact is then classified:

- **Exon(s) hit** → `HIGH` for Deletion / Inversion / Translocation, `MODERATE` for Duplication and
  everything else.
- **Gene overlap but no exon** → `MODIFIER`.
- **No gene overlap** → `LOW`.

`isPathogenic` is `true` when the impact is `HIGH` or `MODERATE`. `populationFrequency` is always
`0` (no population database attached).

## Core Documentation Reference

- Source: [StructuralVariantAnalyzer.cs#L1135](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/StructuralVariantAnalyzer.cs#L1135)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `variants` | array | Yes | Structural variants `{ id, chromosome, start, end, type, length, quality, supportingReads, insertedSequence? }` |
| `genes` | array | Yes | Gene models `{ geneId, chromosome, start, end, exons[] }` |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `annotations` | array | Per-SV `{ svId, affectedGenes[], affectedExons[], functionalImpact, populationFrequency, isPathogenic }` |

## Errors

| Code | Message |
|------|---------|
| 1001 | variants cannot be null |
| 1001 | genes cannot be null |

## Examples

### Example 1: Deletion hitting an exon (HIGH)

**Response:**
```json
{ "annotations": [ { "svId": "sv1", "affectedGenes": ["GENE1"], "affectedExons": ["GENE1:exon1"], "functionalImpact": "HIGH", "isPathogenic": true } ] }
```

### Example 2: No gene overlap (LOW)

**Response:**
```json
{ "annotations": [ { "svId": "sv1", "affectedGenes": [], "affectedExons": [], "functionalImpact": "LOW", "isPathogenic": false } ] }
```

## Performance

- **Time Complexity:** O(v · g) with early chromosome filtering
- **Space Complexity:** O(g)

## See Also

- [genotype_sv](genotype_sv.md) - Genotype an SV from read counts
- [filter_svs](filter_svs.md) - Filter SVs by quality/size
- [merge_overlapping_svs](merge_overlapping_svs.md) - Merge overlapping SV calls
