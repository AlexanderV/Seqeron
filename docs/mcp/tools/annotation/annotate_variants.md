# annotate_variants

Call and annotate all variants between a reference and a query DNA sequence.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Annotation |
| **Tool Name** | `annotate_variants` |
| **Method ID** | `VariantCaller.AnnotateVariants` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Globally aligns the query against the reference, calls every variant (SNP / insertion / deletion),
and annotates each one with:

- **effect** — the protein-level consequence. When `isCodingSequence` is `true` the reference is
  treated as a coding sequence and each SNP's codon is translated to classify the variant as
  `Synonymous`, `Missense`, `Nonsense` or `StopLoss`; indels are reported as `Frameshift`. When
  `isCodingSequence` is `false`, the effect is `Unknown`.
- **mutationType** — for SNPs, `Transition` (purine↔purine or pyrimidine↔pyrimidine) or
  `Transversion` (purine↔pyrimidine); for non-SNPs, `Other`.

Variant `position` is 0-based in the reference.

## Core Documentation Reference

- Source: [VariantCaller.cs#L318](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/VariantCaller.cs#L318)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `reference` | string | Yes | Reference DNA sequence |
| `query` | string | Yes | Query DNA sequence |
| `isCodingSequence` | boolean | No | Interpret reference as a coding sequence and predict effects (default `false`) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `annotated` | array | Per-variant `{ variant { position, referenceAllele, alternateAllele, type, queryPosition }, effect, mutationType }` |

## Errors

| Code | Message |
|------|---------|
| 1001 | Reference cannot be null or empty |
| 1001 | Query cannot be null or empty |

## Examples

### Example 1: Coding missense SNP

Reference `ATGCAT` (codon `ATG` = Met) vs query `ATTCAT` (codon `ATT` = Ile). The single `G>T`
SNP at position 2 is a transversion and a missense change.

**Response:**
```json
{ "annotated": [ { "variant": { "position": 2, "referenceAllele": "G", "alternateAllele": "T", "type": "SNP", "queryPosition": 2 }, "effect": "Missense", "mutationType": "Transversion" } ] }
```

### Example 2: Same SNP, non-coding

**Response:**
```json
{ "annotated": [ { "variant": { "position": 2, "referenceAllele": "G", "alternateAllele": "T", "type": "SNP", "queryPosition": 2 }, "effect": "Unknown", "mutationType": "Transversion" } ] }
```

## Performance

- **Time Complexity:** O(n·m) global alignment plus O(v) annotation
- **Space Complexity:** O(n·m)

## See Also

- [call_variants](call_variants.md) - Call variants without annotation
- [predict_variant_effect](predict_variant_effect.md) - Effect of a single variant in a CDS
- [titv_ratio](titv_ratio.md) - Transition/transversion ratio over a variant set
