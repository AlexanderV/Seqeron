# annotate_variant_on_transcripts

VEP-like annotation of a single variant against one or more transcript models.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Annotation |
| **Tool Name** | `annotate_variant_on_transcripts` |
| **Method ID** | `VariantAnnotator.AnnotateVariant` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

For a single variant, filters the supplied transcripts to those on the same chromosome that
overlap or lie near the variant (5000 bp upstream / 500 bp downstream), then for each relevant
transcript predicts:

- the most severe **consequence** (e.g. `MissenseVariant`, `IntronVariant`, `SpliceAcceptorVariant`,
  `FrameshiftVariant`, `UpstreamGeneVariant`, `FivePrimeUtrVariant`),
- the corresponding **impact** level (`High` / `Moderate` / `Low` / `Modifier`),
- for coding consequences with a CDS: the **codon change**, **amino-acid change**, protein position
  and CDS position, and
- for **missense** variants: heuristic **SIFT** and **PolyPhen** scores.

If no transcript is relevant, a single `IntergenicVariant` / `Modifier` annotation is returned.

## Core Documentation Reference

- Source: [VariantAnnotator.cs#L266](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/VariantAnnotator.cs#L266)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `variant` | object | Yes | `{ chromosome, position, reference, alternate, type, quality?, id? }` (1-based position) |
| `transcripts` | array | Yes | Transcript models `{ transcriptId, geneId, geneName, chromosome, start, end, strand, exons[], codingExons[], cdsStart?, cdsEnd? }` (min 1) |
| `referenceSequence` | string | No | Reference sequence used for codon/AA change prediction |
| `populationFrequencies` | object | No | Map of population frequency labels to values |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `annotations` | array | Per-transcript `{ variant, transcriptId, geneId, geneName, consequence, impact, codonChange?, aminoAcidChange?, proteinPosition?, cdsPosition?, siftScore?, polyphenScore?, ... }` |

## Errors

| Code | Message |
|------|---------|
| 1001 | Variant cannot be null |
| 1001 | Variant chromosome / reference / alternate cannot be empty |
| 1001 | Transcripts cannot be null or empty |
| 1002 | Unknown variant type |

## Examples

### Example 1: Missense variant in a coding exon (Moderate impact)

Variant `chr1:1450 A>G` (SNV) inside the coding exon `[1400,1500]` of a `+`-strand transcript
`[1000,2000]` (CDS `1050..1900`).

**Response:**
```json
{ "annotations": [ { "transcriptId": "T1", "geneId": "G1", "geneName": "Gene1", "consequence": "MissenseVariant", "impact": "Moderate" } ] }
```

### Example 2: Intronic variant (Modifier impact)

Variant `chr1:1250 A>G` (SNV) between exons `[1000,1100]` and `[1400,1500]` of the same transcript.

**Response:**
```json
{ "annotations": [ { "transcriptId": "T1", "geneId": "G1", "geneName": "Gene1", "consequence": "IntronVariant", "impact": "Modifier" } ] }
```

## Performance

- **Time Complexity:** O(t · e) — transcripts times exons, with early chromosome/window filtering
- **Space Complexity:** O(t)

## See Also

- [classify_variant](classify_variant.md) - Classify a (ref, alt) allele pair
- [impact_level](impact_level.md) - Map a consequence to its impact level
- [predict_pathogenicity](predict_pathogenicity.md) - ACMG-like pathogenicity prediction
