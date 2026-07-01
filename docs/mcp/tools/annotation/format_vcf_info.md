# format_vcf_info

Format a variant annotation as a VCF INFO field string.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Annotation |
| **Tool Name** | `format_vcf_info` |
| **Method ID** | `VariantAnnotator.FormatAsVcfInfo` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Serialises a `VariantAnnotation` into a semicolon-delimited VCF INFO string. The four mandatory keys
`GENE`, `TRANSCRIPT`, `CONSEQUENCE`, `IMPACT` are always emitted, in that order. `HGVSP`
(from `aminoAcidChange`), `HGVSC` (from `codonChange`), `SIFT` and `POLYPHEN` are appended only when
their source field is non-null. Numeric SIFT/PolyPhen scores are rendered with three decimals using the
invariant culture (`.` decimal separator, per the VCF v4.x spec).

## Core Documentation Reference

- Source: [VariantAnnotator.cs#L1484](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/VariantAnnotator.cs#L1484)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `annotation` | object | Yes | Variant annotation to serialise |
| `annotation.geneName` | string | Yes | Emitted as `GENE=` |
| `annotation.transcriptId` | string | Yes | Emitted as `TRANSCRIPT=` |
| `annotation.consequence` | string | Yes | `ConsequenceType` enum name |
| `annotation.impact` | string | Yes | `ImpactLevel`: High/Moderate/Low/Modifier |
| `annotation.aminoAcidChange` | string \| null | No | Emitted as `HGVSP=` when present |
| `annotation.codonChange` | string \| null | No | Emitted as `HGVSC=` when present |
| `annotation.siftScore` | number \| null | No | Emitted as `SIFT=` (F3) when present |
| `annotation.polyphenScore` | number \| null | No | Emitted as `POLYPHEN=` (F3) when present |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `info` | string | Semicolon-delimited VCF INFO string |

## Errors

| Code | Message |
|------|---------|
| 1001 | annotation cannot be null |
| 1002 | Unknown consequence type '&lt;value&gt;'. |
| 1002 | Unknown impact level '&lt;value&gt;'. |

## Examples

### Example 1: Full missense annotation

**Response:**
```json
{ "info": "GENE=BRCA1;TRANSCRIPT=ENST1;CONSEQUENCE=MissenseVariant;IMPACT=Moderate;HGVSP=p.Arg170His;HGVSC=c.509G>A;SIFT=0.020;POLYPHEN=0.950" }
```

### Example 2: Only mandatory fields

With no HGVS/SIFT/PolyPhen values supplied:

**Response:**
```json
{ "info": "GENE=BRCA1;TRANSCRIPT=ENST1;CONSEQUENCE=MissenseVariant;IMPACT=Moderate" }
```

## Performance

- **Time Complexity:** O(1) — fixed number of fields
- **Space Complexity:** O(1)

## See Also

- [annotate_variant_on_transcripts](annotate_variant_on_transcripts.md) — produce the annotations to format
- [impact_level](impact_level.md) — impact severity for a consequence
