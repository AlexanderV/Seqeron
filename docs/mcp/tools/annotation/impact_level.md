# impact_level

Map a ConsequenceType to an ImpactLevel (High/Moderate/Low/Modifier).

## Overview

| Property | Value |
|----------|-------|
| **Server** | Annotation |
| **Tool Name** | `impact_level` |
| **Method ID** | `VariantAnnotator.GetImpactLevel` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Maps a VEP-style variant consequence to its Ensembl impact category. **High**: transcript ablation,
splice acceptor/donor, stop gained, frameshift, stop/start lost, transcript amplification.
**Moderate**: inframe insertion/deletion, missense, protein-altering. **Low**: splice region,
incomplete terminal codon, start/stop retained, synonymous. Everything else (e.g. `IntronVariant`)
is a **Modifier**. The consequence name is matched case-insensitively.

## Core Documentation Reference

- Source: [VariantAnnotator.cs#L783](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/VariantAnnotator.cs#L783)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `consequence` | string | Yes | Consequence type name (e.g. MissenseVariant, StopGained, IntronVariant) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `impact` | string | `High` / `Moderate` / `Low` / `Modifier` |

## Errors

| Code | Message |
|------|---------|
| 1001 | Consequence cannot be null or empty |
| 1001 | Unknown consequence type |

## Examples

### Example 1: StopGained is High impact

**Response:**
```json
{ "impact": "High" }
```

### Example 2: MissenseVariant is Moderate impact

**Response:**
```json
{ "impact": "Moderate" }
```

## Performance

- **Time Complexity:** O(1)
- **Space Complexity:** O(1)

## See Also

- [predict_pathogenicity](predict_pathogenicity.md) - ACMG-like pathogenicity prediction
- [annotate_variant_on_transcripts](annotate_variant_on_transcripts.md) - Full consequence annotation
