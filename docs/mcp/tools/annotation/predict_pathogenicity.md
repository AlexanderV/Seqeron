# predict_pathogenicity

ACMG-like pathogenicity prediction combining annotation, frequency, conservation, ClinVar, and
functional evidence.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Annotation |
| **Tool Name** | `predict_pathogenicity` |
| **Method ID** | `VariantAnnotator.PredictPathogenicity` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Applies a simplified ACMG-style scoring scheme. It accumulates pathogenic and benign points from the
annotation's impact (PVS1 for high impact), population frequency (PM2 for rare, BA1/BS1 for common),
computational scores (PP3/BP4 from SIFT/PolyPhen), conservation, ClinVar significance (PP5/BP6), and
functional evidence (PS3 for LOF). The net score maps to a `classification` (Pathogenic /
LikelyPathogenic / UncertainSignificance / LikelyBenign / Benign), with a `confidenceScore`, the
list of `evidenceCriteria` applied, and an `isActionable` flag (true for pathogenic / likely
pathogenic).

## Core Documentation Reference

- Source: [VariantAnnotator.cs#L1004](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/VariantAnnotator.cs#L1004)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `annotation` | object | Yes | Variant annotation (one transcript) |
| `populationFrequency` | number | No | Population allele frequency (0..1) |
| `conservationScore` | number | No | Conservation score (e.g. PhyloP) |
| `inClinvar` | boolean | No | Whether variant is present in ClinVar |
| `clinvarSignificance` | string | No | ClinVar clinical significance string |
| `functionalEvidence` | array | No | Functional evidence labels (e.g. 'LOF confirmed') |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `prediction` | object | `{ classification, confidenceScore, evidenceCriteria, clinicalSignificance, isActionable }` |

## Errors

| Code | Message |
|------|---------|
| 1001 | Annotation cannot be null |

## Examples

### Example 1: High-impact rare variant

A `StopGained` (High impact) variant at AF 0.00001 → likely pathogenic and actionable (PVS1 applied).

**Response:**
```json
{ "prediction": { "classification": "LikelyPathogenic", "isActionable": true } }
```

### Example 2: Common synonymous variant

A `SynonymousVariant` (Low impact) at 10% AF → likely benign (BA1 + BP7), not actionable.

**Response:**
```json
{ "prediction": { "classification": "LikelyBenign", "isActionable": false } }
```

## Performance

- **Time Complexity:** O(functional evidence count)
- **Space Complexity:** O(evidence)

## See Also

- [impact_level](impact_level.md) - Consequence → impact mapping
- [annotate_variant_on_transcripts](annotate_variant_on_transcripts.md) - Produce the annotation
