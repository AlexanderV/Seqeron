# predict_tf_binding_change

Predict transcription-factor binding score changes induced by a single-nucleotide variant.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Annotation |
| **Tool Name** | `predict_tf_binding_change` |
| **Method ID** | `VariantAnnotator.PredictTfBindingChange` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Scores a panel of transcription-factor (TF) IUPAC motifs against the reference sequence context and
against the context after applying the SNV's alternate base, then reports the motifs whose best-match
score changes appreciably. The match score is the best fraction of matching motif positions across all
windows that overlap the variant (IUPAC ambiguity codes A/C/G/T/N/R/Y/W/S are honoured). Only motifs
whose absolute score change exceeds `0.1` are returned; non-SNV variants yield no results.

## Core Documentation Reference

- Source: [VariantAnnotator.cs#L1316](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/VariantAnnotator.cs#L1316)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `variant` | object | Yes | Variant to evaluate; `type` must be `SNV` |
| `motifs` | array | Yes | TF motifs (`tfName`, IUPAC `motif`, `threshold`); at least one |
| `referenceContext` | string | Yes | Reference sequence context centred near the variant |
| `contextOffset` | integer | No | 0-based offset of the variant within the context (default 20) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `changes` | array | One entry per affected TF |
| `changes[].tfName` | string | Transcription-factor name |
| `changes[].refScore` | number | Best motif-match fraction on the reference context (0–1) |
| `changes[].altScore` | number | Best motif-match fraction on the mutated context (0–1) |
| `changes[].scoreDifference` | number | `refScore - altScore` (positive ⇒ binding weakened) |

## Errors

| Code | Message |
|------|---------|
| 1001 | variant cannot be null |
| 1001 | Reference context cannot be null or empty |
| 1001 | Motifs cannot be null or empty |
| 1002 | Variant must be of type SNV |

## Examples

### Example 1: SNV disrupts an AP-1 (TGACTCA) motif

**User Prompt:**
> Does the T>A SNV at the start of this AP-1 site weaken TF binding?

**Expected Tool Call:**
```json
{
  "tool": "predict_tf_binding_change",
  "arguments": {
    "variant": { "chromosome": "chr1", "position": 100, "reference": "T", "alternate": "A", "type": "SNV" },
    "motifs": [ { "tfName": "AP-1", "motif": "TGACTCA", "threshold": 0.9 } ],
    "referenceContext": "NNNNNNNNNNNNNNNNNNNNTGACTCANNNNNNNNNNNNNN",
    "contextOffset": 20
  }
}
```

**Response:**
```json
{
  "changes": [
    { "tfName": "AP-1", "refScore": 1.0, "altScore": 0.8571428571428571, "scoreDifference": 0.14285714285714285 }
  ]
}
```

The reference is a perfect 7/7 match (score 1.0); mutating the first motif base breaks one position,
giving 6/7 ≈ 0.857. The difference 1/7 ≈ 0.143 exceeds the 0.1 reporting threshold.

### Example 2: Silent change is filtered

**User Prompt:**
> The alternate base equals the reference — is anything reported?

**Response:**
```json
{ "changes": [] }
```

With `reference == alternate`, `refScore == altScore`, so the change (0) is below the 0.1 threshold and
nothing is returned.

## Performance

- **Time Complexity:** O(m · w · L) for m motifs, window count w, motif length L
- **Space Complexity:** O(context length)

## See Also

- [annotate_regulatory_elements](annotate_regulatory_elements.md) — regulatory regions overlapping a variant
- [predict_variant_effect](predict_variant_effect.md) — coding consequence of a variant
