# predict_variant_effect

Predict the protein-level effect of a variant in a coding sequence.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Annotation |
| **Tool Name** | `predict_variant_effect` |
| **Method ID** | `VariantCaller.PredictEffect` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Predicts a SNP's protein-level consequence within a coding sequence. The reference and mutant codons
at `variantPosition` are translated (standard genetic code): equal amino acids → `Synonymous`; a new
stop → `Nonsense`; loss of a stop → `StopLoss`; any other change → `Missense`. Insertions/deletions
return `Frameshift`, and out-of-range positions or other cases return `Unknown`.

## Core Documentation Reference

- Source: [VariantCaller.cs#L269](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/VariantCaller.cs#L269)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `variant` | object | Yes | Variant to evaluate |
| `codingSequence` | string | Yes | Coding DNA sequence (CDS) |
| `variantPosition` | integer | Yes | 0-based position of the variant within the CDS |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `effect` | string | `Synonymous` / `Missense` / `Nonsense` / `StopLoss` / `Frameshift` / `Unknown` |

## Errors

| Code | Message |
|------|---------|
| 1001 | Variant cannot be null |
| 1001 | Coding sequence cannot be null or empty |

## Examples

### Example 1: Synonymous

`TTA` (Leu) with A>G at position 2 → `TTG` (still Leu).

**Response:**
```json
{ "effect": "Synonymous" }
```

### Example 2: Nonsense

`TAT` (Tyr) with T>A at position 2 → `TAA` (Stop).

**Response:**
```json
{ "effect": "Nonsense" }
```

## Performance

- **Time Complexity:** O(1)
- **Space Complexity:** O(1)

## See Also

- [classify_mutation](classify_mutation.md) - Transition/transversion class
- [annotate_variants](annotate_variants.md) - Call and annotate all variants
