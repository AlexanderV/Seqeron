# classify_mutation

Classify a SNP as transition, transversion, or other.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Annotation |
| **Tool Name** | `classify_mutation` |
| **Method ID** | `VariantCaller.ClassifyMutation` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Classifies a single-nucleotide variant by base-change class. A substitution between two purines
(A↔G) or two pyrimidines (C↔T) is a **Transition**; a purine↔pyrimidine substitution is a
**Transversion**. Any variant whose `type` is not `SNP` classifies as **Other**.

## Core Documentation Reference

- Source: [VariantCaller.cs#L184](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/VariantCaller.cs#L184)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `variant` | object | Yes | Variant to classify (`referenceAllele`, `alternateAllele`, `type` required) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `mutationType` | string | `Transition`, `Transversion`, or `Other` |

## Errors

| Code | Message |
|------|---------|
| 1001 | Variant cannot be null |
| 1001 | Variant alleles cannot be empty |

## Examples

### Example 1: A>G transition

A purine→purine SNP.

**Response:**
```json
{ "mutationType": "Transition" }
```

### Example 2: A>C transversion

A purine→pyrimidine SNP.

**Response:**
```json
{ "mutationType": "Transversion" }
```

## Performance

- **Time Complexity:** O(1)
- **Space Complexity:** O(1)

## See Also

- [titv_ratio](titv_ratio.md) - Transition/transversion ratio across many variants
- [find_snps](find_snps.md) - SNP detection
