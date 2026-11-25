# compute_mipred_features

The 34 MiPred features V1..V34 from a sequence, its dot-bracket structure, the MFE and a p-value.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Annotation |
| **Tool Name** | `compute_mipred_features` |
| **Method ID** | `MiRnaAnalyzer.ComputeMiPredFeatures` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Computes the MiPred feature vector (Jiang et al. 2007; `microRNAcheck_parallel.pl`): V1..V32 are the frequencies of the 32 triplet elements (centre nucleotide A, U, G, C x local structure `(((`, `((.`, `(.(`, `(..`, `.((`, `.(.`, `..(`, `...`; closing brackets are not counted; positions 2..n-1) divided by the sequence **length**; V33 is the MFE and V34 the p-value, both supplied by the caller. Pure counting: no folding, no shuffling.

## Core Documentation Reference

- Source: [MiRnaAnalyzer.MiPred.cs#L146](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/MiRnaAnalyzer.MiPred.cs#L146)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | Nucleotide sequence (A/C/G/T/U only) |
| `dotBracket` | string | Yes | Structure in ".()" notation, same length as the sequence |
| `freeEnergy` | number | Yes | MFE (kcal/mol) -> V33 |
| `pValue` | number | Yes | Shuffle p-value in [0, 1] -> V34 |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `features` | array | V1..V32 triplet-element frequencies (count / sequence length), V33 = freeEnergy, V34 = pValue |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1001 | dotBracket cannot be null or empty |
| 1001 | pValue must be in [0, 1] |

## Examples

### Example 1: Hand-computed triplets

**Request:**
```json
{"sequence": "GAUCC", "dotBracket": "((.((", "freeEnergy": -12.5, "pValue": 0.25}
```

**Response (abridged):**
```json
{"features": "V2 = V11 = V29 = 0.2, V33 = -12.5, V34 = 0.25, others 0"}
```

## Performance

- **Time Complexity:** O(n), instant

## See Also

- [classify_pre_mirna](classify_pre_mirna.md)
- [classify_pre_mirna_triplet_svm](classify_pre_mirna_triplet_svm.md)
- [calculate_randfold_pvalue](calculate_randfold_pvalue.md)
