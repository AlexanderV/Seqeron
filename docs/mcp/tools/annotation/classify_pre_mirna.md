# classify_pre_mirna

In-repo ridge-logistic classifier of natural pre-miRNA vs dinucleotide-shuffled background (not a published classifier).

## Overview

| Property | Value |
|----------|-------|
| **Server** | Annotation |
| **Tool Name** | `classify_pre_mirna` |
| **Method ID** | `MiRnaAnalyzer.ClassifyPreMiRna` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Scores a candidate hairpin with the bundled **in-repo** ridge-logistic model over five features of the canonical MFE fold (MFE, AMFE, MFEI — Zhang 2006, GC content, paired fraction). The model was trained on 13 public-domain miRBase precursors against 52 dinucleotide-shuffled negatives (Bonnet 2004 convention): it is an in-repo model, **not a published classifier**. For a published classifier with the authors' own trained model use [classify_pre_mirna_triplet_svm](classify_pre_mirna_triplet_svm.md) (Triplet-SVM, Xue et al. 2005). Returns the extracted features, `naturalProbability` P(natural) and `isNatural` (P ≥ `threshold`). `found = false` for an empty sequence.

## Core Documentation Reference

- Source: [MiRnaAnalyzer.cs#L3650](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/MiRnaAnalyzer.cs#L3650)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | Candidate nucleotide sequence (RNA or DNA; T read as U) |
| `threshold` | number | No | Decision threshold on P(natural) in [0, 1] (default 0.5) |
| `minLoopSize` | integer | No | Minimum hairpin loop size of the folder (NNDB minimum 3) (default 3) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `found` | boolean | False for an empty sequence |
| `features` | object | freeEnergy, amfe, mfei, gcContent, pairedFraction, stemBasePairs, loopSize, length, basePairingPropensity |
| `naturalProbability` | number | P(natural) in [0, 1] |
| `isNatural` | boolean | naturalProbability ≥ threshold |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1001 | threshold must be in [0, 1] |
| 1001 | minLoopSize must be >= 0 |

## Examples

### Example 1: hsa-mir-21 is natural

**Request:**
```json
{"sequence": "UGUCGGGUAGCUUAUCAGACUGAUGUUGACUGUUGAAUCUCAUGGCAACACCAGUCGAUGGGCUGUCUGACA"}
```

**Response (abridged):**
```json
{"found": true, "naturalProbability": ">0.95", "isNatural": true}
```

## Performance

- **Time Complexity:** One MFE fold, O(n³)

## See Also

- [classify_pre_mirna_triplet_svm](classify_pre_mirna_triplet_svm.md)
- [calculate_randfold_pvalue](calculate_randfold_pvalue.md)
- [find_pre_mirna_hairpins](find_pre_mirna_hairpins.md)
