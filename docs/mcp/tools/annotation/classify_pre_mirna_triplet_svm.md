# classify_pre_mirna_triplet_svm

Published Triplet-SVM pre-miRNA classifier (Xue et al. 2005) with the authors' trained libsvm model.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Annotation |
| **Tool Name** | `classify_pre_mirna_triplet_svm` |
| **Method ID** | `MiRnaAnalyzer.ClassifyPreMiRnaTripletSvm` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Triplet-SVM (Xue, Li, He, Liu, Li & Zhang 2005, BMC Bioinformatics 6:310) with the authors' own trained libsvm RBF model. The query is folded with the canonical MFE engine (or `dotBracket` is used; the package used RNAfold), stem-loops with at least `minBasePairs` paired bases per arm are retrieved (`2_get_stemloop.pl`), each single stem-loop is coded as 32 triplet-element frequencies (`3_step_triplet_coding_for_queries.pl`) and scored; `decisionValue > 0` ⇒ real pre-miRNA (`isPreMiRna`). Coordinates are 1-based. `found = false` (empty `candidates`) when no stem-loop is retrieved — no stem-loop with ≥ `minBasePairs` pairs per arm (default 22, the package example), a non-ACGU query (dropped like the package), or a multi-branch structure. `minBasePairs` must be ≥ 1 and `minLoopSize` ≥ 0 (values 0–2 are raised to the NNDB minimum 3 by the folder).

## Core Documentation Reference

- Source: [MiRnaAnalyzer.TripletSvm.cs#L59](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/MiRnaAnalyzer.TripletSvm.cs#L59)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | Query nucleotide sequence (RNA or DNA; T read as U) |
| `dotBracket` | string | No | Optional structure in '.()' notation (same length); folded when omitted |
| `minBasePairs` | integer | No | Minimum paired bases per arm of a retrieved stem-loop (≥ 1) (default 22) |
| `minLoopSize` | integer | No | Minimum hairpin loop size of the fold (≥ 0; only when folding) (default 3) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `found` | boolean | True iff at least one stem-loop was scored |
| `candidates[].start` | integer | 1-based start in the query |
| `candidates[].end` | integer | 1-based end in the query |
| `candidates[].sequence` | string | Stem-loop sequence |
| `candidates[].dotBracket` | string | Stem-loop structure |
| `candidates[].features` | array | The 32 triplet-element frequencies |
| `candidates[].decisionValue` | number | libsvm decision value |
| `candidates[].isPreMiRna` | boolean | decisionValue > 0 |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1001 | minBasePairs must be >= 1 |
| 1001 | minLoopSize must be >= 0 |
| 1001 | dotBracket must have the same length as the sequence. |

## Examples

### Example 1: Supplied structure, 4 bp arms

**Request:**
```json
{"sequence": "GGGGAAAACCCC", "dotBracket": "((((....))))", "minBasePairs": 3}
```

**Response (abridged):**
```json
{"found": true, "candidates": [{"start": 1, "end": 12, "dotBracket": "((((....))))"}]}
```

## Performance

- **Time Complexity:** One MFE fold (when folding), O(n³)

## See Also

- [classify_pre_mirna](classify_pre_mirna.md)
- [calculate_randfold_pvalue](calculate_randfold_pvalue.md)
- [find_pre_mirna_hairpins](find_pre_mirna_hairpins.md)
