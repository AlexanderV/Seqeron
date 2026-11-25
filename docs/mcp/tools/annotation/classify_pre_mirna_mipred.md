# classify_pre_mirna_mipred

Published MiPred random-forest pre-miRNA classifier (Jiang et al. 2007): 34 features, authors' 1000-tree forest.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Annotation |
| **Tool Name** | `classify_pre_mirna_mipred` |
| **Method ID** | `MiRnaAnalyzer.ClassifyPreMiRnaMiPred` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Folds the candidate with the canonical Zuker-Stiegler MFE engine (Turner 2004; the MiPred package used RNAfold 1.x), applies the package pre-judgement (MFE <= -15 kcal/mol, >= 18 base pairs, no multiloop, 51-137 nt), computes the 34 features of `microRNAcheck_parallel.pl` (32 triplet-element frequencies, MFE, shuffle p-value) and lets the authors' trained **1000-tree random forest** vote (real pre-miRNA when more than 50% of the trees vote P). Source: Jiang P, Wu H, Wang W, Ma W, Sun X, Lu Z (2007) *Nucleic Acids Res* 35:W339-W344. The p-value feature is statistical: dinucleotide shuffles are drawn from `System.Random(seed)`, so a fixed seed is reproducible and other seeds vary within sampling error (about +-0.01 at p = 0.1 with 999 shuffles). Differences from the package (documented in the library): Turner 2004 fold instead of RNAfold 1.x, uShuffle-exact shuffle instead of randfold's triplet swap. Does not replace [classify_pre_mirna](classify_pre_mirna.md) (in-repo model) or [classify_pre_mirna_triplet_svm](classify_pre_mirna_triplet_svm.md).

## Core Documentation Reference

- Source: [MiRnaAnalyzer.MiPred.cs#L92](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/MiRnaAnalyzer.MiPred.cs#L92)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | Candidate nucleotide sequence (RNA or DNA; characters other than A/C/G/T/U are dropped; T read as U) |
| `seed` | integer | No | Seed of the System.Random driving the p-value shuffles (default 1); fixed seed = reproducible (default 1) |
| `shuffleCount` | integer | No | Dinucleotide shuffles of the p-value feature, 1..10000 (default 999 = package/training value) (default 999) |
| `minLoopSize` | integer | No | Minimum hairpin loop size of the folder (NNDB minimum 3) (default 3) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `found` | boolean | False when nothing remains of the sequence after the clean-up |
| `sequence` | string | Cleaned sequence (upper case RNA) |
| `dotBracket` | string | Canonical MFE structure |
| `freeEnergy` | number | MFE (kcal/mol) |
| `isHairpinLike` | boolean | Package pre-judgement: MFE <= -15, >= 18 pairs, no multiloop, 51-137 nt |
| `basePairs` | integer | Pairs of the structure (0 when it contains a second hairpin) |
| `hasMultiloop` | boolean | A "(" follows the first ")" |
| `length` | integer | Cleaned length |
| `pValue` | number | (1 + #{shuffle dG <= dG})/(shuffleCount + 1); null when not hairpin-like |
| `features` | array | The 34 features V1..V34; null when not hairpin-like |
| `realVotes` | integer | Trees voting real pre-miRNA |
| `treeCount` | integer | 1000 |
| `realProbability` | number | realVotes/treeCount; null when not hairpin-like |
| `isPreMiRna` | boolean | Hairpin-like and realProbability > 0.5 |
| `seed` | integer | Seed used |
| `shuffleCount` | integer | Shuffles used |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1001 | shuffleCount must be in [1, 10000] |
| 1001 | minLoopSize must be >= 0 |

## Examples

### Example 1: hsa-mir-21 hairpin is a real pre-miRNA

**Request:**
```json
{"sequence": "UGUCGGGUAGCUUAUCAGACUGAUGUUGACUGUUGAAUCUCAUGGCAACACCAGUCGAUGGGCUGUCUGACA"}
```

**Response (abridged):**
```json
{"found": true, "isHairpinLike": true, "realVotes": 781, "treeCount": 1000, "realProbability": 0.781, "isPreMiRna": true}
```

## Performance

- **Cost:** one MFE fold plus `shuffleCount` shuffle-folds of a sequence of at most 137 nt, i.e. about 1000 folds for the default 999; measured ~9 s in a Debug build for the 71 nt hsa-mir-21 hairpin. Sequences failing the pre-judgement cost one fold. Lower `shuffleCount` for quick checks (the forest was trained with 999).

## See Also

- [classify_pre_mirna](classify_pre_mirna.md)
- [classify_pre_mirna_triplet_svm](classify_pre_mirna_triplet_svm.md)
- [calculate_randfold_pvalue](calculate_randfold_pvalue.md)
