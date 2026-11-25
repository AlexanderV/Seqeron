# predict_drosha_dicer_cleavage

Predict the Drosha (basal) and Dicer (apical) cleavage sites and the 5p/3p products of a pri-miRNA hairpin.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Annotation |
| **Tool Name** | `predict_drosha_dicer_cleavage` |
| **Method ID** | `MiRnaAnalyzer.PredictDroshaDicerCleavage` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Applies the published measuring rules to the folded stem of a pri-miRNA: Drosha cuts ~11 bp from the basal stem–ssRNA junction (`droshaRule = BasalJunction`, Han 2006) or ~22 bp from the loop–stem junction (`ApicalJunction`, Zeng 2005), leaving a 2-nt 3′ overhang; Dicer then cuts ~22 nt from the 5′ end (`dicerRule = FivePrimeCounting`, default, Park 2011) or from the 3′ end (`ThreePrimeCounting`). Base pairs are counted up the unbranched stem from `basalJunction` (bulged nucleotides do not count). Also reports the CNNC (+16..+18 nt downstream of the 3p cut) and basal UG motifs (Auyeung 2013). `dotBracket` is optional (balanced, same length as the sequence); without it the canonical Turner 2004 MFE fold is used. `predicted = false` when `basalJunction` is not the 5′-arm base of a stem pair or the stem is too short. Coordinates are 0-based, spans inclusive.

## Core Documentation Reference

- Source: [MiRnaAnalyzer.cs#L3111](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/MiRnaAnalyzer.cs#L3111)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | pri-miRNA nucleotide sequence (RNA or DNA) |
| `basalJunction` | integer | Yes | 0-based index of the 5′-arm base of the first stem base pair |
| `dotBracket` | string | No | Optional dot-bracket structure; default = MFE fold |
| `dicerRule` | string | No | `FivePrimeCounting` (default) or `ThreePrimeCounting` |
| `droshaRule` | string | No | `BasalJunction` (default) or `ApicalJunction` |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `predicted` | boolean | False when no prediction is possible |
| `cleavage.droshaCut5Prime / droshaCut3Prime` | integer | Drosha cut positions (pre-miRNA ends) |
| `cleavage.matureStart / matureEnd` | integer | 5p product span |
| `cleavage.starStart / starEnd` | integer | 3p product span |
| `cleavage.matureSequence / starSequence` | string | 5p / 3p products |
| `cleavage.basalJunction / apicalJunction` | integer | Junction pair 5′-arm indices |
| `cleavage.structure` | string | Dot-bracket the stem was measured on |
| `cleavage.threePrimeOverhang` | integer | RNase III 3′ overhang (2) |
| `cleavage.hasCnncMotif / hasBasalUgMotif` | boolean | Auyeung 2013 sequence motifs |
| `cleavage.dicerRule / droshaRule` | string | Rules applied |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1001 | basalJunction {n} is out of range for sequence of length {L} |
| 1001 | Unknown dicerRule / droshaRule (valid values listed) |
| 1001 | Dot-bracket must be balanced and as long as the sequence. |

## Examples

### Example 1: 108-nt model hairpin

5-nt flank, 35-bp stem, 8-nt loop, 35-bp stem, 25-nt flank with `basalJunction = 5`: the 12th pair is (16, 71), so Drosha cuts at 16 and 73; the 5p product is [16..37] and the 3p product [52..73].

**Request:**
```json
{"sequence": "AAAAAGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGAAAAAAAACCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAAAAAAAAAA", "basalJunction": 5, "dotBracket": ".....(((((((((((((((((((((((((((((((((((........)))))))))))))))))))))))))))))))))))........................."}
```

**Response (abridged):**
```json
{"predicted": true, "cleavage": {"droshaCut5Prime": 16, "droshaCut3Prime": 73, "matureStart": 16, "matureEnd": 37, "starStart": 52, "starEnd": 73}}
```

## Performance

- **Time Complexity:** O(n³) when the MFE fold is computed (no `dotBracket`), otherwise O(n)
- **Space Complexity:** O(n²)

## See Also

- [find_pre_mirna_hairpins](find_pre_mirna_hairpins.md)
- [create_mirna](create_mirna.md)
