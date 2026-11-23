# msa_conservation

The per-column conservation scores of Capra and Singh (2007): Jensen-Shannon divergence, Shannon / property entropy, relative entropy, sum of pairs, with gap penalty, sequence weights and a window.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `msa_conservation` |
| **Method ID** | `ConservationScores` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

The conservation score of every column of a protein alignment as in the program of Capra and Singh (2007, Bioinformatics 23:1875): method js_divergence (Jensen-Shannon divergence in bits of the column distribution from the background, 0 for a column that looks like the background, close to 1 for a conserved column of a rare residue; the default method), shannon_entropy (the information score 1 - h with h = sum f ln f / ln(min(21, rows)): 1 conserved, 0 maximal entropy), property_entropy (the same over the groups AVLIMC, FWYH, STNQ, KR, DE, GP and the gap), relative_entropy (Kullback-Leibler divergence from the background / ln 20) or sum_of_pairs (the weighted mean BLOSUM62 score of all pairs of residues of the column). Frequencies carry the pseudo count 1e-7; gapPenalty (default true) multiplies each score by the weighted fraction of non-gap rows; sequenceWeights weights the rows by the position-based weights of Henikoff and Henikoff (1994); windowSize smooths the scores over that many columns on either side (the score of a column 0.5 of itself + 0.5 of the mean of its non-negative neighbours); background blosum62 (default) or swissprot, or a JSON array of 20 / 21 numbers in the order A R N D C Q E G H I L K M F P S T W Y V (-). Rows are upper-cased and read as the program's readers do (B = D, Z = Q, X and . = gap); a letter that is not one of the 20 amino acids is ignored by the entropies and is an error for sum_of_pairs and the weights. VERIFICATION: equal to an independent implementation of the formulas of the paper (every score of the unit tests to 1e-13) and to hand-derived values (the program itself is a Python 2 CGI script that is not run here).

## Core Documentation Reference

- Source: [ConservationScores.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/ConservationScores.cs)
- Algorithm doc: [Multiple_Sequence_Alignment.md](../../../algorithms/Alignment/Multiple_Sequence_Alignment.md)
- Oracle: [conservation_scores_reference.py](../../../../scripts/alignment/conservation_scores_reference.py)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequences` | array<string> | Yes | At least one aligned protein sequence; equal lengths, gaps '-'. |
| `method` | string | No | js_divergence (default), shannon_entropy, property_entropy, relative_entropy or sum_of_pairs. |
| `background` | string | No | Background distribution: blosum62 (default), swissprot, or a JSON array of 20 / 21 numbers. |
| `sequenceWeights` | boolean | No | Weight the rows (Henikoff and Henikoff 1994). (default: false) |
| `gapPenalty` | boolean | No | Multiply every score by the weighted fraction of non-gap rows (default true). |
| `windowSize` | integer | No | Window size on either side (default 0: no smoothing). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `method` | string | The method |
| `windowSize` | integer | The window size |
| `scores` | number[] | The score of every column |
| `sequenceWeights` | number[] | The weights of the rows (empty without sequenceWeights) |
| `columns` | string[] | The columns as text |

## Errors

| Code | Message |
|------|---------|
| 1001 | At least one sequence is required |
| 1001 | Unknown method or background |
| 1001 | Rows of unequal length, a letter outside the 20 amino acids in sum_of_pairs or the weights |

## Examples

### Example 1: Five conserved columns

**User Prompt:**
> Score the conservation of these columns.

**Expected Tool Call:**
```json
{"tool": "msa_conservation", "arguments": {"sequences": ["AAAAA", "AAAAC", "ACDEF", "AAA--", "ARNDC"], "method": "js_divergence"}}
```

**Response:**
```json
{"method": "JsDivergence", "scores": "[...]"}
```

The first column is five A (score 0.79906: the divergence of a conserved alanine column from the BLOSUM62 background, A 7.8 %).

## See Also

- [nucmer_clusters](nucmer_clusters.md) - the anchor-clustering stage of nucmer
