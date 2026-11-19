# sum_of_pairs_score

Sum-of-pairs score of a multiple alignment (linear / affine gaps, simple scores or a substitution matrix).

## Overview

| Property | Value |
|----------|-------|
| **Server** | Alignment |
| **Tool Name** | `sum_of_pairs_score` |
| **Method ID** | `SequenceAligner.SumOfPairsScore` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Sum-of-pairs score of a multiple alignment (equal-length rows, '-' = gap). gapModel 'linear' (default): column sum of pairs - residue pairs match / mismatch, residue against gap gapExtend, gap against gap 0 (the totalScore of the linear aligners); 'affine': the natural affine objective (Altschul 1989) - for every pair of rows the gap-gap columns are removed and every maximal gap run costs gapOpen + length x gapExtend. Scores: match / mismatch / gapOpen / gapExtend (integers; negative = penalty); with matrix (BLOSUM62, PAM250, EDNAFULL or any EMBOSS matrix) the residue pairs are scored by that table (protein / IUPAC rows, residues case-insensitive) and the gap scores may be fractional - the score is then a real number. This is the objective the exact (Carrillo-Lipman) and affine aligners optimise; MUSCLE's own SP score is muscle_msa_command mode 'spscore'.

## Core Documentation Reference

- Source: [Seqeron.Genomics.Alignment/SequenceAligner.MultipleAffine.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.MultipleAffine.cs)
- Algorithm doc: [Multiple_Sequence_Alignment.md](../../../algorithms/Alignment/Multiple_Sequence_Alignment.md)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `alignedRows` | array<string> | Yes | Equal-length aligned rows, '-' = gap. |
| `match` | integer | No | Score for a matching residue pair (ignored with matrix). |
| `mismatch` | integer | No | Score for a mismatching residue pair (ignored with matrix). |
| `gapOpen` | number | No | Gap-open score added once per gap run (affine model only; negative = penalty; integer unless matrix is given). |
| `gapExtend` | number | No | Gap score per gap position (negative = penalty; integer unless matrix is given). |
| `gapModel` | string | No | Gap model: 'linear' (default) or 'affine'. |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `score` | number | Sum of pairs |
| `gapModel` | string | Linear | Affine |
| `scoring` | string | Scoring used |

## Errors

| Code | Message |
|------|---------|
| 1001 | alignedRows must contain at least one row |
| 1001 | All aligned rows must have the same length |
| 1001 | gapModel must be linear or affine |
| 1001 | gapOpen / gapExtend is not an integer (fractional gap scores need a matrix) |

## Examples

### Example 1: Linear

**User Prompt:**
> Sum of pairs of AC-, AC-, ACA (match 2, mismatch -1, gap -2).

**Expected Tool Call:**
```json
{"tool": "sum_of_pairs_score", "arguments": {"alignedRows": ["AC-", "AC-", "ACA"], "match": 2, "mismatch": -1, "gapOpen": -3, "gapExtend": -2}}
```

**Response:**
```json
{"score": 8, "gapModel": "Linear", "scoring": "match 2 mismatch -1"}
```

Pairs (0,1) 4, (0,2) 2, (1,2) 2.

### Example 2: Affine

**User Prompt:**
> The same alignment under the affine objective.

**Expected Tool Call:**
```json
{"tool": "sum_of_pairs_score", "arguments": {"alignedRows": ["AC-", "AC-", "ACA"], "match": 2, "mismatch": -1, "gapOpen": -3, "gapExtend": -2, "gapModel": "affine"}}
```

**Response:**
```json
{"score": 2, "gapModel": "Affine", "scoring": "match 2 mismatch -1"}
```

Each of the two pairs with the gap column costs gapOpen + gapExtend = -5: 4 + (4 - 5) + (4 - 5) = 2.

## See Also

- [multiple_align_exact](multiple_align_exact.md) - optimises this objective
- [muscle_msa_command](muscle_msa_command.md) - MUSCLE SP score (mode spscore)
