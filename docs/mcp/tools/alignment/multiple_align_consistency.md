# multiple_align_consistency

Consistency-based MSA (T-Coffee; Notredame, Higgins & Heringa 2000): a primary library of global (Needleman–Wunsch) and local (Smith–Waterman) pairwise alignments weighted by percent identity, extended by the triplet consistency transformation (direct weight + Σ min(W1, W2) over third sequences), then progressively aligned along a UPGMA guide tree with the extended-library weights as position-specific scores and zero gap penalties (T-Coffee default). The scoring parameters seed the pairwise library (gapModel='affine': Gotoh global and Smith–Waterman–Gotoh local alignments with gapOpen + k·gapExtend) and define totalScore, the sum of pairs of the result (linear column score, or natural affine costs with gapModel='affine'). distance='clustalw' takes the guide-tree distances from ClustalW 2.1's own pairwise stage.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Alignment |
| **Tool Name** | `multiple_align_consistency` |
| **Method ID** | `SequenceAligner.MultipleAlignConsistency` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Consistency-based MSA (T-Coffee; Notredame, Higgins & Heringa 2000): a primary library of global (Needleman–Wunsch) and local (Smith–Waterman) pairwise alignments weighted by percent identity, extended by the triplet consistency transformation (direct weight + Σ min(W1, W2) over third sequences), then progressively aligned along a UPGMA guide tree with the extended-library weights as position-specific scores and zero gap penalties (T-Coffee default). The scoring parameters seed the pairwise library (gapModel='affine': Gotoh global and Smith–Waterman–Gotoh local alignments with gapOpen + k·gapExtend) and define totalScore, the sum of pairs of the result (linear column score, or natural affine costs with gapModel='affine'). distance='clustalw' takes the guide-tree distances from ClustalW 2.1's own pairwise stage.

## Core Documentation Reference

- Source: [Seqeron.Genomics.Alignment/SequenceAligner.cs#L3104](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.cs#L3104)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequences` | array<string> | Yes | DNA sequences to align (at least one). |
| `match` | integer | No | Score for matching base. (default: 1) |
| `mismatch` | integer | No | Score for mismatching base. (default: -1) |
| `gapOpen` | integer | No | Gap-open penalty. Used only with gapModel='affine' (a gap of length k then costs gapOpen + k * gapExtend, Gotoh); ignored by the default linear model. (default: -2) |
| `gapExtend` | integer | No | Gap penalty per gap position (linear model: every gap position; affine model: the per-position part of gapOpen + k * gapExtend). (default: -1) |
| `distance` | string | No | Guide-tree distance: 'linear' (default) = 1 - identities / min(length) of the repository's global alignment under the given scores (linear NW, or Gotoh when gapModel='affine'); 'clustalw' = the same formula on ClustalW 2.1's own pairwise stage (affine local Smith-Waterman-Gotoh alignment with its IUB matrix x100 and gap penalties 15 / 6.66, the 'slow' full-DP mode of clustalw -ALIGN; the match/mismatch/gap parameters are then used only for the profile alignment). (default: linear) |
| `gapModel` | string | No | 'linear' (default): every gap position costs gapExtend and gapOpen is ignored; totalScore is the column sum of pairs (gap-gap columns 0). 'affine': gap runs cost gapOpen + k * gapExtend; the pairwise / profile alignments use affine (Gotoh) dynamic programming and totalScore is the sum of pairs with natural affine gap costs (Altschul 1989: each pair of rows is projected out of the alignment, gap-gap columns skipped, every gap run of the projection pays gapOpen once). For two sequences the affine result is an optimal pairwise alignment; for three or more sequences every method is a heuristic. (default: linear) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `alignedSequences` | array<string> | The aligned rows of the MSA. |
| `consensus` | string | Majority-vote consensus of the aligned rows (gaps vote; a nucleotide beats a gap on ties). |
| `totalScore` | integer | Sum-of-pairs score of the MSA: linear model = column score (match/mismatch from the matrix, residue-gap = gapExtend, gap-gap = 0); gapModel='affine' = natural affine gap costs (every gap run of each pair projection pays gapOpen once). |

## Errors

| Code | Message |
|------|---------|
| 1010 | Invalid DNA in sequences |

## Examples

### Example 1: Three short sequences

**Tool Call:**
```json
{
  "tool": "multiple_align_consistency",
  "arguments": {
    "sequences": [
      "CTT",
      "ACT",
      "CAA"
    ]
  }
}
```

**Response:**
```json
{
  "alignedSequences": [
    "-CTT",
    "AC-T",
    "-CAA"
  ],
  "consensus": "-CAT",
  "totalScore": -3
}
```

## Worked Example

Residue pairs supported by other sequences through the triplet extension (e.g. S1.C-S2.C through S0.C: min(200, 133) = 133) outweigh unsupported pairs; no substitution matrix or gap penalty is used in the progressive DP.

## References

- Algorithm source: [Seqeron.Genomics.Alignment/SequenceAligner.cs#L3104](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.cs#L3104)
- Binding: [AlignmentTools.cs](../../../../src/Seqeron/Mcp/Seqeron.Mcp.Alignment/Tools/AlignmentTools.cs)
