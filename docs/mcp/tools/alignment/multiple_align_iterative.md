# multiple_align_iterative

Progressive MSA followed by MUSCLE-style tree-dependent restricted-partitioning refinement (Edgar 2004, stage 3): each guide-tree edge splits the alignment into two profiles that are re-aligned and kept only if the sum-of-pairs score strictly improves; repeated until convergence or maxIterations passes. The result's sum-of-pairs score is never below multiple_align_progressive's. Default linear gap model (gapOpen ignored); gapModel='affine' refines under the natural affine sum-of-pairs objective with the affine profile DP. distance='clustalw' selects ClustalW 2.1's pairwise stage for the guide tree (see multiple_align_progressive).

## Overview

| Property | Value |
|----------|-------|
| **Server** | Alignment |
| **Tool Name** | `multiple_align_iterative` |
| **Method ID** | `SequenceAligner.MultipleAlignIterative` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Progressive MSA followed by MUSCLE-style tree-dependent restricted-partitioning refinement (Edgar 2004, stage 3): each guide-tree edge splits the alignment into two profiles that are re-aligned and kept only if the sum-of-pairs score strictly improves; repeated until convergence or maxIterations passes. The result's sum-of-pairs score is never below multiple_align_progressive's. Default linear gap model (gapOpen ignored); gapModel='affine' refines under the natural affine sum-of-pairs objective with the affine profile DP. distance='clustalw' selects ClustalW 2.1's pairwise stage for the guide tree (see multiple_align_progressive).

## Core Documentation Reference

- Source: [Seqeron.Genomics.Alignment/SequenceAligner.cs#L2778](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.cs#L2778)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequences` | array<string> | Yes | DNA sequences to align (at least one). |
| `match` | integer | No | Score for matching base. (default: 1) |
| `mismatch` | integer | No | Score for mismatching base. (default: -1) |
| `gapOpen` | integer | No | Gap-open penalty. Used only with gapModel='affine' (a gap of length k then costs gapOpen + k * gapExtend, Gotoh); ignored by the default linear model. (default: -2) |
| `gapExtend` | integer | No | Gap penalty per gap position (linear model: every gap position; affine model: the per-position part of gapOpen + k * gapExtend). (default: -1) |
| `maxIterations` | integer | No | Upper bound on full refinement passes (>= 1). (default: 16) |
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
| 1002 | maxIterations must be >= 1 |

## Examples

### Example 1: Gap relocation

**Tool Call:**
```json
{
  "tool": "multiple_align_iterative",
  "arguments": {
    "sequences": [
      "TCA",
      "TCT",
      "GCAC",
      "TAAC"
    ]
  }
}
```

**Response:**
```json
{
  "alignedSequences": [
    "TCA-",
    "TCT-",
    "GCAC",
    "TAAC"
  ],
  "consensus": "TCAC",
  "totalScore": -3
}
```

## Worked Example

The progressive seed T--CA / T--CT / GCAC- / TAAC- scores -7; refinement moves the gap block to the last column: TCA- / TCT- / GCAC / TAAC, SP -3.

## References

- Algorithm source: [Seqeron.Genomics.Alignment/SequenceAligner.cs#L2778](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.cs#L2778)
- Binding: [AlignmentTools.cs](../../../../src/Seqeron/Mcp/Seqeron.Mcp.Alignment/Tools/AlignmentTools.cs)
