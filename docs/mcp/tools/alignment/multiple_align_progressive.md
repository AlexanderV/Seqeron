# multiple_align_progressive

Guide-tree progressive MSA (Feng & Doolittle 1987 / Clustal-style): the ClustalW pairwise distance of every pairwise alignment, d = 1 − (identical aligned residue pairs) / (length of the shorter sequence) (Thompson, Higgins & Gibson 1994; pairalign.c tmat = (100 − 100·count/MIN(len1, len2)) / 100), gives the distance matrix, a guide tree (UPGMA by default, guideTree='nj' = ClustalW's neighbour-joining tree) orders the merges, and profiles are aligned by profile–profile dynamic programming with sum-of-pairs column scores (once a gap, always a gap). Default linear gap model (gapOpen ignored); gapModel='affine' uses a Gotoh three-state profile DP (occupancy-scaled gapOpen/gapExtend, ClustalW/MUSCLE style) and reports the natural affine sum-of-pairs score. distance='clustalw' takes the guide-tree distances from ClustalW 2.1's own affine pairwise stage instead (percent identity of clustalw -ALIGN). Returns rows, majority consensus and sum-of-pairs score.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Alignment |
| **Tool Name** | `multiple_align_progressive` |
| **Method ID** | `SequenceAligner.MultipleAlignProgressive` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Guide-tree progressive MSA (Feng & Doolittle 1987 / Clustal-style): the ClustalW pairwise distance of every pairwise alignment, d = 1 − (identical aligned residue pairs) / (length of the shorter sequence) (Thompson, Higgins & Gibson 1994; pairalign.c tmat = (100 − 100·count/MIN(len1, len2)) / 100), gives the distance matrix, a guide tree (UPGMA by default, guideTree='nj' = ClustalW's neighbour-joining tree) orders the merges, and profiles are aligned by profile–profile dynamic programming with sum-of-pairs column scores (once a gap, always a gap). Default linear gap model (gapOpen ignored); gapModel='affine' uses a Gotoh three-state profile DP (occupancy-scaled gapOpen/gapExtend, ClustalW/MUSCLE style) and reports the natural affine sum-of-pairs score. distance='clustalw' takes the guide-tree distances from ClustalW 2.1's own affine pairwise stage instead (percent identity of clustalw -ALIGN). Returns rows, majority consensus and sum-of-pairs score.

## Core Documentation Reference

- Source: [Seqeron.Genomics.Alignment/SequenceAligner.cs#L2288](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.cs#L2288)

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
| `guideTree` | string | No | Guide-tree method: 'upgma' (default; Sokal & Michener 1958, the Feng-Doolittle / CLUSTAL-1 choice) or 'nj' (neighbour-joining, Saitou & Nei 1987 - ClustalW 2.1's default tree, T-Coffee/Clustal-style; the unrooted NJ tree is rooted and ordered exactly as clustalw does: Tree::reRoot, left subtree first). With distance='clustalw' and guideTree='nj' the tree equals the .dnd file of clustalw -ALIGN. Both are built by the canonical phylogenetics builders. (default: upgma) |

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

### Example 1: One deletion

**Tool Call:**
```json
{
  "tool": "multiple_align_progressive",
  "arguments": {
    "sequences": [
      "ACGT",
      "ACGT",
      "AGT"
    ]
  }
}
```

**Response:**
```json
{
  "alignedSequences": [
    "ACGT",
    "ACGT",
    "A-GT"
  ],
  "consensus": "ACGT",
  "totalScore": 8
}
```

## Worked Example

ACGT, ACGT merge first; AGT is aligned to that profile with the gap opposite C: SP = +3 (A) - 1 (C,C,-) + 3 (G) + 3 (T) = 8.

## References

- Algorithm source: [Seqeron.Genomics.Alignment/SequenceAligner.cs#L2288](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.cs#L2288)
- Binding: [AlignmentTools.cs](../../../../src/Seqeron/Mcp/Seqeron.Mcp.Alignment/Tools/AlignmentTools.cs)
