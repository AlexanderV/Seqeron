# multiple_align_clustalw

ClustalW 2.1's multiple alignment of DNA (clustalw -TYPE=DNA -ALIGN), ported step by step and verified row-for-row against the installed clustalw 2.1: pairwise distances (distance='clustalw' = the full affine Smith-Waterman-Gotoh stage, 'clustalw-quick' = the fast k-tuple stage of -QUICKTREE with ktuple / window / topDiagonals / pairGap), ClustalW's neighbour-joining guide tree (rooted by Tree::reRoot), sequence weights derived from the tree branch lengths (Thompson, Higgins & Gibson 1994; noWeights = equal weights), delayed divergent sequences (maxDiv, default 30: sequences with at most that percent identity to every other sequence are aligned last), and profile-profile alignment with position-specific gap penalties (positions that already hold gaps get a reduced opening penalty and half the extension penalty, end gaps are free) on the IUB weight matrix shifted to non-negative scores (negative = keep it signed) with transitions weighted by transitionWeight (clustalw -TRANSWEIGHT, default 0.5) and Myers-Miller linear-space alignment. gapOpen / gapExtend are ClustalW's -GAPOPEN 15 / -GAPEXT 6.66 (the DNA residue-specific, hydrophilic and gap-separation penalties of clustalw are protein-only). Rows are returned in input order. match / mismatch / scoreGapOpen / scoreGapExtend define only the reported totalScore (natural affine sum of pairs, Altschul 1989).

## Overview

| Property | Value |
|----------|-------|
| **Server** | Alignment |
| **Tool Name** | `multiple_align_clustalw` |
| **Method ID** | `SequenceAligner.MultipleAlignClustalW` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

ClustalW 2.1's multiple alignment of DNA (clustalw -TYPE=DNA -ALIGN), ported step by step and verified row-for-row against the installed clustalw 2.1: pairwise distances (distance='clustalw' = the full affine Smith-Waterman-Gotoh stage, 'clustalw-quick' = the fast k-tuple stage of -QUICKTREE with ktuple / window / topDiagonals / pairGap), ClustalW's neighbour-joining guide tree (rooted by Tree::reRoot), sequence weights derived from the tree branch lengths (Thompson, Higgins & Gibson 1994; noWeights = equal weights), delayed divergent sequences (maxDiv, default 30: sequences with at most that percent identity to every other sequence are aligned last), and profile-profile alignment with position-specific gap penalties (positions that already hold gaps get a reduced opening penalty and half the extension penalty, end gaps are free) on the IUB weight matrix shifted to non-negative scores (negative = keep it signed) with transitions weighted by transitionWeight (clustalw -TRANSWEIGHT, default 0.5) and Myers-Miller linear-space alignment. gapOpen / gapExtend are ClustalW's -GAPOPEN 15 / -GAPEXT 6.66 (the DNA residue-specific, hydrophilic and gap-separation penalties of clustalw are protein-only). Rows are returned in input order. match / mismatch / scoreGapOpen / scoreGapExtend define only the reported totalScore (natural affine sum of pairs, Altschul 1989).

## Core Documentation Reference

- Source: [Seqeron.Genomics.Alignment/SequenceAligner.ClustalWMultiple.cs#L68](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.ClustalWMultiple.cs#L68)
- Algorithm doc: [Multiple_Sequence_Alignment.md](../../../algorithms/Alignment/Multiple_Sequence_Alignment.md) (section 4.6)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequences` | array<string> | Yes | DNA sequences to align (at least one; A, C, G, T). |
| `gapOpen` | number | No | ClustalW -GAPOPEN: gap-opening penalty of the profile alignment (default 15; a penalty, positive). (default: 15.0) |
| `gapExtend` | number | No | ClustalW -GAPEXT: gap-extension penalty of the profile alignment (default 6.66). (default: 6.66) |
| `distance` | string | No | Pairwise stage that gives the guide-tree distances: 'clustalw' (default; full DP, clustalw -ALIGN) or 'clustalw-quick' (fast k-tuple mode, clustalw -QUICKTREE). (default: clustalw) |
| `pwGapOpen` | number | No | ClustalW -PWGAPOPEN: pairwise gap-opening penalty of the full pairwise stage (default 15). (default: 15.0) |
| `pwGapExtend` | number | No | ClustalW -PWGAPEXT: pairwise gap-extension penalty of the full pairwise stage (default 6.66). (default: 6.66) |
| `transitionWeight` | number | No | ClustalW -TRANSWEIGHT: weight of transitions (A-G, C-T/U) in the profile matrix (default 0.5; 0 = as mismatches, 1 = as matches). (default: 0.5) |
| `maxDiv` | integer | No | ClustalW -MAXDIV (default 30): a sequence whose highest percent identity to any other sequence is at most this value is delayed (aligned last, after the rest). (default: 30) |
| `noWeights` | boolean | No | ClustalW -NOWEIGHTS: no tree-derived sequence weights. (default: false) |
| `negative` | boolean | No | ClustalW -NEGATIVE: keep the substitution matrix signed instead of shifting it to non-negative scores. (default: false) |
| `ktuple` | integer | No | distance='clustalw-quick' only: k-tuple (word) size, 1-4 for DNA (clustalw -KTUPLE, default 2). |
| `window` | integer | No | distance='clustalw-quick' only: diagonals scanned around each top diagonal (clustalw -WINDOW, default 4). |
| `topDiagonals` | integer | No | distance='clustalw-quick' only: number of best diagonals used (clustalw -TOPDIAGS, default 4). |
| `pairGap` | integer | No | distance='clustalw-quick' only: gap penalty of the fragment chain (clustalw -PAIRGAP; default 5, or ktuple + 4 when ktuple is given, as the clustalw command line does). |
| `match` | integer | No | Score for matching base in the reported totalScore only. (default: 1) |
| `mismatch` | integer | No | Score for mismatching base in the reported totalScore only. (default: -1) |
| `scoreGapOpen` | integer | No | Gap-open score of the reported totalScore only (negative = penalty). (default: -2) |
| `scoreGapExtend` | integer | No | Gap-extend score of the reported totalScore only (negative = penalty). (default: -1) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `alignedSequences` | array<string> | The aligned rows of the MSA (input order). |
| `consensus` | string | Majority-vote consensus of the aligned rows (gaps vote; a nucleotide beats a gap on ties). |
| `totalScore` | integer | Natural affine sum-of-pairs score of the alignment under match / mismatch / scoreGapOpen / scoreGapExtend (ClustalW itself reports no such score). |

## Errors

| Code | Message |
|------|---------|
| 1010 | Invalid DNA at sequences[i] |
| 1002 | gapOpen, gapExtend, pwGapOpen, pwGapExtend, transitionWeight must be non-negative; distance must be clustalw or clustalw-quick |

## Examples

### Example 1: Five related 30-mers

**Tool Call:**
```json
{
  "tool": "multiple_align_clustalw",
  "arguments": {
    "sequences": [
      "ACGTACGTTTGACCAGTAGGCATCGATGGA",
      "ACGTACGATTGACCAGTAGCCATCGATGGA",
      "TTGACCAGTAGGCATCGATCCCAAGT",
      "ACGTAGGATTGACCGTAGCCATCGATGGA",
      "GGGTTTACCAGTAGGCATCGATGGAAACCT"
    ]
  }
}
```

**Response:**
```json
{
  "alignedSequences": [
    "ACGTACGTTTGACCAGTAGGCATCGATGGA-----",
    "ACGTACGATTGACCAGTAGCCATCGATGGA-----",
    "--------TTGACCAGTAGGCATCGATCCCAAGT-",
    "ACGTAGGATTGACC-GTAGCCATCGATGGA-----",
    "-----GGGTTTACCAGTAGGCATCGATGGAAACCT"
  ],
  "consensus": "ACGTACGATTGACCAGTAGGCATCGATGGA-----",
  "totalScore": 77
}
```

## Worked Example

The rows are exactly what `clustalw -TYPE=DNA -ALIGN -OUTPUT=FASTA -OUTORDER=INPUT` writes for these five sequences (ClustalW 2.1, run here): the five sequences are aligned along ClustalW's neighbour-joining guide tree with tree-derived weights and position-specific gap penalties; sequences with at most 30 % identity to every other would be delayed and added last.

## References

- Algorithm source: [Seqeron.Genomics.Alignment/SequenceAligner.ClustalWMultiple.cs#L68](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.ClustalWMultiple.cs#L68)
- Binding: [AlignmentTools.cs](../../../../src/Seqeron/Mcp/Seqeron.Mcp.Alignment/Tools/AlignmentTools.cs)
