# clustalw_guide_tree

ClustalW 2.1 neighbour-joining guide-tree file (.dnd) from pairwise distances.

## Overview

| Property | Value |
|----------|-------|
| **Server** | MultipleAlignment |
| **Tool Name** | `clustalw_guide_tree` |
| **Method ID** | `Seqeron.Genomics.ClustalW.ClustalWAligner.ClustalWGuideTreeText` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

The guide-tree file clustalw writes for its neighbour-joining clustering (UnRootedClusterTree::treeFromDistMatrix + ClusterTreeOutput::printPhylipTree / twoWaySplit) of the sequences `names` from their pairwise distances: '(a:d,b:d);' for two sequences (d = half the distance, five significant digits, as the program writes it), else the unrooted tree as three subtrees under the first parenthesis, every node on its own lines, branch lengths with five decimals (the text of the .dnd file). The distances are ClustalW's own (the fractional differences clustalw computes in its pairwise stage; see clustalw_pairwise_align / clustalw_quick_pairwise_align). VERIFICATION: the neighbour-joining and printing routines are the ones multiple_align_clustalw runs; their guide-tree text equals the .dnd files of the installed clustalw 2.1 (committed guide-tree oracle, random DNA sets).

## Core Documentation Reference

- Source: [Seqeron.Genomics.ClustalW/ClustalWAligner.Trees.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.ClustalW/ClustalWAligner.Trees.cs)
- Algorithm doc: [Multiple_Sequence_Alignment.md](../../../algorithms/Alignment/Multiple_Sequence_Alignment.md)
- Oracle: `scripts/alignment/generate_clustalw_guidetree_oracle.py` -> `tests/Seqeron/Seqeron.Genomics.Tests/TestData/Alignment/` (clustalw_guidetree_oracle.tsv)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `names` | array of string | Yes | Sequence names in file order (at least two). |
| `distances` | array of array of number | Yes | ClustalW's pairwise distances: a square symmetric matrix with one row and one column per name (diagonal 0). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `tree` | string | Text of the guide-tree (.dnd) file. |
| `sequences` | integer | Number of sequences. |

## Errors

| Code | Message |
|------|---------|
| 1001 | a sequence argument is null or empty / names and matrix have the wrong shape |
| 1002 | a gap penalty is outside the range of the program / unknown matrix |

## Examples

### Example 1: Four sequences from ClustalW distances

**Tool Call:**
```json
{
  "tool": "clustalw_guide_tree",
  "arguments": {
    "names": [
      "a",
      "b",
      "c",
      "d"
    ],
    "distances": [
      [
        0,
        0.1,
        0.5,
        0.6
      ],
      [
        0.1,
        0,
        0.55,
        0.62
      ],
      [
        0.5,
        0.55,
        0,
        0.2
      ],
      [
        0.6,
        0.62,
        0.2,
        0
      ]
    ]
  }
}
```

**Response:**
```json
{
  "tree": "(\n(\na:0.03250,\nb:0.06750)\n:0.41750,\nc:0.05750,\nd:0.14250);\n",
  "sequences": 4
}
```

## See Also

- [multiple_align_clustalw](multiple_align_clustalw.md)
- [clustalw_pairwise_align](clustalw_pairwise_align.md)
- [clustalw_command](clustalw_command.md)
