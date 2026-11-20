# clustalw_quick_pairwise_align

ClustalW 2.1 fast k-tuple pairwise stage (-QUICKTREE).

## Overview

| Property | Value |
|----------|-------|
| **Server** | Alignment |
| **Tool Name** | `clustalw_quick_pairwise_align` |
| **Method ID** | `SequenceAligner.ClustalWQuickPairwiseAlign` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

The fast pairwise stage of ClustalW 2.1 (clustalw -QUICKTREE -ALIGN: FastPairwiseAlign::pairAlign; Wilbur & Lipman 1983 k-tuple diagonals), an exact port with the score conversion of pairwiseAlign: the best chain of k-tuple fragments (fragmentScore), the percent score fragmentScore / shorter length x 100 (percent = true; the raw fragment score otherwise) and the guide-tree distance (100 - score) / 100. type 'dna' (default; -KTUPLE 2, -WINDOW 4, -TOPDIAGS 4, -PAIRGAP 5) or 'protein' (1 / 5 / 5 / 3). Not an alignment identity. maxSequenceLength = the longest sequence of the whole input set (ClustalW bounds its fragment list by 2 x 2 x this value and reports partialAlignment beyond it; default the longer of the two). Verified against the installed clustalw 2.1.

## Core Documentation Reference

- Source: [Seqeron.Genomics.Alignment/SequenceAligner.ClustalWQuick.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.ClustalWQuick.cs)
- Algorithm doc: [Multiple_Sequence_Alignment.md](../../../algorithms/Alignment/Multiple_Sequence_Alignment.md)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence1` | string | Yes | First sequence: DNA (A, C, G, T/U, IUB codes) or protein (letters A-Z except J). |
| `sequence2` | string | Yes | Second sequence. |
| `type` | string | No | 'dna' (default) or 'protein'. |
| `ktuple` | integer | No | ClustalW -KTUPLE (DNA 1-4, protein 1-2; default 2 for DNA, 1 for protein). |
| `window` | integer | No | ClustalW -WINDOW (default 4 DNA, 5 protein). |
| `topDiagonals` | integer | No | ClustalW -TOPDIAGS (default 4 DNA, 5 protein). |
| `pairGap` | integer | No | ClustalW -PAIRGAP (default 5 DNA, 3 protein). |
| `percent` | boolean | No | ClustalW -SCORE=PERCENT (true, default) or ABSOLUTE (false). (default: true) |
| `maxSequenceLength` | integer | No | Longest sequence of the whole input set (>= the longer of the two); default the longer of the two. |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `fragmentScore` | integer | Best k-tuple fragment chain score |
| `score` | number | Percent (or absolute) score |
| `distance` | number | (100 - score) / 100 |
| `aligned` | boolean | Score above 0.1 (ClustalW aligns the pair) |
| `partialAlignment` | boolean | The "(Partial alignment)" exit |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence is required |
| 1001 | type must be dna or protein |
| 1002 | KTuple / Window / TopDiagonals / PairGap out of range |

## Examples

### Example 1: Hand run

**User Prompt:**
> Fast ClustalW score of two DNA sequences (longest sequence of the set: 27).

**Expected Tool Call:**
```json
{"tool": "clustalw_quick_pairwise_align", "arguments": {"sequence1": "ACGTACGTTTGACCAGTAGGCATCGAT", "sequence2": "ACGTACGATTGACCAGTAGCCATCGAT", "maxSequenceLength": 27}}
```

**Response:**
```json
{"fragmentScore": 25, "score": 92.5925925925926, "distance": 0.07407407407407404, "aligned": true, "partialAlignment": false}
```

`clustalw -QUICKTREE -ALIGN` prints 92.5926 for this pair.

## See Also

- [clustalw_pairwise_align](clustalw_pairwise_align.md)
