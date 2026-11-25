# clustalw_pairwise_align

ClustalW 2.1 full pairwise stage (FullPairwiseAlign) for DNA or protein pairs.

## Overview

| Property | Value |
|----------|-------|
| **Server** | MultipleAlignment |
| **Tool Name** | `clustalw_pairwise_align` |
| **Method ID** | `Seqeron.Genomics.ClustalW.ClustalWAligner.ClustalWPairwiseAlign` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

The full pairwise stage of ClustalW 2.1 (clustalw -ALIGN: FullPairwiseAlign::pairwiseAlign, the step that gives the guide-tree distances), ported step by step: an affine local (Smith-Waterman-Gotoh) alignment under ClustalW's matrix whose co-optimal choice follows the program's Myers-Miller midpoint rule, with the percent identity (identities / the shorter sequence, truncated as the program prints it) and the guide-tree distance 1 - identities / shorter length. type 'dna' (default; -PWDNAMATRIX IUB (default) or CLUSTALW, -PWGAPOPEN 15, -PWGAPEXT 6.66, -TRANSWEIGHT 0.5; IUB nucleotide codes) or 'protein' (-PWMATRIX gonnet (default, GONNET250), blosum (BLOSUM30), pam (PAM350) or id; -PWGAPOPEN 10, -PWGAPEXT 0.1). Verified against the installed clustalw 2.1 on 11 000+ DNA pairs. start / end are 0-based inclusive coordinates of the aligned region (-1 when no positive-scoring alignment exists). The multi-sequence pipeline is multiple_align_clustalw.

## Core Documentation Reference

- Source: [Seqeron.Genomics.ClustalW/ClustalWAligner.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.ClustalW/ClustalWAligner.cs)
- Algorithm doc: [Multiple_Sequence_Alignment.md](../../../algorithms/Alignment/Multiple_Sequence_Alignment.md)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence1` | string | Yes | First sequence: DNA (A, C, G, T/U, IUB codes) or protein (letters A-Z except J), case-insensitive. |
| `sequence2` | string | Yes | Second sequence. |
| `type` | string | No | 'dna' (default) or 'protein' (ClustalW -TYPE). |
| `gapOpen` | number | No | ClustalW -PWGAPOPEN (a penalty, positive; default 15 for DNA, 10 for protein). |
| `gapExtend` | number | No | ClustalW -PWGAPEXT (default 6.66 for DNA, 0.1 for protein). |
| `transitionWeight` | number | No | DNA only, ClustalW -TRANSWEIGHT: weight of transitions in the pairwise matrix (default 0.5). |
| `dnaMatrix` | string | No | type='dna': ClustalW -PWDNAMATRIX built-in, 'iub' (default) or 'clustalw'. |
| `pwMatrix` | string | No | type='protein': ClustalW -PWMATRIX, 'gonnet' (default), 'blosum', 'pam' or 'id'. |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `score` | integer | Local alignment score |
| `identities` | integer | Identical aligned pairs |
| `percentIdentity` | number | 100 x identities / shorter length |
| `distance` | number | Guide-tree distance 1 - identities / shorter length |
| `clustalWScore` | integer | The truncated percent identity clustalw prints |
| `alignedSequence1` | string | Aligned row 1 |
| `alignedSequence2` | string | Aligned row 2 |
| `start1` | integer | 0-based start (-1 when empty) |
| `end1` | integer | 0-based inclusive end |
| `start2` | integer | Start in sequence2 |
| `end2` | integer | End in sequence2 |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence is required |
| 1001 | type must be dna or protein |
| 1001 | A symbol is not a valid ClustalW residue code |

## Examples

### Example 1: DNA pair

**User Prompt:**
> ClustalW pairwise stage of two DNA sequences.

**Expected Tool Call:**
```json
{"tool": "clustalw_pairwise_align", "arguments": {"sequence1": "ACGTACGTTTGACCAGTAGGCATCGAT", "sequence2": "ACGTACGATTGACCAGTAGCCATCGAT"}}
```

**Response:**
```json
{"score": 23200, "identities": 25, "percentIdentity": 92.5925925925926, "distance": 0.07407407407407407, "clustalWScore": 92, "alignedSequence1": "ACGTACGTTTGACCAGTAGGCATCGAT", "alignedSequence2": "ACGTACGATTGACCAGTAGCCATCGAT", "start1": 0, "end1": 26, "start2": 0, "end2": 26}
```

### Example 2: Protein pair (globin fragments)

**User Prompt:**
> ClustalW pairwise score of two globin fragments.

**Expected Tool Call:**
```json
{"tool": "clustalw_pairwise_align", "arguments": {"sequence1": "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHF", "sequence2": "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKV", "type": "protein"}}
```

**Response:**
```json
{"score": 9147, "identities": 19, "percentIdentity": 40.42553191489362, "distance": 0.5957446808510638, "clustalWScore": 40, "alignedSequence1": "LSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHF", "alignedSequence2": "LTPEEKSAVTALWGKV--NVDEVGGEALGRLLVVYPWTQRFFESF", "start1": 2, "end1": 46, "start2": 3, "end2": 45}
```

clustalw 2.1 prints 40 for this pair.

## See Also

- [clustalw_quick_pairwise_align](clustalw_quick_pairwise_align.md)
- [multiple_align_clustalw](multiple_align_clustalw.md) - the whole pipeline
