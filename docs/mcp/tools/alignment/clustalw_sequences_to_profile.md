# clustalw_sequences_to_profile

ClustalW 2.1 alignment of sequences to an existing profile (clustalw -SEQUENCES).

## Overview

| Property | Value |
|----------|-------|
| **Server** | Alignment |
| **Tool Name** | `clustalw_sequences_to_profile` |
| **Method ID** | `SequenceAligner.AlignClustalWSequencesToProfile` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

ClustalW 2.1's sequence-to-profile alignment (clustalw -PROFILE1=a -PROFILE2=b -SEQUENCES, Clustal::sequencesAlignToProfile / MSA::seqsAlignToProfile), ported step by step and verified row-for-row against the installed clustalw 2.1 on 600 DNA and protein sets: the sequences are added to the existing alignment one at a time, the one with the highest identity to the aligned set first. Distances: the identities inside the profile (countid), the full (distance='clustalw') or fast k-tuple (distance='clustalw-quick', ktuple / window / topDiagonals / pairGap) pairwise stage between every other pair of rows (computed on the rows with their gaps; pwGapOpen / pwGapExtend / pwMatrix), then the neighbour-joining tree of all rows (or UPGMA) weights them; each sequence is aligned with weights = tree weight x identity to it and Myers-Miller profile alignment with position-specific gap penalties. Other parameters as multiple_align_clustalw. Rows are returned: profile, then sequences, in input order.

## Core Documentation Reference

- Source: [Seqeron.Genomics.Alignment/SequenceAligner.ClustalWProfile.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.ClustalWProfile.cs)
- Algorithm doc: [Multiple_Sequence_Alignment.md](../../../algorithms/Alignment/Multiple_Sequence_Alignment.md)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `profile` | array<string> | Yes | Rows of the existing alignment (at least one; letters and '-'). |
| `sequences` | array<string> | Yes | The sequences to add (at least one; gaps '-' are kept as in the program). |
| `type` | string | No | 'dna' (default) or 'protein' (ClustalW -TYPE). (default: dna) |
| `gapOpen` | number | No | ClustalW -GAPOPEN (default 15 for DNA, 10 for protein). |
| `gapExtend` | number | No | ClustalW -GAPEXT (default 6.66 for DNA, 0.2 for protein). |
| `distance` | string | No | Pairwise stage: 'clustalw' (default; full DP) or 'clustalw-quick' (-QUICKTREE). (default: clustalw) |
| `pwGapOpen` | number | No | ClustalW -PWGAPOPEN (default 15 for DNA, 10 for protein). |
| `pwGapExtend` | number | No | ClustalW -PWGAPEXT (default 6.66 for DNA, 0.1 for protein). |
| `pwMatrix` | string | No | type='protein': ClustalW -PWMATRIX, 'gonnet' (default), 'blosum', 'pam' or 'id'. |
| `ktuple` | integer | No | distance='clustalw-quick' only: k-tuple (word) size, 1-4 for DNA (clustalw -KTUPLE, default 2), 1-2 for protein (default 1). |
| `window` | integer | No | distance='clustalw-quick' only: diagonals scanned around each top diagonal (clustalw -WINDOW, default 4 for DNA, 5 for protein). |
| `topDiagonals` | integer | No | distance='clustalw-quick' only: number of best diagonals used (clustalw -TOPDIAGS, default 4 for DNA, 5 for protein). |
| `pairGap` | integer | No | distance='clustalw-quick' only: gap penalty of the fragment chain (clustalw -PAIRGAP; DNA default 5, or ktuple + 4 when ktuple is given; protein default 3, or ktuple + 3, as the clustalw command line does). |
| `transitionWeight` | number | No | DNA only, ClustalW -TRANSWEIGHT (default 0.5). (default: 0.5) |
| `noWeights` | boolean | No | ClustalW -NOWEIGHTS: no tree-derived sequence weights. (default: false) |
| `negative` | boolean | No | ClustalW -NEGATIVE: keep the matrix signed. (default: false) |
| `matrix` | string | No | type='protein': ClustalW -MATRIX, 'gonnet' (default), 'blosum', 'pam' or 'id'. |
| `noPreferredGaps` | boolean | No | type='protein': ClustalW -NOPGAP. (default: false) |
| `noHydrophilicGaps` | boolean | No | type='protein': ClustalW -NOHGAP. (default: false) |
| `variableRegionGaps` | boolean | No | type='protein': ClustalW -NOVGAP (sic: switches the variable-region penalties ON). (default: false) |
| `gapDist` | integer | No | type='protein': ClustalW -GAPDIST (default 4; 0 = off). (default: 4) |
| `hydrophilicResidues` | string | No | type='protein': ClustalW -HGAPRESIDUES (default GPSNDQEKR). (default: GPSNDQEKR) |
| `clustering` | string | No | ClustalW -CLUSTERING of the guide tree: 'nj' (default) or 'upgma'. (default: nj) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `alignedSequences` | array<string> | The aligned rows of the MSA (input order). |
| `consensus` | string | Majority-vote consensus of the aligned rows (gaps vote; a nucleotide beats a gap on ties). |
| `totalScore` | integer | Natural affine sum-of-pairs score of the alignment under match / mismatch / scoreGapOpen / scoreGapExtend (ClustalW itself reports no such score). |

## Errors

| Code | Message |
|------|---------|
| 1001 | Both inputs must contain at least one sequence |
| 1010 | Invalid DNA (type=dna) or protein (type=protein) residue in a row; letters and '-' only |
| 1002 | type must be dna or protein; distance must be clustalw or clustalw-quick; clustering must be nj or upgma; matrix / pwMatrix must be gonnet, blosum, pam or id |

## Examples

### Example 1: Two sequences added to a three-row profile

**Tool Call:**
```json
{
  "tool": "clustalw_sequences_to_profile",
  "arguments": {
    "profile": [
      "ACGTACGTTTGACCAGTAGGCATCGATGGA",
      "ACGTACGATTGACCAGTAGCCATCGATGGA",
      "ACGTAGGATTGACCGTAGCC-ATCGATGGA"
    ],
    "sequences": [
      "TTGACCAGTAGGCATCGATCCCAAGT",
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
    "ACGTAGGATTGACCGTAGCC-ATCGATGGA-----",
    "--------TTGACCAGTAGGCATCGATCCCAAGT-",
    "-----GGGTTTACCAGTAGGCATCGATGGAAACCT"
  ],
  "consensus": "ACGTACGATTGACCAGTAGGCATCGATGGA-----",
  "totalScore": 37
}
```

`clustalw -PROFILE1=a -PROFILE2=b -SEQUENCES -TYPE=DNA` prints exactly these rows.

## See Also

- [multiple_align_clustalw](multiple_align_clustalw.md)
