# clustalw_profile_align

ClustalW 2.1 profile-profile alignment of two existing alignments (clustalw -PROFILE).

## Overview

| Property | Value |
|----------|-------|
| **Server** | MultipleAlignment |
| **Tool Name** | `clustalw_profile_align` |
| **Method ID** | `SequenceAligner.AlignClustalWProfiles` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

ClustalW 2.1's profile-profile alignment of two existing alignments (clustalw -PROFILE1=a -PROFILE2=b -PROFILE, Clustal::profileAlign), ported step by step and verified row-for-row against the installed clustalw 2.1 on 600 DNA and protein sets (profiles produced by clustalw itself, option variants included): the rows of both profiles are kept (gaps '-' allowed; a profile may be one sequence), the identities inside each profile weight its sequences by the neighbour-joining tree of the profile (or UPGMA, clustering='upgma'), an initial alignment with equal weights measures the identities between the profiles, and the two profiles are aligned by Myers-Miller profile-profile alignment with position-specific gap penalties (weights = identity to the other profile x tree weight). Reproduced defects of the program: the tree of a profile 2 that does not start at sequence 1 is built from the identities of the FIRST sequences of the alignment, and the identities inside the profiles are truncated to integers. Parameters as multiple_align_clustalw (type 'dna' or 'protein'; gapOpen / gapExtend, transitionWeight, noWeights, negative; protein: matrix, noPreferredGaps, noHydrophilicGaps, variableRegionGaps, gapDist, hydrophilicResidues). Rows are returned: profile 1, then profile 2, in input order.

## Core Documentation Reference

- Source: [Seqeron.Genomics.ClustalW/ClustalWAligner.Profile.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.ClustalW/ClustalWAligner.Profile.cs)
- Algorithm doc: [Multiple_Sequence_Alignment.md](../../../algorithms/Alignment/Multiple_Sequence_Alignment.md)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `profile1` | array<string> | Yes | Rows of the first alignment (at least one; letters and '-'). |
| `profile2` | array<string> | Yes | Rows of the second alignment (at least one; letters and '-'). |
| `type` | string | No | 'dna' (default) or 'protein' (ClustalW -TYPE). |
| `gapOpen` | number | No | ClustalW -GAPOPEN (default 15 for DNA, 10 for protein). |
| `gapExtend` | number | No | ClustalW -GAPEXT (default 6.66 for DNA, 0.2 for protein). |
| `transitionWeight` | number | No | DNA only, ClustalW -TRANSWEIGHT (default 0.5). |
| `noWeights` | boolean | No | ClustalW -NOWEIGHTS: no tree-derived sequence weights. (default: false) |
| `negative` | boolean | No | ClustalW -NEGATIVE: keep the matrix signed (protein: BLOSUM40 / PAM120 / GONNET250 instead of the series). (default: false) |
| `matrix` | string | No | type='protein': ClustalW -MATRIX, 'gonnet' (default), 'blosum', 'pam' or 'id'. |
| `noPreferredGaps` | boolean | No | type='protein': ClustalW -NOPGAP. (default: false) |
| `noHydrophilicGaps` | boolean | No | type='protein': ClustalW -NOHGAP. (default: false) |
| `variableRegionGaps` | boolean | No | type='protein': ClustalW -NOVGAP (sic: switches the variable-region penalties ON). (default: false) |
| `gapDist` | integer | No | type='protein': ClustalW -GAPDIST (default 4; 0 = off). |
| `hydrophilicResidues` | string | No | type='protein': ClustalW -HGAPRESIDUES (default GPSNDQEKR). |
| `clustering` | string | No | ClustalW -CLUSTERING of the profile trees: 'nj' (default) or 'upgma'. |

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

### Example 1: Two profiles: three related 30-mers (one with a gap) and two shorter sequences

**Tool Call:**
```json
{
  "tool": "clustalw_profile_align",
  "arguments": {
    "profile1": [
      "ACGTACGTTTGACCAGTAGGCATCGATGGA",
      "ACGTACGATTGACCAGTAGCCATCGATGGA",
      "ACGTAGGATTGACCGTAGCC-ATCGATGGA"
    ],
    "profile2": [
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
    "-----TTGACCAGTAGGCATCGATCCCAAGT----",
    "-----GGGTTTACCAGTAGGCATCGATGGAAACCT"
  ],
  "consensus": "ACGTACGATTGACCAGTAGCCATCGATGGA-----",
  "totalScore": -62
}
```

`clustalw -PROFILE1=a -PROFILE2=b -PROFILE -TYPE=DNA` prints exactly these rows.

## See Also

- [multiple_align_clustalw](multiple_align_clustalw.md)
