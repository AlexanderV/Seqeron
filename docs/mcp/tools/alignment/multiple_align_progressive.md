# multiple_align_progressive

Guide-tree progressive MSA (Feng & Doolittle 1987 / Clustal-style): the ClustalW pairwise distance of every pairwise alignment, d = 1 − (identical aligned residue pairs) / (length of the shorter sequence) (Thompson, Higgins & Gibson 1994; pairalign.c tmat = (100 − 100·count/MIN(len1, len2)) / 100), gives the distance matrix, a guide tree (UPGMA by default, guideTree='nj' = ClustalW's neighbour-joining tree) orders the merges, and profiles are aligned by profile–profile dynamic programming with sum-of-pairs column scores (once a gap, always a gap). Default linear gap model (gapOpen ignored); gapModel='affine' uses a Gotoh three-state profile DP (occupancy-scaled gapOpen/gapExtend, ClustalW/MUSCLE style) and reports the natural affine sum-of-pairs score. distance='clustalw' takes the guide-tree distances from ClustalW 2.1's own affine pairwise stage instead (percent identity of clustalw -ALIGN). Returns rows, majority consensus and sum-of-pairs score. PROTEIN / MATRIX: type='protein' (or matrix) uses the substitution matrix for the pairwise alignments, the profile column scores and the sum-of-pairs objective (fractional gaps allowed); rows equal an independent Python implementation on 238 sets. With type='protein', distance='clustalw' / 'clustalw-quick' use ClustalW 2.1's PROTEIN pairwise stages (GONNET250 / BLOSUM30 / PAM350 / identity matrix via pwMatrix, gap 10 / 0.1, k-tuple 1; exact ports verified on 4 204 + 3 774 pairs of clustalw -TYPE=PROTEIN -ALIGN / -QUICKTREE).

## Overview

| Property | Value |
|----------|-------|
| **Server** | Alignment |
| **Tool Name** | `multiple_align_progressive` |
| **Method ID** | `SequenceAligner.MultipleAlignProgressive` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Guide-tree progressive MSA (Feng & Doolittle 1987 / Clustal-style): the ClustalW pairwise distance of every pairwise alignment, d = 1 − (identical aligned residue pairs) / (length of the shorter sequence) (Thompson, Higgins & Gibson 1994; pairalign.c tmat = (100 − 100·count/MIN(len1, len2)) / 100), gives the distance matrix, a guide tree (UPGMA by default, guideTree='nj' = ClustalW's neighbour-joining tree) orders the merges, and profiles are aligned by profile–profile dynamic programming with sum-of-pairs column scores (once a gap, always a gap). Default linear gap model (gapOpen ignored); gapModel='affine' uses a Gotoh three-state profile DP (occupancy-scaled gapOpen/gapExtend, ClustalW/MUSCLE style) and reports the natural affine sum-of-pairs score. distance='clustalw' takes the guide-tree distances from ClustalW 2.1's own affine pairwise stage instead (percent identity of clustalw -ALIGN). Returns rows, majority consensus and sum-of-pairs score. PROTEIN / MATRIX: type='protein' (or matrix) uses the substitution matrix for the pairwise alignments, the profile column scores and the sum-of-pairs objective (fractional gaps allowed); rows equal an independent Python implementation on 238 sets. With type='protein', distance='clustalw' / 'clustalw-quick' use ClustalW 2.1's PROTEIN pairwise stages (GONNET250 / BLOSUM30 / PAM350 / identity matrix via pwMatrix, gap 10 / 0.1, k-tuple 1; exact ports verified on 4 204 + 3 774 pairs of clustalw -TYPE=PROTEIN -ALIGN / -QUICKTREE).

## Core Documentation Reference

- Source: [Seqeron.Genomics.Alignment/SequenceAligner.cs#L2275](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.cs#L2275)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequences` | array<string> | Yes | Sequences to align (at least one): DNA (A, C, G, T) by default, amino acids with type='protein'. |
| `match` | integer | No | Score for matching base. (default: 1) |
| `mismatch` | integer | No | Score for mismatching base. (default: -1) |
| `gapOpen` | number | No | Gap-open penalty. Used only with gapModel='affine' (a gap of length k then costs gapOpen + k * gapExtend, Gotoh); ignored by the default linear model. May be fractional only with a substitution matrix (type='protein' or matrix). (default: -2) |
| `gapExtend` | number | No | Gap penalty per gap position (linear model: every gap position; affine model: the per-position part of gapOpen + k * gapExtend). May be fractional (e.g. -0.5) only with a substitution matrix (type='protein' or matrix). (default: -1) |
| `distance` | string | No | Guide-tree distance: 'linear' (default) = 1 - identities / min(length) of the repository's global alignment under the given scores (linear NW, or Gotoh when gapModel='affine'); 'clustalw' = the same formula on ClustalW 2.1's own pairwise stage (affine local Smith-Waterman-Gotoh alignment with its IUB matrix x100 and gap penalties 15 / 6.66, the 'slow' full-DP mode of clustalw -ALIGN; type='protein': its PROTEIN stage with the GONNET250 / BLOSUM30 / PAM350 / identity matrix (pwMatrix) and gap penalties 10 / 0.1 (pwGapOpen / pwGapExtend); the match/mismatch/gap parameters are then used only for the profile alignment); 'clustalw-quick' = ClustalW 2.1's fast k-tuple mode (clustalw -QUICKTREE, FastPairwiseAlign.cpp: best chain of exactly matching k-mers on the top diagonals, d = (100 - 100*fragmentScore/min(length))/100; tuned by ktuple / window / topDiagonals / pairGap). |
| `gapModel` | string | No | 'linear' (default): every gap position costs gapExtend and gapOpen is ignored; totalScore is the column sum of pairs (gap-gap columns 0). 'affine': gap runs cost gapOpen + k * gapExtend; the pairwise / profile alignments use affine (Gotoh) dynamic programming and totalScore is the sum of pairs with natural affine gap costs (Altschul 1989: each pair of rows is projected out of the alignment, gap-gap columns skipped, every gap run of the projection pays gapOpen once). For two sequences the affine result is an optimal pairwise alignment; for three or more sequences every method is a heuristic. |
| `guideTree` | string | No | Guide-tree method: 'upgma' (default; Sokal & Michener 1958, the Feng-Doolittle / CLUSTAL-1 choice) or 'nj' (neighbour-joining, Saitou & Nei 1987 - ClustalW 2.1's default tree, T-Coffee/Clustal-style; the unrooted NJ tree is rooted and ordered exactly as clustalw does: Tree::reRoot, left subtree first). With distance='clustalw' and guideTree='nj' the tree equals the .dnd file of clustalw -ALIGN. UPGMA is built by the canonical phylogenetics builder; NJ uses ClustalW's own neighbour-joining arithmetic (exact ties as in clustalw). |
| `ktuple` | integer | No | distance='clustalw-quick' only: k-tuple (word) size, 1-4 for DNA (clustalw -KTUPLE, default 2), 1-2 for protein (default 1). |
| `window` | integer | No | distance='clustalw-quick' only: diagonals scanned around each top diagonal (clustalw -WINDOW, default 4 for DNA, 5 for protein). |
| `topDiagonals` | integer | No | distance='clustalw-quick' only: number of best diagonals used (clustalw -TOPDIAGS, default 4 for DNA, 5 for protein). |
| `pairGap` | integer | No | distance='clustalw-quick' only: gap penalty of the fragment chain (clustalw -PAIRGAP; DNA default 5, or ktuple + 4 when ktuple is given; protein default 3, or ktuple + 3, as the clustalw command line does). |
| `type` | string | No | Sequence type: 'dna' (default; A, C, G, T scored with match/mismatch) or 'protein' (amino-acid sequences scored with the substitution matrix, default BLOSUM62; match/mismatch are then ignored). The protein path shares the dynamic programs, merge, profile and objective code of the DNA path (no second implementation). |
| `matrix` | string | No | Optional substitution matrix: any built-in protein matrix by EMBOSS or NCBI name (BLOSUM62, PAM250, EBLOSUM45, EPAM120, ...) for type='protein', a nucleotide matrix (EDNAFULL, EDNAMAT, EDNASIMPLE, ENUC.4.4) for type='dna'. When given (or type='protein'), match/mismatch are ignored, residues are scored by the published table (case-insensitive, any symbol of the matrix), gapOpen/gapExtend may be fractional (EMBOSS needle numbers gapopen G / gapextend E = gapOpen -(G - E), gapExtend -E) and the result also carries the real-valued exactTotalScore (totalScore = rounded). Verified against independent Python implementations of every method, Biopython 1.88 pairwise optima and exhaustive enumeration of tiny sets. |
| `pwMatrix` | string | No | type='protein' with distance='clustalw': ClustalW's pairwise matrix (clustalw -PWMATRIX): 'gonnet' (default; GONNET250), 'blosum' (BLOSUM30), 'pam' (PAM350) or 'id' (identity). |
| `pwGapOpen` | number | No | type='protein' with distance='clustalw': ClustalW -PWGAPOPEN (default 10; length-dependent: 2 * matrixAvg * (open + ln min(length)) in matrix units). |
| `pwGapExtend` | number | No | type='protein' with distance='clustalw': ClustalW -PWGAPEXT (default 0.1). |
| `format` | string | No | Optional output format; the alignment written in that format is returned in the result field 'formatted'. Biopython Bio.Align names: fasta, a2m, clustal, phylip, nexus, stockholm (byte for byte as Bio.Align.write; 'emboss:' + name selects the EMBOSS writer of the same name). EMBOSS -aformat names: pair, simple, srspair, srs, markx0, markx1, markx2, markx3, markx10, match, score, sam, tcoffee (two rows) and msf, mega, meganon, nexusnon, paup, paupnon, phylipnon, selex, treecon (any number of rows; EMBOSS prints the run date: now, UTC). Omit for no formatted text. |
| `names` | array of string | No | Optional names of the rows for the formatted output (default seq1, seq2, ...). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `alignedSequences` | array<string> | The aligned rows of the MSA. |
| `consensus` | string | Majority-vote consensus of the aligned rows (gaps vote; a nucleotide beats a gap on ties). |
| `totalScore` | integer | Sum-of-pairs score of the MSA: linear model = column score (match/mismatch from the matrix, residue-gap = gapExtend, gap-gap = 0); gapModel='affine' = natural affine gap costs (every gap run of each pair projection pays gapOpen once). |
| `exactTotalScore` | number | Real-valued sum-of-pairs score (= totalScore for DNA with integer scores; exact, e.g. 49.5, for substitution-matrix runs with fractional gaps, where totalScore is rounded half away from zero). |
| `formatted` | string | No | The alignment written in the requested format (null when format is omitted). |

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

- Algorithm source: [Seqeron.Genomics.Alignment/SequenceAligner.cs#L2275](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.cs#L2275)
- Binding: [AlignmentTools.cs](../../../../src/Seqeron/Mcp/Seqeron.Mcp.Alignment/Tools/AlignmentTools.cs)
