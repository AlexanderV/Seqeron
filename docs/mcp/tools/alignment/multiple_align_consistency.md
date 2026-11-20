# multiple_align_consistency

Consistency-based MSA (T-Coffee; Notredame, Higgins & Heringa 2000): a primary library of global (Needleman–Wunsch) and local (Smith–Waterman) pairwise alignments weighted by percent identity, extended by the triplet consistency transformation (direct weight + Σ min(W1, W2) over third sequences), then progressively aligned along a guide tree (UPGMA by default; guideTree='nj' = neighbour-joining as in the T-Coffee paper) with the extended-library weights as position-specific scores and zero gap penalties (T-Coffee default). The scoring parameters seed the pairwise library (gapModel='affine': Gotoh global and Smith–Waterman–Gotoh local alignments with gapOpen + k·gapExtend) and define totalScore, the sum of pairs of the result (linear column score, or natural affine costs with gapModel='affine'). distance='clustalw' takes the guide-tree distances from ClustalW 2.1's own pairwise stage. consistencyEngine='tcoffee' swaps the progressive stage for the T-Coffee 13.41 program's (library relaxation + linked pair-wise sparse DP; exact for the same library and tree). PROTEIN / MATRIX: type='protein' (or matrix) builds the primary library from the matrix global + local pair alignments (weights = percent identity) and otherwise runs the same pipeline (consistencyEngine='tcoffee' too); rows equal an independent Python implementation, and the T-Coffee progressive stage equals the t_coffee program on 160 protein sets. With type='protein', distance='clustalw' / 'clustalw-quick' use ClustalW 2.1's PROTEIN pairwise stages (GONNET250 / BLOSUM30 / PAM350 / identity matrix via pwMatrix, gap 10 / 0.1, k-tuple 1; exact ports verified on 4 204 + 3 774 pairs of clustalw -TYPE=PROTEIN -ALIGN / -QUICKTREE).

## Overview

| Property | Value |
|----------|-------|
| **Server** | Alignment |
| **Tool Name** | `multiple_align_consistency` |
| **Method ID** | `SequenceAligner.MultipleAlignConsistency` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Consistency-based MSA (T-Coffee; Notredame, Higgins & Heringa 2000): a primary library of global (Needleman–Wunsch) and local (Smith–Waterman) pairwise alignments weighted by percent identity, extended by the triplet consistency transformation (direct weight + Σ min(W1, W2) over third sequences), then progressively aligned along a guide tree (UPGMA by default; guideTree='nj' = neighbour-joining as in the T-Coffee paper) with the extended-library weights as position-specific scores and zero gap penalties (T-Coffee default). The scoring parameters seed the pairwise library (gapModel='affine': Gotoh global and Smith–Waterman–Gotoh local alignments with gapOpen + k·gapExtend) and define totalScore, the sum of pairs of the result (linear column score, or natural affine costs with gapModel='affine'). distance='clustalw' takes the guide-tree distances from ClustalW 2.1's own pairwise stage. consistencyEngine='tcoffee' swaps the progressive stage for the T-Coffee 13.41 program's (library relaxation + linked pair-wise sparse DP; exact for the same library and tree). PROTEIN / MATRIX: type='protein' (or matrix) builds the primary library from the matrix global + local pair alignments (weights = percent identity) and otherwise runs the same pipeline (consistencyEngine='tcoffee' too); rows equal an independent Python implementation, and the T-Coffee progressive stage equals the t_coffee program on 160 protein sets. With type='protein', distance='clustalw' / 'clustalw-quick' use ClustalW 2.1's PROTEIN pairwise stages (GONNET250 / BLOSUM30 / PAM350 / identity matrix via pwMatrix, gap 10 / 0.1, k-tuple 1; exact ports verified on 4 204 + 3 774 pairs of clustalw -TYPE=PROTEIN -ALIGN / -QUICKTREE).

## Core Documentation Reference

- Source: [Seqeron.Genomics.Alignment/SequenceAligner.cs#L3011](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.cs#L3011)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequences` | array<string> | Yes | Sequences to align (at least one): DNA (A, C, G, T) by default, amino acids with type='protein'. |
| `match` | integer | No | Score for matching base (pairwise library and reported score). (default: 1) |
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
| `consistencyEngine` | string | No | Progressive stage: 'paper' (default) = the published formula (extended-library weight summed over every cross pair of two columns, full DP, zero gap penalty); 'tcoffee' = the T-Coffee 13.41 program's own progressive stage (SequenceAligner.MultipleAlignTCoffee: library relaxation, 'linked pair-wise' sparse DP over scored column pairs, weights 10 x percent identity). With the same primary library and guide tree 'tcoffee' reproduces t_coffee's rows exactly (800 explicit-library sets); the library here is this repository's global + local alignments, not the program's slow_pair / lalign_id_pair / proba_pair libraries, and the tree is guideTree (not the program's k-tuple NJ tree), so a plain 't_coffee seqs.fa' run is not reproduced end to end. |
| `type` | string | No | Sequence type: 'dna' (default; A, C, G, T scored with match/mismatch) or 'protein' (amino-acid sequences scored with the substitution matrix, default BLOSUM62; match/mismatch are then ignored). The protein path shares the dynamic programs, merge, profile and objective code of the DNA path (no second implementation). |
| `matrix` | string | No | Optional substitution matrix: any built-in protein matrix by EMBOSS or NCBI name (BLOSUM62, PAM250, EBLOSUM45, EPAM120, ...) for type='protein', a nucleotide matrix (EDNAFULL, EDNAMAT, EDNASIMPLE, ENUC.4.4) for type='dna'. When given (or type='protein'), match/mismatch are ignored, residues are scored by the published table (case-insensitive, any symbol of the matrix), gapOpen/gapExtend may be fractional (EMBOSS needle numbers gapopen G / gapextend E = gapOpen -(G - E), gapExtend -E) and the result also carries the real-valued exactTotalScore (totalScore = rounded). Verified against independent Python implementations of every method, Biopython 1.88 pairwise optima and exhaustive enumeration of tiny sets. |
| `pwMatrix` | string | No | type='protein' with distance='clustalw': ClustalW's pairwise matrix (clustalw -PWMATRIX): 'gonnet' (default; GONNET250), 'blosum' (BLOSUM30), 'pam' (PAM350) or 'id' (identity). |
| `pwGapOpen` | number | No | type='protein' with distance='clustalw': ClustalW -PWGAPOPEN (default 10; length-dependent: 2 * matrixAvg * (open + ln min(length)) in matrix units). |
| `pwGapExtend` | number | No | type='protein' with distance='clustalw': ClustalW -PWGAPEXT (default 0.1). |
| `primaryLocalAlignments` | integer | No | T-Coffee Lalign library: the number of non-intersecting (Waterman-Eggert) local alignments per sequence pair that enter the primary library, each weighted by its own percent identity (lalign_id_pair, -lalign_n_top; default 1 = the single best Smith-Waterman alignment of the published formula). Paper engine and TCoffeeLinked engine. (default: 1) |
| `sequenceWeights` | string | No | Sequence redundancy weights (T-Coffee -seq_weight): 'none' (default) or 'tcoffee' (compute_t_coffee_weight: w = 2 / (1 + sum (10 x identity)^3) rescaled to mean 1; an extended library weight is multiplied by the smaller of the two sequence weights). |
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

- Algorithm source: [Seqeron.Genomics.Alignment/SequenceAligner.cs#L3011](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.cs#L3011)
- Binding: [AlignmentTools.cs](../../../../src/Seqeron/Mcp/Seqeron.Mcp.Alignment/Tools/AlignmentTools.cs)
