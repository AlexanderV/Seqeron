# multiple_align_anchored

Approximate, faster center-star MSA: picks the center by 4-mer count cosine similarity, builds a suffix tree on it and aligns every other sequence with exact-match anchors (MUMmer-style seed, weighted colinear chain, Needleman–Wunsch — or Gotoh affine with gapModel='affine' — on the gaps between anchors). The pairwise alignments are NOT guaranteed optimal, so the sum-of-pairs score can be lower than multiple_align's. Default linear gap model (gapOpen ignored); gapModel='affine' charges gapOpen + k·gapExtend per gap run and reports the natural affine sum-of-pairs score. PROTEIN / MATRIX: type='protein' (or matrix) scores anchors as the sum of matrix(c, c) and fills the gaps between anchors with the matrix aligner (any alphabet); rows equal an independent Python implementation (anchors = suffix-tree matching statistics, colinear chain, Needleman-Wunsch / Gotoh between anchors).

## Overview

| Property | Value |
|----------|-------|
| **Server** | Alignment |
| **Tool Name** | `multiple_align_anchored` |
| **Method ID** | `SequenceAligner.MultipleAlignAnchored` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Approximate, faster center-star MSA: picks the center by 4-mer count cosine similarity, builds a suffix tree on it and aligns every other sequence with exact-match anchors (MUMmer-style seed, weighted colinear chain, Needleman–Wunsch — or Gotoh affine with gapModel='affine' — on the gaps between anchors). The pairwise alignments are NOT guaranteed optimal, so the sum-of-pairs score can be lower than multiple_align's. Default linear gap model (gapOpen ignored); gapModel='affine' charges gapOpen + k·gapExtend per gap run and reports the natural affine sum-of-pairs score. PROTEIN / MATRIX: type='protein' (or matrix) scores anchors as the sum of matrix(c, c) and fills the gaps between anchors with the matrix aligner (any alphabet); rows equal an independent Python implementation (anchors = suffix-tree matching statistics, colinear chain, Needleman-Wunsch / Gotoh between anchors).

## Core Documentation Reference

- Source: [Seqeron.Genomics.Alignment/SequenceAligner.cs#L1821](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.cs#L1821)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequences` | array<string> | Yes | Sequences to align (at least one): DNA (A, C, G, T) by default, amino acids with type='protein'. |
| `match` | integer | No | Score for matching base. (default: 1) |
| `mismatch` | integer | No | Score for mismatching base. (default: -1) |
| `gapOpen` | number | No | Gap-open penalty. Used only with gapModel='affine' (a gap of length k then costs gapOpen + k * gapExtend, Gotoh); ignored by the default linear model. May be fractional only with a substitution matrix (type='protein' or matrix). (default: -2) |
| `gapExtend` | number | No | Gap penalty per gap position (linear model: every gap position; affine model: the per-position part of gapOpen + k * gapExtend). May be fractional (e.g. -0.5) only with a substitution matrix (type='protein' or matrix). (default: -1) |
| `gapModel` | string | No | 'linear' (default): every gap position costs gapExtend and gapOpen is ignored; totalScore is the column sum of pairs (gap-gap columns 0). 'affine': gap runs cost gapOpen + k * gapExtend; the pairwise / profile alignments use affine (Gotoh) dynamic programming and totalScore is the sum of pairs with natural affine gap costs (Altschul 1989: each pair of rows is projected out of the alignment, gap-gap columns skipped, every gap run of the projection pays gapOpen once). For two sequences the affine result is an optimal pairwise alignment; for three or more sequences every method is a heuristic. (default: linear) |
| `type` | string | No | Sequence type: 'dna' (default; A, C, G, T scored with match/mismatch) or 'protein' (amino-acid sequences scored with the substitution matrix, default BLOSUM62; match/mismatch are then ignored). The protein path shares the dynamic programs, merge, profile and objective code of the DNA path (no second implementation). (default: dna) |
| `matrix` | string | No | Optional substitution matrix: any built-in protein matrix by EMBOSS or NCBI name (BLOSUM62, PAM250, EBLOSUM45, EPAM120, ...) for type='protein', a nucleotide matrix (EDNAFULL, EDNAMAT, EDNASIMPLE, ENUC.4.4) for type='dna'. When given (or type='protein'), match/mismatch are ignored, residues are scored by the published table (case-insensitive, any symbol of the matrix), gapOpen/gapExtend may be fractional (EMBOSS needle numbers gapopen G / gapextend E = gapOpen -(G - E), gapExtend -E) and the result also carries the real-valued exactTotalScore (totalScore = rounded). Verified against independent Python implementations of every method, Biopython 1.88 pairwise optima and exhaustive enumeration of tiny sets. |
| `format` | string | No | Optional output format; the alignment written in that format is returned in 'formatted'. Biopython Bio.Align names: fasta, a2m, clustal, phylip, nexus, stockholm (byte for byte as Bio.Align.write; 'emboss:' + name selects the EMBOSS writer of the same name) or EMBOSS -aformat names (pair, simple, srspair, srs, markx0-3, markx10, match, score, sam, tcoffee, msf, mega, meganon, nexusnon, paup, paupnon, phylipnon, selex, treecon). Omit for no formatted text. |
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

### Example 1: Three identical

**Tool Call:**
```json
{
  "tool": "multiple_align_anchored",
  "arguments": {
    "sequences": [
      "ACGTACGT",
      "ACGTACGT",
      "ACGTACGT"
    ]
  }
}
```

**Response:**
```json
{
  "alignedSequences": [
    "ACGTACGT",
    "ACGTACGT",
    "ACGTACGT"
  ],
  "consensus": "ACGTACGT",
  "totalScore": 24
}
```

## Worked Example

Three identical 8-mers align gaplessly; SP = 3 pairs x 8 columns x (+1) = 24.

## References

- Algorithm source: [Seqeron.Genomics.Alignment/SequenceAligner.cs#L1821](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.cs#L1821)
- Binding: [AlignmentTools.cs](../../../../src/Seqeron/Mcp/Seqeron.Mcp.Alignment/Tools/AlignmentTools.cs)
