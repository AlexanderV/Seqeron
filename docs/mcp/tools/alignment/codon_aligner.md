# codon_aligner

Bio.Align.CodonAligner of Biopython 1.88: a protein aligned to its nucleotide sequence over the three translated frames with frame-shift scores (-2, -1, +1, +2); score and every alignment as vertices (refuses the inputs on which Biopython's iterator reads uninitialised memory).

## Overview

| Property | Value |
|----------|-------|
| **Server** | Alignment |
| **Tool Name** | `codon_aligner` |
| **Method ID** | `CodonAligner.Align` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Bio.Align.CodonAligner of Biopython 1.88: aligns a nucleotide sequence to its protein by dynamic programming over the three translated frames of the nucleotide sequence (genetic code by NCBI table number; DNA or RNA, a codon with a letter outside IUPAC is an error, ambiguous codons translate to X). Every residue is aligned to one codon; between two codons the nucleotide sequence may advance by 1, 2, 4 or 5 instead of 3 (a frame shift) at frameshiftMinusTwoScore / MinusOne / PlusOne / PlusTwo (default -3 each, or all four with frameshiftScore); leading and trailing nucleotides are free; residues and translated codons score matchScore / mismatchScore (defaults 1 / 0), 0 against the wildcard letter (default X). Returns the optimal score and every alignment of it as the vertices (target = residue positions, query = nucleotide positions; a frame shift is two vertices at the same residue position, the query going back for -2 / -1 and skipping for +1 / +2). Biopython's path iterator reads memory its traceback never wrote when the optimal alignment skips leading residues (it then lists one-vertex alignments for ever, with a heap overflow): those inputs are refused with an error, as the score() is defined for all inputs (scoreOnly = true). VERIFICATION: score, aligned score and every alignment in the order of the iterator equal to Biopython 1.88 on 1 400 of 1 500 recorded inputs (the other 100 are the refused ones; scripts/alignment/generate_codon_aligner_oracle.py).

## Core Documentation Reference

- Source: [CodonAligner.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/CodonAligner.cs)
- Algorithm: [Pairwise_Output_Formats_And_EMBOSS_Programs.md](../../../algorithms/Alignment/Pairwise_Output_Formats_And_EMBOSS_Programs.md)
- Oracle: [generate_codon_aligner_oracle.py](../../../../scripts/alignment/generate_codon_aligner_oracle.py)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `protein` | string | Yes | The protein sequence (amino-acid letters). |
| `nucleotides` | string | Yes | The nucleotide sequence (DNA or RNA). |
| `tableNumber` | integer | No | NCBI genetic code table number (default 1, the standard code). |
| `matchScore` | number | No | Score of a residue equal to the translated codon (default 1). |
| `mismatchScore` | number | No | Score of a residue different from the translated codon (default 0). |
| `wildcard` | string | No | A single wildcard letter that scores 0 against every letter (default X). |
| `frameshiftScore` | number | No | Score of all four frame shifts (default -3); the four below override it. |
| `frameshiftMinusTwoScore` | number | No | Score of the shift -2 (the next codon overlaps the previous one by two nucleotides). |
| `frameshiftMinusOneScore` | number | No | Score of the shift -1. |
| `frameshiftPlusOneScore` | number | No | Score of the shift +1 (one nucleotide is skipped between two codons). |
| `frameshiftPlusTwoScore` | number | No | Score of the shift +2. |
| `maxAlignments` | integer | No | Return at most this many alignments (default 100). |
| `scoreOnly` | boolean | No | true: only the score (defined for every input). (default: false) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `score` | number | The optimal score |
| `count` | integer | Number of alignments returned |
| `alignments` | array | coordinates: the vertices (target = residues, query = nucleotides) of every alignment |

## Errors

| Code | Message |
|------|---------|
| 1001 | A codon has a letter outside IUPAC (Biopython: TranslationError) |
| 1002 | Unknown genetic code table, or a wildcard that is not a single character |
| 1003 | The reference iterator reads memory its traceback never wrote on this input (InvalidOperationException; scoreOnly = true still gives the score) |

## Examples

### Example 1: A frame shift between codon 7 and 8

**User Prompt:**
> Align this protein to its nucleotide sequence allowing frame shifts.

**Expected Tool Call:**
```json
{"tool": "codon_aligner", "arguments": {"protein": "MLGSNESRVCPKLSPS", "nucleotides": "ATGCTGGGCTCGAACGAGTCCGTGTATGCCCTAAGCTGAGCCCGTCG"}}
```

**Response:**
```json
{"score": 13.0, "count": 1, "alignments": [{"coordinates": [[0, 7, 7, 16], [0, 21, 20, 47]]}]}
```

The query goes back by one nucleotide at the frame shift (21 -> 20): the codon of residue 8 overlaps the previous one by one nucleotide (Bio.Align doctest of CodonAligner.align).

## See Also

- [pairwise_aligner](pairwise_aligner.md)
- [alignment_operation](alignment_operation.md)
