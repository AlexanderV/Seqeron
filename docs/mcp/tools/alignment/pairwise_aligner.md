# pairwise_aligner

Bio.Align.PairwiseAligner of Biopython 1.88 on two sequences: the optimal score, the algorithm Biopython selects (Needleman-Wunsch, Gotoh, Smith-Waterman, FOGSAA) and every alignment of the optimal score (target = sequenceA, query = sequenceB; block coordinates, printed rows, the Bio.Align text). mode global / local / fogsaa (the Fast Optimal Global Sequence Alignment Algorithm: one optimal alignment found by branch and bound, with Biopython's warnings when match < mismatch or a gap score exceeds the mismatch score), matchScore / mismatchScore or a substitution matrix (name or text), the wildcard letter (scores 0 against every letter, scalar scoring only), the presets scoring = blastn / megablast / blastp, the twelve gap scores by name, strand '-' (sequenceA against the reverse complement of sequenceB; the query coordinates decrease). gapScoresJson: a JSON object of Biopython gap-score properties (open_gap_score, extend_gap_score, gap_score, left_gap_score, right_open_gap_score, end_gap_score, target_gap_score, query_extend_gap_score, open_internal_insertion_score ... every name of the property tree; insertion = a gap in the target (sequence A), deletion = a gap in the query); the default of every gap score is -1. Letters are compared case-insensitively. Gap-score functions of (position, length) (gapFunctionsJson: Biopython's Waterman-Smith-Beyer algorithm for gap_score / insertion_score / deletion_score, strands + and -; the functions are named families because a tool cannot receive code) and gap scores above zero (Needleman-Wunsch, Gotoh, Smith-Waterman and Waterman-Smith-Beyer with a positive score, local alignments included) are supported; FOGSAA ignores gap functions, as in Biopython. VERIFICATION: score, algorithm name, warnings and the alignments (as a set; FOGSAA the alignment itself) equal to Biopython 1.88 on 698 / 698 recorded configurations (scripts/alignment/generate_pairwise_aligner_oracle.py) and on 1 501 / 1 501 configurations with gap-score functions or gap scores above zero (scripts/alignment/generate_pairwise_wsb_oracle.py; including the RuntimeError the reference raises on one local alignment).

## Overview

| Property | Value |
|----------|-------|
| **Server** | Alignment |
| **Tool Name** | `pairwise_aligner` |
| **Method ID** | `PairwiseAligner.Align` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Bio.Align.PairwiseAligner of Biopython 1.88 on two sequences: the optimal score, the algorithm Biopython selects (Needleman-Wunsch, Gotoh, Smith-Waterman, FOGSAA) and every alignment of the optimal score (target = sequenceA, query = sequenceB; block coordinates, printed rows, the Bio.Align text). mode global / local / fogsaa (the Fast Optimal Global Sequence Alignment Algorithm: one optimal alignment found by branch and bound, with Biopython's warnings when match < mismatch or a gap score exceeds the mismatch score), matchScore / mismatchScore or a substitution matrix (name or text), the wildcard letter (scores 0 against every letter, scalar scoring only), the presets scoring = blastn / megablast / blastp, the twelve gap scores by name, strand '-' (sequenceA against the reverse complement of sequenceB; the query coordinates decrease). gapScoresJson: a JSON object of Biopython gap-score properties (open_gap_score, extend_gap_score, gap_score, left_gap_score, right_open_gap_score, end_gap_score, target_gap_score, query_extend_gap_score, open_internal_insertion_score ... every name of the property tree; insertion = a gap in the target (sequence A), deletion = a gap in the query); the default of every gap score is -1. Letters are compared case-insensitively. Gap-score functions of (position, length) (gapFunctionsJson: Biopython's Waterman-Smith-Beyer algorithm for gap_score / insertion_score / deletion_score, strands + and -; the functions are named families because a tool cannot receive code) and gap scores above zero (Needleman-Wunsch, Gotoh, Smith-Waterman and Waterman-Smith-Beyer with a positive score, local alignments included) are supported; FOGSAA ignores gap functions, as in Biopython. VERIFICATION: score, algorithm name, warnings and the alignments (as a set; FOGSAA the alignment itself) equal to Biopython 1.88 on 698 / 698 recorded configurations (scripts/alignment/generate_pairwise_aligner_oracle.py) and on 1 501 / 1 501 configurations with gap-score functions or gap scores above zero (scripts/alignment/generate_pairwise_wsb_oracle.py; including the RuntimeError the reference raises on one local alignment).

## Core Documentation Reference

- Source: [PairwiseAligner.Align.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/PairwiseAligner.Align.cs)
- Algorithm: [Pairwise_Output_Formats_And_EMBOSS_Programs.md](../../../algorithms/Alignment/Pairwise_Output_Formats_And_EMBOSS_Programs.md)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequenceA` | string | Yes | The target sequence (any letters). |
| `sequenceB` | string | Yes | The query sequence (any letters). |
| `mode` | string | No | global (default), local or fogsaa. |
| `matchScore` | number | No | Score of an identical pair (default 1; ignored with a matrix). |
| `mismatchScore` | number | No | Score of a different pair (default 0; ignored with a matrix). |
| `matrix` | string | No | Optional substitution matrix: any built-in matrix by EMBOSS, NCBI or Biopython name - BLOSUM62, PAM250, EBLOSUM30 ... EBLOSUM90, EPAM10 ... EPAM500, BLASTP (BLOSUM62 with the residues U = selenocysteine, O = pyrrolysine and J = leucine / isoleucine, as blastp scores them), GONNET1992, JONES, DAYHOFF, MDM78, BENNER6 / 22 / 74, JOHNSON, RISLER, GENETIC, FENG, LEVIN, MCLACHLAN, RAO, STR (protein; the fractional ones are scored exactly) or EDNAFULL, EDNAMAT, EDNASIMPLE, ENUC.4.4, BLASTN, MEGABLAST, HOXD70, TRANS (IUPAC DNA); or THE TEXT of a matrix instead of a name: several lines in the NCBI / EMBOSS / Biopython data-file format (a header line of symbols, then one row per symbol: the symbol and its scores; # lines are comments; scores may be decimals). When given, match/mismatch are ignored, residues are scored by the published table (any protein or IUPAC symbol of the matrix; case-insensitive) and gapOpen/gapExtend may be fractional (the EMBOSS needle/water models: gapOpen = -(G - E), gapExtend = -E for -gapopen G -gapextend E). The result then also carries the real-valued exactScore. |
| `wildcard` | string | No | A single wildcard letter that scores 0 against every letter (scalar scoring only). |
| `scoring` | string | No | Scoring preset: blastn, megablast or blastp (matrix and gap scores of Bio.Align); individual settings below override it. |
| `gapScoresJson` | string | No | Gap scores: a JSON object of Biopython property names to numbers (see the description). |
| `gapFunctionsJson` | string | No | Gap-score functions of (position, length) (Biopython's Waterman-Smith-Beyer algorithm): a JSON object gap_score / insertion_score / deletion_score -> {"kind": lin \| log \| sq \| pos \| step \| tab \| zero \| bonus, "p": [a, b, c]}: lin -(a + b*length), log -(a + b*ln(length + 1)), sq -(a + b*length^2), pos -(a + b*position + c*length), step -a up to c positions and -b beyond, tab -(((position*7 + length*3) mod b)*0.5 + a), zero 0, bonus a*length - b. |
| `strand` | string | No | + (default) or -: align sequenceA with the reverse complement of sequenceB. |
| `maxAlignments` | integer | No | Return at most this many alignments (default 100; the number of optimal alignments can be exponential). |
| `scoreOnly` | boolean | No | true: only the score (no alignments; memory O(n) for large inputs). (default: false) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `score` | number | The optimal score |
| `algorithm` | string | The algorithm Biopython selects |
| `warnings` | array | The FOGSAA warnings |
| `count` | integer | Number of optimal alignments |
| `alignments` | array | alignedSequence1 / alignedSequence2 (rows), coordinates (target, query), score, text |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence is required |
| 1002 | Unknown mode / strand / preset, gap-score property, or an unknown gap-function family (gapFunctionsJson) |
| 1003 | The reference raises RuntimeError on this input (a local alignment whose traceback meets a cell without a move; reproduced with the message of Biopython) |

## Examples

### Example 1: Example

**Input:** `{"sequenceA": "TTACGTTT", "sequenceB": "ACGT", "mode": "fogsaa", "matchScore": 2, "mismatchScore": -1, "gapScoresJson": "{\"gap_score\": -2}"}`

**Response:** score 0.0, one alignment with coordinates [[0, 2, 6, 8], [0, 0, 4, 4]] (--ACGT--)

## See Also

- [global_align](global_align.md)
- [local_align](local_align.md)
