# local_align_top_n

Waterman–Eggert (1987, J Mol Biol 197:723) declumping: the best local alignment, then the best local alignment that shares no ALIGNED RESIDUE PAIR with it, then the best sharing no pair with either, ... up to maxAlignments alignments of score >= minScore (and always > 0) - the definition of EMBOSS matcher -alternatives N (Huang-Miller SIM) and FASTA lalign. An aligned pair is a diagonal column (a_i opposite b_j); gap columns are not pairs, so two reported alignments may cross through a gap column but never align the same two residues. Alignment 1 is exactly local_align (affine = false) / local_align_affine (affine = true); scores are non-increasing; ties between co-optimal alignments are broken, as local_align does (first maximal cell in row-major order, traceback diagonal > up > left, M > X > Y), so on tie-rich inputs a later alignment can differ from matcher's (its score is always the true maximum of the remaining problem); the order of EMBOSS matcher (Huang-Miller SIM) is offered by the matcher_align tool of the EMBOSS pairwise server. affine = false: linear gaps (every gap position costs gapExtend, gapOpen ignored); affine = true: a gap of length k costs gapOpen + k*gapExtend (EMBOSS -gapopen G -gapextend E = gapOpen -(G - E), gapExtend -E). With matrix = BLOSUM62 | PAM250 | EDNAFULL residues are scored by that table, gaps may be fractional and minScore/exactScore are real scores. Returns alignments (same fields as local_align) in the order found; fewer than maxAlignments when nothing positive >= minScore is left (empty list when no pair scores > 0). Memory: up to 25,000,000 cells ((len1+1)*(len2+1)) the original recomputation scheme (matrices kept, only cells whose value can change are recomputed after each report; O(len1*len2) memory); ABOVE that (e.g. two 5-20 kb sequences) the same alignments - identical rows, scores, coordinates, order and ties - are computed in exact LINEAR SPACE (O(len1 + len2) memory per alignment: one run of the Hirschberg / Myers-Miller / Huang-Miller locate-then-solve local program per reported alignment with the aligned pairs of the earlier ones forbidden, the SIM / EMBOSS matcher scheme), so the running time grows with maxAlignments (about one extra local alignment per reported alignment) while the memory does not (20 kb x 20 kb: tens of MB instead of 1.6 GB); the linear-space program needs non-positive gap scores (a positive gap score on such an input is rejected with a clear error).

## Overview

| Property | Value |
|----------|-------|
| **Server** | Alignment |
| **Tool Name** | `local_align_top_n` |
| **Method ID** | `SequenceAligner.LocalAlignTopN` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Waterman–Eggert (1987, J Mol Biol 197:723) declumping: the best local alignment, then the best local alignment that shares no ALIGNED RESIDUE PAIR with it, then the best sharing no pair with either, ... up to maxAlignments alignments of score >= minScore (and always > 0) - the definition of EMBOSS matcher -alternatives N (Huang-Miller SIM) and FASTA lalign. An aligned pair is a diagonal column (a_i opposite b_j); gap columns are not pairs, so two reported alignments may cross through a gap column but never align the same two residues. Alignment 1 is exactly local_align (affine = false) / local_align_affine (affine = true); scores are non-increasing; ties between co-optimal alignments are broken, as local_align does (first maximal cell in row-major order, traceback diagonal > up > left, M > X > Y), so on tie-rich inputs a later alignment can differ from matcher's (its score is always the true maximum of the remaining problem); the order of EMBOSS matcher (Huang-Miller SIM) is offered by the matcher_align tool of the EMBOSS pairwise server. affine = false: linear gaps (every gap position costs gapExtend, gapOpen ignored); affine = true: a gap of length k costs gapOpen + k*gapExtend (EMBOSS -gapopen G -gapextend E = gapOpen -(G - E), gapExtend -E). With matrix = BLOSUM62 | PAM250 | EDNAFULL residues are scored by that table, gaps may be fractional and minScore/exactScore are real scores. Returns alignments (same fields as local_align) in the order found; fewer than maxAlignments when nothing positive >= minScore is left (empty list when no pair scores > 0). Memory: up to 25,000,000 cells ((len1+1)*(len2+1)) the original recomputation scheme (matrices kept, only cells whose value can change are recomputed after each report; O(len1*len2) memory); ABOVE that (e.g. two 5-20 kb sequences) the same alignments - identical rows, scores, coordinates, order and ties - are computed in exact LINEAR SPACE (O(len1 + len2) memory per alignment: one run of the Hirschberg / Myers-Miller / Huang-Miller locate-then-solve local program per reported alignment with the aligned pairs of the earlier ones forbidden, the SIM / EMBOSS matcher scheme), so the running time grows with maxAlignments (about one extra local alignment per reported alignment) while the memory does not (20 kb x 20 kb: tens of MB instead of 1.6 GB); the linear-space program needs non-positive gap scores (a positive gap score on such an input is rejected with a clear error).

## Core Documentation Reference

- Source: [Seqeron.Genomics.Alignment/SequenceAligner.TopN.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.TopN.cs)
- Algorithm doc: [Local_Alignment_Smith_Waterman.md](../../../algorithms/Alignment/Local_Alignment_Smith_Waterman.md) (section 5.5)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence1` | string | Yes | First sequence. |
| `sequence2` | string | Yes | Second sequence. |
| `maxAlignments` | integer | No | Maximum number of alignments to report (>= 1). (default: 10) |
| `minScore` | number | No | Smallest accepted score (inclusive); alignments scoring <= 0 are never reported. Real number (fractional with matrix). (default: 0) |
| `affine` | boolean | No | false = linear gap model (gapExtend per gap position); true = affine (gapOpen + k*gapExtend). (default: false) |
| `match` | integer | No | Score for matching base (ignored with matrix). (default: 1) |
| `mismatch` | integer | No | Score for mismatching base (ignored with matrix). (default: -1) |
| `gapOpen` | number | No | Gap-open score (affine = true only; ignored by the linear model). May be fractional when matrix is given. (default: -2) |
| `gapExtend` | number | No | Gap score charged per gap position (<= 0). May be fractional when matrix is given. (default: -1) |
| `matrix` | string | No | Optional substitution matrix: any built-in matrix by EMBOSS, NCBI or Biopython name - BLOSUM62, PAM250, EBLOSUM30 ... EBLOSUM90, EPAM10 ... EPAM500, BLASTP (BLOSUM62 with the residues U = selenocysteine, O = pyrrolysine and J = leucine / isoleucine, as blastp scores them), GONNET1992, JONES, DAYHOFF, MDM78, BENNER6 / 22 / 74, JOHNSON, RISLER, GENETIC, FENG, LEVIN, MCLACHLAN, RAO, STR (protein; the fractional ones are scored exactly) or EDNAFULL, EDNAMAT, EDNASIMPLE, ENUC.4.4, BLASTN, MEGABLAST, HOXD70, TRANS (IUPAC DNA); or THE TEXT of a matrix instead of a name: several lines in the NCBI / EMBOSS / Biopython data-file format (a header line of symbols, then one row per symbol: the symbol and its scores; # lines are comments; scores may be decimals). When given, match/mismatch are ignored, residues are scored by the published table (any protein or IUPAC symbol of the matrix; case-insensitive) and gapOpen/gapExtend may be fractional (the EMBOSS needle/water models: gapOpen = -(G - E), gapExtend = -E for -gapopen G -gapextend E). The result then also carries the real-valued exactScore. |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `alignments` | array | Local alignments, best first; each has alignedSequence1/2, score, alignmentType, startPosition1/2, endPosition1/2, exactScore (same fields as local_align). |
| `count` | integer | Number of alignments returned. |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence is required (null). An empty string is accepted: the library aligns it to all gaps (the local tools return no alignment: score 0 / empty list) |
| 1002 | maxAlignments must be at least 1 |
| 1003 | A gap score is positive (GapOpen and GapExtend must be <= 0 in every scoring model) or the integer score range of the dynamic program is exceeded (max|step| x (len1 + len2 + 2) must stay below 2^28) |

## Examples

### Example 1

**Tool Call:**
```json
{
  "tool": "local_align_top_n",
  "arguments": {
    "sequence1": "ACGTACTTTACGTACGGGACGTAC",
    "sequence2": "ACGTAC",
    "maxAlignments": 3,
    "match": 2,
    "mismatch": -3,
    "gapExtend": -2
  }
}
```

**Response:**
```json
{
  "alignments": [
    {
      "alignedSequence1": "ACGTAC",
      "alignedSequence2": "ACGTAC",
      "score": 12,
      "alignmentType": "Local",
      "startPosition1": 0,
      "startPosition2": 0,
      "endPosition1": 5,
      "endPosition2": 5,
      "exactScore": 12.0
    },
    {
      "alignedSequence1": "ACGTAC",
      "alignedSequence2": "ACGTAC",
      "score": 12,
      "alignmentType": "Local",
      "startPosition1": 9,
      "startPosition2": 0,
      "endPosition1": 14,
      "endPosition2": 5,
      "exactScore": 12.0
    },
    {
      "alignedSequence1": "ACGTAC",
      "alignedSequence2": "ACGTAC",
      "score": 12,
      "alignmentType": "Local",
      "startPosition1": 18,
      "startPosition2": 0,
      "endPosition1": 23,
      "endPosition2": 5,
      "exactScore": 12.0
    }
  ],
  "count": 3
}
```

## Worked Example

Sequence 1 holds three copies of ACGTAC (indices 0, 9 and 18), sequence 2 one. Smith-Waterman (`local_align`) reports only the first copy (12); Waterman-Eggert declumping reports all three (12, 12, 12) because the copies share no aligned residue pair (every pair of a reported alignment is forbidden for the next one). Oracle: independent full-recompute Waterman-Eggert reference + exhaustive enumeration; EMBOSS matcher -alternatives and Biopython give the same scores.

## References

- Algorithm source: [Seqeron.Genomics.Alignment/SequenceAligner.TopN.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.TopN.cs)
- Binding: [AlignmentTools.cs](../../../../src/Seqeron/Mcp/Seqeron.Mcp.Alignment/Tools/AlignmentTools.cs)
