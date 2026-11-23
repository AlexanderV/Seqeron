# alignment_counts

The counts of an alignment of any number of sequences as Bio.Align's alignment.counts(scoring): aligned letters, identities, mismatches, positives (with a substitution matrix: letters with a positive score), and the gaps split by side (left / internal / right), direction (insertion = a gap in the target, deletion = a gap in the query) and opened / extended; the counts add up over every pair of sequences. Scoring: none (letter equality); wildcard (letters equal to it are neither identities nor mismatches); matrix (name or text; adds the substitution score); or an aligner (matchScore / mismatchScore / matrix / wildcard / gapScoresJson as in pairwise_aligner) which adds the gap score and the total score. gapScoresJson: a JSON object of Biopython gap-score properties (open_gap_score, extend_gap_score, gap_score, left_gap_score, right_open_gap_score, end_gap_score, target_gap_score, query_extend_gap_score, open_internal_insertion_score ... every name of the property tree; insertion = a gap in the target (sequence A), deletion = a gap in the query); the default of every gap score is -1. VERIFICATION: equal to Biopython 1.88 on 260 alignments x 8 scorings (every count and the printed text).

## Overview

| Property | Value |
|----------|-------|
| **Server** | Alignment |
| **Tool Name** | `alignment_counts` |
| **Method ID** | `AlignmentOperations.Counts` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

The counts of an alignment of any number of sequences as Bio.Align's alignment.counts(scoring): aligned letters, identities, mismatches, positives (with a substitution matrix: letters with a positive score), and the gaps split by side (left / internal / right), direction (insertion = a gap in the target, deletion = a gap in the query) and opened / extended; the counts add up over every pair of sequences. Scoring: none (letter equality); wildcard (letters equal to it are neither identities nor mismatches); matrix (name or text; adds the substitution score); or an aligner (matchScore / mismatchScore / matrix / wildcard / gapScoresJson as in pairwise_aligner) which adds the gap score and the total score. gapScoresJson: a JSON object of Biopython gap-score properties (open_gap_score, extend_gap_score, gap_score, left_gap_score, right_open_gap_score, end_gap_score, target_gap_score, query_extend_gap_score, open_internal_insertion_score ... every name of the property tree; insertion = a gap in the target (sequence A), deletion = a gap in the query); the default of every gap score is -1. VERIFICATION: equal to Biopython 1.88 on 260 alignments x 8 scorings (every count and the printed text).

## Core Documentation Reference

- Source: [AlignmentOperations.Counts.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/AlignmentOperations.Counts.cs)
- Algorithm: [Pairwise_Output_Formats_And_EMBOSS_Programs.md](../../../algorithms/Alignment/Pairwise_Output_Formats_And_EMBOSS_Programs.md)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `alignmentJson` | string | Yes | An alignment as JSON: {"sequences": ["ACGT...", ...] (or objects {"id": ..., "length": ..., "parts": {"0": "ACGT..."}}), "coordinates": [[target positions], [query positions], ...]} (one coordinate row per sequence; a position that decreases along the alignment is a sequence on the reverse strand), as parse_coordinate_alignment returns it. |
| `wildcard` | string | No | A single wildcard letter (neither identity nor mismatch). |
| `matrix` | string | No | Optional substitution matrix: any built-in matrix by EMBOSS, NCBI or Biopython name - BLOSUM62, PAM250, EBLOSUM30 ... EBLOSUM90, EPAM10 ... EPAM500, BLASTP (BLOSUM62 with the residues U = selenocysteine, O = pyrrolysine and J = leucine / isoleucine, as blastp scores them), GONNET1992, JONES, DAYHOFF, MDM78, BENNER6 / 22 / 74, JOHNSON, RISLER, GENETIC, FENG, LEVIN, MCLACHLAN, RAO, STR (protein; the fractional ones are scored exactly) or EDNAFULL, EDNAMAT, EDNASIMPLE, ENUC.4.4, BLASTN, MEGABLAST, HOXD70, TRANS (IUPAC DNA); or THE TEXT of a matrix instead of a name: several lines in the NCBI / EMBOSS / Biopython data-file format (a header line of symbols, then one row per symbol: the symbol and its scores; # lines are comments; scores may be decimals). When given, match/mismatch are ignored, residues are scored by the published table (any protein or IUPAC symbol of the matrix; case-insensitive) and gapOpen/gapExtend may be fractional (the EMBOSS needle/water models: gapOpen = -(G - E), gapExtend = -E for -gapopen G -gapextend E). The result then also carries the real-valued exactScore. |
| `matchScore` | number | No | Aligner scoring: score of an identical pair (giving any aligner setting scores the alignment like PairwiseAligner). |
| `mismatchScore` | number | No | Aligner scoring: score of a different pair. |
| `gapScoresJson` | string | No | Aligner scoring: a JSON object of Biopython gap-score properties (see the description). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `aligned`, `identities`, `mismatches`, `positives` | integer | The letter counts |
| `gaps`, `openGaps`, `extendGaps`, `leftGaps`, `internalGaps`, `rightGaps` | integer | The gap counts |
| `openLeftInsertions` ... `extendRightDeletions` | integer | The twelve gap counters |
| `gapScore`, `substitutionScore`, `score` | number | The scores (null without a scoring) |
| `text` | string | The text of Biopython's str(counts) |

## Errors

| Code | Message |
|------|---------|
| 1001 | Alignment JSON is required / not an alignment |
| 1002 | A letter is not in the substitution matrix |

## Examples

### Example 1: Example

**Input:** `{"alignmentJson": "{\"sequences\": [\"TTACGTTT\", \"ACGT\"], \"coordinates\": [[0, 2, 6, 8], [0, 0, 4, 4]]}"}`

**Response:** 4 aligned, 4 identities, 4 gaps (2 left deletions, 2 right deletions)

## See Also

- [alignment_operation](alignment_operation.md)
- [alignment_statistics](alignment_statistics.md)
