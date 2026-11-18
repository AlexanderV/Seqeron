# alignment_statistics

Compute identity / similarity / gap statistics for a pairwise alignment (EMBOSS needle convention: every percentage is over the full alignment length including gap columns; a column is a gap column when either row has '-'; residues compared case-insensitively). Similarity counts columns with a positive substitution score (BLAST positives): with matrix = BLOSUM62 | PAM250 | EDNAFULL the published residue-pair table is used, otherwise the scalar match/mismatch model. Also returns positives and gapOpens (number of gap runs).

## Overview

| Property | Value |
|----------|-------|
| **Server** | Alignment |
| **Tool Name** | `alignment_statistics` |
| **Method ID** | `SequenceAligner.CalculateStatistics` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Compute identity / similarity / gap statistics for a pairwise alignment (EMBOSS needle convention: every percentage is over the full alignment length including gap columns; a column is a gap column when either row has '-'; residues compared case-insensitively). Similarity counts columns with a positive substitution score (BLAST positives): with matrix = BLOSUM62 | PAM250 | EDNAFULL the published residue-pair table is used, otherwise the scalar match/mismatch model. Also returns positives and gapOpens (number of gap runs).

## Core Documentation Reference

- Source: [Seqeron.Genomics.Alignment/SequenceAligner.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.cs)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `alignedSequence1` | string | Yes | Aligned representation of sequence1 (with '-' gaps). |
| `alignedSequence2` | string | Yes | Aligned representation of sequence2 (must equal alignedSequence1.Length). |
| `matrix` | string | No | Optional substitution matrix deciding which non-gap columns are similar: BLOSUM62, PAM250 (protein) or EDNAFULL (IUPAC DNA). Omit to use match/mismatch. |
| `match` | integer | No | Scalar model (used when matrix is omitted): score of an identical pair; the pair is 'similar' iff this is > 0. (default: 1) |
| `mismatch` | integer | No | Scalar model (used when matrix is omitted): score of a non-identical pair; the pair is 'similar' iff this is > 0. (default: -1) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `matches` | integer | Number of identical columns. |
| `mismatches` | integer | Number of mismatched (non-gap, non-identical) columns. |
| `gaps` | integer | Number of columns containing a gap in either sequence. |
| `alignmentLength` | integer | Total alignment length including gap columns. |
| `identity` | number | Percent identity = matches / length * 100. |
| `similarity` | number | Percent similarity = positives / length * 100 (positive-score columns; an identical column whose own score is not positive, e.g. X/X in BLOSUM62, is not similar). |
| `gapPercent` | number | Percent gap columns = gaps / length * 100. |
| `positives` | integer | Number of similar columns (BLAST positives / EMBOSS similarity count). |
| `gapOpens` | integer | Number of gap runs (maximal '-' runs per row; BLAST tabular gapopen, Biopython open_gaps). |

## Errors

| Code | Message |
|------|---------|
| 1001 | Aligned sequence cannot be null or empty |
| 1002 | Aligned sequences must have equal length |

## Examples

### Example 1: One gap column

**Tool Call:**
```json
{
  "tool": "alignment_statistics",
  "arguments": {
    "alignedSequence1": "ACGT-A",
    "alignedSequence2": "ACGTGA"
  }
}
```

**Response:**
```json
{
  "matches": 5,
  "mismatches": 0,
  "gaps": 1,
  "alignmentLength": 6,
  "identity": 83.33333333333334,
  "similarity": 83.33333333333334,
  "gapPercent": 16.666666666666664,
  "positives": 5,
  "gapOpens": 1
}
```

### Example 2: Perfect identity

**Tool Call:**
```json
{
  "tool": "alignment_statistics",
  "arguments": {
    "alignedSequence1": "ACGT",
    "alignedSequence2": "ACGT"
  }
}
```

**Response:**
```json
{
  "matches": 4,
  "mismatches": 0,
  "gaps": 0,
  "alignmentLength": 4,
  "identity": 100,
  "similarity": 100,
  "gapPercent": 0,
  "positives": 4,
  "gapOpens": 0
}
```

### Example 3: BLOSUM62 similarity

**Tool Call:**
```json
{
  "tool": "alignment_statistics",
  "arguments": {
    "alignedSequence1": "LAAK-A",
    "alignedSequence2": "ISTRAA",
    "matrix": "BLOSUM62"
  }
}
```

**Response:**
```json
{
  "matches": 1,
  "mismatches": 4,
  "gaps": 1,
  "alignmentLength": 6,
  "identity": 16.666666666666664,
  "similarity": 66.66666666666666,
  "gapPercent": 16.666666666666664,
  "positives": 4,
  "gapOpens": 1
}
```

## Worked Example

For `ACGT-A` vs `ACGTGA`, five columns are identical, one column is a gap, and the total length is 6. Identity = 5/6 * 100 = 83.33%. With `matrix = BLOSUM62`, L/I (+2), A/S (+1), K/R (+2) and A/A (+4) are positive, A/T (0) is not: 4 positives over 6 columns = 66.67 % similarity. The published EMBOSS needle example (HBA_HUMAN vs HBB_HUMAN, EBLOSUM62) is reproduced exactly: length 149, identity 65, similarity 90, gaps 9, 4 gap runs.

## References

- Algorithm source: [Seqeron.Genomics.Alignment/SequenceAligner.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.cs)
- Binding: [AlignmentTools.cs](../../../../src/Seqeron/Mcp/Seqeron.Mcp.Alignment/Tools/AlignmentTools.cs)
