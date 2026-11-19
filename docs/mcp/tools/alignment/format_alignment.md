# format_alignment

Render a pairwise alignment as three-line blocks in the EMBOSS needle/water (srspair) style: the markup line has '|' for identical columns, ':' for similar columns (positive substitution score), '.' for mismatches and ' ' for gap columns; wrapped to a configurable line width. matrix (BLOSUM62 | PAM250 | EDNAFULL) or the scalar match/mismatch decide which mismatches are similar.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Alignment |
| **Tool Name** | `format_alignment` |
| **Method ID** | `SequenceAligner.FormatAlignment` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Render a pairwise alignment as three-line blocks in the EMBOSS needle/water (srspair) style: the markup line has '|' for identical columns, ':' for similar columns (positive substitution score), '.' for mismatches and ' ' for gap columns; wrapped to a configurable line width. matrix (BLOSUM62 | PAM250 | EDNAFULL) or the scalar match/mismatch decide which mismatches are similar.

## Core Documentation Reference

- Source: [Seqeron.Genomics.Alignment/SequenceAligner.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.cs)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `alignedSequence1` | string | Yes | Aligned representation of sequence1. |
| `alignedSequence2` | string | Yes | Aligned representation of sequence2. |
| `lineWidth` | integer | No | Maximum characters per line (>= 1). (default: 60) |
| `matrix` | string | No | Optional substitution matrix deciding which non-gap columns are similar: any built-in matrix by EMBOSS or NCBI name - BLOSUM62, PAM250, EBLOSUM30 ... EBLOSUM90, EPAM10 ... EPAM500 (protein) or EDNAFULL, EDNAMAT, EDNASIMPLE, ENUC.4.4 (IUPAC DNA). Omit to use match/mismatch. |
| `match` | integer | No | Scalar model (used when matrix is omitted): score of an identical pair. (default: 1) |
| `mismatch` | integer | No | Scalar model (used when matrix is omitted): score of a non-identical pair; > 0 marks mismatches ':' (similar). (default: -1) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `formatted` | string | The formatted multi-line alignment text. |

## Errors

| Code | Message |
|------|---------|
| 1001 | Aligned sequence cannot be null or empty |
| 1007 | lineWidth must be >= 1 |

## Examples

### Example 1: Single block

**Tool Call:**
```json
{
  "tool": "format_alignment",
  "arguments": {
    "alignedSequence1": "ACGT",
    "alignedSequence2": "ACGA",
    "lineWidth": 80
  }
}
```

**Response:**
```json
{
  "formatted": "ACGT\n|||.\nACGA\n\n"
}
```

### Example 2: With a gap

**Tool Call:**
```json
{
  "tool": "format_alignment",
  "arguments": {
    "alignedSequence1": "ACGT",
    "alignedSequence2": "AC-T",
    "lineWidth": 80
  }
}
```

**Response:**
```json
{
  "formatted": "ACGT\n||  \nAC-T\n\n"
}
```

### Example 3: BLOSUM62 similarity marks

**Tool Call:**
```json
{
  "tool": "format_alignment",
  "arguments": {
    "alignedSequence1": "LAAK-A",
    "alignedSequence2": "ISTRAA",
    "lineWidth": 80,
    "matrix": "BLOSUM62"
  }
}
```

**Response:**
```json
{
  "formatted": "LAAK-A\n::.: |\nISTRAA\n\n"
}
```

## Worked Example

`ACGT` vs `ACGA` shares three identical columns (`|||`) and one mismatch (`.`), rendered as a three-line block followed by a blank separator line.

## References

- Algorithm source: [Seqeron.Genomics.Alignment/SequenceAligner.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.cs)
- Binding: [AlignmentTools.cs](../../../../src/Seqeron/Mcp/Seqeron.Mcp.Alignment/Tools/AlignmentTools.cs)
