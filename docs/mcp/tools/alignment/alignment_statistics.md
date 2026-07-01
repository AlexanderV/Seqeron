# alignment_statistics

Compute match/mismatch/gap counts and percent identity, similarity, and gap percent for a pairwise alignment (EMBOSS needle convention: denominator is the full alignment length including gap columns).

## Overview

| Property | Value |
|----------|-------|
| **Server** | Alignment |
| **Tool Name** | `alignment_statistics` |
| **Method ID** | `SequenceAligner.CalculateStatistics` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Compute match/mismatch/gap counts and percent identity, similarity, and gap percent for a pairwise alignment (EMBOSS needle convention: denominator is the full alignment length including gap columns).

## Core Documentation Reference

- Source: [Seqeron.Genomics.Alignment/SequenceAligner.cs#L570](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.cs#L570)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `alignedSequence1` | string | Yes | Aligned representation of sequence1 (with '-' gaps). |
| `alignedSequence2` | string | Yes | Aligned representation of sequence2 (must equal alignedSequence1.Length). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `matches` | integer | Number of identical columns. |
| `mismatches` | integer | Number of mismatched (non-gap, non-identical) columns. |
| `gaps` | integer | Number of columns containing a gap in either sequence. |
| `alignmentLength` | integer | Total alignment length including gap columns. |
| `identity` | number | Percent identity = matches / length * 100. |
| `similarity` | number | Percent similarity = (matches + similar substitutions) / length * 100. |
| `gapPercent` | number | Percent gap columns = gaps / length * 100. |

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
  "gapPercent": 16.666666666666664
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
  "gapPercent": 0
}
```

## Worked Example

For `ACGT-A` vs `ACGTGA`, five columns are identical, one column is a gap, and the total length is 6. Identity = 5/6 * 100 = 83.33%.

## References

- Algorithm source: [Seqeron.Genomics.Alignment/SequenceAligner.cs#L570](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.cs#L570)
- Binding: [AlignmentTools.cs](../../../../src/Seqeron/Mcp/Seqeron.Mcp.Alignment/Tools/AlignmentTools.cs)
