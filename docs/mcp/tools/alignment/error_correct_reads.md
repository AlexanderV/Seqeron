# error_correct_reads

Correct single-base substitution errors using the k-mer spectrum (Musket/Quake two-sided): k-mers with multiplicity below minKmerFrequency are untrusted; a position is corrected only when a unique alternative base makes its covering k-mers trusted.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Alignment |
| **Tool Name** | `error_correct_reads` |
| **Method ID** | `SequenceAssembler.ErrorCorrectReads` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Correct single-base substitution errors using the k-mer spectrum (Musket/Quake two-sided): k-mers with multiplicity below minKmerFrequency are untrusted; a position is corrected only when a unique alternative base makes its covering k-mers trusted.

## Core Documentation Reference

- Source: [Seqeron.Genomics.Alignment/SequenceAssembler.cs#L1104](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAssembler.cs#L1104)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `reads` | array<string> | Yes | Input reads. |
| `kmerSize` | integer | No | K-mer size used to build the frequency table (>= 1). (default: 21) |
| `minKmerFrequency` | integer | No | Minimum k-mer frequency considered trusted. (default: 3) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `corrected` | array<string> | Corrected reads, upper-cased, in input order (length and count preserved). |

## Errors

| Code | Message |
|------|---------|
| 1004 | k-mer size must be at least 1 |

## Examples

### Example 1: Fix one error

**Tool Call:**
```json
{
  "tool": "error_correct_reads",
  "arguments": {
    "reads": [
      "ACGTACGT",
      "ACGTACGT",
      "ACGTACGT",
      "ACGTTCGT"
    ],
    "kmerSize": 4,
    "minKmerFrequency": 2
  }
}
```

**Response:**
```json
{
  "corrected": [
    "ACGTACGT",
    "ACGTACGT",
    "ACGTACGT",
    "ACGTACGT"
  ]
}
```

### Example 2: Empty input

**Tool Call:**
```json
{
  "tool": "error_correct_reads",
  "arguments": {
    "reads": [],
    "kmerSize": 4,
    "minKmerFrequency": 2
  }
}
```

**Response:**
```json
{
  "corrected": []
}
```

## Worked Example

Three copies of `ACGTACGT` make its 4-mers trusted (freq >= 2). The fourth read `ACGTTCGT` has one middle-base error that is uniquely corrected back to `ACGTACGT`.

## References

- Algorithm source: [Seqeron.Genomics.Alignment/SequenceAssembler.cs#L1104](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAssembler.cs#L1104)
- Binding: [AlignmentTools.cs](../../../../src/Seqeron/Mcp/Seqeron.Mcp.Alignment/Tools/AlignmentTools.cs)
