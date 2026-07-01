# format_alignment

Render a human-readable three-line alignment block (EMBOSS srspair legend: '|' identity, ':' similarity, ' ' gap/mismatch), wrapped to a configurable line width.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Alignment |
| **Tool Name** | `format_alignment` |
| **Method ID** | `SequenceAligner.FormatAlignment` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Render a human-readable three-line alignment block (EMBOSS srspair legend: '|' identity, ':' similarity, ' ' gap/mismatch), wrapped to a configurable line width.

## Core Documentation Reference

- Source: [Seqeron.Genomics.Alignment/SequenceAligner.cs#L634](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.cs#L634)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `alignedSequence1` | string | Yes | Aligned representation of sequence1. |
| `alignedSequence2` | string | Yes | Aligned representation of sequence2. |
| `lineWidth` | integer | No | Maximum characters per line (>= 1). (default: 60) |

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
  "formatted": "ACGT\n||| \nACGA\n\n"
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

## Worked Example

`ACGT` vs `ACGA` shares three identical columns (`|||`) and one mismatch (space), rendered as a three-line block followed by a blank separator line.

## References

- Algorithm source: [Seqeron.Genomics.Alignment/SequenceAligner.cs#L634](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.cs#L634)
- Binding: [AlignmentTools.cs](../../../../src/Seqeron/Mcp/Seqeron.Mcp.Alignment/Tools/AlignmentTools.cs)
