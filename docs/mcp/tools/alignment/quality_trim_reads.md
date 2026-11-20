# quality_trim_reads

Trims Phred+33 quality-encoded reads from both ends, dropping bases whose decoded quality is below minQuality. Reads shorter than minLength after trimming are filtered out.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Alignment |
| **Tool Name** | `quality_trim_reads` |
| **Method ID** | `SequenceAssembler.QualityTrimReads` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Trims Phred+33 quality-encoded reads from both ends, dropping bases whose decoded quality is below minQuality. Reads shorter than minLength after trimming are filtered out.

## Core Documentation Reference

- Source: [Seqeron.Genomics.Alignment/SequenceAssembler.cs#L950](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAssembler.cs#L950)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `reads` | array<object> | Yes | Reads, each with sequence and matching Phred+33 quality string. |
| `minQuality` | integer | No | Minimum acceptable quality score (Phred). (default: 20) |
| `minLength` | integer | No | Minimum length after trimming for a read to be kept. (default: 50) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `trimmed` | array<string> | Surviving trimmed sequences, in input order. |

## Errors

| Code | Message |
|------|---------|
| 1011 | sequence and quality must be the same length |
| 1012 | read is null |

## Examples

### Example 1: High quality kept

**Tool Call:**
```json
{
  "tool": "quality_trim_reads",
  "arguments": {
    "reads": [
      {
        "sequence": "ACGTACGTAC",
        "quality": "IIIIIIIIII"
      }
    ],
    "minQuality": 20,
    "minLength": 5
  }
}
```

**Response:**
```json
{
  "trimmed": [
    "ACGTACGTAC"
  ]
}
```

### Example 2: Low quality dropped

**Tool Call:**
```json
{
  "tool": "quality_trim_reads",
  "arguments": {
    "reads": [
      {
        "sequence": "ACGTACGTAC",
        "quality": "!!!!!!!!!!"
      }
    ],
    "minQuality": 20,
    "minLength": 5
  }
}
```

**Response:**
```json
{
  "trimmed": []
}
```

## Worked Example

Phred+33 `I` = 40 (>= cutoff 20) so the read is kept whole; `!` = 0 (< cutoff) so the read is fully trimmed and dropped for being under minLength.

## References

- Algorithm source: [Seqeron.Genomics.Alignment/SequenceAssembler.cs#L950](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAssembler.cs#L950)
- Binding: [AlignmentTools.cs](../../../../src/Seqeron/Mcp/Seqeron.Mcp.Alignment/Tools/AlignmentTools.cs)
