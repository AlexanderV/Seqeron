# sequence_identity

Compute percent identity (fraction in [0,1]) between two equal-length sequences by case-insensitive position-by-position comparison. Returns 0 for unequal-length inputs.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Alignment |
| **Tool Name** | `sequence_identity` |
| **Method ID** | `SequenceAssembler.CalculateIdentity` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Compute percent identity (fraction in [0,1]) between two equal-length sequences by case-insensitive position-by-position comparison. Returns 0 for unequal-length inputs.

## Core Documentation Reference

- Source: [Seqeron.Genomics.Alignment/SequenceAssembler.cs#L247](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAssembler.cs#L247)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence1` | string | Yes | First sequence. |
| `sequence2` | string | Yes | Second sequence (must equal sequence1.Length). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `identity` | number | Fraction of matching positions in [0,1]; 0 if lengths differ. |

## Examples

### Example 1: Three of four match

**Tool Call:**
```json
{
  "tool": "sequence_identity",
  "arguments": {
    "sequence1": "ACGT",
    "sequence2": "ACGA"
  }
}
```

**Response:**
```json
{
  "identity": 0.75
}
```

### Example 2: Identical

**Tool Call:**
```json
{
  "tool": "sequence_identity",
  "arguments": {
    "sequence1": "ACGT",
    "sequence2": "ACGT"
  }
}
```

**Response:**
```json
{
  "identity": 1
}
```

## Worked Example

`ACGT` vs `ACGA` matches 3 of 4 positions, so identity = 0.75.

## References

- Algorithm source: [Seqeron.Genomics.Alignment/SequenceAssembler.cs#L247](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAssembler.cs#L247)
- Binding: [AlignmentTools.cs](../../../../src/Seqeron/Mcp/Seqeron.Mcp.Alignment/Tools/AlignmentTools.cs)
