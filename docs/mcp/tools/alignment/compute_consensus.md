# compute_consensus

Build a consensus from pre-aligned reads (same length; '-'/'.' ignored) by column-wise majority vote (Biopython dumb_consensus, threshold 0.5). Ties and sub-threshold columns emit 'N'.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Alignment |
| **Tool Name** | `compute_consensus` |
| **Method ID** | `SequenceAssembler.ComputeConsensus` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Build a consensus from pre-aligned reads (same length; '-'/'.' ignored) by column-wise majority vote (Biopython dumb_consensus, threshold 0.5). Ties and sub-threshold columns emit 'N'.

## Core Documentation Reference

- Source: [Seqeron.Genomics.Alignment/SequenceAssembler.cs#L861](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAssembler.cs#L861)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `alignedReads` | array<string> | Yes | Pre-aligned reads of identical length. |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `consensus` | string | The consensus sequence (length equals the longest read). |

## Errors

| Code | Message |
|------|---------|
| 1003 | All aligned reads must have the same length |

## Examples

### Example 1: Majority vote

**Tool Call:**
```json
{
  "tool": "compute_consensus",
  "arguments": {
    "alignedReads": [
      "ACGT",
      "ACGT",
      "ACGA"
    ]
  }
}
```

**Response:**
```json
{
  "consensus": "ACGT"
}
```

### Example 2: Tie -> ambiguous

**Tool Call:**
```json
{
  "tool": "compute_consensus",
  "arguments": {
    "alignedReads": [
      "A",
      "C"
    ]
  }
}
```

**Response:**
```json
{
  "consensus": "N"
}
```

## Worked Example

For columns `A`,`C`,`G`,`{T,T,A}`, the fourth column has T at 2/3 (>= 0.5) so the consensus is `ACGT`. A 1:1 tie emits `N`.

## References

- Algorithm source: [Seqeron.Genomics.Alignment/SequenceAssembler.cs#L861](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAssembler.cs#L861)
- Binding: [AlignmentTools.cs](../../../../src/Seqeron/Mcp/Seqeron.Mcp.Alignment/Tools/AlignmentTools.cs)
