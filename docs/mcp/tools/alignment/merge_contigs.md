# merge_contigs

Concatenate two contigs, collapsing the specified suffix(contig1)/prefix(contig2) overlap of overlapLength bases. A non-positive or too-large overlap falls back to plain concatenation.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Alignment |
| **Tool Name** | `merge_contigs` |
| **Method ID** | `SequenceAssembler.MergeContigs` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Concatenate two contigs, collapsing the specified suffix(contig1)/prefix(contig2) overlap of overlapLength bases. A non-positive or too-large overlap falls back to plain concatenation.

## Core Documentation Reference

- Source: [Seqeron.Genomics.Alignment/SequenceAssembler.cs#L629](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAssembler.cs#L629)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `contig1` | string | Yes | First contig (suffix donor). |
| `contig2` | string | Yes | Second contig (prefix donor). |
| `overlapLength` | integer | Yes | Length of the overlap to collapse. |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `merged` | string | The merged superstring. |

## Examples

### Example 1: Collapse 5-bp overlap

**Tool Call:**
```json
{
  "tool": "merge_contigs",
  "arguments": {
    "contig1": "AAAAACCCCC",
    "contig2": "CCCCCGGGGG",
    "overlapLength": 5
  }
}
```

**Response:**
```json
{
  "merged": "AAAAACCCCCGGGGG"
}
```

### Example 2: No overlap -> concat

**Tool Call:**
```json
{
  "tool": "merge_contigs",
  "arguments": {
    "contig1": "AAA",
    "contig2": "GGG",
    "overlapLength": 0
  }
}
```

**Response:**
```json
{
  "merged": "AAAGGG"
}
```

## Worked Example

With overlap 5 the merged length is 10 + 10 - 5 = 15: `AAAAACCCCC` + `GGGGG`.

## References

- Algorithm source: [Seqeron.Genomics.Alignment/SequenceAssembler.cs#L629](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAssembler.cs#L629)
- Binding: [AlignmentTools.cs](../../../../src/Seqeron/Mcp/Seqeron.Mcp.Alignment/Tools/AlignmentTools.cs)
