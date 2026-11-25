# merge_contigs

Merge two contigs over a known suffix(contig1)/prefix(contig2) overlap: contig1 + contig2 without its first overlapLength bases (overlapLength 0 = concatenation; must be in [0, shorter contig length], else an error). By default contig1's bases are kept in the overlap; resolveOverlapConflicts = true resolves overlap mismatches/Ns with the EMBOSS merger rule (N yields to the other base, otherwise the base from the contig with the higher local base-quality window wins; ties keep contig1).

## Overview

| Property | Value |
|----------|-------|
| **Server** | Alignment |
| **Tool Name** | `merge_contigs` |
| **Method ID** | `SequenceAssembler.MergeContigs` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Merge two contigs over a known suffix(contig1)/prefix(contig2) overlap: contig1 + contig2 without its first overlapLength bases (overlapLength 0 = concatenation; must be in [0, shorter contig length], else an error). By default contig1's bases are kept in the overlap; resolveOverlapConflicts = true resolves overlap mismatches/Ns with the EMBOSS merger rule (N yields to the other base, otherwise the base from the contig with the higher local base-quality window wins; ties keep contig1).

## Core Documentation Reference

- Source: [Seqeron.Genomics.Alignment/SequenceAssembler.cs#L950](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAssembler.cs#L950)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `contig1` | string | Yes | First contig (suffix donor). |
| `contig2` | string | Yes | Second contig (prefix donor). |
| `overlapLength` | integer | Yes | Length of the suffix/prefix overlap to collapse, in [0, min(\|contig1\|, \|contig2\|)]. |
| `resolveOverlapConflicts` | boolean | No | Resolve mismatches / N in the overlap by the EMBOSS merger rule instead of keeping contig1's bases (default `false`). |

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

### Example 3: Overlap mismatch resolved by the EMBOSS merger rule

**Tool Call:**
```json
{
  "tool": "merge_contigs",
  "arguments": {
    "contig1": "AAAACCCCGGGG",
    "contig2": "CCCCGTGGTTTT",
    "overlapLength": 8,
    "resolveOverlapConflicts": true
  }
}
```

**Response:**
```json
{
  "merged": "AAAACCCCGTGGTTTT"
}
```

## Errors

| Code | Message |
|------|---------|
| 4001 | overlapLength is negative or exceeds the shorter contig |

## Worked Example

With overlap 5 the merged length is 10 + 10 - 5 = 15: `AAAAACCCCC` + `GGGGG`.

Example 3: the overlap column 6 holds `G` (contig1, 3 bases from its end) vs `T` (contig2, interior). EMBOSS `merger_quality` with window 5 scores contig1 14/11 (its forward window runs off the end) and contig2 18/11, so the contig2 base `T` is used; the default mode would keep `G` (`AAAACCCCGGGGTTTT`).

## References

- Algorithm source: [Seqeron.Genomics.Alignment/SequenceAssembler.cs#L950](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAssembler.cs#L950)
- Binding: [AlignmentTools.cs](../../../../src/Seqeron/Mcp/Seqeron.Mcp.Alignment/Tools/AlignmentTools.cs)
