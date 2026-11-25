# calculate_coverage

Align each read end-to-end without gaps at its minimum-Hamming position on the reference (placed when >= minOverlap bases match; one placement per read) and return per-base coverage depth (samtools depth semantics) as an integer array of length reference.Length. Optionally search the reverse-complement strand and skip multi-mapped reads.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Alignment |
| **Tool Name** | `calculate_coverage` |
| **Method ID** | `SequenceAssembler.CalculateCoverage` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Align each read end-to-end without gaps at its minimum-Hamming position on the reference (placed when >= minOverlap bases match; one placement per read) and return per-base coverage depth (samtools depth semantics) as an integer array of length reference.Length. Optionally search the reverse-complement strand and skip multi-mapped reads.

## Core Documentation Reference

- Source: [Seqeron.Genomics.Alignment/SequenceAssembler.cs#L1369](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAssembler.cs#L1369)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `reference` | string | Yes | Reference sequence. |
| `reads` | array<string> | Yes | Reads to map onto reference. |
| `minOverlap` | integer | No | Minimum matching bases required to map a read, >= 1. (default: 20) |
| `includeReverseComplement` | boolean | No | Also align each read's reverse complement (both strands, as Bowtie 2/BWA do by default). (default: false) |
| `excludeMultiMapped` | boolean | No | Skip reads whose best score is reached at more than one locus (Bowtie 2 MAPQ <= 1; `samtools depth -Q 2`). (default: false) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `coverage` | array<integer> | Per-base coverage depth; coverage[i] is the number of placed reads spanning reference position i. |

## Examples

### Example 1: Two reads

**Tool Call:**
```json
{
  "tool": "calculate_coverage",
  "arguments": {
    "reference": "ACGTACGTAC",
    "reads": [
      "ACGTA",
      "GTACG"
    ],
    "minOverlap": 5
  }
}
```

**Response:**
```json
{
  "coverage": [
    1,
    1,
    2,
    2,
    2,
    1,
    1,
    0,
    0,
    0
  ]
}
```

### Example 2: No reads

**Tool Call:**
```json
{
  "tool": "calculate_coverage",
  "arguments": {
    "reference": "ACGTACGTAC",
    "reads": [],
    "minOverlap": 5
  }
}
```

**Response:**
```json
{
  "coverage": [
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0
  ]
}
```

### Example 3: Reverse-strand read

`TTGACC` is the reverse complement of `GGTCAA` at [2,8).

```json
{ "tool": "calculate_coverage",
  "arguments": { "reference": "ACGGTCAATGCC", "reads": ["TTGACC"], "minOverlap": 6, "includeReverseComplement": true } }
```

Response: `{ "coverage": [0,0,1,1,1,1,1,1,0,0,0,0] }` (forward strand only: all zeros).

## Errors

| Condition | Error |
|-----------|-------|
| `minOverlap` < 1 | ArgumentOutOfRangeException |

## Worked Example

Read `ACGTA` matches at positions 0 and 4 and is counted once at the leftmost (0, covers [0,5)); read `GTACG` places at position 2 (covers [2,7)). Depths sum per position. With `excludeMultiMapped: true` the multi-mapped `ACGTA` is skipped: `[0,0,1,1,1,1,1,0,0,0]`. All values cross-checked with `samtools depth -a` (Evidence ASSEMBLY-COVER-001).

## References

- Algorithm source: [Seqeron.Genomics.Alignment/SequenceAssembler.cs#L1369](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAssembler.cs#L1369)
- Binding: [AlignmentTools.cs](../../../../src/Seqeron/Mcp/Seqeron.Mcp.Alignment/Tools/AlignmentTools.cs)
