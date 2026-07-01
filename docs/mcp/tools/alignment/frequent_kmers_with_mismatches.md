# frequent_kmers_with_mismatches

Find the most-frequent k-mers within sequence allowing up to d mismatches (each window's full Hamming-d neighborhood is tallied). Returns all k-mers tied for the maximum Count_d.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Alignment |
| **Tool Name** | `frequent_kmers_with_mismatches` |
| **Method ID** | `ApproximateMatcher.FindFrequentKmersWithMismatches` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Find the most-frequent k-mers within sequence allowing up to d mismatches (each window's full Hamming-d neighborhood is tallied). Returns all k-mers tied for the maximum Count_d.

## Core Documentation Reference

- Source: [Seqeron.Genomics.Alignment/ApproximateMatcher.cs#L312](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/ApproximateMatcher.cs#L312)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | Sequence to analyze. |
| `k` | integer | Yes | K-mer length (> 0). |
| `d` | integer | Yes | Maximum mismatches in neighborhood (>= 0). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `items` | array<object> | Most-frequent k-mers (ties included), each with kmer and count. |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1008 | k must be > 0 |
| 1009 | d must be >= 0 |

## Examples

### Example 1: d=0 frequent words

**Tool Call:**
```json
{
  "tool": "frequent_kmers_with_mismatches",
  "arguments": {
    "sequence": "AAAAA",
    "k": 2,
    "d": 0
  }
}
```

**Response:**
```json
{
  "items": [
    {
      "kmer": "AA",
      "count": 4
    }
  ]
}
```

### Example 2: Frequent with mismatch

**Tool Call:**
```json
{
  "tool": "frequent_kmers_with_mismatches",
  "arguments": {
    "sequence": "AACAAGCTGATAAACATTTAAAGAG",
    "k": 3,
    "d": 0
  }
}
```

**Response:**
```json
{
  "items": [
    {
      "kmer": "AAA",
      "count": 2
    }
  ]
}
```

## Worked Example

`AAAAA` has four length-2 windows, all `AA`; with d=0 the only frequent k-mer is `AA` with Count_d = 4.

## References

- Algorithm source: [Seqeron.Genomics.Alignment/ApproximateMatcher.cs#L312](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/ApproximateMatcher.cs#L312)
- Binding: [AlignmentTools.cs](../../../../src/Seqeron/Mcp/Seqeron.Mcp.Alignment/Tools/AlignmentTools.cs)
