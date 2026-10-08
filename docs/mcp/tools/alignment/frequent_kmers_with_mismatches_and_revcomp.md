# frequent_kmers_with_mismatches_and_revcomp

Most frequent k-mers with mismatches and reverse complements (ROSALIND BA1J).

## Overview

| Property | Value |
|----------|-------|
| **Server** | Alignment |
| **Tool Name** | `frequent_kmers_with_mismatches_and_revcomp` |
| **Method ID** | `ApproximateMatcher.FindFrequentKmersWithMismatchesAndReverseComplements` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Frequent Words with Mismatches and Reverse Complements (Compeau & Pevzner, Bioinformatics Algorithms ch. 1; ROSALIND BA1J): all DNA k-mers P maximising Count_d(sequence, P) + Count_d(sequence, reverseComplement(P)), ties included. Neighbourhoods are over {A,C,G,T} (windows with other symbols contribute nothing beyond their DNA neighbours, as in BA1I); case-insensitive. A reverse-complement palindrome scores 2·Count_d. The result set is closed under reverse complement. The library leaves the order unspecified; the tool sorts items by k-mer (ordinal).

## Core Documentation Reference

- Source: [Seqeron.Genomics.Alignment/ApproximateMatcher.cs#L967](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/ApproximateMatcher.cs#L967)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | Sequence to analyze. |
| `k` | integer | Yes | K-mer length (> 0). |
| `d` | integer | Yes | Maximum mismatches in neighborhood (>= 0). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `items` | array<object> | `{kmer, count}` with count = Count_d(P) + Count_d(rc(P)), sorted by k-mer |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1007 | k must be > 0 |
| 1008 | d must be >= 0 |

## Examples

### Example 1: ROSALIND BA1J sample

**User Prompt:**
> Most frequent 4-mers with up to 1 mismatch, counting reverse complements, in the BA1J sample.

**Tool Call:**
```json
{
  "tool": "frequent_kmers_with_mismatches_and_revcomp",
  "arguments": {
    "sequence": "ACGTTGCATGTCGCATGATGCATGAGAGCT",
    "k": 4,
    "d": 1
  }
}
```

**Response:**
```json
{
  "items": [
    {
      "kmer": "ACAT",
      "count": 9
    },
    {
      "kmer": "ATGT",
      "count": 9
    }
  ]
}
```

### Example 2: Homopolymer

**Tool Call:**
```json
{
  "tool": "frequent_kmers_with_mismatches_and_revcomp",
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
    },
    {
      "kmer": "TT",
      "count": 4
    }
  ]
}
```

## Worked Example

The published BA1J answer is {ATGT, ACAT} (each other's reverse complement) with combined count 9.

## See Also

- [frequent_kmers_with_mismatches](frequent_kmers_with_mismatches.md) — BA1I — without reverse complements

## References

- Algorithm source: [Seqeron.Genomics.Alignment/ApproximateMatcher.cs#L967](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/ApproximateMatcher.cs#L967)
- Binding: [AlignmentTools.cs](../../../../src/Seqeron/Mcp/Seqeron.Mcp.Alignment/Tools/AlignmentTools.cs)
