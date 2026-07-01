# kmer_positions

Zero-based positions of all (overlapping) occurrences of a k-mer.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `kmer_positions` |
| **Method ID** | `KmerAnalyzer.FindKmerPositions` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Solves the **Pattern Matching Problem** (Rosalind BA1D; Compeau & Pevzner): returns
the ascending 0-based start positions of every occurrence of `kmer` in `sequence`.
Occurrences may overlap and every overlapping start is reported — e.g. `AA` in `AAAA`
yields `0, 1, 2`. Matching is case-insensitive. The result is empty when the k-mer is
longer than the sequence or does not occur.

## Core Documentation Reference

- Source: [KmerAnalyzer.cs#L432](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/KmerAnalyzer.cs#L432)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | Sequence to scan (min length 1) |
| `kmer` | string | Yes | k-mer / pattern to locate (min length 1) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `positions` | array of integer | Ascending 0-based start positions of every overlapping occurrence |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1001 | k-mer cannot be null or empty |

## Examples

### Example 1: Overlapping occurrences (AA in AAAA)

**User Prompt:**
> Where does "AA" occur in "AAAA"?

**Expected Tool Call:**
```json
{
  "tool": "kmer_positions",
  "arguments": { "sequence": "AAAA", "kmer": "AA" }
}
```

**Response:**
```json
{ "positions": [0, 1, 2] }
```
All three overlapping starts are reported.

### Example 2: Repeated trimer (ATG in ATGATG)

**User Prompt:**
> Positions of "ATG" in "ATGATG".

**Expected Tool Call:**
```json
{
  "tool": "kmer_positions",
  "arguments": { "sequence": "ATGATG", "kmer": "ATG" }
}
```

**Response:**
```json
{ "positions": [0, 3] }
```

## Performance

- **Time Complexity:** O(n·k) single forward scan.
- **Space Complexity:** O(number of occurrences).

## See Also

- [count_kmers](count_kmers.md) — occurrence counts
- [find_motif](find_motif.md) — suffix-tree motif search
