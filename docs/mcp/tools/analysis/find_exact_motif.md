# find_exact_motif

Exact-match motif positions in a DNA sequence via suffix tree.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `find_exact_motif` |
| **Method ID** | `MotifFinder.FindExactMotif` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Returns the 0-based start positions of every exact occurrence of a motif in a DNA
sequence, located via the sequence's suffix tree. Matching is case-insensitive and
overlapping occurrences are reported. Positions are enumerated in suffix-tree order
(not necessarily sorted).

## Core Documentation Reference

- Source: [MotifFinder.cs#L24](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/MotifFinder.cs#L24)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | DNA sequence to search (min length 1) |
| `motif` | string | Yes | Motif pattern (min length 1) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `positions` | array of integer | 0-based start positions of every occurrence |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1001 | Invalid DNA sequence |
| 1001 | Motif cannot be null or empty |

## Examples

### Example 1: Repeated motif (ATG in ATGATG)

**User Prompt:**
> Exact positions of "ATG" in "ATGATG".

**Expected Tool Call:**
```json
{
  "tool": "find_exact_motif",
  "arguments": { "sequence": "ATGATG", "motif": "ATG" }
}
```

**Response:**
```json
{ "positions": [0, 3] }
```

### Example 2: Overlapping (AA in AAAA)

**User Prompt:**
> Positions of "AA" in "AAAA".

**Expected Tool Call:**
```json
{
  "tool": "find_exact_motif",
  "arguments": { "sequence": "AAAA", "motif": "AA" }
}
```

**Response:**
```json
{ "positions": [0, 1, 2] }
```

## Performance

- **Time Complexity:** O(n) suffix-tree build + O(m + occ) per query.
- **Space Complexity:** O(n) for the suffix tree.

## See Also

- [find_motif](find_motif.md) — GenomicAnalyzer equivalent
- [find_degenerate_motif](find_degenerate_motif.md) — IUPAC ambiguity search
