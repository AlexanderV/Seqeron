# find_motif

Exact motif occurrences in a DNA sequence via suffix tree.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `find_motif` |
| **Method ID** | `GenomicAnalyzer.FindMotif` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Returns the 0-based start positions of **every exact occurrence** of a motif in a DNA
sequence, located with the sequence's suffix tree. Matching is case-insensitive and
overlapping occurrences are all reported (e.g. `AA` in `AAAA` yields `0, 1, 2`).
Positions are enumerated in suffix-tree DFS order (not necessarily sorted).

## Core Documentation Reference

- Source: [GenomicAnalyzer.cs#L164](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/GenomicAnalyzer.cs#L164)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | DNA sequence (min length 1) |
| `motif` | string | Yes | Motif to locate (min length 1) |

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
> Where does the motif "ATG" occur in "ATGATG"?

**Expected Tool Call:**
```json
{
  "tool": "find_motif",
  "arguments": { "sequence": "ATGATG", "motif": "ATG" }
}
```

**Response:**
```json
{ "positions": [0, 3] }
```

### Example 2: Overlapping occurrences (AA in AAAA)

**User Prompt:**
> Positions of "AA" in "AAAA".

**Expected Tool Call:**
```json
{
  "tool": "find_motif",
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

- [find_exact_motif](find_exact_motif.md) — MotifFinder exact search
- [find_known_motifs](find_known_motifs.md) — search a set of motifs
- [kmer_positions](kmer_positions.md) — k-mer occurrence positions
