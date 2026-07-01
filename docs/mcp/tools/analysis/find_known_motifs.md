# find_known_motifs

Search a DNA sequence for a set of known motifs simultaneously.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `find_known_motifs` |
| **Method ID** | `GenomicAnalyzer.FindKnownMotifs` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Solves the **exact set-matching problem**: given a set of motifs, returns for each
motif that occurs its 0-based start positions (sorted ascending) in the sequence.
Matching is case-insensitive (keys are upper-cased) and overlapping occurrences are
all reported (`AA` in `AAAA` → 0, 1, 2). Empty/whitespace motifs are skipped and
motifs with no occurrence are omitted from the result.

## Core Documentation Reference

- Source: [GenomicAnalyzer.cs#L218](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/GenomicAnalyzer.cs#L218)
- Algorithm: `docs/algorithms/Motif_Analysis/Known_Motif_Search.md`

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | DNA sequence (min length 1) |
| `motifs` | array of string | Yes | Set of motifs to search |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `matches` | object | Map from upper-cased motif → ascending 0-based positions |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1001 | Invalid DNA sequence |
| 1001 | Motifs cannot be null |

## Examples

### Example 1: Two hits, one miss

**User Prompt:**
> Search "ATGATGCC" for ATG, CC and GG.

**Expected Tool Call:**
```json
{
  "tool": "find_known_motifs",
  "arguments": { "sequence": "ATGATGCC", "motifs": ["ATG", "CC", "GG"] }
}
```

**Response:**
```json
{ "matches": { "ATG": [0, 3], "CC": [6] } }
```
GG does not occur, so it is omitted.

### Example 2: Overlapping motif

**User Prompt:**
> Find "AA" in "AAAA".

**Expected Tool Call:**
```json
{
  "tool": "find_known_motifs",
  "arguments": { "sequence": "AAAA", "motifs": ["AA"] }
}
```

**Response:**
```json
{ "matches": { "AA": [0, 1, 2] } }
```

## Performance

- **Time Complexity:** O(n) suffix-tree build + O(|m| + occ) per motif.
- **Space Complexity:** O(n) for the suffix tree.

## See Also

- [find_motif](find_motif.md) — single-motif search
- [find_regulatory_elements](find_regulatory_elements.md) — built-in regulatory motif set
