# find_known_motifs

Search a DNA sequence for a set of known motifs simultaneously (exact, or IUPAC / both strands).

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

By default the search is exact and forward-strand only (IUPAC symbols are literal and match
nothing). `iupac = true` interprets IUPAC ambiguity codes (`TATAWAW`); `bothStrands = true`
also reports, in `reverseMatches`, the windows that match the reverse complement of each motif
(forward-strand start coordinates; a palindromic site such as `GAATTC` appears in both maps).
Verified against EMBOSS fuzznuc 6.6.0 (`-pmismatch 0 -complement Y`) on 2448 sequence/pattern pairs.

## Core Documentation Reference

- Source: [GenomicAnalyzer.cs#L257](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/GenomicAnalyzer.cs#L257)
- Algorithm: `docs/algorithms/Motif_Analysis/Known_Motif_Search.md`

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | DNA sequence (min length 1) |
| `motifs` | array of string | Yes | Set of motifs to search |
| `bothStrands` | boolean | No | Also search the reverse strand (default false) |
| `iupac` | boolean | No | Interpret IUPAC codes in the motifs (default false) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `matches` | object | Map from upper-cased motif → ascending 0-based forward-strand positions |
| `reverseMatches` | object | Map from upper-cased motif → ascending forward-strand window starts matching the motif's reverse complement (empty unless `bothStrands`) |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1001 | Invalid DNA sequence |
| 1001 | Motifs cannot be null |
| 1001 | Invalid IUPAC code in motif (`iupac = true` only) |

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

### Example 3: Both strands, IUPAC (EMBOSS fuzznuc)

```json
{ "tool": "find_known_motifs",
  "arguments": { "sequence": "GAATTCAAAGAATTCTTTATATAAAGTTCAAGAAC", "motifs": ["TATAWAW"], "bothStrands": true, "iupac": true } }
```
**Response:** `{ "matches": { "TATAWAW": [17] }, "reverseMatches": { "TATAWAW": [16] } }`
(fuzznuc: 18-24 `+`, 17-23 `-`, 1-based).

## Performance

- **Time Complexity:** fewer than 128 distinct queries (motifs × strands): one KMP scan per query, O(k·n); from 128 queries: O(n) suffix-tree build + O(|m| + occ·log occ) per query.
- **Space Complexity:** O(occ) for scans; O(n) for the suffix tree (heavy: 1.7 GB at 10 Mbp).

## See Also

- [find_motif](find_motif.md) — single-motif search
- [find_regulatory_elements](find_regulatory_elements.md) — built-in regulatory motif set
