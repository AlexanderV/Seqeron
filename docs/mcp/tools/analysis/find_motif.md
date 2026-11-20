# find_motif

Exact motif occurrences in a DNA sequence (ascending, forward strand).

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
sequence (forward strand), in **ascending** order. Matching is case-insensitive and
overlapping occurrences are all reported (e.g. `AA` in `AAAA` yields `0, 1, 2`).
IUPAC symbols are literal and match nothing; use `find_known_motifs` with `iupac` /
`bothStrands` for degenerate motifs and the reverse strand. Same engine and result as
`find_exact_motif` (`MotifFinder.FindExactMotif`); the suffix tree is built per call, so for
one-off queries on very large sequences prefer `kmer_positions` (KMP, no index).

**IUPAC input.** The sequence may contain IUPAC codes (N): a symbol other than A/C/G/T matches nothing, so no window containing one is an occurrence.

## Core Documentation Reference

- Source: [GenomicAnalyzer.cs#L195](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/GenomicAnalyzer.cs#L195)

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

- **Time Complexity:** O(n) suffix-tree build + O(m + occ·log occ) per query (positions sorted).
- **Space Complexity:** O(n) for the suffix tree (measured 0.9 s / 171 MB at 1 Mbp, 23 s / 1.7 GB at 10 Mbp).

## See Also

- [find_exact_motif](find_exact_motif.md) — MotifFinder exact search
- [find_known_motifs](find_known_motifs.md) — search a set of motifs
- [kmer_positions](kmer_positions.md) — k-mer occurrence positions
