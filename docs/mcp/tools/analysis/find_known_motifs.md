# find_known_motifs

Search a DNA sequence for a set of known motifs simultaneously (exact, or IUPAC / both strands / up to k mismatches).

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

`mismatches = k > 0` is fuzznuc `-pmismatch k`: a window is reported when at most k of its positions differ from the motif
(Hamming distance, substitutions only, no gaps; with `iupac` a position agrees when the base is in the code's base set,
without it every motif symbol is literal and a degenerate symbol such as `N` costs one mismatch). The reverse strand is
searched with the reverse-complement motif and the same bound. The result then also carries `matchMismatches` /
`reverseMismatches`: the Hamming distance of each reported position (parallel arrays). A bound at least the motif length
makes every window a hit. Verified against fuzznuc 6.6.0 `-pmismatch 1..3` (positions and `*mismatch` counts) on 6576
sequence/pattern/k rows.

**IUPAC input.** The sequence may contain IUPAC codes (N): such a symbol matches no motif symbol (not even N) and costs a mismatch when `mismatches` > 0.

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
| `mismatches` | integer | No | Maximum mismatching positions per hit, Hamming (default 0 = exact; fuzznuc `-pmismatch`) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `matches` | object | Map from upper-cased motif → ascending 0-based forward-strand positions |
| `reverseMatches` | object | Map from upper-cased motif → ascending forward-strand window starts matching the motif's reverse complement (empty unless `bothStrands`) |
| `matchMismatches` | object | Only with `mismatches` > 0: motif → Hamming distance of each `matches` position (parallel array) |
| `reverseMismatches` | object | Only with `mismatches` > 0: motif → Hamming distance of each `reverseMatches` position |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1001 | Invalid DNA sequence |
| 1001 | Motifs cannot be null |
| 1001 | Invalid IUPAC code in motif (`iupac = true` only) |
| 1001 | Mismatches cannot be negative |

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

### Example 4: Mismatches (EMBOSS fuzznuc `-pmismatch 1`)

```json
{ "tool": "find_known_motifs",
  "arguments": { "sequence": "GAATTCAAGATTCGGATTTTGAATTTT", "motifs": ["GAATTC"], "bothStrands": true, "mismatches": 1 } }
```
**Response:** `{ "matches": { "GAATTC": [0, 20] }, "reverseMatches": { "GAATTC": [0, 20] }, "matchMismatches": { "GAATTC": [0, 1] }, "reverseMismatches": { "GAATTC": [0, 1] } }`
(fuzznuc: `+` 1-6, 21-26 `*mismatch 1`; `-` 1-6, 21-26 `*mismatch 1`).

## Performance

- **Time Complexity:** `mismatches > 0`: one Hamming scan per query, O(n·m) per motif and strand (no index). Exact: fewer than 128 distinct queries (motifs × strands): one KMP scan per query, O(k·n); from 128 queries: O(n) suffix-tree build + O(|m| + occ·log occ) per query.
- **Space Complexity:** O(occ) for scans; O(n) for the suffix tree (heavy: 1.7 GB at 10 Mbp).

## See Also

- [find_motif](find_motif.md) — single-motif search
- [find_regulatory_elements](find_regulatory_elements.md) — built-in regulatory motif set
