# find_repeats

All repeated substrings of length ≥ `minLength` in a DNA sequence.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `find_repeats` |
| **Method ID** | `GenomicAnalyzer.FindRepeats` |
| **Version** | 1.1.0 |
| **Stability** | Stable |

## Description

Finds **every distinct substring** of length ≥ `minLength` that occurs at least twice
in a DNA sequence, reporting each with all its 0-based occurrence positions (ascending).
Uses a single traversal of the sequence's suffix tree: an internal node with path-label
length `d` and parent depth `p` carries the repeats of lengths `p+1 … d` (CMU 15-451
Lecture #10 §2.1; Gusfield 1997 §7.1). Overlapping occurrences are counted (e.g. `AA` in
`AAAA` occurs at 0, 1, 2). Results are ordered by the lexicographically smallest suffix
that starts with the repeat, then by length. The list includes every prefix and suffix of
a longer repeat (`ACGTACGT` → `A`, `AC`, `ACG`, `ACGT`, `C`, `CG`, … `T`).

With `maximalOnly: true` only the **maximal repeats** are returned (Gusfield 1997 §7.12):
substrings with ≥ 2 occurrences that are *left-diverse* (not all preceded by the same
character; the start of the text counts as a distinct symbol) and *right-diverse* (not all
followed by the same character; the end of the text likewise), i.e. those that cannot be
extended on either side without losing an occurrence (`ACGTACGT` → `ACGT` only). At most
n results; linear time. For maximal repeated *pairs* use `find_direct_repeats`; for
supermaximal repeats `find_supermaximal_repeats`.

## Core Documentation Reference

- Source: [GenomicAnalyzer.cs#L65](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/GenomicAnalyzer.cs#L63)
- Algorithm: `docs/algorithms/Repeat_Analysis/Repeat_Detection.md`

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | DNA sequence (min length 1) |
| `minLength` | integer | Yes | Minimum repeat length (≥ 1) |
| `maximalOnly` | boolean | No | Only Gusfield maximal repeats (default `false`) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `items` | array | Repeats: `{ sequence, positions[], length, count }` |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1001 | Invalid DNA sequence |
| 1002 | Minimum repeat length must be at least 1 |

## Examples

### Example 1: Homopolymer repeats (AAAA, minLength 2)

**User Prompt:**
> Find repeats of length ≥ 2 in "AAAA".

**Expected Tool Call:**
```json
{
  "tool": "find_repeats",
  "arguments": { "sequence": "AAAA", "minLength": 2 }
}
```

**Response:**
```json
{ "items": [ { "sequence": "AA", "positions": [0, 1, 2], "length": 2, "count": 3 }, { "sequence": "AAA", "positions": [0, 1], "length": 3, "count": 2 } ] }
```

### Example 2: Tandem ATG (ATGATG, minLength 3)

**User Prompt:**
> Repeats of length ≥ 3 in "ATGATG".

**Expected Tool Call:**
```json
{
  "tool": "find_repeats",
  "arguments": { "sequence": "ATGATG", "minLength": 3 }
}
```

**Response:**
```json
{ "items": [ { "sequence": "ATG", "positions": [0, 3], "length": 3, "count": 2 } ] }
```
Only ATG recurs; TGA and GAT occur once each.

## Performance

- **Time Complexity:** O(n + output): one suffix-tree traversal plus, per internal node, one occurrence query; the output itself can be quadratic (a homopolymer of length n has n − 1 repeats with up to n positions each). 1 Mbp random DNA with `minLength` 20: ≈ 4.7 s for `FindRepeats` + `FindMaximalRepeats` including the tree build.
- **Space Complexity:** O(n) for the suffix tree (≈ 170 bytes per input character).

## See Also

- [find_tandem_repeats](find_tandem_repeats.md) — consecutive repeat units
- [find_direct_repeats](find_direct_repeats.md) — spaced identical copies
- [find_motif](find_motif.md) — locate a known motif
