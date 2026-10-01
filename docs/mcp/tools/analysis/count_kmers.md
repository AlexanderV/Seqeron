# count_kmers

Count every k-mer occurrence in a sequence.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `count_kmers` |
| **Method ID** | `KmerAnalyzer.CountKmers` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Counts every overlapping k-mer (substring of length `k`) in a sequence and returns a
map of k-mer → occurrence count. Counting is case-insensitive. The counts sum to
`L − k + 1` (Wikipedia — K-mer). When `k` exceeds the sequence length the map is empty.

By default every symbol (including `N` and other IUPAC codes) is counted literally. Two optional
Jellyfish modes delegate to `KmerAnalyzer.CountKmers(sequence, k, KmerCountingOptions)`:

- `acgtOnly = true` — skip every window that contains a symbol other than A/C/G/T after case
  folding (Jellyfish `mer_iterator` resets its window on such a base). Counts sum to the number of
  all-ACGT windows.
- `canonical = true` — key each k-mer by the lexicographically smaller of itself and its reverse
  complement (`jellyfish count -C`, Marçais & Kingsford 2011). Implies `acgtOnly`, because the
  canonical form is defined only over ACGT. Cross-checked against Jellyfish 2.3.1 `count -C` + `dump -c`.

## Core Documentation Reference

- Source: [KmerAnalyzer.cs#L20](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/KmerAnalyzer.cs#L20)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | Sequence to analyze (min length 1) |
| `k` | integer | Yes | k-mer length (> 0) |
| `canonical` | boolean | No | Canonical counting (`jellyfish count -C`); implies `acgtOnly`. Default `false` |
| `acgtOnly` | boolean | No | Skip windows containing a non-ACGT symbol (Jellyfish convention). Default `false` |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `counts` | object | Map of k-mer → occurrence count |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1003 | k must be positive |

## Examples

### Example 1: Homopolymer

**User Prompt:**
> Count 2-mers in "AAAA".

**Expected Tool Call:**
```json
{ "tool": "count_kmers", "arguments": { "sequence": "AAAA", "k": 2 } }
```

**Response:**
```json
{ "counts": { "AA": 3 } }
```
Overlapping windows AA, AA, AA ⇒ AA=3 (= L−k+1 = 3).

### Example 2: Mixed 3-mers

**Input:** `{ "sequence": "ATGATG", "k": 3 }`
→ **Response:** `{ "counts": { "ATG": 2, "TGA": 1, "GAT": 1 } }`

### Example 3: Canonical counting (Jellyfish -C)

**Input:** `{ "sequence": "ATGATG", "k": 3, "canonical": true }`
→ **Response:** `{ "counts": { "ATG": 2, "ATC": 1, "TCA": 1 } }`
(GAT → ATC, TGA → TCA; equals `jellyfish count -m 3 -C` + `dump -c`.)

### Example 4: ACGT-only (N resets the window)

**Input:** `{ "sequence": "ACGTNACGT", "k": 4, "acgtOnly": true }`
→ **Response:** `{ "counts": { "ACGT": 2 } }` (literal mode would also return CGTN, GTNA, TNAC, NACG).

## Performance

- **Time Complexity:** O(L·k) (string keys; canonical adds O(D·k) for the fold over D distinct k-mers). **Space Complexity:** O(distinct k-mers).

## See Also

- [analyze_kmers](analyze_kmers.md)
- [kmer_frequencies](kmer_frequencies.md)
- [most_frequent_kmers](most_frequent_kmers.md)
