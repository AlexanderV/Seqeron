# most_frequent_kmers

All k-mers tied for the maximum occurrence count.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `most_frequent_kmers` |
| **Method ID** | `KmerAnalyzer.FindMostFrequentKmers` |
| **Version** | 1.1.0 |
| **Stability** | Stable |

## Description

Solves the **Frequent Words Problem** (Rosalind BA1B; Compeau & Pevzner): after
counting every overlapping k-mer, it returns all k-mers whose count equals the
maximum. There may be a single winner or several tied k-mers. Counting is
case-insensitive; order of the returned list is unspecified. When `k` exceeds the
sequence length the result is empty.

Optional Jellyfish modes (`KmerAnalyzer.FindMostFrequentKmers(sequence, k, KmerCountingOptions)`): `acgtOnly`
skips every window containing a non-ACGT symbol; `canonical` keys each k-mer by min(k-mer, reverse complement), so
a k-mer and its reverse complement are one entry whose count pools both strands — the arg-max of
`jellyfish count -C -m k` + `jellyfish dump -c` (implies `acgtOnly`). This is the exact-match reverse-complement-aware
frequent-words question; mismatches (Rosalind BA1J) are the Alignment server's
`frequent_kmers_with_mismatches_and_revcomp`.

## Core Documentation Reference

- Source: [KmerAnalyzer.cs#L612](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/KmerAnalyzer.cs#L612)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | Sequence to analyze (min length 1) |
| `k` | integer | Yes | k-mer length (> 0) |
| `canonical` | boolean | No | Canonical k-mers (jellyfish count -C); implies `acgtOnly`. Default false |
| `acgtOnly` | boolean | No | Skip windows containing a non-ACGT symbol. Default false |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `kmers` | array of string | All k-mers tied for the maximum occurrence count |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1003 | k must be positive |

## Examples

### Example 1: Single winner (ATGATG, k=3)

**User Prompt:**
> Which 3-mers are most frequent in "ATGATG"?

**Expected Tool Call:**
```json
{
  "tool": "most_frequent_kmers",
  "arguments": { "sequence": "ATGATG", "k": 3 }
}
```

**Response:**
```json
{ "kmers": ["ATG"] }
```
ATG occurs twice; TGA and GAT once each.

### Example 2: Tie (Rosalind BA1B sample, k=4)

**User Prompt:**
> Most frequent 4-mers of "ACGTTGCATGTCGCATGATGCATGAGAGCT".

**Expected Tool Call:**
```json
{
  "tool": "most_frequent_kmers",
  "arguments": { "sequence": "ACGTTGCATGTCGCATGATGCATGAGAGCT", "k": 4 }
}
```

**Response:**
```json
{ "kmers": ["CATG", "GCAT"] }
```
Both CATG and GCAT occur 3 times (the documented Rosalind BA1B answer).

### Example 3: Canonical (jellyfish count -C), BA1B sample, k=4

**Expected Tool Call:**
```json
{
  "tool": "most_frequent_kmers",
  "arguments": { "sequence": "ACGTTGCATGTCGCATGATGCATGAGAGCT", "k": 4, "canonical": true }
}
```

**Response:**
```json
{ "kmers": ["ATGC"] }
```
`jellyfish count -m 4 -C` + `dump -c` (Jellyfish 2.3.1): ATGC 4 (ATGC + GCAT), CATG 3 (palindrome), others ≤ 2.

## Performance

- **Time Complexity:** O(n) to build the k-mer table.
- **Space Complexity:** O(distinct k-mers).

## See Also

- [count_kmers](count_kmers.md) — full k-mer → count map
- [kmers_with_min_count](kmers_with_min_count.md) — count-thresholded k-mers
- [find_clumps](find_clumps.md) — locally over-represented k-mers
