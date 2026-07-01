# generate_all_kmers

Enumerate the entire k-mer space for an alphabet.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `generate_all_kmers` |
| **Method ID** | `KmerAnalyzer.GenerateAllKmers` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Generates **all possible k-mers** over a given alphabet: the k-fold Cartesian product.
The number produced is `alphabet.Length^k` (4^k for the default DNA alphabet). When
the alphabet is sorted the k-mers are emitted in lexicographic order with the
rightmost position advancing fastest (odometer ordering), so the default `"ACGT"`
yields `AAA, AAC, …, TTT`.

## Core Documentation Reference

- Source: [KmerAnalyzer.cs#L299](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/KmerAnalyzer.cs#L299)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `k` | integer | Yes | k-mer length (> 0) |
| `alphabet` | string | No | Alphabet (default `"ACGT"`, non-empty) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `kmers` | array of string | All `alphabet.Length^k` distinct k-mers |

## Errors

| Code | Message |
|------|---------|
| 1003 | k must be positive |
| 1001 | Alphabet cannot be null or empty |

## Examples

### Example 1: DNA monomers (k=1)

**User Prompt:**
> List all 1-mers of the DNA alphabet.

**Expected Tool Call:**
```json
{
  "tool": "generate_all_kmers",
  "arguments": { "k": 1 }
}
```

**Response:**
```json
{ "kmers": ["A", "C", "G", "T"] }
```

### Example 2: Binary alphabet 2-mers (k=2, alphabet "AT")

**User Prompt:**
> All 2-mers over the alphabet "AT".

**Expected Tool Call:**
```json
{
  "tool": "generate_all_kmers",
  "arguments": { "k": 2, "alphabet": "AT" }
}
```

**Response:**
```json
{ "kmers": ["AA", "AT", "TA", "TT"] }
```
2² = 4 k-mers in odometer order.

## Performance

- **Time Complexity:** O(alphabet.Length^k) — exponential in k.
- **Space Complexity:** O(k) recursion depth (streamed).

## See Also

- [count_kmers](count_kmers.md) — observed k-mer counts
- [kmer_frequencies](kmer_frequencies.md) — normalized frequencies
