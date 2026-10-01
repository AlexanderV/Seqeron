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

The tool returns the whole space as one array, so it refuses results larger than
**1,048,576 k-mers** (= 4^10, `AnalysisTools.MaxGeneratedKmers`; DNA k ≤ 10, protein k ≤ 4). The size
`alphabet.Length^k` is checked with overflow-safe arithmetic before anything is enumerated: without the cap, DNA
k ≥ 16 exceeds the maximum .NET array length and crashes the server. For larger spaces, call the streaming
library method `KmerAnalyzer.GenerateAllKmers` directly.

## Core Documentation Reference

- Source: [KmerAnalyzer.cs#L1941](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/KmerAnalyzer.cs#L1941)

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
| 1003 | alphabet.Length^k = {a}^{k} exceeds the maximum of 1,048,576 k-mers per call |

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
- **Space Complexity:** O(alphabet.Length^k) for the returned array (at most 1,048,576 k-mers); the library enumerator itself streams in O(k).

## See Also

- [count_kmers](count_kmers.md) — observed k-mer counts
- [kmer_frequencies](kmer_frequencies.md) — normalized frequencies
