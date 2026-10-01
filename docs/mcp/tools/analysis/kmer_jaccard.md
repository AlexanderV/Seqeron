# kmer_jaccard

Exact k-mer Jaccard similarity of two sequences and the Mash distance derived from it.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `kmer_jaccard` |
| **Method ID** | `KmerAnalyzer.JaccardSimilarity` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Delegates to `KmerAnalyzer.JaccardSimilarity` and `KmerAnalyzer.MashDistance` (same options). Decomposes each sequence into its **set** of distinct k-mers and returns the Jaccard index
J = |A ∩ B| / |A ∪ B| (Jaccard 1901/1912) as a fraction in [0, 1], plus the Mash distance
D = −(1/k)·ln(2J/(1+J)) (Ondov et al. 2016, Genome Biol 17:132, eq. 4) with Mash's boundary rules
(identical sets → 0, nothing shared → 1, capped at 1). The value is exact (no MinHash sketch): it equals
`mash dist` whenever the sketch size is at least the union size.

- Default: literal k-mers, case-insensitive (every symbol forms k-mers).
- `canonical = true`: k-mers keyed by min(w, revcomp(w)) and windows containing a non-ACGT base skipped —
  the k-mers of `mash dist` (default) and sourmash DNA MinHash.
- `acgtOnly = true`: non-ACGT windows skipped, strands kept apart (`mash sketch -n`).
- Both k-mer sets empty (sequences shorter than k): `jaccard` = 0 (undefined 0/0), `mashDistance` = 0 (Mash).

## Core Documentation Reference

- Source: [KmerAnalyzer.cs#L1136](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/KmerAnalyzer.cs#L1136)
- Algorithm: [K-mer_Euclidean_Distance.md](../../../algorithms/K-mer/K-mer_Euclidean_Distance.md) §2.7

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `seq1` | string | Yes | First sequence (min length 1) |
| `seq2` | string | Yes | Second sequence (min length 1) |
| `k` | integer | Yes | k-mer length (> 0) |
| `canonical` | boolean | No | Canonical (strand-collapsed) k-mers, implies `acgtOnly` (default false) |
| `acgtOnly` | boolean | No | Skip k-mers containing a non-ACGT symbol (default false) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `jaccard` | number | Jaccard index in [0, 1] |
| `mashDistance` | number | Mash distance in [0, 1] |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1003 | k must be positive |

## Examples

### Example 1: Mash k-mers (canonical)

**User Prompt:**
> Mash distance between ACGTTGCAACGGT and ACGTAGCATCGGTA with k=2.

**Expected Tool Call:**
```json
{
  "tool": "kmer_jaccard",
  "arguments": { "seq1": "ACGTTGCAACGGT", "seq2": "ACGTAGCATCGGTA", "k": 2, "canonical": true }
}
```

**Response:**
```json
{ "jaccard": 0.5, "mashDistance": 0.20273255405408222 }
```
Mash 2.3 `mash dist -k 2 -s 100000`: 0.202733, shared-hashes 5/10; sourmash `MinHash(scaled=1)` Jaccard 0.5.

### Example 2: Literal k-mers

```json
{ "tool": "kmer_jaccard", "arguments": { "seq1": "ATGTGTG", "seq2": "CATGTG", "k": 3 } }
```
**Response:** `{ "jaccard": 0.75, "mashDistance": 0.05138355994241945 }` — sets {ATG, GTG, TGT} and {ATG, CAT, GTG, TGT}.

## Performance

- **Time Complexity:** O(n·k) to build the two k-mer sets.
- **Space Complexity:** O(distinct k-mers).

## See Also

- [kmer_distance](kmer_distance.md) — count/frequency vector metrics
- [count_kmers](count_kmers.md) — the k-mer counting modes
