# kmer_frequencies

Normalized k-mer frequencies for a sequence.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `kmer_frequencies` |
| **Method ID** | `KmerAnalyzer.GetKmerFrequencies` |
| **Version** | 1.1.0 |
| **Stability** | Stable |

## Description

Returns the **normalized k-mer frequency vector**: each distinct k-mer's overlapping
occurrence count divided by the total number of k-mers (`L − k + 1`). Every value is
in `[0, 1]` and the values sum to 1. This is the composition vector used for
alignment-free comparison (see `kmer_distance`). Counting is case-insensitive; when
`k` exceeds the sequence length the result is empty.

Optional modes (`KmerAnalyzer.GetKmerFrequencies(sequence, k, KmerCountingOptions)`): `acgtOnly` skips every window
containing a non-ACGT symbol — a kPAL profile (`kpal/klib.py` splits the sequence on non-ACGT symbols) divided by
its total; `canonical` uses Jellyfish `count -C` counts (implies `acgtOnly`). The denominator is always the number
of counted windows, so the values sum to 1.

## Core Documentation Reference

- Source: [KmerAnalyzer.cs#L660](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/KmerAnalyzer.cs#L660)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | Sequence to analyze (min length 1) |
| `k` | integer | Yes | k-mer length (> 0) |
| `canonical` | boolean | No | Canonical k-mers (jellyfish count -C); implies `acgtOnly`. Default false |
| `acgtOnly` | boolean | No | Skip windows containing a non-ACGT symbol (kPAL / Jellyfish). Default false |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `frequencies` | object | Map from k-mer to its normalized frequency in `[0,1]` |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1003 | k must be positive |

## Examples

### Example 1: Monomer frequencies (ATGC, k=1)

**User Prompt:**
> Base frequencies of "ATGC".

**Expected Tool Call:**
```json
{
  "tool": "kmer_frequencies",
  "arguments": { "sequence": "ATGC", "k": 1 }
}
```

**Response:**
```json
{ "frequencies": { "A": 0.25, "T": 0.25, "G": 0.25, "C": 0.25 } }
```
Each of the 4 bases occurs once out of 4 ⇒ 0.25.

### Example 2: Homopolymer (AAAA, k=2)

**User Prompt:**
> 2-mer frequencies of "AAAA".

**Expected Tool Call:**
```json
{
  "tool": "kmer_frequencies",
  "arguments": { "sequence": "AAAA", "k": 2 }
}
```

**Response:**
```json
{ "frequencies": { "AA": 1.0 } }
```
AA is the only 2-mer (3/3 = 1.0).

### Example 3: Canonical profile of a sequence with N / IUPAC symbols (k=3)

**Expected Tool Call:**
```json
{
  "tool": "kmer_frequencies",
  "arguments": { "sequence": "ACGTNACGTAAcgtRTT", "k": 3, "canonical": true }
}
```

**Response:**
```json
{ "frequencies": { "ACG": 0.6666666666666666, "AAC": 0.1111111111111111, "GTA": 0.1111111111111111, "TAA": 0.1111111111111111 } }
```
`jellyfish count -m 3 -C` + `dump -c`: ACG 6, AAC 1, GTA 1, TAA 1 over the 9 all-ACGT windows.

## Performance

- **Time Complexity:** O(n) to build the k-mer table.
- **Space Complexity:** O(distinct k-mers).

## See Also

- [count_kmers](count_kmers.md) — raw counts
- [kmer_distance](kmer_distance.md) — Euclidean distance over frequency vectors
- [kmer_spectrum](kmer_spectrum.md) — frequency-of-frequencies
