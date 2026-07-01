# kmer_spectrum

Frequency-of-frequencies (k-mer spectrum) for a sequence.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `kmer_spectrum` |
| **Method ID** | `KmerAnalyzer.GetKmerSpectrum` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Computes the **k-mer spectrum** (frequency-of-frequencies): after counting every
overlapping k-mer, it reports, for each occurrence count `c`, how many *distinct*
k-mers occur exactly `c` times. This is the classic sequencing-coverage histogram
used to distinguish erroneous (low-count) from genomic (high-count) k-mers.
Counting is case-insensitive. When `k` exceeds the sequence length the spectrum is
empty.

## Core Documentation Reference

- Source: [KmerAnalyzer.cs#L136](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/KmerAnalyzer.cs#L136)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | Sequence to analyze (min length 1) |
| `k` | integer | Yes | k-mer length (> 0) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `spectrum` | object | Map from occurrence count → number of distinct k-mers reaching that count |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1003 | k must be positive |

## Examples

### Example 1: Monomer spectrum (GTAGAGCTGT, k=1)

**User Prompt:**
> Give me the k-mer spectrum of "GTAGAGCTGT" for k=1.

**Expected Tool Call:**
```json
{
  "tool": "kmer_spectrum",
  "arguments": { "sequence": "GTAGAGCTGT", "k": 1 }
}
```

**Response:**
```json
{ "spectrum": { "4": 1, "3": 1, "2": 1, "1": 1 } }
```
Monomer counts are G=4, T=3, A=2, C=1 — four distinct counts, one k-mer each.

### Example 2: Repeated trimer (ATGATG, k=3)

**User Prompt:**
> k-mer spectrum of "ATGATG" for k=3.

**Expected Tool Call:**
```json
{
  "tool": "kmer_spectrum",
  "arguments": { "sequence": "ATGATG", "k": 3 }
}
```

**Response:**
```json
{ "spectrum": { "2": 1, "1": 2 } }
```
Trimers: ATG×2, TGA×1, GAT×1 ⇒ one k-mer at count 2, two k-mers at count 1.

## Performance

- **Time Complexity:** O(n) to build the k-mer table.
- **Space Complexity:** O(distinct k-mers).

## See Also

- [count_kmers](count_kmers.md) — raw k-mer → count map
- [analyze_kmers](analyze_kmers.md) — aggregate statistics
- [kmer_frequencies](kmer_frequencies.md) — normalized frequencies
