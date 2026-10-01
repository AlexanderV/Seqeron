# spaced_word_distance

Multiple-pattern spaced-word distance between two sequences (Leimeister et al. 2014).

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `spaced_word_distance` |
| **Method ID** | `KmerAnalyzer.SpacedWordDistance` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

A pattern P ∈ {0,1}^ℓ (starting and ending with `1`) selects the symbols at its `1` (match) positions of each
length-ℓ window; the selected k symbols (k = weight = number of `1`s) form the window's **spaced word**
(`KmerAnalyzer.CountSpacedWords`, one word per window, case-insensitive, literal symbols). For a set of patterns
of equal weight the distance is the **average of the per-pattern distances** of the spaced-word vectors
(Leimeister, Boden, Horwege, Lindner & Morgenstern 2014, Bioinformatics 30:1991):

d_P(S₁, S₂) = (1/m) Σᵢ d(N_{Pᵢ}(S₁), N_{Pᵢ}(S₂)).

The per-pattern metric is any word-vector metric of [`kmer_distance`](kmer_distance.md) except the
background-adjusted ones: `euclidean` (default; relative frequencies, the paper's Euclidean distance),
`jensen_shannon` (JS divergence base 2 of the relative frequencies, the paper's JS distance), `euclidean_counts`,
`squared_euclidean_counts`, `manhattan`, `chebyshev`, `canberra`, `cosine`, `d2`. The all-`1` pattern of length k
gives the contiguous k-mer distance.

Reference program `spaced` 1.2.0 (`-r -f patterns`): `-d JS` equals `jensen_shannon`; `-d EU` equals
`euclidean_counts` (the program takes the Euclidean distance of raw counts, not frequencies). Conventions not
reproduced: `spaced` drops words with a non-ACGT symbol at a match position (keeping their windows in the JS
denominator), and its default both-strand mode (no `-r`) compares forward words of one sequence with both strands
of the other. Pattern-set generation (`spaced`/rasbhari random optimisation) is not provided; pass the patterns.

## Core Documentation Reference

- Source: [KmerAnalyzer.cs#L1731](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/KmerAnalyzer.cs#L1731)
- Algorithm: [K-mer_Euclidean_Distance.md](../../../algorithms/K-mer/K-mer_Euclidean_Distance.md) §7.5

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `seq1` | string | Yes | First sequence (min length 1) |
| `seq2` | string | Yes | Second sequence (min length 1) |
| `patterns` | string[] | Yes | One or more patterns over {0,1}, each starting and ending with `1`, all of the same weight |
| `metric` | string | No | `euclidean` (default), `jensen_shannon` (`js`), `euclidean_counts`, `squared_euclidean_counts`, `manhattan`, `chebyshev`, `canberra`, `cosine`, `d2` |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `distance` | number | Mean of the per-pattern metric values |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1002 | At least one pattern is required. |
| 1002 | All patterns must have the same weight (number of '1' positions). |
| 1002 | Pattern must be a non-empty string over {0,1} that starts and ends with '1'. |
| 1002 | D2*/D2S need each sequence's background model; use KmerDistance(string, string, int, metric) or BackgroundAdjustedD2. |
| 1002 | metric must be one of: euclidean, squared_euclidean_counts, manhattan, chebyshev, canberra, cosine, d2, d2star, d2shepherd, jensen_shannon, euclidean_counts |

## Examples

### Example 1: Jensen–Shannon over three weight-4 patterns (= spaced -d JS)

**Expected Tool Call:**
```json
{
  "tool": "spaced_word_distance",
  "arguments": {
    "seq1": "AGGTAAGGTGGTTGAGATCTGGACTTTTGACGCCTGGAGCCCGCAGTGCTCCTCGAAAAGTAGCCATGCCTTGGGCTGCT",
    "seq2": "CAAAGGCCCTACCTTCTTATAGTCCTTTCAACATACAAGTATAGTTGGAAGTTCTAAGTTCAGTTTAATC",
    "patterns": ["11011", "10111", "11101"],
    "metric": "jensen_shannon"
  }
}
```

**Response:**
```json
{ "distance": 0.8163228541607376 }
```
`spaced -r -d JS -f patterns` prints 0.816322854161; scipy `jensenshannon(base=2)²` averaged over the patterns gives the same value.

### Example 2: Paper's Euclidean distance on relative frequencies (default metric)

**Expected Tool Call:**
```json
{
  "tool": "spaced_word_distance",
  "arguments": {
    "seq1": "AGGTAAGGTGGTTGAGATCTGGACTTTTGACGCCTGGAGCCCGCAGTGCTCCTCGAAAAGTAGCCATGCCTTGGGCTGCT",
    "seq2": "CAAAGGCCCTACCTTCTTATAGTCCTTTCAACATACAAGTATAGTTGGAAGTTCTAAGTTCAGTTTAATC",
    "patterns": ["11011", "10111", "11101"]
  }
}
```

**Response:**
```json
{ "distance": 0.17567404832368613 }
```
With `"metric": "euclidean_counts"` the value is 12.40897581662776 (`spaced -r -d EU` prints 12.4089758166).

## Performance

- **Time Complexity:** O(m · n · k) for m patterns of weight k. **Space Complexity:** O(distinct spaced words per pattern).

## See Also

- [kmer_distance](kmer_distance.md) — contiguous k-mer word-vector distances
- [kmer_d2_statistics](kmer_d2_statistics.md) — background-adjusted D2* / D2S
