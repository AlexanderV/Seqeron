# spaced_word_distance

Multiple-pattern spaced-word distance between two sequences (Leimeister et al. 2014).

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `spaced_word_distance` |
| **Method ID** | `KmerAnalyzer.SpacedWordDistance` |
| **Version** | 1.1.0 |
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
`euclidean_counts` (the program takes the Euclidean distance of raw counts, not frequencies). Two further `spaced`
conventions are opt-in (`KmerAnalyzer.SpacedWordDistance(seq1, seq2, patterns, metric, KmerCountingOptions, bothStrands)`):

- `acgtOnly: true` — a window with a non-ACGT symbol at a **match** position gives no word (`spaced` stores every
  non-ACGT letter as N and drops words reading N; don't-care positions are ignored). Frequencies stay
  count ÷ (L − ℓ + 1), the number of windows, as in `spaced`, so they sum to less than 1 when words are dropped.
- `bothStrands: true` — `spaced`'s default mode (no `-r`): the first sequence of the input (`seq1`) is counted on
  its forward strand plus its reverse-complement strand (total 2·W₁), the second (`seq2`) on its forward strand
  only. The value depends on the argument order (a convention of the tool; the paper does not define it).

With both options on, every `spaced` 1.2.0 value checked (3 pairs incl. N/IUPAC/lower case × 3 pattern sets × JS/EU,
with and without `-r`, 36 runs) is reproduced to the 12 printed digits. Pattern-set generation (`spaced`/rasbhari
random optimisation) is not provided; pass the patterns.

## Core Documentation Reference

- Source: [KmerAnalyzer.cs#L1883](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/KmerAnalyzer.cs#L1883)
- Algorithm: [K-mer_Euclidean_Distance.md](../../../algorithms/K-mer/K-mer_Euclidean_Distance.md) §7.5, §7.7

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `seq1` | string | Yes | First sequence (min length 1) |
| `seq2` | string | Yes | Second sequence (min length 1) |
| `patterns` | string[] | Yes | One or more patterns over {0,1}, each starting and ending with `1`, all of the same weight |
| `metric` | string | No | `euclidean` (default), `jensen_shannon` (`js`), `euclidean_counts`, `squared_euclidean_counts`, `manhattan`, `chebyshev`, `canberra`, `cosine`, `d2` |
| `acgtOnly` | boolean | No | Drop words with a non-ACGT symbol at a match position (`spaced` rule). Default false |
| `bothStrands` | boolean | No | `spaced` default mode: `seq1` on both strands vs `seq2` forward (order-dependent). Default false (= `spaced -r`) |

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

### Example 3: `spaced` default mode on sequences with N / IUPAC symbols

**Expected Tool Call:**
```json
{
  "tool": "spaced_word_distance",
  "arguments": {
    "seq1": "AGGTAAGGTGNGTTGAGATctggacTTTTGACGCCTRGAGCCCGCAGTGCTCCTCGAAAAGTAGCNNATGCCTTGGGCTGCT",
    "seq2": "CAAAGGCCCTACCTTCTTATAGTCCTTYCAACATACAAGTAtagttgGAAGTTCTAAGTTCAGNTTAATC",
    "patterns": ["11011", "10111", "11101"],
    "metric": "jensen_shannon",
    "acgtOnly": true,
    "bothStrands": true
  }
}
```

**Response:**
```json
{ "distance": 0.6447798019943918 }
```
`spaced -t 1 -f patterns -d JS` (no `-r`) on the FASTA file with `seq1` first prints 0.644779801994; with `-r`
(`bothStrands: false`) 0.706137367852; swapping the records prints 0.632325490347.

## Performance

- **Time Complexity:** O(m · n · k) for m patterns of weight k. **Space Complexity:** O(distinct spaced words per pattern).

## See Also

- [kmer_distance](kmer_distance.md) — contiguous k-mer word-vector distances
- [kmer_d2_statistics](kmer_d2_statistics.md) — background-adjusted D2* / D2S
