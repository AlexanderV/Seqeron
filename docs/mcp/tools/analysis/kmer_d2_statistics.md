# kmer_d2_statistics

Background-adjusted word-match statistics D2* and D2S of two DNA sequences, with the Markov orders and BIC values behind them.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `kmer_d2_statistics` |
| **Method ID** | `KmerAnalyzer.BackgroundAdjustedD2` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Returns the full result of `KmerAnalyzer.BackgroundAdjustedD2` (plus `KmerAnalyzer.MarkovOrderBics` over orders 0 … `AutoMarkovOrderLimit(k)`), which
[`kmer_distance`](kmer_distance.md) reduces to one dissimilarity:

- **D2*** = Σ X̃Ỹ/√(E_X E_Y) and **D2S** = Σ X̃Ỹ/√(X̃² + Ỹ²) (Reinert, Chew, Sun & Waterman 2009; Wan et al. 2010),
  over all 4^k words, with background-centred counts X̃ = X − E_X, Ỹ = Y − E_Y;
- the dissimilarities **d2*** = ½(1 − D2*/√(Σ X̃²/E_X · Σ Ỹ²/E_Y)) and
  **d2S** = ½(1 − D2S/√(Σ X̃²/√(X̃²+Ỹ²) · Σ Ỹ²/√(X̃²+Ỹ²))) in [0, 1] (Song et al. 2014; CAFE `D2star` / `D2shepp`);
- the background Markov orders used for each sequence (`markovOrder`, or the BIC choice when `markovOrder` = −1);
- each sequence's BIC(r) = −2 ln L̂_r + 3·4^r·ln N_r for r = 0 … min(k − 1, 10) (Schwarz 1978; the criterion of
  `markovOrder` = −1, CAFE `-M -1`); index = order.

E_X(w) = n̄·p̂_X(w) under the order-r Markov chain fitted to that sequence by maximum likelihood. Windows with a
non-ACGT symbol are skipped (case-insensitive); k ≤ 12. A dissimilarity is NaN when its normaliser is 0.
`bothStrands = true` is CAFE's `-R` mode: counts X(w) + X(RC(w)), background ½(p̂(w) + p̂(RC(w))) (expected count
n̄·(p̂(w) + p̂(RC(w)))), all 4^k words; the chain and the BIC are computed on the given strand.

## Core Documentation Reference

- Source: [KmerAnalyzer.cs#L898](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/KmerAnalyzer.cs#L898)
- Algorithm: [K-mer_Euclidean_Distance.md](../../../algorithms/K-mer/K-mer_Euclidean_Distance.md) §2.9, §7.4, §7.5

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `seq1` | string | Yes | First DNA sequence (min length 1; needs at least one ACGT k-mer window) |
| `seq2` | string | Yes | Second DNA sequence |
| `k` | integer | Yes | Word length, 1 ≤ k ≤ 12 |
| `markovOrder` | integer | No | Background Markov order r, 0 ≤ r < k (default 0), or −1 = each sequence's order chosen by BIC |
| `bothStrands` | boolean | No | CAFE `-R` both-strand mode (default false) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `d2Star` | number | Raw D2* statistic |
| `d2Shepherd` | number | Raw D2S statistic |
| `d2StarDistance` | number | d2* dissimilarity in [0, 1] (NaN when undefined) |
| `d2ShepherdDistance` | number | d2S dissimilarity in [0, 1] (NaN when undefined) |
| `markovOrder1` | integer | Background order used for `seq1` |
| `markovOrder2` | integer | Background order used for `seq2` |
| `bic1` | number[] | BIC of `seq1` for orders 0 … min(k − 1, 10) |
| `bic2` | number[] | BIC of `seq2` for orders 0 … min(k − 1, 10) |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1002 | Sequence has no k-mer window over A/C/G/T. |
| 1003 | K must be in [1, 12]. |
| 1003 | Markov order must be in [0, k) or -1 (BIC). |

## Examples

### Example 1: Order-0 background, single strand

**Expected Tool Call:**
```json
{
  "tool": "kmer_d2_statistics",
  "arguments": {
    "seq1": "AGGTAAGGTGGTTGAGATCTGGACTTTTGACGCCTGGAGCCCGCAGTGCTCCTCGAAAAGTAGCCATGCCTTGGGCTGCT",
    "seq2": "CAAAGGCCCTACCTTCTTATAGTCCTTTCAACATACAAGTATAGTTGGAAGTTCTAAGTTCAGTTTAATC",
    "k": 3
  }
}
```

**Response:**
```json
{
  "d2Star": 7.136981012975184, "d2Shepherd": -0.6884013403313263,
  "d2StarDistance": 0.44457941706964565, "d2ShepherdDistance": 0.5084031847096179,
  "markovOrder1": 0, "markovOrder2": 0,
  "bic1": [231.88984329117494, 255.74056400499072, 370.90681674887685],
  "bic2": [199.825085479026, 217.12856241056033, 322.30742386881116]
}
```

### Example 2: Both strands (CAFE `-R`), order 1

**Expected Tool Call:**
```json
{
  "tool": "kmer_d2_statistics",
  "arguments": {
    "seq1": "AGGTAAGGTGGTTGAGATCTGGACTTTTGACGCCTGGAGCCCGCAGTGCTCCTCGAAAAGTAGCCATGCCTTGGGCTGCT",
    "seq2": "CAAAGGCCCTACCTTCTTATAGTCCTTTCAACATACAAGTATAGTTGGAAGTTCTAAGTTCAGTTTAATC",
    "k": 3, "markovOrder": 1, "bothStrands": true
  }
}
```

**Response (statistics):**
```json
{ "d2Star": -5.014593950349383, "d2Shepherd": -1.3145456673597726, "d2StarDistance": 0.5716659331844273, "d2ShepherdDistance": 0.5152519131905425, "markovOrder1": 1, "markovOrder2": 1 }
```
Reference values: Python replica of the formulas with a per-sequence maximum-likelihood background (its formula
engine reproduces the CAFE binary, single strand and `-R`, to 6 digits); BIC replica to 1e-9.

## Performance

- **Time Complexity:** Θ(4^k · k) for the statistics, O(n) per BIC order. **Space Complexity:** O(distinct k-mers + distinct (r+1)-mers) (sparse Markov tables above r + 1 = 8).

## See Also

- [kmer_distance](kmer_distance.md) — `d2star` / `d2shepherd` as a single distance
- [spaced_word_distance](spaced_word_distance.md) — spaced-word frequency distance
