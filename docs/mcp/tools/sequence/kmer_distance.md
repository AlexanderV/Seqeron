# kmer_distance

Calculate k-mer based distance between two sequences.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Sequence |
| **Tool Name** | `kmer_distance` |
| **Method ID** | `KmerAnalyzer.KmerDistance` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Computes the k-mer distance between two sequences using Euclidean distance of k-mer frequencies. This alignment-free method is useful for comparing sequences of different lengths and provides a quick measure of sequence similarity. Lower values indicate more similar sequences; identical sequences have a distance of 0.

The optional `metric` / `markovOrder` parameters are the same as on the Analysis server's [`kmer_distance`](../analysis/kmer_distance.md) and call the same library overload `KmerAnalyzer.KmerDistance(seq1, seq2, k, metric, markovOrder, bothStrands)`:

| `metric` | Definition | Vector |
|----------|------------|--------|
| `euclidean` (default) | √Σ(f₁−f₂)² | relative frequencies |
| `squared_euclidean_counts` | Σ(c₁−c₂)² — Blaisdell (1986) d_E | raw counts |
| `manhattan` / `chebyshev` / `canberra` | L1 / L∞ / Canberra | relative frequencies |
| `cosine` | 1 − c₁·c₂/(‖c₁‖‖c₂‖) | counts |
| `d2` | Σ c₁·c₂ (a similarity) | raw counts |
| `d2star` | d2* = ½(1 − D2*/√(Σ X̃²/E_X · Σ Ỹ²/E_Y)), D2* = Σ X̃Ỹ/√(E_X E_Y) — Reinert et al. (2009), Song et al. (2014) | background-centred counts X̃ = X − E_X (ACGT words) |
| `d2shepherd` (alias `d2s`) | d2S = ½(1 − D2S/√(Σ X̃²/√(X̃²+Ỹ²) · Σ Ỹ²/√(X̃²+Ỹ²))), D2S = Σ X̃Ỹ/√(X̃²+Ỹ²) | background-centred counts |
| `jensen_shannon` (alias `js`) | ½Σ f₁ log₂(f₁/m) + ½Σ f₂ log₂(f₂/m), m = ½(f₁+f₂) — Jensen–Shannon divergence (Lin 1991), = scipy `jensenshannon(p, q, base=2)²` | relative frequencies |
| `euclidean_counts` | √Σ(c₁−c₂)² (the `spaced` program's `-d EU` value) | raw counts |

For `d2star` / `d2shepherd` the expected counts E_X = n̄·p̂_X(w) come from an order-`markovOrder` Markov chain
fitted to each sequence (maximum likelihood on its ACGT r- and (r+1)-mer counts; order 0 = letter frequencies),
the sums run over all 4^k words, windows with a non-ACGT symbol are skipped, and k ≤ 12
(`KmerAnalyzer.BackgroundAdjustedD2`). Formulas as CAFE (Lu et al. 2017) `D2star` / `D2shepp`.

`bothStrands = true` (d2star / d2shepherd only) is CAFE's `-R` mode: each word's count is X(w) + X(RC(w)) and the
background probability is ½(p̂(w) + p̂(RC(w))), so E = n̄·(p̂(w) + p̂(RC(w))), summed over all 4^k words; the Markov
chain (and the BIC order) is fitted to the given strand. Cross-checked: a replica that reproduces the CAFE `-R`
binary on 20 runs to 6 digits gives the Example 4 pair at k=3, r=0 d2* = 0.3838876158581438 here. The raw statistics,
orders and BIC values are returned by [`kmer_d2_statistics`](../analysis/kmer_d2_statistics.md).

## Core Documentation Reference

- Source: [KmerAnalyzer.cs#L772](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/KmerAnalyzer.cs#L772)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence1` | string | Yes | First sequence (min length: 1) |
| `sequence2` | string | Yes | Second sequence (min length: 1) |
| `k` | integer | No | K-mer length (default: 3, minimum: 1) |
| `metric` | string | No | `euclidean` (default), `squared_euclidean_counts`, `manhattan`, `chebyshev`, `canberra`, `cosine`, `d2`, `d2star`, `d2shepherd` (`d2s`), `jensen_shannon` (`js`), `euclidean_counts` |
| `markovOrder` | integer | No | Background Markov order r, 0 ≤ r < k, for `d2star`/`d2shepherd` (default 0), or −1 = each sequence's order chosen by BIC; must be 0 for the other metrics |
| `bothStrands` | boolean | No | `d2star`/`d2shepherd` only: CAFE `-R` both-strand counts and background (default false) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `distance` | number | K-mer distance (0 = identical k-mer composition) |
| `k` | integer | K-mer length used |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1003 | K must be at least 1 |
| 1002 | metric must be one of: euclidean, squared_euclidean_counts, manhattan, chebyshev, canberra, cosine, d2, d2star, d2shepherd, jensen_shannon, euclidean_counts |
| 1002 | markovOrder applies only to the background-adjusted metrics D2Star and D2Shepherd. |
| 1002 | bothStrands applies only to the background-adjusted metrics D2Star and D2Shepherd; count canonical k-mers (KmerCountingOptions.Canonical) for the other metrics. |

## Examples

### Example 1: Compare identical sequences

**User Prompt:**
> What's the k-mer distance between "ATGCATGC" and itself?

**Expected Tool Call:**
```json
{
  "tool": "kmer_distance",
  "arguments": {
    "sequence1": "ATGCATGC",
    "sequence2": "ATGCATGC",
    "k": 3
  }
}
```

**Response:**
```json
{
  "distance": 0,
  "k": 3
}
```

### Example 2: Compare different sequences

**User Prompt:**
> Compare "ATGCATGC" and "CCCCCCCC" using k-mer distance

**Expected Tool Call:**
```json
{
  "tool": "kmer_distance",
  "arguments": {
    "sequence1": "ATGCATGC",
    "sequence2": "CCCCCCCC",
    "k": 3
  }
}
```

**Response:**
```json
{
  "distance": 1.130388330520878,
  "k": 3
}
```

Euclidean distance of the 3-mer frequency vectors: ATGCATGC → ATG 2/6, TGC 2/6, GCA 1/6, CAT 1/6; CCCCCCCC → CCC 6/6; √(2·(1/3)² + 2·(1/6)² + 1²) = 1.130388330520878 (`scipy.spatial.distance.euclidean`, identical).

### Example 4: Both-strand d2* (CAFE `-R`)

**Expected Tool Call:**
```json
{
  "tool": "kmer_distance",
  "arguments": { "sequence1": "AGGTAAGGTGGTTGAGATCTGGACTTTTGACGCCTGGAGCCCGCAGTGCTCCTCGAAAAGTAGCCATGCCTTGGGCTGCT", "sequence2": "CAAAGGCCCTACCTTCTTATAGTCCTTTCAACATACAAGTATAGTTGGAAGTTCTAAGTTCAGTTTAATC", "k": 3, "metric": "d2star", "bothStrands": true }
}
```

**Response:**
```json
{ "distance": 0.3838876158581438, "k": 3 }
```
Single strand (`bothStrands` false) gives 0.44457941706964565. Reference: Python replica (sequence MLE background)
whose formula engine reproduces the CAFE `-R` binary (K-mer_Euclidean_Distance.md §7.5).

## Performance

- **Time Complexity:** O(n + m) where n, m are sequence lengths
- **Space Complexity:** O(4^k) for k-mer frequency storage

## See Also

- [kmer_count](kmer_count.md) - Count k-mer frequencies
- [kmer_analyze](kmer_analyze.md) - Comprehensive k-mer statistics
- [edit_distance](../core/edit_distance.md) - Levenshtein distance
