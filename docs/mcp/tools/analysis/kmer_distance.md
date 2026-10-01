# kmer_distance

Alignment-free Euclidean distance between two sequences' k-mer composition.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `kmer_distance` |
| **Method ID** | `KmerAnalyzer.KmerDistance` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Computes the **alignment-free k-mer distance**: the Euclidean distance between the
two sequences' normalized k-mer frequency vectors, taken over the union of k-mers
occurring in either sequence (a k-mer absent from one sequence contributes a 0
component). Identical sequences yield 0; sharing no k-mers yields the largest
distance. This is the frequency (relative-count) word-composition distance of
Zielezinski et al. (2017) / Vinga & Almeida (2003). Counting is case-insensitive.

The optional `metric` selects another word-vector metric (`KmerAnalyzer.KmerDistance(seq1, seq2, k, metric, markovOrder, bothStrands)`):

| `metric` | Definition | Vector |
|----------|------------|--------|
| `euclidean` (default) | √Σ(f₁−f₂)² | relative frequencies |
| `squared_euclidean_counts` | Σ(c₁−c₂)² — Blaisdell (1986) d_E | raw counts |
| `manhattan` | Σ\|f₁−f₂\| | relative frequencies |
| `chebyshev` | max\|f₁−f₂\| | relative frequencies |
| `canberra` | Σ\|f₁−f₂\|/(f₁+f₂) | relative frequencies |
| `cosine` | 1 − c₁·c₂/(‖c₁‖‖c₂‖) (zero vector → 1) | counts (scale-invariant) |
| `d2` | Σ c₁·c₂ — D2 statistic (a similarity, not a distance) | raw counts |
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
| `seq1` | string | Yes | First sequence (min length 1) |
| `seq2` | string | Yes | Second sequence (min length 1) |
| `k` | integer | Yes | k-mer length (> 0) |
| `metric` | string | No | `euclidean` (default), `squared_euclidean_counts`, `manhattan`, `chebyshev`, `canberra`, `cosine`, `d2`, `d2star`, `d2shepherd` (`d2s`), `jensen_shannon` (`js`), `euclidean_counts` |
| `markovOrder` | integer | No | Background Markov order r, 0 ≤ r < k, for `d2star`/`d2shepherd` (default 0), or −1 = each sequence's order chosen by BIC; must be 0 for the other metrics |
| `bothStrands` | boolean | No | `d2star`/`d2shepherd` only: CAFE `-R` both-strand counts and background (default false) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `distance` | number | Non-negative Euclidean distance between the two frequency vectors |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1002 | metric must be one of: euclidean, squared_euclidean_counts, manhattan, chebyshev, canberra, cosine, d2, d2star, d2shepherd, jensen_shannon, euclidean_counts |
| 1002 | markovOrder applies only to the background-adjusted metrics D2Star and D2Shepherd. |
| 1002 | bothStrands applies only to the background-adjusted metrics D2Star and D2Shepherd; count canonical k-mers (KmerCountingOptions.Canonical) for the other metrics. |
| 1003 | k must be positive |

## Examples

### Example 1: Identical sequences (distance = 0)

**User Prompt:**
> k-mer distance between "ACGTACGT" and itself for k=2.

**Expected Tool Call:**
```json
{
  "tool": "kmer_distance",
  "arguments": { "seq1": "ACGTACGT", "seq2": "ACGTACGT", "k": 2 }
}
```

**Response:**
```json
{ "distance": 0.0 }
```

### Example 2: Disjoint compositions (distance = √2)

**User Prompt:**
> k-mer distance between "AAAA" and "TTTT" for k=1.

**Expected Tool Call:**
```json
{
  "tool": "kmer_distance",
  "arguments": { "seq1": "AAAA", "seq2": "TTTT", "k": 1 }
}
```

**Response:**
```json
{ "distance": 1.4142135623730951 }
```
Frequency vectors {A:1} and {T:1} are orthogonal ⇒ √(1² + 1²) = √2.

### Example 3: Blaisdell squared Euclidean on counts (= 3)

**Expected Tool Call:**
```json
{
  "tool": "kmer_distance",
  "arguments": { "seq1": "ATGTGTG", "seq2": "CATGTG", "k": 3, "metric": "squared_euclidean_counts" }
}
```

**Response:**
```json
{ "distance": 3.0 }
```
Counts (1,0,2,2) vs (1,1,1,1) over {ATG, CAT, GTG, TGT} (Zielezinski et al. 2017 Fig. 1). scipy/alfpy give the
other metrics as manhattan 0.6, chebyshev 0.25, canberra 1.5726495726495728, cosine 0.16666666666666663, d2 5.

### Example 4: Both-strand d2* (CAFE `-R`)

**Expected Tool Call:**
```json
{
  "tool": "kmer_distance",
  "arguments": { "seq1": "AGGTAAGGTGGTTGAGATCTGGACTTTTGACGCCTGGAGCCCGCAGTGCTCCTCGAAAAGTAGCCATGCCTTGGGCTGCT", "seq2": "CAAAGGCCCTACCTTCTTATAGTCCTTTCAACATACAAGTATAGTTGGAAGTTCTAAGTTCAGTTTAATC", "k": 3, "metric": "d2star", "bothStrands": true }
}
```

**Response:**
```json
{ "distance": 0.3838876158581438 }
```
Single strand (`bothStrands` false) gives 0.44457941706964565. Reference: Python replica (sequence MLE background)
whose formula engine reproduces the CAFE `-R` binary (K-mer_Euclidean_Distance.md §7.5).

## Performance

- **Time Complexity:** O(n) to build the two k-mer tables.
- **Space Complexity:** O(distinct k-mers).

## See Also

- [kmer_frequencies](kmer_frequencies.md) — the underlying frequency vectors
- [count_kmers](count_kmers.md) — raw counts
- [kmer_jaccard](kmer_jaccard.md) — k-mer set Jaccard index and Mash distance
