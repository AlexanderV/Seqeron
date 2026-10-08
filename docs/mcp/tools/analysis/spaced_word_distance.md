# spaced_word_distance

Multiple-pattern spaced-word distance between two sequences (Leimeister et al. 2014), and the `spaced -d EV` evolutionary distance (Morgenstern et al. 2015).

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `spaced_word_distance` |
| **Method ID** | `KmerAnalyzer.SpacedWordDistance` |
| **Version** | 1.2.0 |
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

- `acgtOnly: true` — the spaced-faithful mode. Each sequence is first read as `spaced`'s FASTA reader does: every
  character that is not a letter (gap `-`, `*`, digits, blanks) is **deleted**, so `ACG-TACGT` is windowed as
  `ACGTACGT`; letters are upper-cased and every non-ACGT letter becomes N. A window with N at a **match** position
  gives no word (don't-care positions are ignored). Frequencies stay count ÷ (L − ℓ + 1), the number of windows of the
  read sequence, as in `spaced`, so they sum to less than 1 when words are dropped. The default (literal) mode windows
  the string as given.
- `bothStrands: true` — `spaced`'s default mode (no `-r`): the first sequence of the input (`seq1`) is counted on
  its forward strand plus its reverse-complement strand (total 2·W₁), the second (`seq2`) on its forward strand
  only. The value depends on the argument order (a convention of the tool; the paper does not define it).

`metric: "ev"` (alias `evolutionary`) is `spaced -d EV`, the evolutionary distance of Morgenstern, Zhu, Horwege &
Leimeister (2015, Algorithms Mol Biol 10:5), computed exactly as `spaced` 1.2.0 `sort.h` (EV branch). It is not averaged
per pattern: N = Σ over patterns and words of min(c₂(w), c₁(w)) spaced-word matches (seq2 forward vs seq1, both strands
with `bothStrands`), m = min(L₁, L₂) − ℓ + 1, M = max(L₁, L₂) − ℓ + 1, q = Σ_a f₁(a)·f₂(a) (strand-averaged frequencies
with `bothStrands`), V = N/(|P|·m) − s·M·q^w (s = 2 with `bothStrands`, else 1); V ≥ 0 → p = V^(1/w),
d = −¾·ln(4p/3 − 1/3); V < 0 → 1.2 (spaced's saturation value); p < ¼ gives NaN as in `spaced`. Always read in the
spaced-faithful mode; all patterns must have the same length, and both sequences at least ℓ letters.

With `acgtOnly` (and `bothStrands` for the default mode), every `spaced` 1.2.0 value checked is reproduced to the 12
printed digits: WP8's 36 JS/EU runs, and WP9's 248 EV/JS/EU runs on 11 pairs (N, IUPAC, lower case, gaps, `*`,
digits; with and without `-r`; up to 7 pattern sets). Pattern-set generation (`spaced`/rasbhari
random optimisation) is not provided; pass the patterns.

## Core Documentation Reference

- Source: [KmerAnalyzer.cs#L2340](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/KmerAnalyzer.cs#L2340)
- Algorithm: [K-mer_Euclidean_Distance.md](../../../algorithms/K-mer/K-mer_Euclidean_Distance.md) §7.5, §7.7

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `seq1` | string | Yes | First sequence (min length 1) |
| `seq2` | string | Yes | Second sequence (min length 1) |
| `patterns` | string[] | Yes | One or more patterns over {0,1}, each starting and ending with `1`, all of the same weight |
| `metric` | string | No | `euclidean` (default), `jensen_shannon` (`js`), `euclidean_counts`, `squared_euclidean_counts`, `manhattan`, `chebyshev`, `canberra`, `cosine`, `d2`; or `ev` (`evolutionary`) = `spaced -d EV` |
| `acgtOnly` | boolean | No | `spaced` reader and word rule: delete non-letters, drop words with a non-ACGT symbol at a match position. Default false (`ev` always applies it) |
| `bothStrands` | boolean | No | `spaced` default mode: `seq1` on both strands vs `seq2` forward (order-dependent). Default false (= `spaced -r`) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `distance` | number | Mean of the per-pattern metric values, or the `ev` distance (NaN when p < ¼, as `spaced`) |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1002 | At least one pattern is required. |
| 1002 | All patterns must have the same weight (number of '1' positions). |
| 1002 | Pattern must be a non-empty string over {0,1} that starts and ends with '1'. |
| 1002 | D2*/D2S need each sequence's background model; use KmerDistance(string, string, int, metric) or BackgroundAdjustedD2. |
| 1002 | metric must be one of: euclidean, squared_euclidean_counts, manhattan, chebyshev, canberra, cosine, d2, d2star, d2shepherd, jensen_shannon, euclidean_counts, ev |
| 1002 | The spaced EV distance needs all patterns of the same length (spaced keeps one weight and one don't-care count). |
| 1002 | The spaced EV distance needs both sequences to have at least as many letters as the pattern length. |

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

### Example 4: `spaced -d EV` evolutionary distance

```json
{ "tool": "spaced_word_distance", "arguments": {
  "seq1": "AGGTAAGGTGNGTTGAGATctggacTTTTGACGCCTRGAGCCCGCAGTGCTCCTCGAAAAGTAGCNNATGCCTTGGGCTGCT",
  "seq2": "CAAAGGCCCTACCTTCTTATAGTCCTTYCAACATACAAGTAtagttgGAAGTTCTAAGTTCAGNTTAATC",
  "patterns": ["1101011", "1011101", "1110011"], "metric": "ev", "bothStrands": true } }
```
**Response:** `{ "distance": 0.767483449137… }` — `spaced -t 1 -f patterns -d EV` (no `-r`) prints 0.767483449137; with
`-r` it prints `-nan` (p < ¼), returned as NaN.

### Example 5: gapped input in the spaced-faithful mode

```json
{ "tool": "spaced_word_distance", "arguments": { "seq1": "ACG-TACGT", "seq2": "ACGTTACGA",
  "patterns": ["1011", "1101"], "metric": "jensen_shannon", "acgtOnly": true } }
```
**Response:** `{ "distance": 0.57013316426… }` — `spaced -r -d JS` prints 0.57013316426 (the gap is deleted by the reader);
`-d EV` gives 0.427666596474. Keeping the gap as an N window (the pre-WP9 behaviour) gave 0.25.

## Performance

- **Time Complexity:** O(m · n · k) for m patterns of weight k. **Space Complexity:** O(distinct spaced words per pattern).

## See Also

- [kmer_distance](kmer_distance.md) — contiguous k-mer word-vector distances
- [kmer_d2_statistics](kmer_d2_statistics.md) — background-adjusted D2* / D2S
