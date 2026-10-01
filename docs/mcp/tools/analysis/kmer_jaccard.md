# kmer_jaccard

k-mer Jaccard similarity of two sequences (exact, or from Mash MinHash sketches), the Mash distance derived from it, the Mash p-value in sketch mode, and the exact containment indices.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `kmer_jaccard` |
| **Method ID** | `KmerAnalyzer.JaccardSimilarity` |
| **Version** | 1.1.0 |
| **Stability** | Stable |

## Description

Delegates to `KmerAnalyzer.JaccardSimilarity`, `KmerAnalyzer.MashDistance` and `KmerAnalyzer.ContainmentIndex` (same options); with `sketchSize > 0` to `KmerAnalyzer.CreateMinHashSketch` + `KmerAnalyzer.CompareMinHashSketches`. Decomposes each sequence into its **set** of distinct k-mers and returns the Jaccard index
J = |A ∩ B| / |A ∪ B| (Jaccard 1901/1912) as a fraction in [0, 1], plus the Mash distance
D = −(1/k)·ln(2J/(1+J)) (Ondov et al. 2016, Genome Biol 17:132, eq. 4) with Mash's boundary rules
(identical sets → 0, nothing shared → 1, capped at 1). The value is exact (no MinHash sketch): it equals
`mash dist` whenever the sketch size is at least the union size.

- Default: literal k-mers, case-insensitive (every symbol forms k-mers).
- `canonical = true`: k-mers keyed by min(w, revcomp(w)) and windows containing a non-ACGT base skipped —
  the k-mers of `mash dist` (default) and sourmash DNA MinHash.
- `acgtOnly = true`: non-ACGT windows skipped, strands kept apart (`mash sketch -n`).
- Both k-mer sets empty (sequences shorter than k): `jaccard` = 0 (undefined 0/0), `mashDistance` = 0 (Mash).
- Containment (always exact): `containmentSeq1InSeq2` = |A ∩ B| / |A|, `containmentSeq2InSeq1` = |A ∩ B| / |B|
  (Koslicki & Zabeti 2019; sourmash `compare --containment`, exact for `scaled=1`); 0 when the first set is empty.
- `sketchSize = s > 0`: `jaccard` and `mashDistance` are estimated as `mash dist -k k -s s` does (Ondov et al. 2016;
  Mash 2.3 `Sketch.cpp`/`CommandDistance.cpp`): bottom-s sketch of the MurmurHash3_x64_128 (seed 42) hashes of the
  k-mers (64-bit h1 for k ≥ 17, its low 32 bits for k ≤ 16), merged until s union hashes; `sharedHashes`/`sketchDenominator`
  are Mash's "x/s" column and `pValue` the Mash binomial p-value (sequence lengths = the input lengths). Non-ACGT
  windows are always skipped; `canonical = false` gives `mash dist -n`. k must be ≤ 32.

## Core Documentation Reference

- Source: [KmerAnalyzer.cs#L1293](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/KmerAnalyzer.cs#L1293)
- Algorithm: [K-mer_Euclidean_Distance.md](../../../algorithms/K-mer/K-mer_Euclidean_Distance.md) §2.7, §2.10

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `seq1` | string | Yes | First sequence (min length 1) |
| `seq2` | string | Yes | Second sequence (min length 1) |
| `k` | integer | Yes | k-mer length (> 0) |
| `canonical` | boolean | No | Canonical (strand-collapsed) k-mers, implies `acgtOnly` (default false) |
| `acgtOnly` | boolean | No | Skip k-mers containing a non-ACGT symbol (default false) |
| `sketchSize` | integer | No | 0 (default) = exact sets; s > 0 = Mash MinHash sketch size (Mash default 1000) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `jaccard` | number | Jaccard index in [0, 1] |
| `mashDistance` | number | Mash distance in [0, 1] |
| `containmentSeq1InSeq2` | number | Exact containment |A ∩ B| / |A| |
| `containmentSeq2InSeq1` | number | Exact containment |A ∩ B| / |B| |
| `sharedHashes` | integer \| null | Sketch mode: shared hashes x (Mash "x/s"); null when exact |
| `sketchDenominator` | integer \| null | Sketch mode: union hashes compared (s unless both sketches are smaller); null when exact |
| `pValue` | number \| null | Sketch mode: Mash p-value; null when exact |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1003 | k must be positive |
| 1003 | sketchSize must be >= 0 |
| 1003 | k must be <= 32 when sketchSize > 0 |

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
{ "jaccard": 0.5, "mashDistance": 0.20273255405408222, "containmentSeq1InSeq2": 0.8333333333333334, "containmentSeq2InSeq1": 0.5555555555555556, "sharedHashes": null, "sketchDenominator": null, "pValue": null }
```
Mash 2.3 `mash dist -k 2 -s 100000`: 0.202733, shared-hashes 5/10; sourmash `MinHash(scaled=1)` Jaccard 0.5,
`contained_by` 5/6 and 5/9.

### Example 2: Literal k-mers

```json
{ "tool": "kmer_jaccard", "arguments": { "seq1": "ATGTGTG", "seq2": "CATGTG", "k": 3 } }
```
**Response:** `{ "jaccard": 0.75, "mashDistance": 0.05138355994241945, "containmentSeq1InSeq2": 1, "containmentSeq2InSeq1": 0.75, "sharedHashes": null, "sketchDenominator": null, "pValue": null }` — sets {ATG, GTG, TGT} and {ATG, CAT, GTG, TGT}.

### Example 3: Mash sketch (x/s and p-value)

```json
{ "tool": "kmer_jaccard", "arguments": {
  "seq1": "AGGTAAGGTGGTTGAGATCTGGACTTTTGACGCCTGGAGCCCGCAGTGCTCCTCGAAAAGTAGCCATGCCTTGGGCTGCT",
  "seq2": "CAAAGGCCCTACCTTCTTATAGTCCTTTCAACATACAAGTATAGTTGGAAGTTCTAAGTTCAGTTTAATC",
  "k": 4, "canonical": true, "sketchSize": 10 } }
```
**Response:** `jaccard` 0.4, `mashDistance` 0.1399036…, `sharedHashes` 4, `sketchDenominator` 10, `pValue` 0.02917…,
`containmentSeq1InSeq2` 0.35, `containmentSeq2InSeq1` 0.42857142857142855. Mash 2.3 `mash dist -k 4 -s 10`:
`0.139904  0.0291702  4/10`; sourmash `contained_by` (scaled=1) 0.35 / 0.428571.

## Performance

- **Time Complexity:** O(n·k) to build the two k-mer sets (sketch mode: plus O(d log d) to sort the d distinct hashes).
- **Space Complexity:** O(distinct k-mers).

## See Also

- [kmer_distance](kmer_distance.md) — count/frequency vector metrics
- [count_kmers](count_kmers.md) — the k-mer counting modes
