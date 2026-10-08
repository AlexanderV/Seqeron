# kmer_jaccard

k-mer Jaccard similarity of two sequences (exact, or from Mash MinHash sketches), the Mash distance derived from it, the Mash p-value in sketch mode, and the exact containment indices; or sourmash FracMinHash (`scaled`) estimates, optionally with the abundance-weighted angular similarity and weighted containments.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `kmer_jaccard` |
| **Method ID** | `KmerAnalyzer.JaccardSimilarity` |
| **Version** | 1.4.0 |
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
- `scaled = S` in 1..4294967295 (sourmash's u32 `scaled`; mutually exclusive with `sketchSize`; requires
  `canonical = true`, since sourmash DNA hashing is always canonical — Rust `signature.rs` `SeqToHashes` hashes
  `min(kmer, krc)` — so `canonical = false` is rejected instead of being ignored; `acgtOnly` is implied, as sourmash
  `add_sequence(force=True)` skips every k-mer with a non-ACGT symbol and raises without `force`): sourmash 4.9.4 FracMinHash, delegating to
  `KmerAnalyzer.CreateFracMinHashSketch` + `KmerAnalyzer.CompareFracMinHashSketches` (Irber et al. 2022; Hera, Pierce-Ward
  & Koslicki 2023). Canonical k-mers (lexicographically smaller of k-mer and reverse complement, non-ACGT windows
  skipped = `add_sequence(force=True)`), MurmurHash3_x64_128 h1 with seed 42, every hash ≤ max_hash kept, where
  max_hash = `(u64::MAX as f64 / S) as u64` (Rust `max_hash_for_scaled`; S = 1 keeps all). `jaccard` = shared / max(1, union)
  (`MinHash.jaccard`); the containments are `contained_by` = shared / (|A|·(1 − (1 − 1/S)^(|A|·S))) clamped to [0, 1];
  `maxContainment` = `max_containment` (min(|A|, |B|) in the denominator); `sharedHashes` = |A ∩ B|,
  `sketchDenominator` = |A ∪ B|; `mashDistance` = the Mash formula applied to the estimated J; `pValue` null.
- `trackAbundance = true` (needs `scaled > 0`): both sketches record each hash's (canonical) k-mer count (sourmash
  `track_abundance=True`) and `angularSimilarity` = sourmash `angular_similarity` (= `similarity(ignore_abundance=False)`;
  Rust `KmerMinHash::angular_similarity`): cos = Σ_{h∈A∩B} a_h·b_h / (‖a‖·‖b‖) over the abundance vectors (capped at 1),
  similarity = 1 − 2·acos(cos)/π (0 when a sketch is empty). `weightedContainmentSeq1InSeq2` /
  `weightedContainmentSeq2InSeq1` = sourmash `contained_by_weighted` = Σ_{h∈A∩B} a_h / Σ_{h∈A} a_h (the
  abundance-weighted fraction of A's k-mers shared with B; not bias-corrected). The other fields are unchanged by
  abundance tracking.

## Core Documentation Reference

- Source: [KmerAnalyzer.cs#L1406](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/KmerAnalyzer.cs#L1406)
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
| `scaled` | integer | No | 0 (default) = no FracMinHash; 1..4294967295 = sourmash scaled factor (sourmash default 1000); requires `canonical = true`; not with `sketchSize` |
| `trackAbundance` | boolean | No | With `scaled > 0`: track abundances and return `angularSimilarity` and the weighted containments (default false) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `jaccard` | number | Jaccard index in [0, 1] |
| `mashDistance` | number | Mash distance in [0, 1] |
| `containmentSeq1InSeq2` | number | Exact containment |A ∩ B| / |A| (scaled mode: sourmash `contained_by`) |
| `containmentSeq2InSeq1` | number | Exact containment |A ∩ B| / |B| (scaled mode: sourmash `contained_by`) |
| `sharedHashes` | integer \| null | Sketch mode: shared hashes x (Mash "x/s"); scaled mode: |A ∩ B|; null when exact |
| `sketchDenominator` | integer \| null | Sketch mode: union hashes compared (s unless both sketches are smaller); scaled mode: |A ∪ B|; null when exact |
| `pValue` | number \| null | Sketch mode: Mash p-value; null otherwise |
| `maxContainment` | number \| null | Scaled mode: sourmash `max_containment`; null otherwise |
| `angularSimilarity` | number \| null | `trackAbundance`: sourmash `angular_similarity`; null otherwise |
| `weightedContainmentSeq1InSeq2` | number \| null | `trackAbundance`: sourmash `seq1.contained_by_weighted(seq2)`; null otherwise |
| `weightedContainmentSeq2InSeq1` | number \| null | `trackAbundance`: sourmash `seq2.contained_by_weighted(seq1)`; null otherwise |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1003 | k must be positive |
| 1003 | sketchSize must be >= 0 |
| 1003 | k must be <= 32 when sketchSize > 0 |
| 1003 | scaled must be in 0..4294967295 |
| 1003 | sketchSize and scaled are mutually exclusive (Mash bottom-s vs sourmash FracMinHash) |
| 1003 | scaled > 0 requires canonical=true (sourmash DNA FracMinHash hashes canonical k-mers only) |
| 1003 | trackAbundance requires scaled > 0 (sourmash FracMinHash) |

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

### Example 4: sourmash FracMinHash (scaled)

```json
{ "tool": "kmer_jaccard", "arguments": { "seq1": "AGGTAAGGTGGTTGAGATCTGGACTTTTGACGCCTGGAGCCCGCAGTGCTCCTCGAAAAGTAGCCATGCCTTGGGCTGCT",
  "seq2": "CAAAGGCCCTACCTTCTTATAGTCCTTTCAACATACAAGTATAGTTGGAAGTTCTAAGTTCAGTTTAATC", "k": 4, "canonical": true, "scaled": 3,
  "trackAbundance": true } }
```
**Response:** `jaccard` 0.3103448275862069, `sharedHashes` 9, `sketchDenominator` 29, `containmentSeq1InSeq2`
0.4285714285748822, `containmentSeq2InSeq1` 0.52941176525941, `maxContainment` 0.52941176525941, `mashDistance`
0.18680360045755526, `angularSimilarity` 0.3123095603640216, `weightedContainmentSeq1InSeq2` 0.4642857142857143,
`weightedContainmentSeq2InSeq1` 0.5238095238095238, `pValue` null. sourmash 4.9.4
`MinHash(n=0, ksize=4, scaled=3, track_abundance=True)`: len 21 / 17, `count_common` 9, `jaccard` 0.3103448275862069,
`contained_by` 0.4285714285748822 / 0.52941176525941, `max_containment` 0.52941176525941, `angular_similarity`
0.3123095603640216, `contained_by_weighted` 0.4642857142857143 / 0.5238095238095238 (scaled=1: 0.2363801370444173,
0.37662337662337664 / 0.417910447761194).

## Performance

- **Time Complexity:** O(n·k) to build the two k-mer sets (sketch / scaled mode: plus O(d log d) to sort the d kept hashes).
- **Space Complexity:** O(distinct k-mers).

## See Also

- [kmer_distance](kmer_distance.md) — count/frequency vector metrics
- [count_kmers](count_kmers.md) — the k-mer counting modes
