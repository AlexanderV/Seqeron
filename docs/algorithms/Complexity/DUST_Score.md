# DUST Score

| Field | Value |
|-------|-------|
| Algorithm Group | Complexity |
| Test Unit ID | SEQ-COMPLEX-DUST-001 |
| Related Projects | Seqeron.Genomics.Analysis |
| Implementation Status | Complete (score + SDUST masking + dustmasker linker/soft mask + longdust k-mer generalisation) |
| Last Reviewed | 2026-09-30 |

## 1. Overview

The DUST score is a heuristic measure of nucleotide-sequence low complexity introduced for masking repetitive DNA in BLAST [1]. It counts how often each overlapping triplet (3-mer) recurs in a sequence and combines those counts into a single score: highly repetitive sequences (e.g. homopolymers, simple satellites) score high, while sequences with all-distinct triplets score 0. It is a deterministic heuristic, not a probabilistic model; this unit implements the per-sequence/per-window complexity score, the quantity DUST/SDUST thresholds to decide masking.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

Low-complexity regions (homopolymers, tandem repeats, simple sequence repeats) produce spurious alignments and inflate database-search noise. DUST ("Decreasing Uncertainty in Sequence Tags") masks such regions by scoring fixed-size windows and masking windows whose score exceeds a threshold [1]. The score function is unchanged between the original and the symmetric (SDUST) implementation; only the masking rule differs [1].

### 2.2 Core Model

For a sequence `x` of length `L`, let `c_t` be the number of occurrences of triplet `t` among the `ℓ = L − 2` overlapping triplets. The complexity score is [1][3][4]:

```
S(x) = ( Σ_t  c_t·(c_t − 1)/2 ) / (ℓ − 1)
```

The numerator `Σ_t c_t(c_t−1)/2` counts pairs of identical triplets (the reference accumulates it as `r += c[t]++` [3]). The normaliser is `ℓ − 1`: NCBI dustmasker tests `10·r > thresholds_[ℓ−1]` with `thresholds_[i] = i·level` [4], and lh3/sdust tests `new_r·10 > T·new_l` with `new_l = kdq_size(w) − i − 1 = ℓ − 1` [3]. Numerically confirmed with the compiled sdust binary (`-t 20`): an isolated 7-A run (ℓ = 5, r = 10, 10/4 = 2.5) is masked while a 6-A run (6/3 = 2.0) is not; `(AC)×6` (20/9) is masked while `ACACACACACA` (16/8 = 2.0) is not. A divisor of `ℓ` would mask none of them. (The longdust README/preprint [2] writes the SDUST score with `ℓ(x)`; this does not match the reference code and is not used.)

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | `S(x) ≥ 0` | Each `c(c−1)/2 ≥ 0` and the divisor `ℓ − 1 > 0` for `L ≥ 4`; `ℓ ≤ 1` ⇒ 0. |
| INV-02 | All-distinct triplets ⇒ `S(x) = 0` | Every `c_t = 1` ⇒ `c(c−1)/2 = 0` [2]. |
| INV-03 | Homopolymer of length `L ≥ 4` ⇒ `S = (L−2)/2` | One triplet repeated `ℓ = L−2` times: `(L−2)(L−3)/2 / (L−3)`. |
| INV-04 | Higher `S` ⇒ lower complexity | Repeated triplets increase `Σ c(c−1)/2` [1][2]. |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| sequence | `DnaSequence` / `string` | required | Sequence to score | Upper-cased for the `string` overload; null `DnaSequence` throws |
| wordSize | `int` | 3 | Word (k-mer) size | must be 3 (kept for compatibility); DUST is defined for triplets only [1][3][4] |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| score | `double` | DUST score `Σ c(c−1)/2 / (ℓ − 1)`, `ℓ = L − 2`; ≥ 0; higher ⇒ lower complexity |

### 3.3 Preconditions and Validation

`string` input is upper-cased (T↔U not performed; DNA alphabet assumed). Null `DnaSequence` ⇒ `ArgumentNullException`; null/empty `string` ⇒ 0. `wordSize ≠ 3` ⇒ `ArgumentOutOfRangeException` (B04 F34: DUST is defined for triplets only — Morgulis 2006 [1], NCBI symdust `triplet_type` [4], lh3/sdust `SD_WLEN = 3` [3], longdust README "It [SDUST] hardcodes k=3" [2], and sdust itself marks `SD_WLEN != 3` as untested (`TODO: is this right for SD_WLEN!=3?`); the former extrapolation to other k had no source; the sourced k-mer generalisation is longdust, §4.3). Fewer than two words (`ℓ ≤ 1`) ⇒ 0 (no word pair exists; defined-output convention, not a source value). `MaskLowComplexity` / `FindLowComplexityIntervals`: `windowSize < 3`, a negative/NaN/infinite `threshold` or `linker` outside 1–32 ⇒ `ArgumentOutOfRangeException`; null ⇒ `ArgumentNullException`.

## 4. Algorithm

### 4.1 High-Level Steps

1. If fewer than two triplets exist (`ℓ = L − 2 ≤ 1`), return 0.
2. Tally the count `c_t` of each overlapping word with the canonical `KmerAnalyzer.CountKmers`.
3. Sum `c_t·(c_t − 1)/2` over all distinct words.
4. Divide by `ℓ − 1` and return.

**SDUST masking (`MaskLowComplexity`)** — a port of lh3/sdust `sdust_core` [3]: slide a window of `W` bases (≤ W − 2 triplets); maintain the window score `rw` and the score `rv` of the longest window suffix in which no triplet occurs more than `2·threshold` times; when `rw > threshold·L_suffix`, extend leftwards from that suffix to find every *perfect interval* (score > threshold and ≥ the best score of any perfect interval it contains); intervals leaving the window are emitted and merged when overlapping/adjacent (dustmasker `linker`, below). Non-ACGT characters split the scan: every maximal ACGT run is scanned as an independent sequence (sdust's stated contract "N effectively breaks input into pieces of independent sequences" [3]; upstream `sdust_core` resets only `l` and `t` and leaks the triplet window across the break — e.g. `sdust -w 64 -t 20` reports `35 72` on a 53-bp input — so the window is reset here, B04 F36).

**dustmasker linker and soft masking (B04 F35)** — NCBI symdust `save_masked_regions` [4] merges a new interval into the previous one when `prev.last + linker ≥ next.first` (closed coordinates), i.e. when fewer than `linker` unmasked bases separate them; `linker` = 1 (dustmasker `DEFAULT_LINKER`; symdust accepts 1–32 and silently substitutes 1 for any other value — `-linker 0` and `-linker 50` both give the linker-1 output — so values outside 1–32 are rejected here) is exactly sdust's overlap/adjacency merge, so the default output is unchanged. `softMask` writes masked bases in lower case and the rest in upper case (dustmasker `-outfmt fasta`). `FindLowComplexityIntervals` returns the merged intervals as 0-based half-open `[start, end)` (sdust output; dustmasker `-outfmt interval` prints `start - end−1`).

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

- **Word size:** k = 3 (triplets), hardcoded in DUST/SDUST [1][3]; the `wordSize` parameter is kept for compatibility and must be 3.
- **Mask threshold:** 2.0 (reference default level `T = 20`, with `rw·10 > L·T` ⇔ score > 2.0) [3]; used by `MaskLowComplexity`.
- **Window size:** 64 bases (default for windowed masking) [1][2].

### 4.3 Longdust — sourced k-mer generalisation (`CalculateLongdustScore`, `FindLongdustRegions`)

Li & Li 2025 [2] generalise DUST to arbitrary k and long windows. For a string `x` with `ℓ(x) = |x| − k + 1` k-mer positions and k-mer counts `c_x(t)`:

```
S_L(x) = Σ_t log c_x(t)!  −  f(ℓ(x)/4^k),   f(λ) = 4^k · e^{−λ} Σ_n log(n!) λ^n / n!
```

`f` is the expected `Σ log c!` of a random sequence (Poisson k-mer counts), so `S_L` is a composite log-likelihood ratio; `x` is low-complexity when `S_L(x) − T·ℓ(x) > 0` (T = 0.6). With GC correction (`-g`) λ is split per GC class (`ld_cal_f2`). `FindLongdustRegions` is a port of lh3/longdust 1.4-r97 `ld_dust1`/`ld_dust2` (MIT): per position a backward scan (X-drop, `max_test` pre-filter, window-count bound) then forward re-scans find the leftmost start of a good interval ending at the position; overlapping hits merge; the default result is the union over both strands. Defaults = longdust `ld_opt_init`: k = 7, w = 5000, T = 0.6, X-drop 50, min start count 3. Non-ACGT symbols make overlapping k-mers ambiguous (they add −T).

Cross-check (compiled lh3/longdust 1.4-r97, commit 9491215): 1 500 random repeat-rich inputs (5 bp–12 kb, 4.5 Mb total, 410 with N/IUPAC, k ∈ 3–8, w ∈ 50–5000, T ∈ 0.3–1.0, X-drop 0/10/50/200, `-b 2/3`, `-f`, `-a`, `-g 0.3/0.41/0.6`; 12 818 reference intervals) → **0 mismatches**; 121 homopolymer edge cases around the window size + one 2 Mb sequence (ASan build) → 0 mismatches; `S_L` vs the reference's own `ld_cal_f`/`ld_cal_f2` on 3 000 inputs (k 1–10) → 3 000 bit-identical, and vs an independent lgamma/Poisson evaluation → max relative difference 8.7·10⁻⁸.

### 4.4 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `CalculateDustScore` | O(L) | O(min(L, 64)) | One pass via `KmerAnalyzer.CountKmers` |
| `MaskLowComplexity` / `FindLowComplexityIntervals` | O(L·W) typical, O(L·W³) worst | O(W + 64) | As lh3/sdust |
| `CalculateLongdustScore` | O(L + 4^k-table) | O(min(L, 4^k)) | f table O(ℓ·series) |
| `FindLongdustRegions` | O(L·w) typical (guaranteed with `approximate`) | O(4^k + w) | As lh3/longdust |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [SequenceComplexity.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/SequenceComplexity.cs)

- `SequenceComplexity.CalculateDustScore(DnaSequence, int)`: scores a `DnaSequence` (word size must be 3).
- `SequenceComplexity.CalculateDustScore(string, int)`: scores a raw string (upper-cased).
- `SequenceComplexity.MaskLowComplexity(DnaSequence, …)`: symmetric DUST (SDUST) masking of perfect intervals; default W = 64, threshold 2.0 (level 20); overload with `linker` + `softMask`.
- `SequenceComplexity.MaskLowComplexity(string, …)`: same on a raw string that may contain N/IUPAC (split into independent ACGT runs), with `linker` and `softMask`.
- `SequenceComplexity.FindLowComplexityIntervals(string|DnaSequence, …)`: the SDUST intervals as `[start, end)`.
- `SequenceComplexity.CalculateLongdustScore(string, k, gc)` / `FindLongdustRegions(string, …)`: longdust (§4.3).

### 5.2 Current Behavior

The score is computed exactly as the formula in §2.2: numerator `Σ c(c−1)/2`, divisor `ℓ − 1`. Word counting delegates to `KmerAnalyzer.CountKmers`. History: the original code divided by `words − 1` (correct); a 2026-06 change "corrected" it to the word count from a misreading of the longdust restatement [2]; the 2026-09 review restored `ℓ − 1` after confirming it against the NCBI symdust and lh3/sdust sources and the compiled sdust binary. `MaskLowComplexity` was previously a fixed 64-bp window scan that masked whole windows (and skipped sequences shorter than the window); it is now the real SDUST algorithm, verified identical to the lh3/sdust binary on 3,000 random repeat-rich sequences (W ∈ {3,5,8,16,30,64,100}, T ∈ {10,12,15,20,25,30}) and a 1-Mb sequence.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Numerator `Σ_t c_t·(c_t − 1)/2` over overlapping triplets [1][3].
- Normalization by `ℓ − 1` [1][3][4].
- Default triplet word size, window 64 and mask threshold 2.0 (level 20) [3][4].
- Symmetric perfect-interval masking (SDUST) [1][3].

- dustmasker `linker` interval merge and `-outfmt fasta` soft masking [4] — identical to dustmasker 2.12.0 on 1 500 ACGT inputs (W 8–64, level 2–30, linker 1–32; 261 outputs changed by the linker): 0 interval and 0 soft-mask mismatches.
- k-mer generalisation: longdust [2] (§4.3), not an extrapolation of the DUST formula.

**Intentionally simplified:** none. (The former `wordSize ≠ 3` extrapolation was removed, B04 F34.)

**Known reference differences:** dustmasker's own SDUST core differs from lh3/sdust for very small windows at high levels (10/1 000 cases, all W = 8, level 30: symdust's `thresholds_` table has only W − 3 entries); dustmasker additionally treats IUPAC codes as bases and only skips long N runs (`GetDustMasks_SkipNs`), so on N-containing input the reference here is sdust per ACGT run.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | `ℓ ≤ 1` ⇒ 0 | Assumption | No score defined by sources (ℓ − 1 = 0) | accepted | Defined-output convention |
| 2 | `wordSize ≠ 3` rejected | Source rule | DUST is triplet-only | resolved (F34) | k-mer generalisation = longdust |
| 3 | Window reset at non-ACGT | Deviation from upstream sdust code, conforms to its stated contract | Intervals after an N are no longer shifted/out of range | resolved (F36) | 3 000 inputs (1 037 with N): 0 mismatches vs sdust per ACGT run; upstream whole-input output differs on 274 |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Null `DnaSequence` | `ArgumentNullException` | Validation contract |
| Null/empty `string` | 0 | No words ⇒ minimal complexity |
| `ℓ ≤ 1` (L ≤ 3 for triplets) | 0 | No word pair exists (convention) |
| All-distinct triplets | 0 | INV-02 |
| Homopolymer length L ≥ 4 | `(L−2)/2` | INV-03 |

### 6.2 Limitations

Heuristic, not a probabilistic significance test. Only triplet scoring (k = 3) is source-defined for DUST; longer k use longdust. Masked intervals equal lh3/sdust output (per ACGT run for inputs with N).

## 7. Examples and Related Material

### 7.1 Worked Example

**Numerical walk-through:** `AAAAAA` (L = 6). Triplets at positions 0–3 are all `AAA` ⇒ `c_AAA = 4`, ℓ = 4. Numerator = `4·3/2 = 6`. Divisor = `ℓ − 1 = 3`. Score = `6 / 3 = 2.0`.

`ACGTACGT` (L = 8, ℓ = 6): `ACG=2, CGT=2, GTA=1, TAC=1` ⇒ numerator = `1 + 1 = 2`, divisor = 5, score = `0.4`.

Masking: `ACGTTGCAGTCATGCGATC` + `A×7` + `TGCATCGGATCCTAGGCTAA` ⇒ sdust/`MaskLowComplexity` mask [19, 26).

Linker (dustmasker 2.12.0, `-window 64 -level 20`): the 96-bp `ACGTGCATGC A×16 GCTAGCATCGACTGCAG (CA)×8 GATCGATCGTACGGTGCATGAC A×13 CT` gives `[10,26) [43,59) [81,94)` for linker 1–17 (17 unmasked bases between the first two), `[10,59) [81,94)` for linker 18–19 and `[10,94)` for linker 32.

N input: `ACGTNNAAAAAAAAAAAANACGTACACACACACACACANNGGGCCCTAGGTCA` ⇒ `[6,18) [23,38)` (sdust per ACGT run; upstream whole-input sdust: `8 20 / 21 34 / 35 72`).

### 7.3 Related Tests, Evidence, or Documents

- Tests: [SequenceComplexity_CalculateDustScore_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Analysis/SequenceComplexity_CalculateDustScore_Tests.cs) — covers `INV-01`–`INV-04`
- Evidence: [SEQ-COMPLEX-DUST-001-Evidence.md](../../../docs/Evidence/SEQ-COMPLEX-DUST-001-Evidence.md)
- Related algorithms: [K-mer_Entropy](./K-mer_Entropy.md)

## 8. References

1. Morgulis A, Gertz EM, Schäffer AA, Agarwala R. 2006. A fast and symmetric DUST implementation to mask low-complexity DNA sequences. Journal of Computational Biology 13(5):1028–1040. https://doi.org/10.1089/cmb.2006.13.1028
2. Li H, Li B. 2025. Finding low-complexity DNA sequences with longdust. arXiv:2509.07357. https://arxiv.org/abs/2509.07357 — opened via github.com/lh3/longdust (`README.md`, `tex/longdust.tex`, `longdust.c` 1.4-r97, MIT).
3. Li H. sdust — Symmetric DUST reference C implementation. https://raw.githubusercontent.com/lh3/sdust/master/sdust.c
4. NCBI C++ Toolkit — dustmasker symmetric DUST (`src/algo/dustmask/symdust.cpp`, `include/algo/dustmask/symdust.hpp`). https://raw.githubusercontent.com/ncbi/ncbi-cxx-toolkit-public/master/src/algo/dustmask/symdust.cpp
5. NCBI BLAST+ 2.12.0 `dustmasker` (Debian package `ncbi-blast+`) and `src/app/dustmask/dust_mask_app.cpp` (`GetDustMasks_SkipNs`).
