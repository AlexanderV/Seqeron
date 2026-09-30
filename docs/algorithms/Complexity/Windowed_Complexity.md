# Windowed Sequence Complexity

| Field | Value |
|-------|-------|
| Algorithm Group | Complexity |
| Test Unit ID | SEQ-COMPLEX-WINDOW-001 |
| Related Projects | Seqeron.Genomics.Analysis |
| Implementation Status | Production |
| Last Reviewed | 2026-09-30 |

## 1. Overview

Windowed sequence complexity produces a *complexity profile*: it slides a fixed-size window along a DNA sequence and, for each window, reports two complexity metrics — Shannon entropy of the per-base distribution (bits) and linguistic complexity (summation form) — together with the window's coordinates. It is the per-position view of the scalar complexity metrics, used to locate low-complexity stretches (simple repeats, homopolymers) along a longer sequence [1][2]. The computation is exact: each window's metrics are the deterministic Shannon-entropy and linguistic-complexity values for that substring.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

Complexity profiles plot a complexity measure as a function of position along a genomic sequence; low-complexity regions appear as troughs in the profile [2]. The two metrics combined here are the classical per-symbol Shannon entropy [3] and linguistic complexity, a vocabulary-richness measure for nucleotide sequences [1][4].

### 2.2 Core Model

For a window `W` of length `w`:

Shannon entropy over the four DNA bases (Shannon 1948 [3]):

$$ H(W) = - \sum_{b \in \{A,C,G,T\}} p_b \log_2 p_b $$

where `p_b` is the frequency of base `b` in the window; `0·log₂0` is taken as 0 [3]. `H` ranges from 0 (a homopolymer / deterministic distribution) to `log₂4 = 2.0` bits (uniform distribution) [3].

Linguistic complexity (summation form, as used by the repository's linguistic-complexity unit) [1][4]:

$$ LC(W) = \frac{\sum_{i=1}^{m} V_i}{\sum_{i=1}^{m} V_{max,i}}, \qquad V_{max,i} = \min(4^{i}, w - i + 1) $$

where `V_i` is the number of distinct length-`i` subwords observed in the window, `V_{max,i}` is the maximum possible, and `m = min(lcMaxWordLength, w)` is the word-length cap (default `lcMaxWordLength = 6`) [4]. `LC ∈ (0, 1]` [1].

The window enumeration emits a `ComplexityPoint` for every window fully contained in the sequence, advancing by `stepSize`:

$$ \text{starts } i \in \{0, s, 2s, \dots\} \text{ with } i + w \le L, \qquad \#\text{windows} = \left\lfloor \tfrac{L - w}{s} \right\rfloor + 1 \;\; (L \ge w) $$

### 2.3 BBDuk k-mer entropy (low-entropy masking)

BBTools BBDuk (`entropy=c entropymask=t`, `tracker/EntropyTracker.java` [5]) scores a window of `w` bases by the
entropy of its `W_k = w − k + 1` overlapping k-mers, normalised by `ln W_k`:

$$ e(W) = \frac{-\sum_j p_j \ln p_j}{\ln W_k}, \qquad p_j = c_j / W_k $$

(defaults `entropyk = 5`, `entropywindow = 50`; the normaliser is `ln W_k` even when `4^k < W_k`). A window fails
when `e < c` (single-precision comparison); only windows with no undefined base (`ns() < 1`; A/C/G/T/U, either case,
are defined) are tested; the masked set is the union of failing windows. With `k = 1`, `e = H_bits / log₂ w`, so the
per-base Shannon scan of §5.2a is BBDuk with `entropyk=1 entropy=t/log₂w`.

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | Number of points = ⌊(L−w)/s⌋+1 for L≥w, else 0 | Loop emits one point per start `i` with `i+w ≤ L`, stepping by `s` |
| INV-02 | Per point: WindowStart=i, WindowEnd=i+w−1, Position=i+⌊w/2⌋ (0-based, end inclusive) | `ComplexityPoint` construction |
| INV-03 | 0 ≤ ShannonEntropy ≤ log₂4 = 2.0 | Shannon entropy bounds for a 4-symbol alphabet [3] |
| INV-04 | 0 < LinguisticComplexity ≤ 1 for DNA windows | LC range (0,1) [1] |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `sequence` | `DnaSequence` or `string` | required | Sequence to profile (string: any symbols, upper-cased) | Null ⇒ `ArgumentNullException` |
| `windowSize` | `int` | `64` | Window length `w` | `< 1` ⇒ `ArgumentOutOfRangeException` |
| `stepSize` | `int` | `10` | Window advance `s` | `< 1` ⇒ `ArgumentOutOfRangeException` |
| `lcMaxWordLength` | `int` | `6` | Per-window LC word-length cap `m` (≥ `w` ⇒ all-length LC) | `< 1` ⇒ `ArgumentOutOfRangeException` |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `Position` | `int` | Window center, `WindowStart + windowSize/2` (0-based) |
| `ShannonEntropy` | `double` | Per-base Shannon entropy of the window (bits) |
| `LinguisticComplexity` | `double` | Summation-form linguistic complexity of the window |
| `WindowStart` | `int` | 0-based inclusive start index |
| `WindowEnd` | `int` | 0-based inclusive end index (`WindowStart + windowSize − 1`) |

### 3.3 Preconditions and Validation

0-based indexing; `WindowEnd` inclusive. Input is upper-cased by `DnaSequence`. Null sequence ⇒ `ArgumentNullException`; `windowSize < 1` or `stepSize < 1` ⇒ `ArgumentOutOfRangeException`. When `L < windowSize` the profile is empty (no partial trailing window is emitted). The `string` overload upper-cases the input and skips (emits no point for) every window containing a symbol other than A/C/G/T/U, the BBDuk rule of scoring only windows with `ns() < 1` [5]; on A/C/G/T input it equals the `DnaSequence` overload. The result is a lazily-evaluated `IEnumerable<ComplexityPoint>`.

## 4. Algorithm

### 4.1 High-Level Steps

1. Validate inputs and normalize the sequence to upper case (via `DnaSequence`).
2. For each start `i = 0, s, 2s, …` while `i + w ≤ L`, extract the window substring.
3. Compute the window's Shannon entropy `H` over the 4 bases.
4. Compute the window's linguistic complexity with `maxWordLength = min(lcMaxWordLength, w)` (default 6).
5. Yield a `ComplexityPoint(Position=i+w/2, H, LC, WindowStart=i, WindowEnd=i+w−1)`.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

- Linguistic-complexity word-length cap `lcMaxWordLength` (default 6). Gabrielian & Bolshoy (1999) bound the word length of the windowed LC "not in the range of 2 to N−1 but only up to W" for efficiency [4]; W is a free parameter there (universalmotif `sequence_complexity` uses 7 [6]), and 6 is this library's historical default, kept for backward compatibility.
- Shannon `0·log₂0 = 0` convention [3].

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `CalculateWindowedComplexity` | O((L/s) · w²) | O(distinct subwords per window) | One pass per window (≈L/s windows); each window's LC enumerates subwords of lengths 1..min(6,w) over the w-length window |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [SequenceComplexity.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/SequenceComplexity.cs)

- `SequenceComplexity.CalculateWindowedComplexity(DnaSequence | string, int windowSize, int stepSize, int lcMaxWordLength = 6)`: sliding-window driver returning one `ComplexityPoint` per fully-contained (string: fully-defined) window.
- `SequenceComplexity.FindLowComplexityRegions(DnaSequence | string, int windowSize, double entropyThreshold)`: per-base Shannon window-union scan (§5.2a).
- `SequenceComplexity.FindLowEntropyRegionsBbduk(string, double entropyCutoff, int windowSize = 50, int k = 5)`: BBDuk `maskLowEntropy` port (§5.2b).
- `SequenceComplexity.CalculateShannonEntropy(...)`: per-window Shannon metric (reused internally).
- `SequenceComplexity.CalculateLinguisticComplexity(...)`: per-window LC metric (reused internally).

### 5.2 Current Behavior

The driver delegates per-window metrics to the existing `CalculateShannonEntropyCore` and `CalculateLinguisticComplexityCore` helpers, so window values match the standalone scalar metrics exactly. Windows are non-overlapping when `stepSize ≥ windowSize` and overlapping otherwise. A suffix tree was **not** used: this is a single left-to-right scan that computes scoring-based (entropy/LC) metrics over each window rather than locating exact-match occurrences, so the suffix-tree occurrence API does not fit; per-window LC subword enumeration is bounded by the small word-length cap (≤6).

### 5.2a Low-complexity regions (`FindLowComplexityRegions`)

`SequenceComplexity.FindLowComplexityRegions(DnaSequence | string, windowSize = 64, entropyThreshold = 1.0)` scans every
step-1 window with the same per-base Shannon kernel (1-mer entropy in bits, not normalised) and flags windows with
`H < entropyThreshold` (strict). A region is a maximal run of positions covered by flagged windows — the union of
flagged windows, the bit-mask reporting rule of BBTools BBDuk `maskLowEntropy` (each failing window sets
`[leftPos, rightPos]`) [5]. The window *statistic* is not BBDuk's default (BBDuk uses 5-mer entropy normalised by
`ln W_k`, §2.3); this method equals BBDuk with `entropyk=1 entropy=t/log₂w` (1 600 random cases vs `bbduk.sh` 40.02,
0 mismatches). `entropyThreshold` must be finite and ≥ 0 (else `ArgumentOutOfRangeException`). The `string` overload
never flags a window containing a non-A/C/G/T/U symbol (BBDuk `ns() < 1`); a region may still span such a symbol when
flagged windows on both sides overlap it.
Overlapping or abutting flagged windows therefore merge; regions are disjoint and ascending; `End` is inclusive and equals
the end of the last flagged window of the run; `MinEntropy` is the minimum window entropy in the run. Validation
(null / `windowSize < 1`) is eager. Worked example: ATGC×20 + A×64 + ATGC×20, w=20, threshold 0.5 → flagged window
starts 79..126 → one region 79..145 (length 67, MinEntropy 0).

### 5.2b BBDuk-faithful low-entropy masking (`FindLowEntropyRegionsBbduk`)

Port of BBMap 40.02 `BBDuk.maskLowEntropy` + `tracker/EntropyTracker` (FAST mode: incremental running sum of the
precomputed `p·ln p` table, undefined bases encoded as A in k-mers but counted by `ns`, `float` comparison). Returns
the regions (inclusive `End`) that `bbduk.sh entropy=c entropymask=t entropywindow=w entropyk=k` masks; `MinEntropy`
is the lowest normalised entropy of the failing windows. Reads shorter than `w` are never masked. Validation mirrors the
EntropyTracker assertions: `1 ≤ k ≤ 15`, `k < w`, `0 ≤ c ≤ 1`. Cross-check vs compiled-release `bbduk.sh` 40.02:
4 500 random cases (1–400 bp, low-complexity segments, N/IUPAC/U, lower case; w ∈ {k+1, k+2, 20, 25, 50, 64, 100},
k ∈ {1..6, 8, 10, 12}, 13 cutoffs) + 20 reads of 20–60 kb → **0 mismatches**. Example: 30 bp + A×40 + 30 bp,
defaults, `entropy=0.5` → masked 11..88 (MinEntropy 0.26749653).

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Per-window Shannon entropy `H = -Σ p log₂ p` in bits, with `0·log₂0 = 0` [3].
- Per-window linguistic complexity in summation form with `V_{max,i} = min(4^i, w−i+1)` [1][4].
- Sliding-window complexity profile over fully-contained windows [2].

**Intentionally simplified:** none. The per-window LC word-length cap is a parameter (`lcMaxWordLength`, default 6; `≥ w` gives the full Troyanskaya LC of each window) [4].

**Not implemented:**

- Suffix-tree-based linear-time profile of Troyanskaya et al. (2002); **users should rely on:** the direct per-window enumeration here, which is exact for the bounded word lengths used.

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| L < windowSize | Empty profile (0 points) | No window is fully contained; partial windows are not emitted |
| L = windowSize | Exactly 1 point (start 0) | A single fully-contained window |
| Homopolymer window | ShannonEntropy = 0 | Deterministic base distribution [3] |
| Uniform window (`ACGT…`) | ShannonEntropy = 2.0 | Uniform 4-base distribution [3] |
| Null DnaSequence | `ArgumentNullException` | Explicit guard |
| windowSize < 1 / stepSize < 1 | `ArgumentOutOfRangeException` | Explicit guard |

### 6.2 Limitations

DNA-oriented: the Shannon metric counts only A/C/G/T (U = T) and the LC denominator assumes a 4-letter alphabet. The `DnaSequence` input admits only A/C/G/T and the `string` overload skips windows with any other symbol, so every window has alphabet size 4 and LC ∈ (0, 1] (the raw-string LC extends the alphabet by any other symbol present, see Linguistic_Complexity.md). The profile reports only windows fully inside the sequence, so the final `(L − w) mod s` bases at the 3′ end may not be covered by any window.

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
var seq = new DnaSequence("ACGTACGTAAAAAAAAACGTACGT");
var profile = SequenceComplexity.CalculateWindowedComplexity(seq, windowSize: 8, stepSize: 8).ToList();
// profile.Count == 3; profile[0].WindowStart == 0, WindowEnd == 7, Position == 4
```

**Numerical walk-through:** window `ACGTACGT` (w=8): bases A=C=G=T=2 ⇒ H = log₂4 = 2.0. Distinct subwords by length 1..6 = 4,4,4,4,4,3 (sum 23); maxima min(4^i,8−i+1) = 4,7,6,5,4,3 (sum 29) ⇒ LC = 23/29 = 0.7931034482758621. Window `AAAAAAAA`: H = 0; distinct = 1 per length (sum 6) ⇒ LC = 6/29 = 0.20689655172413793.

### 7.3 Related Tests, Evidence, or Documents

- Tests: [SequenceComplexity_CalculateWindowedComplexity_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Analysis/SequenceComplexity_CalculateWindowedComplexity_Tests.cs) — covers `INV-01`..`INV-04`
- Evidence: [SEQ-COMPLEX-WINDOW-001-Evidence.md](../../../docs/Evidence/SEQ-COMPLEX-WINDOW-001-Evidence.md)
- Related algorithms: [Linguistic_Complexity](../Sequence_Composition/Linguistic_Complexity.md), [K-mer_Entropy](./K-mer_Entropy.md)

## 8. References

1. Wikipedia. Linguistic sequence complexity. https://en.wikipedia.org/wiki/Linguistic_sequence_complexity (accessed 2026-06-14)
2. Troyanskaya, O.G., Arbell, O., Koren, Y., Landau, G.M., Bolshoy, A. 2002. Sequence complexity profiles of prokaryotic genomic sequences: a fast algorithm for calculating linguistic complexity. Bioinformatics 18(5):679–688. https://doi.org/10.1093/bioinformatics/18.5.679
3. Shannon, C.E. 1948. A Mathematical Theory of Communication. Bell System Technical Journal 27(3):379–423. https://doi.org/10.1002/j.1538-7305.1948.tb01338.x
4. Gabrielian, A., Bolshoy, A. 1999. Sequence complexity and DNA curvature. Computers & Chemistry 23(3–4):263–274. https://doi.org/10.1016/S0097-8485(99)00007-8
5. Bushnell, B. BBTools — BBDuk `maskLowEntropy` / `EntropyTracker` (`calcEntropyFast`, `add`, `passes`). BBMap 40.02 release (`bbtools.jar`: `jgi/BBDuk.java`, `tracker/EntropyTracker.java`, `dna/AminoAcid.java`; sourceforge.net/projects/bbmap) and https://raw.githubusercontent.com/BioInfoTools/BBMap/master/current/structures/EntropyTracker.java (accessed 2026-09-30)
6. Tremblay, B.J.M. universalmotif `sequence_complexity` (R/sequence_complexity.R, `trifonov.max.word.size = 7`). https://raw.githubusercontent.com/bjmt/universalmotif/master/R/sequence_complexity.R (accessed 2026-09-30)
