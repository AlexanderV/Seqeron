# Windowed Sequence Complexity

| Field | Value |
|-------|-------|
| Algorithm Group | Complexity |
| Test Unit ID | SEQ-COMPLEX-WINDOW-001 |
| Related Projects | Seqeron.Genomics.Analysis |
| Implementation Status | Production |
| Last Reviewed | 2026-06-14 |

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

where `V_i` is the number of distinct length-`i` subwords observed in the window, `V_{max,i}` is the maximum possible, and `m = min(6, w)` is the word-length cap [4]. `LC ∈ (0, 1]` [1].

The window enumeration emits a `ComplexityPoint` for every window fully contained in the sequence, advancing by `stepSize`:

$$ \text{starts } i \in \{0, s, 2s, \dots\} \text{ with } i + w \le L, \qquad \#\text{windows} = \left\lfloor \tfrac{L - w}{s} \right\rfloor + 1 \;\; (L \ge w) $$

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
| `sequence` | `DnaSequence` | required | Sequence to profile | Null ⇒ `ArgumentNullException` |
| `windowSize` | `int` | `64` | Window length `w` | `< 1` ⇒ `ArgumentOutOfRangeException` |
| `stepSize` | `int` | `10` | Window advance `s` | `< 1` ⇒ `ArgumentOutOfRangeException` |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `Position` | `int` | Window center, `WindowStart + windowSize/2` (0-based) |
| `ShannonEntropy` | `double` | Per-base Shannon entropy of the window (bits) |
| `LinguisticComplexity` | `double` | Summation-form linguistic complexity of the window |
| `WindowStart` | `int` | 0-based inclusive start index |
| `WindowEnd` | `int` | 0-based inclusive end index (`WindowStart + windowSize − 1`) |

### 3.3 Preconditions and Validation

0-based indexing; `WindowEnd` inclusive. Input is upper-cased by `DnaSequence`. Null sequence ⇒ `ArgumentNullException`; `windowSize < 1` or `stepSize < 1` ⇒ `ArgumentOutOfRangeException`. When `L < windowSize` the profile is empty (no partial trailing window is emitted). The result is a lazily-evaluated `IEnumerable<ComplexityPoint>`.

## 4. Algorithm

### 4.1 High-Level Steps

1. Validate inputs and normalize the sequence to upper case (via `DnaSequence`).
2. For each start `i = 0, s, 2s, …` while `i + w ≤ L`, extract the window substring.
3. Compute the window's Shannon entropy `H` over the 4 bases.
4. Compute the window's linguistic complexity with `maxWordLength = min(6, w)`.
5. Yield a `ComplexityPoint(Position=i+w/2, H, LC, WindowStart=i, WindowEnd=i+w−1)`.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

- Linguistic-complexity word-length cap `WindowLcMaxWordLength = 6` per Gabrielian & Bolshoy (1999) efficiency limitation of vocabulary assessment to a bounded W [4].
- Shannon `0·log₂0 = 0` convention [3].

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `CalculateWindowedComplexity` | O((L/s) · w²) | O(distinct subwords per window) | One pass per window (≈L/s windows); each window's LC enumerates subwords of lengths 1..min(6,w) over the w-length window |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [SequenceComplexity.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/SequenceComplexity.cs)

- `SequenceComplexity.CalculateWindowedComplexity(DnaSequence, int windowSize, int stepSize)`: sliding-window driver returning one `ComplexityPoint` per fully-contained window.
- `SequenceComplexity.CalculateShannonEntropy(...)`: per-window Shannon metric (reused internally).
- `SequenceComplexity.CalculateLinguisticComplexity(...)`: per-window LC metric (reused internally).

### 5.2 Current Behavior

The driver delegates per-window metrics to the existing `CalculateShannonEntropyCore` and `CalculateLinguisticComplexityCore` helpers, so window values match the standalone scalar metrics exactly. Windows are non-overlapping when `stepSize ≥ windowSize` and overlapping otherwise. A suffix tree was **not** used: this is a single left-to-right scan that computes scoring-based (entropy/LC) metrics over each window rather than locating exact-match occurrences, so the suffix-tree occurrence API does not fit; per-window LC subword enumeration is bounded by the small word-length cap (≤6).

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Per-window Shannon entropy `H = -Σ p log₂ p` in bits, with `0·log₂0 = 0` [3].
- Per-window linguistic complexity in summation form with `V_{max,i} = min(4^i, w−i+1)` [1][4].
- Sliding-window complexity profile over fully-contained windows [2].

**Intentionally simplified:**

- Linguistic-complexity word length is capped at `min(6, windowSize)`; **consequence:** distinct subwords longer than 6 are not counted in the per-window LC, matching the repository's LC unit and the bounded-vocabulary efficiency choice [4].

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

DNA-oriented: the Shannon metric counts only A/C/G/T and the LC denominator assumes a 4-letter alphabet. Non-ACGT symbols in a window are ignored by the entropy count and may push LC outside the usual `(0,1]` interpretation. The profile reports only windows fully inside the sequence, so the final `(L − w) mod s` bases at the 3′ end may not be covered by any window.

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

- Tests: [SequenceComplexity_CalculateWindowedComplexity_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/SequenceComplexity_CalculateWindowedComplexity_Tests.cs) — covers `INV-01`..`INV-04`
- Evidence: [SEQ-COMPLEX-WINDOW-001-Evidence.md](../../../docs/Evidence/SEQ-COMPLEX-WINDOW-001-Evidence.md)
- Related algorithms: [Linguistic_Complexity](../Sequence_Composition/Linguistic_Complexity.md), [K-mer_Entropy](./K-mer_Entropy.md)

## 8. References

1. Wikipedia. Linguistic sequence complexity. https://en.wikipedia.org/wiki/Linguistic_sequence_complexity (accessed 2026-06-14)
2. Troyanskaya, O.G., Arbell, O., Koren, Y., Landau, G.M., Bolshoy, A. 2002. Sequence complexity profiles of prokaryotic genomic sequences: a fast algorithm for calculating linguistic complexity. Bioinformatics 18(5):679–688. https://doi.org/10.1093/bioinformatics/18.5.679
3. Shannon, C.E. 1948. A Mathematical Theory of Communication. Bell System Technical Journal 27(3):379–423. https://doi.org/10.1002/j.1538-7305.1948.tb01338.x
4. Gabrielian, A., Bolshoy, A. 1999. Sequence complexity and DNA curvature. Computers & Chemistry 23(3–4):263–274. https://doi.org/10.1016/S0097-8485(99)00007-8
