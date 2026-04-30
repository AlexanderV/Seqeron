# Position Weight Matrix (PWM)

| Field | Value |
|-------|-------|
| Algorithm Group | Pattern Matching |
| Test Unit ID | PAT-PWM-001 |
| Related Projects | N/A |
| Implementation Status | Simplified |
| Last Reviewed | 2026-04-30 |

## 1. Overview

A Position Weight Matrix (PWM), also called a position-specific scoring matrix, models a conserved motif by storing per-position nucleotide scores. In this repository, PWMs are constructed from aligned DNA sequences using log-odds scores with pseudocount smoothing and scanned across target sequences using a simple threshold test.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

PWMs are a standard representation for transcription factor binding sites and other conserved sequence motifs. They are typically derived from aligned sequences, smoothed with pseudocounts, and scored against candidate windows by summing per-position contributions. Sources preserved from the original document: Wikipedia (Position weight matrix), Kel et al. (2003), Nishida et al. (2008), Rosalind CONS, Stormo (2000).

### 2.2 Core Model

The original document's construction equations apply directly to the implementation:

$$
PFM_{k,j} = \sum_{i=1}^{N} \mathbf{1}(X_{i,j} = k)
$$

$$
PPM_{k,j} = \frac{PFM_{k,j} + p}{N + |\Sigma| \cdot p}
$$

$$
PWM_{k,j} = \log_2\left(\frac{PPM_{k,j}}{b_k}\right)
$$

Where `p` is the pseudocount, `|Σ| = 4` for DNA, and the current implementation fixes `b_k = 0.25` for all nucleotides.

The sequence score for a window of length `L` is:

$$
Score(S) = \sum_{j=1}^{L} PWM_{S_j, j}
$$

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | `Length` equals the aligned training-sequence length | `CreatePwm(...)` validates equal sequence lengths before construction |
| INV-02 | `Consensus` uses the highest-scoring base in each PWM column | `PositionWeightMatrix.GenerateConsensus()` picks the per-column maximum |
| INV-03 | `MaxScore` and `MinScore` are sums of per-column extrema | The `PositionWeightMatrix` properties aggregate columnwise maxima and minima |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `sequences` | `IEnumerable<string>` | required | Aligned DNA training sequences used to build the PWM | Must be non-null, non-empty, equal length, and contain only `A/C/G/T` |
| `pseudocount` | `double` | `0.25` | Smoothing parameter used during PWM construction | Applied uniformly to all four bases |
| `sequence` | `DnaSequence` | required | Sequence scanned with an existing PWM | Null input throws `ArgumentNullException` |
| `pwm` | `PositionWeightMatrix` | required | Matrix used for scoring windows | Null input throws `ArgumentNullException` |
| `threshold` | `double` | `0.0` | Minimum score required for a reported match | Match condition is `score >= threshold` |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `Matrix` | `double[,]` | Log-odds PWM with 4 rows (`A,C,G,T`) and `Length` columns |
| `Length` | `int` | Motif width |
| `Consensus` | `string` | Highest-scoring base per position |
| `MaxScore` | `double` | Sum of per-column maxima |
| `MinScore` | `double` | Sum of per-column minima |
| `MotifMatch` | `MotifMatch` | Scan result with `Position`, `MatchedSequence`, `Pattern = pwm.Consensus`, and `Score` |

### 3.3 Preconditions and Validation

`CreatePwm(...)` throws `ArgumentNullException` when `sequences` is null and `ArgumentException` when the collection is empty, when lengths differ, or when any character is outside `A/C/G/T`. `ScanWithPwm(...)` throws `ArgumentNullException` for null `sequence` or `pwm` and returns no matches when the target sequence is shorter than the PWM length.

## 4. Algorithm

### 4.1 High-Level Steps

1. Uppercase all training sequences.
2. Validate that at least one sequence is present, all lengths match, and all characters are `A/C/G/T`.
3. Count nucleotide occurrences at each position.
4. Apply pseudocount smoothing and convert the resulting frequencies to log-odds scores against a uniform background of `0.25`.
5. Build a `PositionWeightMatrix` and derive its consensus from the maximum score in each column.
6. For scanning, slide the PWM across the sequence, sum the per-position scores, and yield matches whose score is at least the threshold.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

The PWM matrix layout is fixed:

```text
Matrix[4, Length]
Row 0 = A
Row 1 = C
Row 2 = G
Row 3 = T
```

`ScanWithPwm(...)` marks a window invalid if any scanned character maps to a negative base index, although `DnaSequence` input normally limits the alphabet to `A/C/G/T`.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `CreatePwm` | `O(N × L)` | `O(L)` auxiliary beyond the `4 × L` matrix | `N` aligned sequences of length `L` |
| `ScanWithPwm` | `O(S × L)` | `O(1)` auxiliary per window | `S` is target sequence length |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [MotifFinder.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/MotifFinder.cs)

- `MotifFinder.CreatePwm(IEnumerable<string>, double)`: Builds a DNA PWM.
- `MotifFinder.ScanWithPwm(DnaSequence, PositionWeightMatrix, double)`: Scores each sequence window against the PWM.
- `PositionWeightMatrix`: Holds `Matrix`, `Length`, `Consensus`, `MaxScore`, and `MinScore`.

### 5.2 Current Behavior

`CreatePwm(...)` uppercases all training sequences, uses a default pseudocount of `0.25`, and always computes log-odds scores against a uniform background frequency of `0.25`. `ScanWithPwm(...)` reports matches whose score is greater than or equal to the threshold and uses the PWM consensus as the `Pattern` field in returned `MotifMatch` values.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- PWM construction from aligned sequences.
- Pseudocount smoothing before log-odds conversion.
- Window scoring by summing per-position PWM values.

**Intentionally simplified:**

- The background distribution is fixed at `0.25` for each DNA base; **consequence:** callers cannot model GC-biased or otherwise nonuniform backgrounds.
- The implementation is DNA-specific with four matrix rows; **consequence:** it does not directly support protein alphabets or other symbol sets.

**Not implemented:**

- Insertions, deletions, or profile-HMM style state transitions; **users should rely on:** other motif models if gap-aware scoring is required.

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Empty training set | Throws `ArgumentException` | At least one aligned sequence is required |
| Unequal training-sequence lengths | Throws `ArgumentException` | PWM columns require alignment |
| Non-ACGT training character | Throws `ArgumentException` | Source performs strict validation |
| Single training sequence | Produces a valid PWM | The implementation permits `Count = 1` |
| Sequence shorter than PWM | Returns no matches | The scan loop does not execute |

### 6.2 Limitations

The current implementation assumes a uniform DNA background and does not expose alternative alphabets or richer probabilistic motif models. It also returns raw PWM scores without converting them to calibrated probabilities or p-values.

## 7. Examples and Related Material

### 7.1 Worked Example

**Numerical / biological walk-through (optional):**

The original document highlights these related motif representations and alternatives:

- Consensus sequence: one best character per column.
- IUPAC degenerate matching: ambiguity-code pattern matching without weighted scores.
- Hidden Markov Models: extension with insertion and deletion probabilities.

## 8. References

1. Wikipedia contributors. "Position weight matrix." *Wikipedia, The Free Encyclopedia*. https://en.wikipedia.org/wiki/Position_weight_matrix
2. Kel, A.E. et al. (2003). "MATCH: A tool for searching transcription factor binding sites." *Nucleic Acids Research* 31(13):3576-3579.
3. Nishida, K.; Frith, M.C.; Nakai, K. (2008). "Pseudocounts for transcription factor binding sites." *Nucleic Acids Research* 37(3):939-944.
4. Rosalind. "Consensus and Profile." https://rosalind.info/problems/cons/
5. Stormo, G.D. (2000). "DNA binding sites: representation and discovery." *Bioinformatics* review article.
