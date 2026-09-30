# Position Weight Matrix (PWM)

| Field | Value |
|-------|-------|
| Algorithm Group | Pattern Matching |
| Test Unit ID | PAT-PWM-001 |
| Related Projects | N/A |
| Implementation Status | Complete |
| Last Reviewed | 2026-09-29 |

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

Where `p` is the pseudocount added to **every cell** (total `4p` per column; identical to Biopython `counts.normalize(pseudocounts=p)`), `|Σ| = 4` for DNA, and `b_k` is the background probability: `0.25` for `CreatePwm(sequences, p)`, or a caller-supplied distribution for `CreatePwm(sequences, p, background)` (normalised to sum 1, as Biopython `log_odds(background=...)`). Logarithm base is 2 (bits), as in Biopython `log_odds` and Wasserman & Sandelin (2004).

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
| `pseudocount` | `double` | `0.25` | Smoothing parameter added to each cell | Finite and `>= 0` (else `ArgumentOutOfRangeException`); `0` gives `-inf` for unseen bases |
| `background` | `IReadOnlyList<double>` | uniform | Background probabilities A,C,G,T (overload) | 4 finite values `> 0`; normalised to sum 1 |
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

`CreatePwm(...)` throws `ArgumentNullException` when `sequences` (or `background`) is null, `ArgumentException` when the collection is empty, contains a null element, lengths differ, or any character is outside `A/C/G/T`, and `ArgumentOutOfRangeException` for a negative/non-finite pseudocount or a non-positive/non-finite background value. The `PositionWeightMatrix(double[,], int)` constructor requires a non-null `4 × length` matrix. `ScanWithPwm(...)` scans the forward strand only; windows containing a non-ACGT symbol are skipped (Biopython `calculate` yields NaN for them, so `search` never reports them); results are in ascending position order with `score >= threshold` (Biopython `search` semantics). For the reverse strand, scan with `pwm.ReverseComplement()` (Biopython `reverse_complement`); positions are forward-strand window starts (Biopython's `search(both=True)` reports the same hit as `position - len(seq)`). `ScanWithPwm(...)` throws `ArgumentNullException` for null `sequence` or `pwm` and returns no matches when the target sequence is shorter than the PWM length.

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
- `MotifFinder.CreatePwm(IEnumerable<string>, double, IReadOnlyList<double>)`: Builds a DNA PWM against a non-uniform background.
- `PositionWeightMatrix`: Holds `Matrix`, `Length`, `Consensus` (first maximum in A,C,G,T order — Biopython `consensus` tie rule), `MaxScore`, `MinScore` (sums of column extrema — Biopython `pssm.max`/`pssm.min`) and `ReverseComplement()`.

Additive members ([MotifFinder.PwmScoring.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/MotifFinder.PwmScoring.cs), 2026-09 B05 follow-up):

- `MotifFinder.ScanWithPwmBothStrands(DnaSequence, PositionWeightMatrix, double)` → `PwmStrandMatch(Position, BiopythonPosition, Strand, MatchedSequence, Pattern, Score)` — Biopython `pssm.search(seq, threshold, both=True)`: minus strand scored with `ReverseComplement()` over the same forward windows; `Position` = forward window start on both strands, `BiopythonPosition` = start on '+', start − n on '−'; ascending start, '+' before '−' on ties (Biopython uses NumPy's unstable `argsort`, so its tie order is implementation-defined — observed: `(0,+),(−12,−),(−7,−),(5,+)` for a palindromic PWM); minus `MatchedSequence` = the site read 5'→3' on the minus strand.
- `MotifFinder.CalculatePwmScores(string | DnaSequence, PositionWeightMatrix)` — Biopython `calculate`: one score per window, NaN for windows with non-ACGT symbols, case-insensitive; empty array when n < m (Biopython raises). Scores are double (Biopython float32).
- `PositionWeightMatrix.FromCounts(double[,] counts, double pseudocount = 0 | IReadOnlyList<double> pseudocounts, background?)` — Biopython `counts.normalize(pseudocounts).log_odds(background)` with per-column totals (JASPAR matrices may have unequal column sums); `MotifFinder.JasparPseudocounts(counts, background?)` — Biopython `motifs.jaspar.calculate_pseudocounts` (√N̄·q[b]). `CreatePwm` now uses the same private log-odds kernel (bit-identical results).
- `PositionWeightMatrix.Mean(background?)`, `Std(background?)` — Biopython `pssm.mean` / `pssm.std`.
- `PositionWeightMatrix.ScoreDistribution(background?, precision = 1000)` → `PwmScoreDistribution` with `ThresholdFpr`, `ThresholdFnr`, `ThresholdBalanced(rateProportion[, out fpr])`, `ThresholdPatser()` — line-by-line port of Biopython `Bio.motifs.thresholds.ScoreDistribution` (Dojer 2008), including CPython float floor-division for grid indices, so thresholds equal Biopython to ~1e-12. Non-finite matrices (pseudocount 0) throw (`InvalidOperationException`; Biopython would produce an infinite grid step). The FPR/FNR loops stop at the grid ends instead of running past them (Biopython would wrap to negative indices / raise).
- The single window-scoring kernel `ScorePwmWindow` is shared by `ScanWithPwm`, `ScanWithPwmBothStrands` and `CalculatePwmScores`.

### 5.2 Current Behavior

`CreatePwm(...)` uppercases all training sequences and uses a default pseudocount of `0.25`; the two-argument overload computes log-odds against a uniform background (0.25), the background overload against a normalised arbitrary background. `ScanWithPwm(...)` reports matches whose score is greater than or equal to the threshold and uses the PWM consensus as the `Pattern` field in returned `MotifMatch` values.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- PWM construction from aligned sequences.
- Pseudocount smoothing before log-odds conversion.
- Window scoring by summing per-position PWM values.

**Scope choices (not simplifications):**

- Pseudocounts are a scalar added to every cell (Biopython float `pseudocounts`); background-proportional pseudocount distribution (Wasserman & Sandelin 2004 `sqrt(N)·b_k`) can be emulated only for a uniform background.
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

DNA alphabet only. Score calibration is the Biopython discretised score distribution (`ScoreDistribution`: FPR/FNR/balanced/patser thresholds); exact p-values per score (e.g. TFM-Pvalue / FIMO's exact DP) are not provided.

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
6. Wasserman, W.W.; Sandelin, A. (2004). "Applied bioinformatics for the identification of regulatory elements." *Nat Rev Genet* 5:276-287. doi:10.1038/nrg1315.
7. Biopython 1.88 `Bio.motifs.matrix` (`normalize`, `log_odds`, `max`/`min`, `consensus`, `reverse_complement`, `calculate`, `search`) — reference implementation used to lock the test values.
8. Biopython 1.88 source (installed package): `Bio/motifs/thresholds.py` (`ScoreDistribution`, N. Dojer 2008), `Bio/motifs/jaspar/__init__.py` (`calculate_pseudocounts`), `Bio/motifs/matrix.py` (`search(both=True)`, `mean`, `std`, `distribution`) — B05 follow-up additions.
9. Hertz G.Z., Stormo G.D. 1999. Identifying DNA and protein patterns with statistically significant alignments of multiple sequences. Bioinformatics 15:563-577 (patser threshold, via Biopython `threshold_patser`).
