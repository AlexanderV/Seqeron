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

Additive members (audit group D, [MotifFinder.PwmPValue.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/MotifFinder.PwmPValue.cs), [MotifFinder.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/MotifFinder.cs)):

- `MotifFinder.PwmScorePValue(pwm, score, background?)` → `PwmPValueResult(Score, PValue, PValueLowerBound, PValueUpperBound, IsExact, Granularity)` — **exact** P(S ≥ score) for an i.i.d. background word, S being the window score exactly as `CalculatePwmScores` computes it (left-to-right double sum). Method: Touzet & Varré (2007) successive refinement (TFM-Pvalue): for g = 10, 100, … the matrix is rounded down, M = ⌊g·W⌋, so g·S − M ∈ [0, E], E = Σ_j max_b (g·W − M) over **all** columns (TFM-Pvalue's `computesIntegerMatrix` starts this sum at column 1); a sparse DP over integer partial scores with look-ahead pruning (TFM-Pvalue `fastPvalue`) gives the mass certainly above (M ≥ ⌈gα⌉ + 2) and the undecided band (gα − E − 2 ≤ M < ⌈gα⌉ + 2; ±2 absorbs floating-point rounding). The band words are then enumerated exactly — a layered DP over (integer partial score, exact double partial window score) states, i.e. branch-and-bound with memoisation (prefixes with equal state have identical completions), pruned by the exact sets of reachable integer suffix scores of the last columns (meet-in-the-middle) — and their exact scores compared with α as soon as the band fits the budget (≤ 2²¹ states per layer), otherwise g is refined. This also resolves scores that a word attains exactly — TFM-Pvalue itself does not converge there and returns 0 for the consensus score. Not resolved by g·Σ max|W| ≤ 2⁴⁰ → certified bounds, `IsExact = false`, `PValue` = upper (conservative) bound. Finite matrices only (pseudocount > 0), as TFM-Pvalue / `ScoreDistribution`.
- `MotifFinder.PwmScoreThresholdForPValue(pwm, pValue, background?)` — exact inverse (TFM-Pvalue `pv2sc`): the smallest word score t with P(S ≥ t) ≤ pValue and its exact p-value (the largest achievable p-value ≤ pValue). The integer distribution at scale g (restricted, as TFM-Pvalue's PvalueToScore loop, to the window certified at the previous scale) brackets t between the last bucket whose tail exceeds pValue and the next one; the bracket (± E) is enumerated with exact scores and accumulated from the top. On an 8-column random matrix TFM-Pvalue returned P = 9.923731321091039e-05 for p = 10⁻⁴ where a larger achievable p-value 9.93067885403891e-05 ≤ 10⁻⁴ exists (confirmed by enumerating all 4^L words; also 2 more cases) — its stop rule (bucket gap > E) does not guarantee the largest p-value; here it is. TFM-Pvalue returns the same p-value but a rounded score ((α − offset)/g, e.g. 7.11 vs 7.150252343998389 for the Wikipedia PWM at 10⁻³). pValue ≥ 1 → minimum word score, p = 1; P(S ≥ max) > pValue → `Score = +∞`, p = 0.
- `MotifFinder.CreatePwm(sequences, IReadOnlyList<double> pseudocounts, background?)` — per-base pseudocounts (Biopython `motifs.create(seqs).counts.normalize(pseudocounts={...}).log_odds(background)`) and `MotifFinder.CreatePwmWithJasparPseudocounts(sequences, background?)` — JASPAR pseudocounts √N·q[b] against the same background (Biopython `jaspar.calculate_pseudocounts` with `m.background`); both build the shared count matrix and delegate to the `FromCounts` kernel.

Generic-alphabet members (audit group D, part 2, [MotifFinder.AlphabetPwm.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/MotifFinder.AlphabetPwm.cs)):

- `MotifFinder.CreateAlphabetPwm(sequences, string alphabet, double pseudocount = 0 | IReadOnlyList<double> pseudocounts, background?, ignoreUnknownSymbols = false)` → `AlphabetPositionWeightMatrix` — Biopython `motifs.create(instances, alphabet).counts.normalize(pseudocounts).log_odds(background)` for any alphabet (e.g. `ACDEFGHIKLMNPQRSTVWY`): K = |alphabet| rows in alphabet order. Symbols are matched case-insensitively (alphabet symbols must be distinct ignoring case). Symbols outside the alphabet throw, or with `ignoreUnknownSymbols` are not counted exactly as Biopython (`alignment.frequencies` keeps only alphabet letters, so that column's total is smaller); a column with zero total and zero pseudocounts throws (Biopython: `ZeroDivisionError`). Background: one finite, positive value per symbol, normalised (Biopython also accepts 0 → ±∞/NaN cells; not supported).
- `AlphabetPositionWeightMatrix`: `Alphabet`, `GetMatrix()`, indexer `[symbol, position]`, `Consensus` / `Anticonsensus` (strict improvement from ∓∞ → first symbol in alphabet order on ties, Biopython), `MaxScore` / `MinScore` (Σ column max/min, Biopython `max`/`min`), `Mean(background?)` / `Std(background?)` (Biopython `mean`/`std`, NaN/−∞ cells skipped), `FromCounts(alphabet, counts, pseudocount(s), background?)`.
- `MotifFinder.CalculateAlphabetPwmScores(string, AlphabetPositionWeightMatrix)` — Σ_j W[s[i+j], j] per window, NaN for a window holding a symbol outside the alphabet (the `_pwm.c` rule; Biopython 1.88 `calculate` raises "Use only with DNA motifs" for any other alphabet, so for protein the oracle is the per-letter sum of Biopython's PSSM); `MotifFinder.ScanWithAlphabetPwm(string, pwm, threshold)` — Biopython `search(seq, threshold, both=False)` (ascending positions, NaN never matches).
- Shared core: `LogOddsFromCounts` (row count from the matrix), `NormalizeBackground(background, size, …)`, `PwmMean` / `PwmStd` and the generic `ScorePwmWindow<TRowIndex>` (struct row mapper: `AcgtRowIndex` for the DNA PWM, a case-folded table for any alphabet) are used by both `PositionWeightMatrix` and `AlphabetPositionWeightMatrix`; the DNA results are unchanged (same arithmetic, same order).
- Cross-check: 400 random cases (protein, random sub-alphabets K = 2…20, DNA; L 1–12, N 1–25; scalar / per-symbol pseudocounts; random backgrounds; instances with X/− under `ignoreUnknownSymbols`; scanned sequences with lower case, X, B, Z, *, −) vs Biopython 1.88: matrices, max/min, mean/std to ≤ 1.1e-14 relative, consensus / anticonsensus / NaN windows / hits identical; for alphabet ACGT (92 cases) the scores equal Biopython `calculate` in float32 exactly and the hits equal `search(both=False)` (Biopython's own `search` fails when n = m — NumPy 0-d array — so those were compared via `calculate`).

### 5.2 Current Behavior

`CreatePwm(...)` uppercases all training sequences and uses a default pseudocount of `0.25`; the two-argument overload computes log-odds against a uniform background (0.25), the background overload against a normalised arbitrary background. `ScanWithPwm(...)` reports matches whose score is greater than or equal to the threshold and uses the PWM consensus as the `Pattern` field in returned `MotifMatch` values.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- PWM construction from aligned sequences.
- Pseudocount smoothing before log-odds conversion.
- Window scoring by summing per-position PWM values.

**Scope choices (not simplifications):**

- Pseudocounts: a scalar added to every cell (Biopython float `pseudocounts`), per-base pseudocounts (`CreatePwm(sequences, pseudocounts[], background)`), or the background-proportional JASPAR pseudocounts √N·b_k of Wasserman & Sandelin 2004 for any background (`CreatePwmWithJasparPseudocounts`, `JasparPseudocounts` + `FromCounts`).
- `PositionWeightMatrix` is the DNA PWM (four rows A, C, G, T, reverse complement, both-strand scan, p-values); arbitrary alphabets (protein, RNA, gapped DNA) use `AlphabetPositionWeightMatrix` (§5.1, generic-alphabet members), which shares the log-odds, mean/std and window-scoring kernels.

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

DNA alphabet for `PositionWeightMatrix` (any alphabet: `AlphabetPositionWeightMatrix` — forward scan only; reverse complement, both-strand search and p-values are DNA-specific). Score calibration: exact p-values and p-value thresholds (`PwmScorePValue` / `PwmScoreThresholdForPValue`, Touzet & Varré 2007) and the Biopython discretised score distribution (`ScoreDistribution`: FPR/FNR/balanced/patser thresholds — a grid approximation; e.g. Wikipedia PWM `threshold_fpr(0.01)` = 4.028388324862519 vs the exact threshold 4.028050165603465). Both assume an i.i.d. (order-0) background; Markov backgrounds are not modelled. Exponential worst case: the p-value problem is NP-hard (Touzet & Varré 2007), so for long, degenerate matrices at moderate thresholds the exact value may be unreachable within the budget and certified bounds are returned (`IsExact = false`).

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
7. Biopython 1.88 `Bio.motifs.matrix` (`normalize`, `log_odds`, `max`/`min`, `consensus`, `anticonsensus`, `reverse_complement`, `calculate`, `search`, `mean`, `std`) and `Bio.motifs.__init__` (`create(instances, alphabet)`, `Motif.__init__` counting via `Alignment.frequencies`), `Bio/motifs/_pwm.c` (NaN for symbols outside ACGT) — reference implementation used to lock the test values (installed package source opened).
8. Biopython 1.88 source (installed package): `Bio/motifs/thresholds.py` (`ScoreDistribution`, N. Dojer 2008), `Bio/motifs/jaspar/__init__.py` (`calculate_pseudocounts`), `Bio/motifs/matrix.py` (`search(both=True)`, `mean`, `std`, `distribution`) — B05 follow-up additions.
9. Hertz G.Z., Stormo G.D. 1999. Identifying DNA and protein patterns with statistically significant alignments of multiple sequences. Bioinformatics 15:563-577 (patser threshold, via Biopython `threshold_patser`).
10. Touzet H., Varré J.-S. (2007). "Efficient and accurate P-value computation for Position Weight Matrices." *Algorithms Mol Biol* 2:15. doi:10.1186/1748-7188-2-15. Reference C++ source opened: CRAN `TFMPvalue` 1.0.0 `src/Matrix.cpp`, `src/Matrix.h`, `src/TFMpvalue.cpp` (raw.githubusercontent.com/cran/TFMPvalue); pytfmpval 0.2.1 sdist `tfmp.py` (score2pval / pval2score driver loops).
