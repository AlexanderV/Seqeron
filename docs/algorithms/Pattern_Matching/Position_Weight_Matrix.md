# Position Weight Matrix (PWM)

| Field | Value |
|-------|-------|
| Algorithm Group | Pattern Matching |
| Test Unit ID | PAT-PWM-001 |
| Related Projects | N/A |
| Implementation Status | Complete |
| Last Reviewed | 2026-10-01 |

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

P-value options, Markov backgrounds and K rows (audit round 2, group G1, B05 F32; [MotifFinder.PwmPValue.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/MotifFinder.PwmPValue.cs)):

- `PwmPValueOptions` (record; `Default`, `Exact`): `InitialGranularity` (0.1), `MaxGranularity` (null), `DecreaseFactor` (10) — the TFM-Pvalue driver parameters (`testScoreToPvalue` / `testPvalueToScore`: `initialGranularity`, `maxGranularity`, `decrgr`; granularity = rounding step 1/g, `PwmPValueResult.Granularity` reports g) — plus `MaxStates` (2²¹ DP states per column), `MaxSuffixSet` (2²⁰) and `Exhaustive`: when the refinement leaves the value undecided, the band is enumerated once more without a state limit (forward: at the finest scale whose integer DP fitted; inverse: the certified bracket lowScore < t, words with integer score ≥ ⌈g·highScore⌉ + pad + 1 counted as a block), so the result is always exact. Overloads `PwmScorePValue(pwm, score, background, options)` / `PwmScoreThresholdForPValue(pwm, p, background, options)`; null/`Default` = the 3-argument overloads bit for bit (712 random cases: all scores, p-values and bounds identical). The 20-column case that is not exact by default (bounds 0.058402542608746444 … 0.05840680768687889) is exact with `MaxStates = 2²³` or `Exact`: 0.058404839602189895 = meet-in-the-middle enumeration of all 4²⁰ words.
- `MotifFinder.PwmMarkovScorePValue(pwm, score, OligoBackgroundModel, options?)` / `PwmMarkovScoreThresholdForPValue` — the same exact engine under an order-m Markov background: the word probability is RSAT `segment_proba` of an `OligoBackgroundModel.MarkovFromOligoFrequencies` table (F21), P(w) = P(w[0..m−1])·∏ P(w_c | w[c−m..c−1]) (prefix marginals for the first m letters, also for L < m). The DP state packs (integer partial score, last min(j, m) letters) — as MACRO-APE's dinucleotide-background distribution keeps one score map per previous letter (`DiPWMScoresGenerator.recalc_score_hash`) and AhoPro's Markov-background DP (Boeva et al. 2007) — and the band enumeration keeps (state, exact partial window score). Words of probability 0 (ψ = 0 tables) are excluded, also from the extreme scores; for a table where some context has no successor the word measure has total mass < 1 and the look-ahead multiplies settled prefixes by the exact completion mass (P(S ≥ min) = that total). `Equiprobable` / `Bernoulli` give the i.i.d. engine (identical results); input-estimated models (`BernoulliFromInput`, `MarkovFromInput`) and `Lexicon` throw — they need sequences, RSAT's input Markov estimate is not a normalised word distribution (m-mer and (m+1)-mer frequencies from different window sets) and the lexicon is a segmentation frequency, so build the chain with `MarkovFromOligoFrequencies` from the (m+1)-mer counts. Order ≤ 10. Scores stay the log-odds window scores of `CalculatePwmScores` (the matrix is scored as built; only the null distribution changes), so no Markov-aware scoring is needed.
- `MotifFinder.AlphabetPwmScorePValue(AlphabetPositionWeightMatrix, score, background?, options?)` / `AlphabetPwmScoreThresholdForPValue` — the same engine with K = |alphabet| rows (one shared engine: the DNA PWM is K = 4), i.i.d. background in alphabet order. `AlphabetPositionWeightMatrix.ScoreDistribution(background?, precision)` — Biopython `pssm.distribution` for any alphabet (its `ScoreDistribution` iterates `pssm[:, position].items()`, the PSSM's own alphabet); the grid DP is the DNA one with K rows.
- Inverse fix found by the Markov oracle: the bracket's upper end highScore (P(S ≥ highScore) ≤ p) does not bound t from above when no word scores in [highScore, t); the scan now widens the window upwards when every decided score is infeasible, and decides every band score when no word reaches the block above (y = +∞). Previously such cases fell back to certified bounds; exact results are unchanged (712-case comparison: values identical, 11 resolved one scale earlier).
- Cross-check: Markov — 3,399 random cases (order 0–3, L 1–8, tables with zeros, ψ ∈ {0, 0.01, 0.3}, strand-insensitive tables; forward + inverse) = exhaustive enumeration weighting every word by its Markov probability (scores identical, p ≤ 4.3e-15 relative; order 0 incl. zero residues: 1,500 cases ≤ 1.6e-13); MACRO-APE 3.0.6 (`ru.autosome.ape.di.FindPvalue --from-mono -d 16 -b <16 dinucleotide frequencies>`, dyadic matrices so its ceil discretisation is exact, circular-sequence dinucleotide tables so its second-letter initial distribution equals RSAT's first-letter one): 199 thresholds / 40 matrices ≤ 3.2e-15. K rows — 289 random cases (protein, RNA, gapped DNA, sub-alphabets; L 1–7) with Biopython 1.88's own PSSM: `threshold_fpr/fnr/balanced/patser` 1,156/1,156 bit-identical; exact p-values and thresholds = enumeration of all K^L words; protein example (20⁶ words) locked in tests.

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
- `PositionWeightMatrix` is the DNA PWM (four rows A, C, G, T, reverse complement, both-strand scan, Markov-background p-values); arbitrary alphabets (protein, RNA, gapped DNA) use `AlphabetPositionWeightMatrix` (§5.1, generic-alphabet members), which shares the log-odds, mean/std, window-scoring, score-distribution and exact p-value kernels.

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

DNA alphabet for `PositionWeightMatrix` (any alphabet: `AlphabetPositionWeightMatrix` — forward scan only, as Biopython; reverse complement and both-strand search are DNA notions). Score calibration for both: exact p-values and p-value thresholds (`PwmScorePValue` / `PwmScoreThresholdForPValue`, `AlphabetPwmScorePValue` / `AlphabetPwmScoreThresholdForPValue`, Touzet & Varré 2007) and the Biopython discretised score distribution (`ScoreDistribution`: FPR/FNR/balanced/patser thresholds — a grid approximation; e.g. Wikipedia PWM `threshold_fpr(0.01)` = 4.028388324862519 vs the exact threshold 4.028050165603465). DNA p-values also under an order-m Markov background (`PwmMarkovScorePValue`, RSAT frequency tables, m ≤ 10); the Biopython grid distribution is i.i.d. by definition (Biopython has no Markov variant). Exponential worst case: the p-value problem is NP-hard (Touzet & Varré 2007), so for long, degenerate matrices at moderate thresholds the exact value may be unreachable within the default budget and certified bounds are returned (`IsExact = false`); a larger `PwmPValueOptions.MaxStates` or `PwmPValueOptions.Exact` resolves it at the cost of time and memory (e.g. the 20-column inverse at p = 0.05: exact in ~28 s / 3 GB).

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
11. Boeva V., Clément J., Régnier M., Roytberg M.A., Makeev V.J. (2007). "Exact p-value calculation for heterotypic clusters of regulatory motifs and its application in computational annotation of cis-regulatory modules." *Algorithms Mol Biol* 2:13. doi:10.1186/1748-7188-2-13 (AhoPro; Markov-background p-values — publisher/PMC pages proxy-blocked, bibliographic record via WebSearch).
12. Vorontsov I.E., Kulakovskiy I.V., Makeev V.J. (2013). "Jaccard index based similarity measure to compare transcription factor binding site models." *Algorithms Mol Biol* 8:23 (MACRO-APE). Source opened: raw.githubusercontent.com/autosome-ru/macro-perfectos-ape master `src/main/java/ru/autosome/commons/backgroundModel/di/DiBackground.java`, `ape/calculation/ScoringModelDistributions/DiPWMScoresGenerator.java`, `commons/motifModel/di/DiPWM.java` (`fromPWM`), `commons/model/Discretizer.java`, `ape/calculation/findPvalue/FindPvalueAPE.java`; manual LaTeX (VorontsovIE/macro-perfectos-ape-manual `sections/formats/dibackground-format.tex`, `sections/discretization-strategy.tex`); `releases/ape-3.0.6.jar` run as the oracle.
