# Validation Report: PAT-PWM-001 — Position Weight Matrix (PWM/PSSM)

- **Validated:** 2026-06-12   **Area:** Pattern Matching
- **Canonical method(s):** `MotifFinder.CreatePwm(IEnumerable<string>, double pseudocount = 0.25)`, `MotifFinder.ScanWithPwm(DnaSequence, PositionWeightMatrix, double threshold = 0.0)`; supporting type `PositionWeightMatrix` (Matrix, Length, Consensus, MaxScore, MinScore).
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS

---

## Stage A — Description

### Sources opened & what they confirm

1. **Wikipedia: Position weight matrix** (https://en.wikipedia.org/wiki/Position_weight_matrix) — fetched.
   - Construction process PFM → PPM → PWM confirmed verbatim: PFM = "counting the occurrences of each nucleotide at each position"; PPM = "dividing that former nucleotide count at each position by the number of sequences"; PWM elements = "calculated as log odds".
   - Log-odds formula verbatim: `M_{k,j} = log₂( M_{k,j} / b_k )`, log base **2**, background `b_k = 1/|k| = 0.25 for nucleotides`.
   - Pseudocounts: "Pseudocounts (or *Laplace estimators*) are often applied when calculating PPMs if based on a small dataset, in order to avoid matrix entries having a value of 0." The article does **not** mandate a specific pseudocount value (neither +1 nor any other) — it only recommends applying one.
   - Scoring verbatim: score "can be calculated by **adding (rather than multiplying) the relevant values at each position** in the PWM."
   - 10-sequence worked example and its PPM table confirmed exactly as quoted in the spec.

2. **Rosalind: Consensus and Profile (CONS)** (https://rosalind.info/problems/cons/) — fetched.
   - Profile matrix = 4×n count matrix; consensus = "taking the most common symbol at each position".
   - **Tie-breaking:** "there may be more than one most common symbol, leading to multiple possible consensus strings" — returning **any one** is acceptable. So a deterministic tie-break (e.g. A<C<G<T) is valid.
   - Sample dataset (7 seqs, length 8) and expected consensus `ATGCAACT` confirmed.

3. **Nishida et al. (2008/2009), Pseudocounts for transcription factor binding sites**, NAR 37(3):939–944 — the pseudocount is applied to the **counts** (PFM), not to probabilities, via the Bayesian/Dirichlet form `(count + p)/(N + Σp)`; recommended values are entropy-dependent (commonly ~0.8 for TFBS). The spec's default `p = 0.25` is a legitimate (Jeffreys-like) choice within this family.

### Formula check

The spec (Contract C5) and implementation use the Bayesian-pseudocount PPM combined with the Wikipedia log-odds:

```
freq(b,j) = (count(b,j) + p) / (N + 4p)          # PPM with per-base pseudocount p
PWM(b,j)  = log2( freq(b,j) / 0.25 )              # log-odds, background 0.25
score(seq) = Σ_j PWM(seq[j], j)                   # additive scoring
consensus[j] = argmax_b count(b,j)                # == argmax_b PWM(b,j) (log2 monotonic)
```

- Pseudocount: applied to PFM counts (correct per Nishida/Wikipedia); default `p = 0.25`. ✅
- Log base: **2** (`Math.Log2`). ✅
- Background: **0.25** uniform, hardcoded. ✅
- Scoring: additive sum of per-position log-odds. ✅
- Consensus: most-frequent base per column. Because `log2` is monotonic and the column normaliser `(N+4p)` and background `0.25` are constant within a column, `argmax` over log-odds equals `argmax` over counts — so deriving consensus from the log-odds matrix is equivalent to Rosalind's count-based rule. ✅

### Edge-case semantics

- Empty alignment → ArgumentException; null → ArgumentNullException; unequal lengths → ArgumentException (Wikipedia requires equal length). ✅
- `p = 0` with an unseen base → `log2(0)` = −∞ (Wikipedia notes this risk). ✅ defined, not guarded (by design — pseudocount is the guard).
- All-same base column with default p → log-odds ≈ +1.898 (N=10) / +2.0 (p=0). ✅
- Consensus ties: deterministic A<C<G<T (allowed by Rosalind). ✅

### Independent cross-check (hand-computed, Python)

| Quantity | Formula input | Computed | Spec/test expected |
|---|---|---|---|
| Wikipedia pos-4 G (all-G, N=10) | log2((10.25/11)/0.25) | **1.898120** | `Math.Log2(10.25/11/0.25)` ✓ |
| Wikipedia pos-4 unseen | log2((0.25/11)/0.25) | **−3.459432** = log2(1/11) | matches ✓ |
| "AG" A@0 (N=1,p=.25) | log2((1.25/2)/0.25) | **1.321928** = log2(2.5) | matches ✓ |
| "AG" C@0 | log2((0.25/2)/0.25) | **−1.0** = log2(0.5) | matches ✓ |
| Uniform AAAA×3 (N=3) A | log2((3.25/4)/0.25) | **1.700440** = log2(3.25) | matches ✓ |
| Uniform unseen | log2((0.25/4)/0.25) | **−2.0** | matches ✓ |
| Max-entropy ACGT (N=4) | log2((1.25/5)/0.25) | **0.0** | matches ✓ |
| Zero-pc "A": A / C | log2(1/0.25) / log2(0/0.25) | **2.0 / −∞** | matches ✓ |
| Rosalind CONS consensus | argmax counts | **ATGCAACT** | matches ✓ |

Note: the test comment labels the Wikipedia G-value "≈1.89849", but both code and test compute `Math.Log2(10.25/11.0/0.25)` = 1.898120 identically, so the assertion is exact and correct; only the human-readable comment rounds differently. No defect.

### Findings / divergences

None material. The only nuance is that Wikipedia leaves the pseudocount value open; the spec/code pick `p = 0.25` (Bayesian Dirichlet form), which is sourced (Nishida family) and documented. Stage A: **PASS**.

---

## Stage B — Implementation

### Code path reviewed

`src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/MotifFinder.cs`:
- `CreatePwm` (L213–273): null/empty/length validation (L215–223), strict A/C/G/T validation (L226–237), PFM count (L243–259), PPM+log-odds with pseudocount/background (L262–270) — matches the validated formula exactly.
- `PositionWeightMatrix.GenerateConsensus` (L639–662): per-column argmax over log-odds with first-max (A<C<G<T) tie-break — equivalent to Rosalind count-based consensus.
- `MaxScore`/`MinScore` (L667–700): per-column max/min sums.
- `ScanWithPwm` (L282–328): additive sum `score += pwm.Matrix[baseIndex, j]` (L315), `score >= threshold` inclusion (L319, `>=` semantics confirmed), windows `i = 0..len-motifLen` (no matches if sequence shorter), non-ACGT window characters cause the window to be skipped (`valid = false`).

### Formula realised correctly?

Yes. `freq = (matrix[b,i] + pseudocount) / (count + 4*pseudocount)` (L267) then `Math.Log2(freq / 0.25)` (L268). Background hardcoded 0.25 (L262). Scoring additive (L315). Index convention Matrix[base, position] is used consistently in both construction and scanning. ✅

### Cross-verification table recomputed vs code

All nine hand-computed values above are asserted by the test suite against the live code (`Pwm_LogOddsFormula_NumericalVerification`, `CreatePwm_WikipediaExample_*`, `CreatePwm_UniformSequence_ExactLogOddsValues`, `CreatePwm_MaximumEntropy_AllLogOddsZero`, `CreatePwm_ZeroPseudocount_ProducesNegativeInfinity`, `CreatePwm_RosalindCONS_TestCase`). All pass.

### Variant/delegate consistency

`ScanWithPwm` consumes the same `Matrix` produced by `CreatePwm`; `Consensus`/`MaxScore`/`MinScore` are derived properties of that matrix — internally consistent. The standalone `GenerateConsensus(IEnumerable<string>)` (IUPAC degenerate consensus, L339) is a *different* feature (>25% threshold, IUPAC codes) and is not the PWM consensus; not in scope for PAT-PWM-001 and not relied on by it.

### Test quality audit

31 canonical + smoke/property/snapshot tests; the `~Pwm` filter executed **45 tests, all pass**. Assertions check exact sourced numeric values (`.Within(1e-10)`), edge cases (−∞ at p=0, threshold `>=` boundary, overlapping matches, short sequence, non-ACGT rejection), and invariants. Tests are deterministic and substantive (not no-throw tautologies).

### Findings / defects

None. Stage B: **PASS**.

---

## Verdict & follow-ups

- **Stage A: PASS** — formulas, log base (2), background (0.25), additive scoring, consensus rule, and tie-break all match Wikipedia + Rosalind; pseudocount form sourced to the Nishida family with documented default p=0.25.
- **Stage B: PASS** — implementation realises the validated formula exactly; all 9 hand-computed cross-checks reproduced by passing tests.
- **State: CLEAN** — no defects found; no code changes required.
- **Tests:** `~Pwm` filter → 45 passed / 0 failed. Full `Seqeron.Genomics.Tests` suite → 4461 passed / 0 failed.
