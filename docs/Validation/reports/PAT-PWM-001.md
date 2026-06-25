# Validation Report: PAT-PWM-001 — Position Weight Matrix (PWM / PSSM)

- **Validated:** 2026-06-24   **Area:** Pattern Matching
- **Canonical method(s):** `MotifFinder.CreatePwm(IEnumerable<string>, double pseudocount = 0.25)`, `MotifFinder.ScanWithPwm(DnaSequence, PositionWeightMatrix, double threshold = 0.0)`; supporting type `PositionWeightMatrix` (`Matrix`, `Length`, `Consensus`, `MaxScore`, `MinScore`).
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS

This is an independent re-validation in a fresh context. The PWM model, formulas, log base, background, scoring, and consensus rule were re-derived from the authoritative sources and the Wikipedia worked example was recomputed cell-by-cell, rather than relying on the prior archived report.

---

## Stage A — Description

### Sources opened & what they confirm

1. **Wikipedia: Position weight matrix** (https://en.wikipedia.org/wiki/Position_weight_matrix) — fetched live.
   - PPM table for the 10-sequence / 9-position example confirmed **verbatim** (all 36 cells; see cross-check below).
   - Log-odds formula verbatim: `M_{k,j} = log₂( M_{k,j} / b_k )`, **log base 2**, background `b_k = 1/|k| = 0.25` for nucleotides.
   - Pseudocounts: "Pseudocounts (or Laplace estimators) are often applied when calculating PPMs if based on a small dataset, in order to avoid matrix entries having a value of 0." The article does **not** prescribe a specific pseudocount value.
   - Scoring verbatim: score "can be calculated by **adding (rather than multiplying) the relevant values at each position** in the PWM."

2. **Rosalind: Consensus and Profile (CONS)** — profile matrix = 4×n count matrix; consensus = "most common symbol at each position"; ties may yield multiple valid consensus strings (a deterministic tie-break is acceptable). Sample dataset (7 seqs len 8) → consensus `ATGCAACT`.

3. **Stormo (2000), Bioinformatics review; Nishida et al. (2008) NAR 37(3):939–944** — pseudocount is the Bayesian/Dirichlet correction applied to the **counts** (PFM), `(count + p)/(N + Σp)`; recommended values are entropy-dependent. The spec/code default `p = 0.25` (Jeffreys-like) is a legitimate member of this family and is documented.

### Formula check

```
freq(b,j)   = (count(b,j) + p) / (N + 4p)     # PPM with per-base pseudocount p
PWM(b,j)    = log2( freq(b,j) / 0.25 )         # log-odds, uniform background 0.25
score(seq)  = Σ_j PWM(seq[j], j)               # additive scoring
consensus[j]= argmax_b count(b,j)              # == argmax_b PWM(b,j) (log2 monotonic)
```

- Pseudocount applied to PFM counts (correct per Nishida/Stormo); default `p = 0.25`. ✅
- Log base **2**, background **0.25** uniform. ✅
- Scoring additive (sum of per-position log-odds). ✅
- Consensus = most-frequent base per column. Because `log2` is monotonic and the column normaliser `(N+4p)` and background `0.25` are constant within a column, `argmax` over the log-odds matrix equals `argmax` over counts — equivalent to Rosalind's count rule. ✅

### Edge-case semantics

- Empty → `ArgumentException`; null → `ArgumentNullException`; unequal lengths → `ArgumentException` (Wikipedia requires equal length). ✅
- `p = 0` with an unseen base → `log2(0) = −∞` (Wikipedia notes this risk; the pseudocount is the guard, by design). ✅
- All-same column, default `p` (N=10) → `log2(10.25/11/0.25) ≈ +1.898`; with `p=0` → `+2.0`. ✅
- Consensus ties → deterministic A<C<G<T (allowed by Rosalind). ✅
- Sequence shorter than PWM → no windows → empty result. ✅

### Independent cross-check (recomputed this session, Python)

PPM recomputed from the 10 sequences matches the Wikipedia table **exactly** (all 36 cells), e.g. row A = `[0.3,0.6,0.1,0,0,0.6,0.7,0.2,0.1]`, G@pos4 = `1.0`, T@pos5 = `1.0`.

| Quantity | Formula | Computed | Source/test |
|---|---|---|---|
| Consensus (10-seq Wikipedia ex.) | argmax counts | **TAGGTAAGT** | M18 ✓ |
| pos-4 G (all-G, N=10, p=.25) | log2(10.25/11/0.25) | **1.8981204** | M18 `expectedG4` ✓ |
| pos-4 unseen (N=10) | log2(0.25/11/0.25)=log2(1/11) | **−3.4594316** | ✓ |
| zero-pc all-G@pos4 / unseen | log2(10/10/0.25) / log2(0) | **2.0 / −∞** | M20 ✓ |
| uniform AAA (N=3) A / unseen | log2(3.25/4/0.25) / log2(0.25/4/0.25) | **1.7004397 / −2.0** | E1 ✓ |
| max-entropy ACGT (N=4) | log2(1.25/5/0.25) | **0.0** | E2 ✓ |
| Rosalind CONS consensus | argmax counts | **ATGCAACT** | M14 ✓ |

### Findings / divergences

None material. Wikipedia leaves the pseudocount value open; the spec/code choose `p = 0.25` (Bayesian Dirichlet form), sourced to the Nishida/Stormo family and documented (Contract C5). **Stage A: PASS.**

---

## Stage B — Implementation

### Code path reviewed

`src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/MotifFinder.cs`:
- `CreatePwm` (L213–273): null (L215) / empty (L218) / equal-length (L222) validation; strict A/C/G/T validation (L226–237, `ToUpperInvariant` at L217 → case-insensitive); PFM count (L243–259); PPM+log-odds `freq=(count+p)/(N+4p)` (L267) then `Math.Log2(freq/0.25)` (L268), background hardcoded 0.25 (L262).
- `PositionWeightMatrix.GenerateConsensus` (L774–797): per-column argmax over log-odds, first-max (A<C<G<T) tie-break — equivalent to count argmax.
- `MaxScore` (L802–816) / `MinScore` (L821–835): per-column max/min sums.
- `ScanWithPwm` (L282–328): windows `i = 0..len-motifLen` (L293, empty if too short); additive `score += pwm.Matrix[baseIndex, j]` (L315); `score >= threshold` inclusion (L319, `>=` semantics); non-ACGT in a window → `valid=false`, window skipped (L309–311). `Matrix[base, position]` index convention used consistently in construction and scanning.

### Formula realised correctly?

Yes — the construction (L262–270) and scanning (L315) realise the validated formula exactly, with log base 2 and background 0.25. No approximation. ✅

### Cross-verification table recomputed vs code

All values in the Stage-A table are asserted against the live code by the test suite (`Pwm_LogOddsFormula_NumericalVerification`, `CreatePwm_WikipediaExample_ConsensusAndScoresMatchSource` asserting `Consensus == "TAGGTAAGT"` and `Math.Log2(10.25/11.0/0.25)`, `CreatePwm_UniformSequence_ExactLogOddsValues`, `CreatePwm_MaximumEntropy_AllLogOddsZero`, `CreatePwm_ZeroPseudocount_ProducesNegativeInfinity`, `CreatePwm_RosalindCONS_TestCase` asserting `"ATGCAACT"`). All pass.

### Variant/delegate consistency

`ScanWithPwm` consumes the same `Matrix` produced by `CreatePwm`; `Consensus`/`MaxScore`/`MinScore` are derived properties of that matrix → internally consistent. The standalone `GenerateConsensus(IEnumerable<string>)` / `CreateConsensusFromAlignment` are separate IUPAC-degenerate features, not the PWM consensus and not relied on by this unit; out of scope.

### Test quality audit

The `~Pwm` filter ran **119 tests, all passing** (the prior report's 45 has grown via the fuzzing/algebraic/combinatorial campaigns). Assertions check exact sourced numeric values (`Math.Log2(...)`, `.Within(1e-10)`), named consensus strings, the `−∞` zero-pseudocount edge, the `>=` threshold boundary, overlapping matches, short-sequence empty result, and non-ACGT rejection. Tests are deterministic and substantive (not no-throw tautologies).

### Findings / defects

None. **Stage B: PASS.**

---

## Verdict & follow-ups

- **Stage A: PASS** — formulas, log base (2), background (0.25), additive scoring, consensus rule and tie-break match Wikipedia + Rosalind; pseudocount form sourced to Stormo/Nishida with documented default `p = 0.25`. Wikipedia PPM recomputed cell-by-cell with exact agreement.
- **Stage B: PASS** — implementation realises the validated formula exactly; all hand-computed cross-checks are reproduced by passing tests.
- **State: CLEAN** — no defects found; no code changes required.
- **Tests:** `~Pwm` filter → 119 passed / 0 failed. (Build: succeeded, 0 warnings.)
