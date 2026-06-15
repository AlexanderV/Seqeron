# Validation Report: SEQ-COMPLEX-WINDOW-001 — Windowed Sequence Complexity

- **Validated:** 2026-06-16   **Area:** Complexity
- **Canonical method(s):** `SequenceComplexity.CalculateWindowedComplexity(DnaSequence, int windowSize = 64, int stepSize = 10)` → `IEnumerable<ComplexityPoint>`
- **Stage A verdict:** PASS-WITH-NOTES
- **Stage B verdict:** PASS

## Stage A — Description

### Sources opened this session (retrieved live)
- **Wikipedia, "Entropy (information theory)"** — WebFetch https://en.wikipedia.org/wiki/Entropy_(information_theory)
  Confirmed: `H(X) = −Σ p(x) log_b p(x)`; base 2 ⇒ bits; maximum `log_b(n)` for `n` equiprobable outcomes ("uncertainty is maximal when all possible events are equiprobable"); deterministic distribution ⇒ `Η₁(1) = 0`; convention `0·log_b 0 = 0` (limit `p→0⁺ p log p = 0`). For DNA (n=4) max = log₂4 = 2.0 bits.
- **Wikipedia, "Linguistic sequence complexity"** — WebFetch https://en.wikipedia.org/wiki/Linguistic_sequence_complexity
  Confirmed: vocabulary usage `Uᵢ` = ratio of actual vocabulary size to maximal possible vocabulary size of a fragment of that length; `V_max` at level i = "either 4^i or N−i+1, whichever is smaller"; range "0<C<1". **The headline overall complexity is the PRODUCT `C = U₁U₂…U_w`** (per-word-length usages multiplied).
- **WebSearch** ("Troyanskaya 2002 … vocabulary usage U_i …" and "linguistic complexity DNA sum-of vocabularies ratio …") surfacing Gabrielian & Bolshoy (1999) / Troyanskaya et al. (2002) / Bolshoy "DNA sequence analysis linguistic tools" (PubMed 15130826) / the lzcomposer review. Confirmed **two distinct published forms** exist:
  1. **Product form** `C = ∏ Uᵢ` (Trifonov 1990; Wikipedia headline).
  2. **Summation / ratio-of-sums form** ("alphabet-capacity L-gram method … the sum of the observed range x_i from 1..L divided by the sum of the expected E for this sequence length") — i.e. `LC = (Σ Vᵢ)/(Σ V_max,i)`.

### Formula check
This unit is a **sliding-window driver** that, per fully-contained window, reports:
- Shannon entropy `H = −Σ_{b∈ACGT} p_b log₂ p_b` (bits) — **matches Shannon 1948 / Wikipedia exactly** (base 2, `0·log0=0`, max log₂4=2.0).
- Linguistic complexity in the **summation form** `LC = (Σ Vᵢ)/(Σ min(4^i, w−i+1))`, `m = min(6, w)`.
- Window geometry: starts `i = 0, s, 2s, …` with `i+w ≤ L`; count `= ⌊(L−w)/s⌋+1` for `L≥w` else 0. This is the standard sliding-window complexity-profile enumeration (Troyanskaya 2002).

### Edge-case semantics
- Uniform window ⇒ H = 2.0 (sourced, Shannon max for n=4).
- Homopolymer window ⇒ H = 0.0 (sourced, deterministic distribution).
- `L < w` ⇒ empty profile (no partial trailing window) — repository contract; consistent with "windows fully contained in the sequence".
- Null ⇒ ArgumentNullException; `windowSize<1`/`stepSize<1` ⇒ ArgumentOutOfRangeException — repository contract.

### Independent cross-check (numbers computed from scratch this session, Python)
| Window | H (computed) | Σ Vᵢ | Σ V_max | LC (computed) |
|--------|------|------|---------|---------------|
| `ACGTACGT` (w=8, m=6) | **2.0** | 23 | 29 | **23/29 = 0.7931034482758621** |
| `AAAAAAAA` (w=8, m=6) | **0.0** | 6 | 29 | **6/29 = 0.20689655172413793** |

Per-length distinct subwords for `ACGTACGT`: 4,4,4,4,4,3 (sum 23); `V_max,i = min(4^i,8−i+1)` = 4,7,6,5,4,3 (sum 29). All independently reproduced.
Window counts: `⌊(24−8)/8⌋+1 = 3`; `⌊(24−8)/4⌋+1 = 5`; `L=8,w=8 ⇒ 1`; `L=5,w=8 ⇒ 0`. All confirmed.

### Findings / divergences (NOTE)
The per-window LC uses the **summation form**, whereas Wikipedia's *headline* linguistic-complexity definition is the **product** `∏ Uᵢ`. The summation form is, however, an independently **published** variant (the alphabet-capacity / L-gram "sum of observed ÷ sum of expected" method) — confirmed by the web search above — and is the exact metric of the already-validated repository unit **SEQ-COMPLEX-001** (ledger #5: Stage A 🟡, CLEAN), to which this driver delegates verbatim. The divergence is therefore a documented, sourced choice rather than an error → **PASS-WITH-NOTES** (inherited from SEQ-COMPLEX-001). The `Position = WindowStart + windowSize/2` center label is a repository convention not mandated by any source; it does not affect any complexity value (documented assumption).

## Stage B — Implementation

### Code path reviewed
`src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/SequenceComplexity.cs`
- `CalculateWindowedComplexity` (L211–221): eager guards `ThrowIfNull` / `windowSize<1` / `stepSize<1`, then delegates to the core.
- `CalculateWindowedComplexityCore` (L223–241): `for (i=0; i+windowSize<=seq.Length; i+=stepSize)`; per window calls `CalculateShannonEntropyCore` and `CalculateLinguisticComplexityCore(window, Math.Min(6, windowSize))`; yields `ComplexityPoint(Position=i+windowSize/2, H, LC, WindowStart=i, WindowEnd=i+windowSize−1)`.
- `CalculateShannonEntropyCore` (L93–120): `−Σ (count/total)·log₂(count/total)` over A/C/G/T; skips count==0 (so `0·log0=0`). Matches Shannon 1948.
- `CalculateLinguisticComplexityCore` (L39–66): distinct-subword `HashSet` per length 1..min(m,len), `min(4^i, len−i+1)` denominator, ratio of sums. Matches summation-form LC.
- `WindowLcMaxWordLength = 6` (L199).

### Formula realised correctly? (evidence)
Yes. The enumeration condition `i + windowSize <= seq.Length` realises "fully-contained windows", giving exactly `⌊(L−w)/s⌋+1` points (verified by tracing the four count cases). Entropy and LC cores reproduce the externally-confirmed formulas; the `min(6, windowSize)` cap is applied (verified by the w=4 invariant case). Guards are **eager** (the public method `return`s the iterator rather than being an iterator itself), so the exception tests are valid regardless of enumeration.

### Cross-verification table recomputed vs code
| Case | Source value | Code (test) value | Match |
|------|--------------|-------------------|-------|
| `ACGTACGT` H | 2.0 | 2.0 (M3) | ✓ |
| `AAAAAAAA` H | 0.0 | 0.0 (M4) | ✓ |
| `ACGTACGT` LC | 23/29 | 23/29 (M5) | ✓ |
| `AAAAAAAA` LC | 6/29 | 6/29 (M6) | ✓ |
| count L24/w8/s8 | 3 | 3 (M1) | ✓ |
| starts/ends/centers | 0,8,16 / 7,15,23 / 4,12,20 | same (M2) | ✓ |
| count L24/w8/s4 | 5; starts 0,4,8,12,16 | same (S1) | ✓ |
| L=w=8 | 1 window, 0..7 | same (S2) | ✓ |
| L<w | empty | empty (M7) | ✓ |
| LC max attainable | 1.0 (w=1 / `ACGT` w=4) | within C1 bound | ✓ |

### Variant/delegate consistency
The driver delegates to the same `*Core` helpers used by the standalone `CalculateShannonEntropy` (SEQ-COMPLEX-*) and `CalculateLinguisticComplexity` (SEQ-COMPLEX-001), so per-window values equal the scalar metrics for those substrings by construction. No separate re-implementation to drift.

### Test quality audit (canonical file `SequenceComplexity_CalculateWindowedComplexity_Tests.cs`, 13 tests)
- **Sourced, exact values:** M1 (count 3), M2 (exact start/end/center arrays), M3 (2.0), M4 (0.0), M5 (`23.0/29.0`), M6 (`6.0/29.0`), S1 (count 5, starts 0,4,8,12,16), S2 (1 window, 0..7) — all assert exact values traced to the external Shannon/LC formulas, **not** code echoes. A deliberately-wrong driver (e.g. off-by-one window end, wrong LC denominator, partial-window emission) would fail these.
- **All public-method branches / Stage-A edges covered:** non-overlap (M1/M2), overlap s<w (S1), exact-fit L=w (S2), L<w empty (M7), null (M8), windowSize<1 (M9), stepSize<1 (M10), `min(6,w)` cap (C1 w=4 case), homopolymer & uniform metric extremes (M3–M6).
- **C1 bounds invariant** uses ranges (`0≤H≤log₂4`, `0<LC≤1`) — this is the *correct* use of a range (a genuine multi-window invariant where individual exact values are not separately known), not green-washing. The `≤1.0+1e-10` tolerance correctly accommodates the attainable LC=1.0 (verified w=1 and `ACGT` w=4 ⇒ exactly 1.0).
- Exception tests force enumeration with `.ToList()` (harmless given eager guards; correct).
- No duplicate windowed tests remain in `SequenceComplexityTests.cs` (only a pointer comment) — consolidation was performed as the TestSpec records.
- **No weakened assertions, no skips, no widened tolerances, no expected-value-to-match-output edits.**

### Findings / defects
**None.** No code defect, no test defect, no coverage gap. The algorithm faithfully realises the (Stage-A-validated, NOTE-qualified) description.

### Test-quality gate result
**PASS.** Assertions are exact sourced values where exact values exist; ranges used only for a genuine invariant; every public-method branch and Stage-A edge/error case is exercised; full unfiltered suite green.

## Verdict & follow-ups
- **Stage A: PASS-WITH-NOTES** — formulas/geometry correct and externally confirmed; the single NOTE is the summation-vs-product LC form, an independently published variant inherited from the already-validated SEQ-COMPLEX-001.
- **Stage B: PASS** — code realises the validated formulas exactly; tests are honest, exact, and complete.
- **End-state: ✅ CLEAN** — no defect found; nothing to fix.
- **Build/test:** `dotnet build` 0 errors (4 pre-existing warnings in unrelated `ApproximateMatcher_EditDistance_Tests.cs`, untouched); full unfiltered `dotnet test` = **6598 passed, 0 failed**.
- No code or test changes were made this session (none warranted). Working tree changes are limited to validation docs.
