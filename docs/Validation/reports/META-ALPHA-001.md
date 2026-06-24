# Validation Report: META-ALPHA-001 — Alpha Diversity Indices

- **Validated:** 2026-06-24   **Area:** Metagenomics
- **Canonical method(s):** `MetagenomicsAnalyzer.CalculateAlphaDiversity(IReadOnlyDictionary<string,double>)`
  (`src/Seqeron/Algorithms/Seqeron.Genomics.Metagenomics/MetagenomicsAnalyzer.cs:480`),
  with private helpers `CalculateShannonIndex` (:510), `CalculateSimpsonIndex` (:524), `CalculateChao1` (:545).
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS

## Stage A — Description

### Sources opened & what they confirm
- **Wikipedia — Diversity index** (fetched 2026-06-24) confirms:
  - Shannon: `H' = −Σ_{i=1}^{R} p_i ln(p_i)`; log base is free, ecology commonly uses natural log (base e, "nats").
  - Simpson: `λ = Σ_{i=1}^{R} p_i²` (probability two random draws are the same type).
  - Gini–Simpson: `1 − λ` (probability they differ) — NOT reported by this API.
  - Inverse Simpson: `1/λ = ²D` (effective number of equally-abundant types / true diversity of order 2).
- **Chao1** (WebSearch 2026-06-24, cross-checked vs MetricGate, scikit-bio, Statology, mothur): classic estimator
  `Ŝ_Chao1 = S_obs + f₁²/(2·f₂)`, f₁ = singletons (count = 1), f₂ = doubletons (count = 2); when f₂ = 0 the
  bias-corrected substitute `S_obs + f₁(f₁−1)/2` is used. (A separate always-on bias-corrected form uses
  `2(f₂+1)` in the denominator; the spec deliberately selects the classic 1984 branch — internally consistent.)
- **Pielou (1966):** `J = H'/ln(S)`, with J = 0 by convention when S ≤ 1 (ln 1 = 0 → undefined).

### Formula check (per spec vs sources)
| Index | Spec formula | Source | Match |
|-------|--------------|--------|-------|
| Shannon | `H' = −Σ pᵢ ln pᵢ`, natural log | Shannon 1948 / Wikipedia | ✅ |
| Simpson | `λ = Σ pᵢ²` (raw D, NOT 1−D) | Simpson 1949 / Wikipedia | ✅ |
| Inverse Simpson | `1/λ` | Hill 1973 / Wikipedia | ✅ |
| Pielou | `J = H'/ln S`; 0 for S ≤ 1 | Pielou 1966 | ✅ |
| Chao1 | `S_obs + f₁²/(2f₂)`; `S_obs + f₁(f₁−1)/2` if f₂=0 | Chao 1984 | ✅ |

**Log base:** spec Decision #2 selects natural log; consistent with the oracle (two equal → ln 2, four equal → ln 4),
which hold only for base e. **Simpson variant:** raw concentration index λ = Σpᵢ² ("Simpson") and 1/λ
("InverseSimpson") — Gini–Simpson (1−λ) is deliberately NOT reported, consistent throughout.

### Edge-case semantics check
Empty/null → all 0; single species → Shannon 0, Simpson 1, InvSimpson 1, Pielou 0 (S ≤ 1 guard); zero abundances
filtered (ln 0 undefined); abundances normalized before computation; Chao1 = S_obs for proportional (non-integer)
data and when f₁ = 0. All defined and sourced.

### Independent cross-check (hand computation, fresh from prior report)
Counts [4,3,2,1] → total 10, p = [0.4,0.3,0.2,0.1]:
- Shannon = 1.27985 nats ✅; Simpson λ = 0.30 ✅; Inverse Simpson = 3.33333 ✅; Pielou J = H/ln4 = 0.92322 ✅.
- f₁ = 1 (one count = 1), f₂ = 1 (one count = 2), S_obs = 4 → Chao1 = 4 + 1²/(2·1) = **4.5** ✅ (verifies the f₂>0 branch off the test fixtures).
- Gini–Simpson would be 0.70 (correctly NOT exposed by this API).

### Findings / divergences
None. The protocol "FAIL traps" are all absent: log base is natural log (Math.Log), Simpson is raw D (claimed D,
returns D), Chao1 main branch uses `2·f₂` (not the erroneous `f₂+1`), singleton/doubleton counting is exact.

## Stage B — Implementation

### Code path reviewed
`MetagenomicsAnalyzer.cs:480–563` (line numbers shifted from prior report due to upstream additions; logic unchanged).

### Formula realised correctly? (evidence)
- **Shannon** (:510–522): `pi = p/sum`; `h -= pi*Math.Log(pi)` over positive values. `Math.Log` = natural log. ✅
- **Simpson** (:524–536): `d += pi*pi` → raw `λ = Σ pᵢ²`. ✅ (NOT 1−λ.)
- **InverseSimpson** (:493): `simpson > 0 ? 1/simpson : 0`. ✅
- **Pielou** (:499): `observedSpecies > 1 ? shannon/Math.Log(observedSpecies) : 0`. ✅ S ≤ 1 guard present.
- **Chao1** (:545–563): integer-data detection (tol 1e-9) → non-integer returns S_obs; counts f₁ (≈1), f₂ (≈2);
  f₁=0 → S_obs; f₂>0 → `S_obs + f₁²/(2·f₂)`; f₂=0 → `S_obs + f₁(f₁−1)/2`. ✅ Exactly the classic 1984 form.
- **ObservedSpecies** (:486): count of positive abundances. ✅ Empty/null/all-zero → AlphaDiversity(0,…,0) (:482–489). ✅

### Cross-verification table recomputed vs code
| Input | Shannon | Simpson | InvSimpson | Pielou | Chao1 | ObsSp | Match |
|-------|---------|---------|-----------|--------|-------|-------|-------|
| counts [4,3,2,1] | 1.27985 | 0.30 | 3.33333 | 0.92322 | 4.5 | 4 | ✅ (hand) |
| {0.5,0.5} | ln2 | 0.5 | 2.0 | 1.0 | 2 | 2 | ✅ (test) |
| {0.25×4} | ln4 | 0.25 | 4.0 | 1.0 | 4 | 4 | ✅ (test) |
| {50,30,1,1,2} | — | — | — | — | 7 | 5 | ✅ (test) |
| {100,1,1,1} | — | — | — | — | 7 | 4 | ✅ (test, bias-corrected) |
| {50,30,20} | — | — | — | — | 3 | 3 | ✅ (test, f₁=0) |

### Variant/delegate consistency
`CalculateBiodiversityMetrics` (:451) reuses the same Shannon/Simpson helpers — consistent. No `*Fast` variant.

### Test quality audit
`MetagenomicsAnalyzer_AlphaDiversity_Tests.cs` — 28 executed cases. Assertions check **exact** sourced values with
tight tolerances (Ln2/Ln4 constants, 1e-10), cover M1–M22, S1–S4, C1–C2, edge cases and invariants. Not tautological;
Chao1 singleton/doubleton/bias-corrected/proportional branches all exercised. Filtered run: 28 passed / 0 failed.

### Findings / defects
None.

## Verdict & follow-ups
- Stage A: PASS. Stage B: PASS. State: **CLEAN** — all formulas, log base (natural), Simpson variant (raw λ = Σpᵢ²),
  and Chao1 variant (classic Chao 1984 with f₂=0 bias-corrected special case) match authoritative sources and fresh
  hand computations. No code or test changes required. Implementation unchanged since prior validation (only line
  numbers shifted from upstream metagenomics additions).
- Alpha-diversity test filter: 28 passed / 0 failed. Project builds clean (0 warnings, 0 errors).
