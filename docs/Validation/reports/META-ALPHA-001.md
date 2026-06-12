# Validation Report: META-ALPHA-001 — Alpha Diversity Indices

- **Validated:** 2026-06-12   **Area:** Metagenomics
- **Canonical method(s):** `MetagenomicsAnalyzer.CalculateAlphaDiversity(IReadOnlyDictionary<string,double>)`
  (`src/Seqeron/Algorithms/Seqeron.Genomics.Metagenomics/MetagenomicsAnalyzer.cs:297`),
  with private helpers `CalculateShannonIndex` (:327), `CalculateSimpsonIndex` (:341), `CalculateChao1` (:362).
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS

## Stage A — Description

### Sources opened & what they confirm
- **Wikipedia — Diversity index** (https://en.wikipedia.org/wiki/Diversity_index), fetched 2026-06-12, confirms verbatim:
  - Shannon: `H' = −Σ_{i=1}^{R} p_i ln(p_i)`. Log base: Shannon discussed bases 2, 10, and e; all are popular. Ecology commonly uses natural log (e, "nats").
  - Simpson: `λ = Σ_{i=1}^{R} p_i²`.
  - Gini–Simpson: `1 − λ`.
  - Inverse Simpson: `1/λ = ²D` (effective number of species / true diversity of order 2).
- **Chao1** (WebSearch, cross-checked against MetricGate, scikit-bio, Statology, mothur): classic estimator
  `Ŝ_Chao1 = S_obs + f₁²/(2·f₂)` with f₁ = singletons (count exactly 1), f₂ = doubletons (count exactly 2);
  bias-corrected substitute `S_obs + f₁(f₁−1)/2` when f₂ = 0. (A separate "bias-corrected" general form uses the
  `f₂+1` denominator always; the spec deliberately selects the classic 1984 form — see Findings.)
- Pielou (1966): `J = H'/ln(S)`, with J = 0 by convention when S ≤ 1 (ln(1)=0 → undefined).

### Formula check (per spec vs sources)
| Index | Spec/Evidence formula | Source | Match |
|-------|-----------------------|--------|-------|
| Shannon | `H' = −Σ pᵢ ln(pᵢ)`, natural log | Shannon 1948 / Wikipedia | ✅ |
| Simpson | `λ = Σ pᵢ²` (raw D, NOT 1−D) | Simpson 1949 / Wikipedia | ✅ |
| Inverse Simpson | `1/λ` | Hill 1973 / Wikipedia | ✅ |
| Pielou | `J = H'/ln(S)`; 0 for S≤1 | Pielou 1966 | ✅ |
| Chao1 | `S_obs + f₁²/(2f₂)`; `S_obs + f₁(f₁−1)/2` if f₂=0 | Chao 1984 | ✅ |

**Log base:** spec/Evidence (Decision #2) explicitly select natural log. Internally consistent with the test
oracle — e.g. two equal species must give Shannon = ln(2), four equal species = ln(4), which hold only for base e.

**Simpson variant:** spec reports the *raw concentration index* `λ = Σ pᵢ²` ("Simpson"), and `1/λ` ("InverseSimpson").
It does NOT report Gini-Simpson (1−λ). This is consistent throughout (single species → Simpson = 1.0, four equal → 0.25).

### Edge-case semantics check
Empty/null → all 0; single species → Shannon 0, Simpson 1, InvSimpson 1, Pielou 0 (S≤1 guard); zero abundances
filtered (ln(0) undefined); abundances normalized before computation; Chao1 = S_obs for proportional (non-integer)
data and when f₁=0. All defined and sourced.

### Independent cross-check (hand computation)
Worked example counts [10,20,30,40] → total 100, p=[0.1,0.2,0.3,0.4]:
- Shannon = −(0.1·ln0.1 + 0.2·ln0.2 + 0.3·ln0.3 + 0.4·ln0.4) = 0.230259+0.321888+0.361192+0.366516 = **1.279854 nats** ✅
- Simpson λ = 0.01+0.04+0.09+0.16 = **0.30** ✅; Inverse Simpson = 1/0.30 = **3.3333** ✅
- Gini-Simpson would be 0.70 (not reported by this API — correctly so per spec).
- These are integer counts but f₁=0 (no species with count 1) → Chao1 = S_obs = **4** ✅
Chao1 worked examples: {50,30,1,1,2} → S_obs=5, f₁=2, f₂=1 → 5 + 4/2 = **7** ✅;
{100,1,1,1} → S_obs=4, f₁=3, f₂=0 → 4 + (3·2)/2 = **7** (bias-corrected) ✅.

### Findings / divergences
None. The "FAIL" traps flagged in the protocol are all absent: log base correctly natural log (not mislabeled);
Simpson correctly reported as raw D (claimed D, returns D); Chao1 main branch uses `2·f₂` (NOT the erroneous
`f₂+1`); singleton/doubleton counting is exact (count==1, count==2).

## Stage B — Implementation

### Code path reviewed
`MetagenomicsAnalyzer.cs:297–380`.

### Formula realised correctly? (evidence)
- **Shannon** (:327–339): normalizes `pi = p/sum`, accumulates `h -= pi*Math.Log(pi)` over positive values.
  `Math.Log` = natural log. ✅ Matches `−Σ pᵢ ln pᵢ`.
- **Simpson** (:341–353): `d += pi*pi` → returns raw `λ = Σ pᵢ²`. ✅ (NOT 1−λ, matching the reported variant.)
- **InverseSimpson** (:310): `simpson > 0 ? 1/simpson : 0`. ✅
- **Pielou** (:316): `observedSpecies > 1 ? shannon/Math.Log(observedSpecies) : 0`. ✅ ln(1)=0 guard present.
- **Chao1** (:362–380): detects integer data (tolerance 1e-9); non-integer → S_obs. Counts f₁ (≈1), f₂ (≈2).
  f₁=0 → S_obs. f₂>0 → `S_obs + f₁²/(2·f₂)`. f₂=0 → `S_obs + f₁(f₁−1)/2`. ✅ Exactly the classic 1984 form.
- **ObservedSpecies** (:303): count of positive abundances. ✅
- Empty/null and all-zero-after-filter → AlphaDiversity(0,0,0,0,0,0) (:299–306). ✅

### Cross-verification table recomputed vs code
| Input | Shannon | Simpson | InvSimpson | Pielou | Chao1 | ObsSp | Code match |
|-------|---------|---------|-----------|--------|-------|-------|------------|
| [10,20,30,40] | 1.279854 | 0.30 | 3.3333 | 0.92327 | 4 | 4 | ✅ |
| single {1.0} | 0 | 1.0 | 1.0 | 0 | 1 | 1 | ✅ (test) |
| {0.5,0.5} | ln2=0.69315 | 0.5 | 2.0 | 1.0 | 2 | 2 | ✅ (test) |
| {0.25×4} | ln4=1.38629 | 0.25 | 4.0 | 1.0 | 4 | 4 | ✅ (test) |
| {50,30,1,1,2} | — | — | — | — | 7 | 5 | ✅ (test) |
| {100,1,1,1} | — | — | — | — | 7 | 4 | ✅ (test) |

### Variant/delegate consistency
A second helper (line 268, `CalculateBiodiversityMetrics`) reuses the same `CalculateShannonIndex`/`CalculateSimpsonIndex`
— consistent. No `*Fast` variant for alpha diversity.

### Test quality audit
`tests/.../MetagenomicsAnalyzer_AlphaDiversity_Tests.cs` — 39 executed cases (20 [Test] + TestCase expansions).
Assertions check **exact** sourced values with tight tolerances (Ln2/Ln4 constants, 1e-10), cover M1–M22, S1–S4,
C1–C2, edge cases, and invariants. Not tautological. Chao1 singleton/doubleton/bias-corrected branches all exercised.

### Findings / defects
None.

## Verdict & follow-ups
- Stage A: PASS. Stage B: PASS. State: **CLEAN** — no defects; all formulas, log base (natural), Simpson variant
  (raw λ = Σpᵢ²), and Chao1 variant (classic Chao 1984 with f₂=0 bias-corrected special case) match authoritative
  sources and hand computations. No code or test changes required.
- Full suite: 4484 passed / 0 failed (baseline preserved). AlphaDiversity filter: 39 passed.
