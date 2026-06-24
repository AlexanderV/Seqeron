# Validation Report: POP-HW-001 — Hardy–Weinberg Equilibrium Test

- **Validated:** 2026-06-24   **Area:** Population Genetics
- **Canonical method(s):** `PopulationGeneticsAnalyzer.TestHardyWeinberg(variantId, observedAA, observedAa, observedaa, significanceLevel = 0.05)`
  - Source: `src/Seqeron/Algorithms/Seqeron.Genomics.Population/PopulationGeneticsAnalyzer.cs:424`
  - Result type: `HardyWeinbergResult` (readonly record struct, line 52)
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS

## Stage A — Description

### Sources opened & what they confirm
- **Wikipedia: Hardy–Weinberg principle** (fetched 2026-06-24) confirms verbatim:
  - Expected genotype frequencies: `f(AA) = p²`, `f(Aa) = 2pq`, `f(aa) = q²`.
  - χ² goodness-of-fit: `χ² = Σ (O − E)² / E`.
  - Degrees of freedom: "1 degree of freedom (degrees of freedom for test for Hardy–Weinberg
    proportions are # genotypes − # alleles)" → biallelic = 3 − 2 = **1 df**.
  - "The 5% significance level for 1 degree of freedom is 3.84."
  - Ford scarlet tiger moth worked example: observed (1469, 138, 5), N=1612, p=0.954, q=0.046;
    expected (1467.4, 141.2, 3.4); **χ² = 0.83** (does not reject HWE).
- **Hardy (1908), Weinberg (1908), Emigh (1980)** are the cited primary references for the
  p²/2pq/q² law and HWE chi-square testing; consistent with the Wikipedia synthesis.

### Formula check
| Quantity | Source formula | Spec / Code |
|----------|----------------|-------------|
| Allele freq | p = (2·n_AA + n_Aa) / 2N, q = 1−p | matches |
| Expected counts | E = freq × N: p²N, 2pqN, q²N | matches |
| Statistic | χ² = Σ (O−E)²/E over 3 genotype classes (E in denominator) | matches |
| df | #genotypes − #alleles = 3 − 2 = **1** (NOT 2) | matches ✓ |
| Critical value | 3.841 at α=0.05, df=1 | matches |

### Edge-case semantics check
- **N=0**: defined as χ²=0, PValue=1, InEquilibrium=true (no data ⇒ no evidence against H₀).
  Defensible per hypothesis-testing framework; sourced in spec.
- **E(genotype)=0** (monomorphic): skip that term in the χ² sum (Wikipedia formula computed
  only for E>0) — avoids div-by-zero.
- **Fixed allele (p=1 or q=1)**: observed = expected ⇒ χ²=0 ⇒ InEquilibrium=true.
- **Large deviation** (excess/deficit het): χ² ≫ 3.841 ⇒ InEquilibrium=false.

### Independent cross-check (numbers — hand-computed this session)
Ford moth: p = (2·1469 + 138)/(2·1612) = 3076/3224 = **0.954094**, q = **0.045906**.
- E(AA) = 0.954094²·1612 = **1467.40**, E(Aa) = 2·0.954094·0.045906·1612 = **141.21**,
  E(aa) = 0.045906²·1612 = **3.397**.
- χ² = (1469−1467.40)²/1467.40 + (138−141.21)²/141.21 + (5−3.397)²/3.397
      = 0.001744 + 0.072975 + 0.756 ≈ **0.831** → Wikipedia's rounded 0.83. ✓ (matches to ~1e-3)

Second hand example (AA=30, Aa=50, aa=20, N=100): p=0.55, q=0.45; E=(30.25, 49.5, 20.25);
χ² = 0.00207 + 0.00505 + 0.00309 = **0.0102** ≪ 3.841 ⇒ in HWE. ✓

### Findings / divergences
None. Description is mathematically and biologically correct; the classic df=2 mistake is avoided.

## Stage B — Implementation

### Code path reviewed
`PopulationGeneticsAnalyzer.cs:424–471` (TestHardyWeinberg) + `:473–485` (ChiSquareCDF) +
`:491–531` (RegularizedGammaP series / continued fraction). Line numbers shifted vs the prior
report (398→424) due to upstream commits; logic unchanged.

### Formula realised correctly? (evidence)
- Allele freqs (`:439–440`): `p = (2*observedAA + observedAa)/(2*n)`, `q = 1-p`. ✓
- Expected counts (`:443–445`): `p*p*n`, `2*p*q*n`, `q*q*n` = p²N, 2pqN, q²N. ✓
- χ² sum (`:448–455`): `Σ (O−E)²/E`, each term guarded by `expected > 0` (E=0 terms skipped). ✓
- df=1: PValue = `1 - ChiSquareCDF(chiSquare, 1)` (`:458`); CDF = `RegularizedGammaP(df/2, x/2)`
  (`:484`), the exact χ²(df) CDF. ✓
- InEquilibrium = `pValue >= significanceLevel`, default 0.05 (`:429, :470`). ✓
- N=0 guard (`:433–436`): returns χ²=0, PValue=1, InEquilibrium=true. ✓

χ² denominator uses Expected (not Observed) — correct. No df=2, no wrong-df p-value.

### Cross-verification table recomputed vs code
| Case | Input (AA,Aa,aa) | Expected χ² | Test assertion |
|------|------------------|-------------|----------------|
| Ford moth | 1469,138,5 | 0.8309 | asserts 0.8309 Within(0.01) ✓ |
| Perfect HWE | 25,50,25 | 0 | asserts 0 ✓ |
| Excess het | 10,80,10 | 36.0 | asserts 36.0 Within(0.01), >3.841 ✓ |
| All het | 0,100,0 | 100.0 | asserts 100.0 ✓ |
| Deficit het | 45,10,45 | 64.0 | asserts 64.0 ✓ |
| Fixed major | 100,0,0 | 0 | E=(100,0,0), χ²=0 ✓ |
| Borderline | 46,49,5 | 3.169 | asserts ∈(2.706,3.841) ✓ |
| N=0 | 0,0,0 | 0 | guard returns p=1, InEq=true ✓ |

### Variant/delegate consistency
Single canonical static method; no `*Fast` / instance variant for HWE. Consistent.

### Test quality audit
`PopulationGeneticsAnalyzer_HardyWeinberg_Tests.cs` — **60 test runs**, all passing (grew from
47 since the prior report). Critical-value constant `ChiSquareCriticalValue_Alpha05_Df1 = 3.841`
(line 20). Tolerances tight (χ² Within(0.01), expected counts Within(0.1)). Invariant tests cover
χ²≥0, PValue∈[0,1], expected counts sum to N, InEquilibrium⟺(PValue≥significance), symmetry under
allele relabeling, extreme/rare frequencies. Plus FsCheck property tests (χ²≥0, PValue∈[0,1]) and
a snapshot regression guard. Assertions check exact sourced values, not tautologies.

### Findings / defects
None.

## Verdict & follow-ups
- **Stage A: PASS** — p²/2pq/q² expected counts, χ²=Σ(O−E)²/E, df=1 (3 genotypes − 2 alleles),
  and 3.841 critical value confirmed against Wikipedia (Hardy 1908 / Weinberg 1908 / Emigh 1980
  lineage); Ford moth example hand-recomputed to ~1e-3.
- **Stage B: PASS** — code faithfully implements expected counts, χ²=Σ(O−E)²/E with E-in-
  denominator and E>0 guard, df=1 p-value via exact regularized incomplete gamma, correct
  N=0 / monomorphic edge guards. Tests lock the sourced values with tight tolerances.
- **State: CLEAN** — no defects. Build succeeds; HardyWeinberg filter = 60/60 passed. No code
  changes (only line-number drift vs prior report).
