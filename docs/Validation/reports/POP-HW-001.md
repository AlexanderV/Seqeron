# Validation Report: POP-HW-001 — Hardy–Weinberg Equilibrium Test

- **Validated:** 2026-06-12   **Area:** Population Genetics
- **Canonical method(s):** `PopulationGeneticsAnalyzer.TestHardyWeinberg(variantId, observedAA, observedAa, observedaa, significanceLevel = 0.05)`
  - Source: `src/Seqeron/Algorithms/Seqeron.Genomics.Population/PopulationGeneticsAnalyzer.cs:398`
  - Result type: `HardyWeinbergResult` (record struct, line 52)
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS

## Stage A — Description

### Sources opened & what they confirm
- **Wikipedia: Hardy–Weinberg principle** (fetched 2026-06-12) confirms verbatim:
  - Expected genotype frequencies: `f(AA) = p²`, `f(aa) = q²`, `f(Aa) = 2pq`.
  - χ² goodness-of-fit: `χ² = Σ (O − E)² / E`.
  - Degrees of freedom: "degrees of freedom for test for Hardy–Weinberg proportions are # genotypes − # alleles" → biallelic = 3 − 2 = **1 df**.
  - "The 5% significance level for 1 degree of freedom is 3.84."
  - Ford scarlet tiger moth worked example: observed (1469, 138, 5), N=1612, p=0.954, q=0.046; expected (1467.4, 141.2, 3.4); **χ² = 0.83** (does not reject HWE).
- **Hardy (1908), Weinberg (1908), Emigh (1980)** are the cited primary references for the p²/2pq/q² law and HWE chi-square testing; they are consistent with the Wikipedia synthesis used.

### Formula check
| Quantity | Source formula | Spec / Evidence |
|----------|----------------|-----------------|
| Allele freq | p = (2·n_AA + n_Aa) / 2N, q = 1−p | matches Evidence §2.2 |
| Expected counts | E = freq × N: p²N, 2pqN, q²N | matches Evidence §2.3 |
| Statistic | χ² = Σ (O−E)²/E over 3 genotype classes | matches Evidence §2.4 |
| df | #genotypes − #alleles = 3 − 2 = 1 | matches (df=1, NOT 2) ✓ |
| Critical value | 3.841 at α=0.05, df=1 | matches |

### Edge-case semantics check
- **N=0**: defined as χ²=0, PValue=1, InEquilibrium=true (no evidence against H₀). Sourced rationale in spec; defensible per hypothesis-testing framework.
- **E(genotype)=0** (monomorphic): skip that term in the χ² sum (Wikipedia formula computed only for E>0). Avoids div-by-zero.
- **Fixed allele (p=1 or q=1)**: observed = expected ⇒ χ²=0 ⇒ InEquilibrium=true.
- **Large deviation** (excess/deficit heterozygotes): χ² ≫ 3.841 ⇒ InEquilibrium=false.

### Independent cross-check (numbers)
Prompt worked example — observed AA=30, Aa=50, aa=20, N=100:
- p = (2·30 + 50)/200 = 110/200 = **0.55**, q = **0.45**
- E(AA)=0.55²·100=**30.25**, E(Aa)=2·0.55·0.45·100=**49.5**, E(aa)=0.45²·100=**20.25**
- χ² = (30−30.25)²/30.25 + (50−49.5)²/49.5 + (20−20.25)²/20.25
      = 0.002066 + 0.005051 + 0.003086 = **0.01020**
- df=1; χ² ≪ 3.841 ⇒ NOT significant, in HWE. ✓ Confirmed.

Ford moth example (independent recompute): p=3076/3224=0.95409, q=0.04591; E=(1467.40, 141.21, 3.40); χ² = 0.00174 + 0.07307 + 0.7560 = **0.8309** → Wikipedia's rounded 0.83. ✓

### Findings / divergences
None. Description is mathematically and biologically correct; df is correctly 1 (the classic df=2 mistake is avoided).

## Stage B — Implementation

### Code path reviewed
`PopulationGeneticsAnalyzer.cs:398-459` (TestHardyWeinberg + ChiSquareCDF), plus RegularizedGammaP series/continued-fraction (`:465-520`).

### Formula realised correctly? (evidence)
- Allele freqs (`:413-414`): `p = (2*observedAA + observedAa)/(2*n)`, `q = 1-p`. ✓
- Expected counts (`:417-419`): `p*p*n`, `2*p*q*n`, `q*q*n` = p²N, 2pqN, q²N. ✓
- χ² sum (`:422-429`): `Σ (O−E)²/E`, each term guarded by `expected > 0` (E=0 terms skipped). ✓
- df=1: PValue = `1 - ChiSquareCDF(chiSquare, 1)` (`:432`); CDF = `RegularizedGammaP(df/2, x/2)` (`:458`), the exact χ²(df) CDF. ✓
- InEquilibrium = `pValue >= significanceLevel`, default 0.05 (`:403, :444`). ✓
- N=0 guard (`:407-409`): returns χ²=0, PValue=1, InEquilibrium=true. ✓

The χ² denominator uses Expected (not Observed) — correct. No df=2, no wrong-df p-value.

### Cross-verification table recomputed vs code
| Case | Input (AA,Aa,aa) | Expected χ² | Code result | df/PValue |
|------|------------------|-------------|-------------|-----------|
| Worked example | 30,50,20 | 0.0102 | matches (manual trace) | df=1, ~0.92 |
| Ford moth | 1469,138,5 | 0.8309 | test asserts 0.8309 Within(0.01) ✓ | p>0.05 |
| Perfect HWE | 25,50,25 | 0 | test asserts 0 ✓ | p=1 |
| Excess het | 10,80,10 | 36.0 | test asserts 36.0 ✓ | p<0.05 |
| All het | 0,100,0 | 100.0 | test asserts 100.0 ✓ | p<0.05 |
| Deficit het | 45,10,45 | 64.0 | test asserts 64.0 ✓ | p<0.05 |
| Fixed major | 100,0,0 | 0 | E=(100,0,0), χ²=0 ✓ | InEq=true |
| Borderline | 46,49,5 | 3.169 | asserts ∈(2.706,3.841) ✓ | α-dependent |
| N=0 | 0,0,0 | 0 | guard returns p=1 ✓ | InEq=true |

### Variant/delegate consistency
Single canonical static method; no `*Fast` / instance variant for HWE. Consistent.

### Test quality audit
`PopulationGeneticsAnalyzer_HardyWeinberg_Tests.cs` — 47 test runs, all passing. Critical value constant `ChiSquareCriticalValue_Alpha05_Df1 = 3.841` (line 20). Tolerances tight (χ² Within(0.01), expected counts Within(0.1)). Invariant tests cover χ²≥0, PValue∈[0,1], expected counts sum to N, InEquilibrium⟺(PValue≥significance), symmetry, extreme/rare frequencies. Plus 2 FsCheck property tests and 1 snapshot test. Assertions check exact sourced values, not tautologies.

### Findings / defects
None.

## Verdict & follow-ups
- **Stage A: PASS** — formulas, df=1, and critical value confirmed against Wikipedia (Hardy 1908 / Weinberg 1908 / Emigh 1980 lineage); worked example hand-computed.
- **Stage B: PASS** — code faithfully implements p²/2pq/q² expected counts, χ²=Σ(O−E)²/E with E-in-denominator, df=1 p-value via exact regularized incomplete gamma, correct edge-case guards.
- **State: CLEAN** — no defects. Build succeeds; HardyWeinberg filter = 47/47 passed; full `Seqeron.Genomics.Tests` = 4484/4484 passed (baseline preserved). No code changes.
