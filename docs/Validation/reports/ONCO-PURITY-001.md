# Validation Report: ONCO-PURITY-001 — Tumor Purity Estimation

- **Validated:** 2026-06-16   **Area:** Oncology
- **Canonical method(s):** `OncologyAnalyzer.EstimatePurityFromVaf(double)`, `OncologyAnalyzer.EstimatePurityFromVAF(IEnumerable<VariantObservation>)`, `OncologyAnalyzer.EstimatePurity(IEnumerable<PurityVariant>)` (private helper `EstimatePurityFromAlleleSpecificVaf`)
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS

## Stage A — Description

### Sources opened & what they confirm (all retrieved this session)

1. **CNAqc vignette** (caravagnalab.github.io/CNAqc) — WebFetched. Expected-VAF of a clonal mutation, verbatim:
   `v_m(c) = mπc / [2(1−π) + π(n_A+n_B)]`. Variables: m = multiplicity, π = purity, c = cancer-cell fraction,
   n_A/n_B = allele-specific copy numbers. For copy-neutral diploid het (m=1, n_A+n_B=2) it reduces to
   `v = πc/2`; with c=1 (clonal), **v = π/2 ⇒ ρ = 2·VAF**.
2. **ABSOLUTE / Carter et al. 2012, Nat Biotechnol 30:413** — WebFetched the PMC full text (PMC4383288). Expected
   allelic fraction, verbatim: `F = α·s_q / D_s` with `D_s = α·q(x) + 2(1−α)`. α = purity, q = local tumour copy
   number, s_q = multiplicity, and **2(1−α) = contribution of contaminating normal diploid cells**. This is
   *algebraically identical* to the CNAqc relation (α↔π, s_q↔m, q↔n_tot) — an **independent peer-reviewed
   confirmation** of both the model and the 2(1−π) normal-diploid term.
3. **TPES, Bioinformatics 2019** — WebFetched. "The VAF distribution of a set of clonal monoallelic SNVs from
   pure tumor samples … should be centered in 0.5"; impure samples scale linearly ⇒ purity = 2·VAF for the
   copy-neutral clonal het case. Independent confirmation of the closed form.
4. **WebSearch corroboration** — All-FIT / TPES / biostars: "for a heterozygous mutation in a diploid cancer
   genome, in the absence of CNAs and normal contamination, CCF is twice the VAF"; general framework
   `VAF = p·CCF·n/N`.

### Formula check
- Forward model `v = m·π / [2(1−π) + π·n_tot]` — matches CNAqc (source 1) and ABSOLUTE (source 2) verbatim.
- Inversion `π = 2v / [m + v(2 − n_tot)]` — hand-derived and confirmed: from `v[2(1−π)+π·n_tot]=mπ` →
  `2v = π(m + v(2 − n_tot))`. Exact algebraic inverse (INV-3). ✔
- Diploid het special case `ρ = 2·VAF` (m=1, n_tot=2) — confirmed by three independent sources. ✔

### Edge-case semantics check
- VAF=0 ⇒ purity 0; VAF=0.5 (diploid) ⇒ purity 1.0; VAF>0.5 (diploid) ⇒ ρ>1 impossible → reject. All sourced
  to the closed form and the 0 ≤ purity ≤ 1 invariant.
- Empty input ⇒ purity undefined (CNAqc corner case); null ⇒ guard.
- Allele-specific domain: m≥1, n_tot≥1; (VAF,m,n_tot) yielding ρ∉[0,1] or non-positive denominator rejected
  (formula domain).

### Independent cross-check (numbers, all from sourced formulas)
- 2:1 segment (n_tot=3) at π=1: m=1 ⇒ v=1/[0+3]=**1/3**; m=2 ⇒ v=2/3=**2/3** (CNAqc "peaks at 33% and 66%"). ✔
- Diploid het, purity 60% ⇒ VAF 30%; band 55–65% ⇒ 27.5–32.5% (CNAqc worked example). ✔
- ABSOLUTE D_s = αq + 2(1−α): diploid het q=2, α=0.6, s_q=1 ⇒ F = 0.6/(1.2+0.8) = 0.6/2 = 0.30. ✔

### Findings / divergences
None. The description is biologically and mathematically correct and faithful to two independent peer-reviewed
models (CNAqc and ABSOLUTE). The VAF-only estimator's fixed m=1/n_tot=2 scope is an explicitly stated,
correctly-sourced modeling choice, not an error.

## Stage B — Implementation

### Code path reviewed
- `OncologyAnalyzer.cs:480-497` `EstimatePurityFromVaf`: `ρ = HeterozygousDiploidPurityFactor(2.0)·vaf`;
  rejects NaN/vaf∉[0,1]; rejects resulting ρ>1.
- `OncologyAnalyzer.cs:451-468` `EstimatePurityFromVAF`: per-variant VAF via `CalculateVAF`, ρ=2·VAF, median.
- `OncologyAnalyzer.cs:515-575` `EstimatePurity` + `EstimatePurityFromAlleleSpecificVaf`:
  `π = 2·v / [m + v·(2 − n_tot)]` with `NormalDiploidCopyNumber = 2.0`; guards VAF∉[0,1], m<1, n_tot<1,
  denominator≤0, purity∉[0,1]; median aggregation.
- `Median` (`:578-585`): sort copy, lower-mid average for even counts. Deterministic, non-mutating.

### Formula realised correctly?
Yes — verbatim. `2.0*vaf` and `2.0*v/(m + v*(2 − n_tot))` are exactly the sourced closed form and its inverse.

### Cross-verification table recomputed vs code
| Case | Input | Formula value | Code (test) | Match |
|------|-------|---------------|-------------|-------|
| M1/M5 | VAF 0.30 / 0.275 / 0.50 / 0.00 | 0.60 / 0.55 / 1.0 / 0.0 | same | ✔ |
| M4 | VAFs {0.10,0.15,0.30} | median{0.20,0.30,0.60}=0.30 | 0.30 | ✔ |
| S1 | VAF 0.02 | 0.04 | 0.04 | ✔ |
| M6 | v0.30,m1,n2 | 0.60 | 0.60 | ✔ |
| M7 | v1/3,m1,n3 | 1.0 | 1.0 | ✔ |
| M8 | v2/3,m2,n3 | 1.0 | 1.0 | ✔ |
| M9 | v1/6,m1,n4 | (1/3)/(2/3)=0.5 | 0.5 | ✔ |
| S2 | mixed | median{0.60,1.0,0.5}=0.60 | 0.60 | ✔ |
| M10 | VAF 0.6 diploid; (0.6,1,2) | ρ=1.2>1 → throw | throws | ✔ |
| M10b (new) | (0.9,1,4) | denom −0.8≤0 → throw | throws | ✔ |
| M11–M14 | empty/null/n0/m0/vaf∉[0,1] | throw | throws | ✔ |

### Variant/delegate consistency
`EstimatePurity` at (m=1,n_tot=2) reduces to ρ=2·VAF, agreeing with `EstimatePurityFromVaf`/`EstimatePurityFromVAF`
(M6 vs M1). Collection overloads delegate to the single-variant helpers + median. Consistent.

### Numerical robustness
Closed-form arithmetic; no overflow on the stated ranges. Denominator-≤0 and purity∉[0,1] explicitly guarded.
NaN VAF rejected. Median uses a copy (no mutation, order-independent — C1).

### Test quality audit (HARD gate)
- **Sourced, not code echoes:** every expected value is hand-derived from the CNAqc/ABSOLUTE formulas (above), not
  read off the implementation. A wrong-constant implementation would fail M1/M6/M7/M8/M9.
- **No green-washing:** all assertions are exact `Is.EqualTo(…).Within(1e-10)` (1e-12 for determinism). No
  Greater/AtLeast/Contains/range where an exact value is known; no widened tolerance; no skip/ignore.
- **Cover all logic:** all three public methods exercised; both VAF-only overloads; allele-specific diploid,
  2:1 peaks, and a general (m1,n4) inversion; median for both overloads; determinism; and the error branches
  null/empty/VAF-out-of-range/VAF>0.5/m<1/n_tot<1/purity>1.
- **Gap found & fixed:** the documented "non-positive denominator" rejection (spec §3.3) was the one Stage-A
  error branch with no test — M10's allele-specific case hits the distinct purity>1 branch. Added
  `EstimatePurity_VafUnreachableForCopyState_NonPositiveDenominator_Throws` (vaf 0.9, m1, n_tot4) with the input
  sourced from the forward VAF range v∈[0,0.25] for that copy state. Fixture 21→22 tests.
- **Honest green:** full unfiltered suite `Failed: 0, Passed: 6661` (one pre-existing unrelated benchmark skip);
  `dotnet build` 0 errors. The 4 NUnit2007 warnings are pre-existing in unrelated `ApproximateMatcher_EditDistance_Tests.cs`;
  the changed test file is warning-free.

### Findings / defects
No algorithm or description defect. One test-coverage gap (non-positive-denominator branch), fixed this session
with a source-derived test. Logged as FINDINGS_REGISTER §A51.

## Verdict & follow-ups
- **Stage A: PASS** — model confirmed verbatim against two independent peer-reviewed sources (CNAqc, ABSOLUTE)
  plus TPES; all worked numbers reproduced.
- **Stage B: PASS** — code realises the closed form and its exact inverse; all cross-check values match; tests
  are exact and sourced. Coverage gap closed.
- **End-state: CLEAN.** No code or description change; one source-derived edge-case test added. Build + full
  unfiltered suite green (6661 passed, 0 failed).
