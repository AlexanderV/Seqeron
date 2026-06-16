# Validation Report: ONCO-VAF-001 — Variant Allele Frequency Analysis

- **Validated:** 2026-06-16   **Area:** Oncology
- **Canonical method(s):** `OncologyAnalyzer.CalculateVAF(int,int)`, `OncologyAnalyzer.CalculateVAFConfidenceInterval(int,int,double)`, `OncologyAnalyzer.AdjustVAFForPurity(double,double,double)`
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS

## Stage A — Description

### Sources opened & what they confirm

1. **Wikipedia — Binomial proportion confidence interval** (retrieved this session via WebFetch).
   Confirms the Wilson score interval verbatim:
   - center = (p̂ + z²/2n) / (1 + z²/n)
   - half-width = (z / (1 + z²/n)) · √( z²/4n² + p̂(1−p̂)/n )
   - interval = center ± half-width
   - The interval is contained in [0,1] and has **non-zero width** at p̂=0 or 1 (unlike Wald).
   - z = 1.96 for 95% (0.975 standard-normal quantile).
   These match the algorithm doc §2.2 and the implementation exactly.

2. **Tarabichi et al. 2017, PMC5538405** (WebFetch). Confirms the diploid-heterozygous clonal
   relationship: expected VAF = purity × 0.5, with the verbatim worked example "80% tumor cells …
   expected proportion of reads 0.5 × 0.8 = 0.4." This anchors INV-04 and M13.

3. **CNAqc, Genome Biology 2024 (10.1186/s13059-024-03170-5)** (WebSearch retrieval of the formula).
   Confirms the expected-VAF relation v = (m·ρ·π) / (2(1−π) + π·n_tot) with the normal contribution
   fixed at 2(1−π). For a clonal mutation (ρ=1) this is v = m·π/(2(1−π)+π·n_tot); the documented
   inversion adjusted = vaf·(2(1−π)+π·n_tot)/π recovers m·CCF. The paper's 2:1 (n_tot=3) peaks at
   33%/66% are consistent with the relation.

4. **GATK Mutect2 FAQ** (via Evidence; definition cross-checked). Empirical allele fraction =
   alt-supporting reads / total reads (AD-based), distinct from Mutect2's Bayesian `AF`. Matches
   `CalculateVAF`.

### Formula check
All three formulas reproduced character-for-character from the cited sources. The empirical VAF,
the Wilson center/margin, and the purity/ploidy inversion are all standard and correctly stated.

### Edge-case semantics check
- totalReads = 0 → VAF defined as 0 (no coverage); CI undefined (throws). Sourced (GATK / Wilson n>0).
- altReads > totalReads (alignment artifact, VAF>1) → invalid input. Sourced (GATK).
- p̂=0 → Wilson lower = 0, upper > 0; p̂=1 → upper = 1, lower < 1. Sourced (Wikipedia, no overshoot).
- purity = 0 → correction divides by π → undefined/throw. Sourced (CNAqc denominator).
All defined and sourced; none "implementation-defined".

### Independent cross-check (numbers)
Wilson intervals (z=1.96) were recomputed **two independent ways** this session:
(a) the closed form, and (b) by solving the Wilson score equation as a quadratic in p directly
(`(p̂−p)² = z²p(1−p)/n`) — a different algebraic route that does not share the closed-form blind spot.
Both routes agree to 10 decimals with the Evidence/TestSpec tables:

| alt/total | center | lower | upper |
|-----------|--------|-------|-------|
| 25/100 | 0.2592487019 | 0.1754509400 | 0.3430464637 |
| 50/100 | 0.5000000000 | 0.4038298286 | 0.5961701714 |
| 5/20   | 0.2902825314 | 0.1118600528 | 0.4687050100 |
| 0/10   | 0.1387700844 | 0.0000000000 | 0.2775401688 |
| 10/10  | 0.8612299156 | 0.7224598312 | 1.0000000000 |

`statsmodels.proportion_confint(method='wilson')` was also run; it reproduces the structure and
differs only in the 5th–6th decimal because it uses the exact quantile z=1.959964 rather than the
source-cited z=1.96 (e.g. 25/100 lower 0.1754521 vs 0.1754509) — confirming the formula is identical
and 1.96 is the right constant for the documented assumption.

Purity adjustment recomputed: (0.4,0.8,2)→1.0, (0.35,0.7,2)→1.0, (0.2,0.5,2)→0.8, (0.3,0.5,4)→1.8 —
all match.

### Findings / divergences
None. Description is mathematically and biologically correct and fully sourced.

## Stage B — Implementation

### Code path reviewed
- `OncologyAnalyzer.cs:344` `CalculateVAF` → delegates to private `CalculateVaf` (`:617`).
- `OncologyAnalyzer.cs:367-395` `CalculateVAFConfidenceInterval` (Wilson, clamp to [0,1]).
- `OncologyAnalyzer.cs:414-435` `AdjustVAFForPurity`.
- `OncologyAnalyzer.cs:595-611` `ZScoreFor` (z=1.96 only for 0.95).
- Constants `ZScore95 = 1.96` (`:58`), `NormalDiploidCopyNumber = 2.0` (`:68`).

### Formula realised correctly?
Yes. The Wilson `center`/`margin` lines are an exact transcription of the validated formula. The
purity correction computes `2(1−π) + π·ploidy` then `vaf·avg/π`. `CalculateVaf` returns
`totalReads==0 ? 0 : (double)alt/total` with the validation (negative, alt>total) throwing
`ArgumentOutOfRangeException`. All matches Stage A.

### Cross-verification table recomputed vs code
The full test suite executes the exact-value cases above against the actual code and passes;
hand/independent recomputation matches the code's outputs to within 1e-9.

### Variant/delegate consistency
`CalculateVAF` and the somatic-calling path share the same private `CalculateVaf`, so the empirical
ratio is consistent across the unit. No `*Fast`/instance variants for these methods.

### Test quality audit (HARD gate)
- **Sourced expectations, not code echoes:** All exact assertions (`Is.EqualTo(x).Within(1e-9/1e-10)`)
  use values derived from the external sources (Wilson 1927 formula, CNAqc/Tarabichi), independently
  reproduced this session. A wrong implementation (Wald interval, wrong z, alt/(alt+ref) mistakes)
  would fail at the 1e-9 tolerance. **PASS.**
- **No green-washing:** No weakened assertions, no widened tolerances, no skipped/ignored tests, no
  expected values tuned to actual output. M9's midpoint check is an additional symmetry assertion on
  top of exact Lower/Upper bounds, not a replacement for them. **PASS.**
- **Cover all logic:** All three public methods + all overload defaults exercised. Branches covered:
  valid ratios (M1–M4), zero coverage VAF (M5) and CI throw (S1), alt>total throw (M6), negative
  counts (M7), exact Wilson for 4 (alt,total) pairs (M8–M10), no-overshoot at p̂=0 (M11) and p̂=1
  (M12), confidence validation 0/1/1.5/0.90 (S1), purity adjustment M13–M15, zero-purity throw (M16),
  vaf/purity/ploidy validation (S2/S3), INV-04 round-trip over 4 purities, and an INV-01/02/03
  property sweep over all alt≤total for n=1..50. **PASS.**
- **Honest green:** Full unfiltered suite = **Failed: 0, Passed: 6621**; `dotnet build` 0 errors. The
  4 build warnings are pre-existing NUnit analyzer notes in unrelated test files
  (ApproximateMatcher, etc.), not in the VAF code/tests touched here. **PASS.**

### Findings / defects
None. No code or test change was required.

## Verdict & follow-ups

- **Stage A:** PASS — formulas and edge-case semantics independently confirmed against Wikipedia
  (Wilson 1927), Tarabichi 2017, CNAqc 2024, and GATK; all expected numbers traced to retrieved
  sources and reproduced two independent ways.
- **Stage B:** PASS — implementation faithfully realises the validated formulas; tests assert exact
  sourced values, cover every method/branch/edge case, and the full suite is honestly green.
- **Test-quality gate:** PASS.
- **End-state:** CLEAN — no defect found; algorithm fully functional.
- **Follow-ups:** none.
