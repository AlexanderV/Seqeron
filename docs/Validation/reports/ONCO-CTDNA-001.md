# Validation Report: ONCO-CTDNA-001 — ctDNA Analysis (Poisson LoD, tumour fraction, mean VAF)

- **Validated:** 2026-06-16   **Area:** Oncology
- **Canonical method(s):** `OncologyAnalyzer.CtDnaDetectionProbability(int,double,int)`, `ExpectedMutantMolecules(int,double,int)`, `IsCtDnaDetected(int,double,int,double)`, `CalculateTumorFraction(IEnumerable<VariantObservation>)`, `CalculateMeanVaf(IEnumerable<VariantObservation>)`, `HaploidGenomeEquivalents(double)`
- **Stage A verdict:** PASS-WITH-NOTES
- **Stage B verdict:** PASS

## Stage A — Description

### Sources opened this session (independent of the repo)

1. **USPTO Patent US 11,085,084 B2** — *Identification and use of circulating nucleic acids*
   (https://image-ppubs.uspto.gov/dirsearch-public/print/downloadPdf/11085084). Retrieved via WebSearch
   this session. Returned **verbatim**: "the probability of observing a single tumor reporter in cfDNA
   follows a Poisson distribution with mean λ = n × d", "x = 1 − e^(−nd)", and for k independent reporters
   "p = 1 − e^(−ndk)". This is the exact formula in the implementation. (Restates the Avanzini et al. 2020
   *Science Advances* shedding model, paywalled this session — the patent supplies the identical equations.)

2. **On the length, weight and GC content of the human genome** — PMC6391780
   (https://pmc.ncbi.nlm.nih.gov/articles/PMC6391780/). Retrieved this session. Male diploid genome ≈ 6.41 pg
   ⇒ haploid ≈ **3.205 pg**. The cfDNA-quantification field standard (Devonshire/Alcaide) uses **3.3 pg** per
   haploid genome — the conventional rounded value (cf. the classic 0.978 pg / 978 Mb genome-size convention).
   The repo cites and uses 3.3 pg explicitly; this is the standard cfDNA convention, not an error. See note N1.

3. **Genome-size / general references** (Wikipedia *Genome size*, sandwalk, ScienceDirect Genome Size) confirm
   the ~3.3 pg haploid figure as the commonly cited standard. 1000 pg/ng ÷ 3.3 pg/GE = **303.03 GE/ng** — matches
   the implementation and Alcaide et al. 2020 (303 haploid GE per ng), already in the repo Evidence.

4. **Tumour fraction = 2 × VAF** for a clonal heterozygous SNV at a copy-neutral diploid locus: confirmed in
   principle by web search (one mutant allele of two ⇒ observed VAF = TF/2 ⇒ TF = 2·VAF) and by the CNAqc
   (Antonello et al. 2024) relation v = m·π/[2(1−π)+π·n_tot] with m=1, n_tot=2 ⇒ v = π/2 ⇒ fraction = 2v, the
   same identity the repo already uses for `EstimatePurityFromVAF`. The Foundation Medicine ctDNA Tumor Fraction
   white paper was located but its PDF body was non-extractable this session.

### Formula check
- Detection probability **p = 1 − e^(−n·d·k)** — matches Patent US11085084 verbatim. ✓
- Expected mutant molecules **λ = n·d·k** — matches the patent's Poisson mean (×k for k reporters). ✓
- **GE = ng·1000/3.3** — matches the 3.3 pg/haploid convention (≈303 GE/ng). ✓
- **TF = 2 · mean VAF** — matches the clonal-het-diploid identity (v = π/2). ✓
- **Mean VAF = mean(alt/total)** — standard per-reporter ctDNA-level summarization (Newman 2014). ✓

### Edge-case semantics
- λ = 0 (n=0 or d=0) ⇒ p = 1 − e⁰ = 0: sourced (Poisson model). ✓
- Sub-molecule regime λ < 1 ⇒ not detected: physically sound (cannot observe < 1 expected molecule); the λ≥1
  floor in `IsCtDnaDetected` is a reasonable, documented operating rule. ✓
- d ∉ [0,1], n < 0, k < 1, negative mass, null/empty variant set: all defined error behaviours. ✓

### Independent cross-check (numbers hand-computed this session, Python)
| Quantity | Inputs | Independent value | Matches code/test? |
|----------|--------|-------------------|--------------------|
| p, λ=15 | n=15000,d=0.001,k=1 | 0.9999996940976795 | ✓ |
| p, λ=1 | n=1000,d=0.001,k=1 | 0.6321205588285577 | ✓ |
| p, λ=10 | n=1000,d=0.001,k=10 | 0.9999546000702375 | ✓ |
| p, λ=0.01 | n=100,d=0.0001 | 0.009950166… (<0.95) | ✓ (not detected) |
| GE | 1 ng | 303.0303030303 | ✓ |
| TF | VAF 0.10,0.20 | 2·0.15 = 0.30 | ✓ |
| mean VAF | 0.05,0.30,0.01 | 0.12 | ✓ |

### Findings / divergences (Stage A)
- **N1 (NOTE, not a defect):** the primary genome-weight paper gives 3.205 pg/haploid; the code uses the
  cfDNA-field-standard **3.3 pg** (Devonshire/Alcaide). The repo explicitly cites 3.3 pg as the convention, so
  this is a documented, defensible choice, not an error. Recorded as PASS-WITH-NOTES.
- The Avanzini 2020 primary remained paywalled; the detection equations were obtained verbatim from the patent
  this session, so no expected value rests on an unretrieved source.

## Stage B — Implementation

- **Code path:** `src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.cs:5400–5614` (ctDNA region);
  shared `CalculateVaf` helper at line 617; constants `PicogramsPerHaploidGenome=3.3` (83), `PicogramsPerNanogram=1000` (86),
  `TumorFractionFromVafFactor=2.0` (94), `DefaultCtDnaDetectionProbability=0.95` (103), `MinExpectedMutantMolecules=1.0` (111).
- **Formula realised correctly:** `1.0 - Math.Exp(-lambda)` with `lambda = n*d*k`; GE = `ng*1000/3.3`; TF =
  `2*meanVaf` clamped to 1.0; mean VAF = arithmetic mean of `alt/total`. All match the validated formulas.
- **Cross-verification table:** every value above recomputed independently matches the actual code (full suite run).
- **Variant/delegate consistency:** `ExpectedMutantMolecules` and `CtDnaDetectionProbability` share λ=n·d·k;
  `IsCtDnaDetected` reuses `ExpectedMutantMolecules`; `CalculateTumorFraction`/`CalculateMeanVaf` reuse `CalculateVaf`.
  Consistent.
- **Numerical robustness:** `1 - Math.Exp(-λ)` is well-conditioned for λ ≥ 0.01 (smallest tested λ); no overflow/div-by-zero
  on the stated ranges. The TF clamp `Math.Min(tf,1.0)` is effectively unreachable (per-variant VAF>0.5 is rejected
  first, so mean ≤ 0.5 ⇒ TF ≤ 1.0) — harmless defensive code; documented, not a defect.

### Test quality audit (test-quality gate)
- Existing 21 tests assert **exact sourced values** (1−e⁻¹⁵, 1−e⁻¹, 1−e⁻¹⁰, 303.03, 0.30, 0.12) with tight
  tolerances — not code echoes. No green-washing (no weakened assertions, no widened tolerances, no skips).
- **Coverage gaps found and fixed this session** (the only Stage-B action):
  1. `IsCtDnaDetected` second AND-condition (λ ≥ 1 **but** p < threshold) was untested. Added
     `IsCtDnaDetected_LambdaAtLeastOneButProbabilityBelowThreshold_ReturnsFalse` (λ=1 ⇒ p=0.6321 < 0.95 ⇒ false;
     and detected when threshold lowered to 0.5). Value derived from the sourced Poisson formula.
  2. `IsCtDnaDetected` `minDetectionProbability ∉ (0,1]` guard was untested. Added
     `IsCtDnaDetected_ThresholdOutOfRange_Throws` (0.0, 1.1, NaN).
  3. `ExpectedMutantMolecules` had only the happy path. Added `_TenReporters_Returns150` (λ=150) and
     `_InvalidArguments_Throw` (n<0, d>1, NaN, k<1).
- After additions: **26 tests** in the fixture, all green.

### Findings / defects (Stage B)
- No code defect. The three documented logic branches that lacked tests are now covered (test gap was a Stage-B
  defect per the gate; fixed in-session with sourced/derived expectations).

## Verdict & follow-ups
- **Stage A:** PASS-WITH-NOTES (N1: 3.3 pg is a documented standard convention vs 3.205 pg primary; intentional).
- **Stage B:** PASS (formulas faithful; test coverage gaps fixed in-session).
- **End-state:** ✅ CLEAN — algorithm fully functional; full unfiltered suite green (Failed: 0, Passed: 6667).
- **Test-quality gate:** PASS (exact sourced values, no green-washing, all public methods + all branches/error
  cases now exercised, honest green on the FULL suite).

---

## Review 2026-09 (campaign `review-2026-09`, batch B25) — 2026-09-28

- **Stage A:** PASS-WITH-NOTES (unchanged N1: 3.3 pg/haploid genome convention). **Stage B:** FAIL → FIXED. **State:** ✅ FIXED.
- **Code now at:** `OncologyAnalyzer.CtdnaMrdChip.cs` (ctDNA region, lines ~1–240); constants remain in `OncologyAnalyzer.cs`.

### Stage A — sources opened this session
1. **WebSearch** (Avanzini 2020 / US 11,085,084): snippet returns verbatim "x = 1 − e^(−nd)" and "p = 1 − (1 − x)^k … p = 1 − e^(−ndk)", n = sequenced genome equivalents, d = detection limit. Confirms the formula.
2. **raw.githubusercontent.com/reiterlab/ctdna** (the Avanzini et al. 2020 reference code; `ctdna/settings.py`, `ctdna/detection.py`): `DIPLOID_GE_WEIGHT_ng = 0.0066` ("average yield of a nucleated cell is 6.6 pg") ⇒ 3.3 pg per haploid genome equivalent — independently confirms the 3.3 pg / ≈303 GE·ng⁻¹ conversion used by `HaploidGenomeEquivalents`. The reference code's full detection model (binomial/NB sampling of mutant fragments against a sequencing-error background with a p-value threshold) is a different, richer model than the closed-form Poisson p = 1 − e^(−ndk) the unit exposes; the unit documents the closed form only (§5.3), which is exactly the patent formula — no mis-attribution.
3. **WebSearch** (CAPP-Seq mean VAF): "VAFs … are averaged at each time point to yield the mean VAF" — confirms the arithmetic mean-of-reporters summarisation of `CalculateMeanVaf`.

### Stage B — findings
- **F1 (numerical robustness, fixed).** `CtDnaDetectionProbability` / `IsCtDnaDetected` evaluated `1.0 − Math.Exp(−λ)`, which cancels catastrophically as λ → 0. Repro (reference = `-numpy.expm1(-λ)` = 60-digit `Decimal` evaluation): λ = 1e-12 → old 9.999778782798785e-13 vs ref 9.999999999995e-13 (2.2e-5 rel. err.); λ = 1e-9 → old 9.999999717180685e-10 vs 9.999999995e-10 (2.8e-8); λ = 1e-300 → old **0** (violates INV-03 strict monotonicity) vs 1e-300. Fix: private `PoissonProbabilityAtLeastOne(λ)` — for λ ≤ ln 2 Kahan's expm1 identity (1 − u)·λ/(−ln u), u = e^(−λ) (Goldberg 1991, ACM Comput. Surv. 23(1)); for λ > ln 2 the direct 1 − u (no cancellation). The split matters: a pure Kahan quotient returned 0.99981 at λ = 744.3 and 1.00034 at λ = 744 (ln of a subnormal u) — caught by the FsCheck property `DetectMRD_MatchesIndependentPanelOracle` during this review, and locked by `DetectionProbability_LambdaWithSubnormalExp_ReturnsOne`. New values match the reference to ≤ 1 ulp on λ ∈ {1e-12, 1e-9, 1e-6, 1e-4, 0.01, 1, 3, 10, 15, 720, 744}. Values for λ ≥ 0.01 unchanged (≤ 1 ulp).
- **Internal duplication removed.** `CtDnaDetectionProbability` re-implemented the three argument guards of `ExpectedMutantMolecules` and `IsCtDnaDetected` re-implemented `1 − e^(−λ)`: both now call `ExpectedMutantMolecules` / `PoissonProbabilityAtLeastOne`. `CalculateTumorFraction` and `CalculateMeanVaf` duplicated the mean-of-VAF loop: both now call private `MeanReporterVaf` (per-variant VAF ≤ 0.5 check enabled for tumour fraction). The unreachable `Math.Min(TF, 1)` clamp was removed (every VAF ≤ 0.5 ⇒ TF ≤ 1); behaviour unchanged.
- **StatisticsHelper:** contains only `NormalCDF`/`Erf`; no canonical mean/Poisson helper to delegate to. No MCP wrapper exposes these methods (grep of `src/Seqeron/Mcp/**`).
- **Note:** a reporter with zero coverage contributes VAF 0 to the mean (convention of shared `CalculateVaf`); documented in the algorithm doc.
- **Cross-batch:** `CalculateTumorFraction`'s 2·v step duplicates `EstimatePurityFromVaf` (SomaticCalling.cs, B22). Delegating would leave `TumorFractionFromVafFactor` (OncologyAnalyzer.cs) unused, which fails the build (S1144/CA1823) and removing it is not an additive change — deferred to the Phase-2 dedup pass.

### Tests added
`DetectionProbability_SmallLambda_FullRelativePrecision` (5 cases), `DetectionProbability_SubnormalLambda_ReturnsLambdaNotZero`, `IsCtDnaDetected_ThresholdEqualToReturnedProbability_IsDetected` (λ = 3 ⇒ 0.950212931632136), `DetectionProbability_LambdaWithSubnormalExp_ReturnsOne` (λ = 720, 744 ⇒ 1.0). Fixed the stale comment value in `DetectionProbability_WorkedExample_15Molecules` (1 − e^(−15) = 0.9999996940976795).
