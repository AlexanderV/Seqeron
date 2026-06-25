# Validation Report: ONCO-MRD-001 — Minimal/Molecular Residual Disease (MRD) Detection

- **Validated:** 2026-06-24   **Area:** Oncology
- **Canonical method(s):** `OncologyAnalyzer.DetectMRD`, `TrackVariantsOverTime`, `IsVariantDetected` (legacy ≥2-of-N panel rule + Poisson LoD, ONCO-CTDNA-001 reuse); **INVAR/INVAR2 extension under review:** `EstimateInvarSignal`, `IntegratedMutantAlleleFractionV2`, `EstimateInvarSignalWithSize` (+`FragmentSizeProfile`/`InvarMolecule`), `SuppressOutlierLoci`, `EstimateLocusBackground`, `PassesBothStrandsFilter`
- **Stage A verdict:** PASS-WITH-NOTES
- **Stage B verdict:** PASS

This session re-validates the INVAR/INVAR2 extension added by commits `28b1199a` (background-subtracted AF-weighted GLRT + IMAFv2) and `01e21a8a` (fragment-size weighting, outlier suppression, locus-noise filtering, control-derived background estimation). The legacy panel-rule path was already validated in the prior report (2026-06-16); it is unchanged and re-confirmed.

## Stage A — Description

### Sources opened this session
- **INVAR2 reference implementation, `R/shared/detectionFunctions.R`** (nrlab-CRUK/INVAR2, the Rosenfeld-lab pipeline implementing Wan et al. 2020) — fetched the verbatim R bodies of `calc_log_likelihood`, `estimate_p_EM`, `calc_likelihood_ratio`, `calc_log_likelihood_with_RL`, `estimate_p_EM_with_RL`, `calc_likelihood_ratio_with_RL`.
- **Wan et al. 2020, Sci Transl Med 12(548):eaaz8084 (INVAR/IMAF)** — integration of variant reads across many patient-specific loci; depth-weighted, background-subtracted IMAF; tumour-derived cfDNA is shorter (size weighting); sensitivity to ~ppm.
- **Reinert 2019 / Signatera white paper / PMC9265001** — legacy ≥2-of-16 positivity rule and Poisson LoD `p = 1 − e^(−nfm)` (unchanged; carried from prior report).
- TestSpec `tests/TestSpecs/ONCO-MRD-001.md` and `docs/Evidence/ONCO-MRD-001-Evidence.md` (the cited verbatim formulas and synthetic datasets).

### Formula check (INVAR2 source vs description vs code)
| Quantity | INVAR2 source (verbatim) | Description / code | Status |
|----------|--------------------------|--------------------|--------|
| mixture `q` | `q = AF(1−e)p + (1−AF)e·p + e(1−p)` | identical | ✅ |
| `g` | `g = AF(1−e) + (1−AF)e` | identical | ✅ |
| logL | `Σ[lchoose(R,M)+M·log q+(R−M)·log(1−q)] / length(R)` | identical (normalised by locus count) | ✅ |
| EM E-step | `Z0=(1−g)p/((1−g)p+(1−e)(1−p))`, `Z1=gp/(gp+e(1−p))` | identical | ✅ |
| EM M-step | `p=Σ(M·Z1+(R−M)·Z0)/ΣR`, init 0.01, 200 iters | identical | ✅ |
| LR | `logL(p̂)−logL(0)` | identical | ✅ |
| IMAFv2 | `pmax(0, MEAN_AF−BACKGROUND_AF)` then `weighted.mean(·, TOTAL_DP)` | identical | ✅ |
| zero-bg floor | `ifelse(bg>0, bg, 1/BACKGROUND_DP)` | floored to `1/depth` | ✅ |
| informative filter | `filter(TUMOUR_AF>0)` | identical | ✅ |
| with-RL `L0/L1` | `L0=(1−e)P0(1−p)+(1−g)P1·p`, `L1=e·P0(1−p)+g·P1·p`; `logL=Σ[M·log L1+(R−M)·log L0]/length(R)` | identical | ✅ |
| with-RL EM | `Z0=(1−g)P1·p/(…+(1−e)P0(1−p))`, `Z1=g·P1·p/(…+e·P0(1−p))` | identical | ✅ |
| outlier | `P_THRESHOLD=α/n_loci`; `P_ESTIMATE=max(EM, weighted.mean(AF,TUMOUR_AF))`; `binom.test(x,DP,p,"greater")`; OUTLIER ⟺ `p.value ≤ P_THRESHOLD` | identical (tail ≤ threshold flags outlier) | ✅ |
| locus-noise/bg | `BACKGROUND_AF=Σ(ALT_F+ALT_R)/ΣDP`; `LOCUS_NOISE.PASS=(N_signal/N)<prop ∧ bg<maxBg` | identical | ✅ |
| both-strands | `ALT_F>0 ∧ ALT_R>0 ∨ AF==0` | identical | ✅ |

Every formula in the Evidence/TestSpec matches the INVAR2 source line-for-line.

### Edge-case semantics
- Pure-background sample ⇒ `p̂≈0`, `LR≈0`, not detected (INV-7). Confirmed.
- No informative locus (all tumour AF=0) ⇒ `ArgumentException` (INVAR `filter(TUMOUR_AF>0)` empties the table). Confirmed.
- Zero background ⇒ floored to `1/depth` so logs are finite (INVAR `doMain`). Confirmed.
- Flat size profile (P1==P0) ⇒ size factor cancels ⇒ with-RL LR == no-size LR (INV-11). Confirmed numerically.
- Outlier `x≤0` ⇒ binomial tail 1, never an outlier. Confirmed.
- Both-strands: `AF==0` passes vacuously; single-strand-only fails. Confirmed.

### Independent cross-check (numbers — Python port of the INVAR2 equations, computed this session)
| Case | Expected (TestSpec/Evidence) | Python reference | Match |
|------|------------------------------|------------------|-------|
| GLRT inj=0 (alt=1) | p̂≈3.3e-5, LR≈0 | p̂=3.276e-5, LR=−8.5e-5 | ✅ |
| GLRT inj=0.01 (alt=5) | p̂≈0.01002, LR≈4.06 | p̂=0.01002, LR=4.0552 | ✅ |
| GLRT inj=0.05 (alt=21) | p̂≈0.0501, LR≈44.14 | p̂=0.05010, LR=44.137 | ✅ |
| size-weighted (M17) | p̂=0.12042621132507245, LR=0.19691792427890276 | identical to 17 sig figs | ✅ |
| binom tail P(X≥1\|1000,0.001) | 0.6323045752290356 | 0.6323045752289 | ✅ |
| binom tail P(X≥50\|1000,0.001) | 3.7264670792676e-66 | 3.7264670792654e-66 | ✅ |
| control bg (M21) | 40/20000 = 0.002 | 0.002 (by hand) | ✅ |

### Findings / divergences (Stage A → PASS-WITH-NOTES)
1. **KDE size-smoothing now implemented opt-in (2026-06-24 update; was Assumption #2).** INVAR2's `estimate_real_length_probability` smooths the per-length size histogram with a weighted Gaussian KDE (`density()`, `adjust=0.03`) and integrates over each integer bin. This is now provided opt-in as `FragmentSizeProfile.FromKernelDensity` — a Gaussian KDE (Silverman 1986 eq. 2.2a) with Silverman's-rule bandwidth (R `bw.nrd0`, `adjust` multiplier) integrated analytically over each integer bin via `Φ(z)=½[1+erf(z/√2)]` and renormalised over the support. The default `FragmentSizeProfile` constructor (discrete `COUNT/TOTAL`) is unchanged. The mixture/EM/GLRT-with-RL equations are unchanged; the KDE only changes how `P0`/`P1` are estimated. No longer a residual.
2. **Per-locus background `e` is caller-supplied in `EstimateInvarSignal`** (the no-size GLRT), with `EstimateLocusBackground` providing the control-derived estimate as a separate step rather than an integrated pipeline. This matches INVAR2's separation of `createLociErrorRateTable` (parse stage) from the detection GLRT. This is the sole remaining honest residual (the optional, separately-invoked background re-estimation), disclosed in TestSpec §7 / Assumption register.
3. INVAR2 returns the LR normalised per-locus (`/length(R)`); the code reproduces this exactly. The absolute LR scale therefore depends on locus count (a property of the source, preserved faithfully).

## Stage B — Implementation

### Code path reviewed
`src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.cs`:
- `IntegratedMutantAlleleFractionV2` L8807–8832; `EstimateInvarSignal` L8857–8937; `EstimateCtDnaFractionEm` L8946–8970; `InvarLogLikelihood` L8978–8998; `LogChoose` L9004–9012.
- `EstimateInvarSignalWithSize` L9177–9241; `EstimateCtDnaFractionEmWithRl` L9249–9274; `InvarLogLikelihoodWithRl` L9282–9304; `FragmentSizeProfile` L9084–9150.
- `SuppressOutlierLoci` L9340–9426; `BinomialUpperTail` L9433–9465; `EstimateLocusBackground` L9517–9565; `PassesBothStrandsFilter` L9575–9587.

### Formula realised correctly? (evidence)
- `InvarLogLikelihood` computes `q = AF(1−e)p+(1−AF)e·p+e(1−p)`, clamps to `(ε,1−ε)`, sums `lchoose+M·log q+(R−M)·log(1−q)` then divides by locus count — exact INVAR2 `calc_log_likelihood`.
- `EstimateCtDnaFractionEm` uses `g`, the exact Z0/Z1 E-step and the `ΣR`-normalised M-step, init 0.01, 200 iters — exact `estimate_p_EM`.
- LR = `logL(p̂)−logL(0)` — exact `calc_likelihood_ratio`.
- `LogChoose` via `LogGamma` = R's `lchoose`; verified against the Python `math.lgamma` reference (binomial tails match to 13 sig figs).
- With-RL: `L0=(1−e)P0(1−p)+(1−g)P1·p`, `L1=e·P0(1−p)+g·P1·p`, EM with the P1/P0 factors — exact `*_with_RL`; reproduces the pinned `0.19691792427890276`/`0.12042621132507245` to full double precision.
- `SuppressOutlierLoci`: null `P_ESTIMATE=max(EM, weighted.mean(VAF,TUMOUR_AF))` restricted to loci with `VAF≤afThreshold ∧ alt≤maxMutantReads ∧ AF>0`; Bonferroni `α/n`; one-sided upper-tail binomial; flags outlier when `tail ≤ threshold` — exact `repolish`.
- `EstimateLocusBackground`/`PassesBothStrandsFilter`: pooled `Σ(ALT_F+ALT_R)/ΣDP`, recurrence+bg AND rule, both-strands rule — exact `createLociErrorRateTable`.

### Cross-verification table recomputed vs code
The 57-test class pins exact INVAR2-derived values (LR 4.06/44.14, size LR `0.19691792427890276`, binomial tails `0.6323045752290356`/`3.7264670792676273e-66`, bg `5e-5`/`0.0125`/`0.002`). All independently reproduced this session by the Python port (table above) and all pass under the actual code.

### Variant/delegate consistency
- The no-size and with-RL GLRTs reduce to one another when P1==P0 (SW2) — verified the with-RL path on a flat profile equals `EstimateInvarSignal` on the equivalent one-molecule-per-locus loci (`Within(1e-9)`).
- `SuppressOutlierLoci` reuses `EstimateCtDnaFractionEm` and `BinomialUpperTail`; `EstimateInvarSignal` reuses `IntegratedMutantAlleleFractionV2`. Shared primitives, consistent.

### Numerical robustness
- `q` clamped to `(double.Epsilon, 1−double.Epsilon)`; with-RL `L0/L1` floored to `double.Epsilon`; zero background floored to `1/depth`; `BinomialUpperTail` uses log-space `exp(lchoose+i·logp+(n−i)·log(1−p))` summation (no overflow at n=1000); `alt` clamped to `[0,total]`; denominators guarded against zero. No div-by-zero or overflow on the stated ranges.

### Test quality audit (HARD gate)
57 tests, all deterministic, asserting **exact** sourced values via `Within` tolerances tight enough to be real (1e-9 on the pinned size/binomial/background values; ≤0.05 on the LR references derived independently). Coverage spans all 25 MUST cases (M1–M25), SHOULD (S1–S6) and COULD validation cases (C1–C21). No `Greater`-only green-washing on the load-bearing assertions (directional INV-9/INV-10 checks are paired with exact pinned magnitudes). Edge/error paths (null, empty, out-of-range AF/bg/threshold/α, no-informative-locus) all covered.

### Findings / defects
- **No code defect.** Every formula realises the INVAR2 source exactly; all pinned values reproduce to full precision against an independent Python port.

## Verdict & follow-ups
- **Stage A: PASS-WITH-NOTES** — INVAR/INVAR2 maths confirmed line-for-line against the nrlab-CRUK reference implementation; two documented by-design residuals (KDE size-smoothing → empirical proportion; caller-supplied vs separately-estimated background) and the per-locus LR normalisation (a faithful property of the source).
- **Stage B: PASS** — implementation faithful; all 57 tests pass; independent cross-checks match to full double precision.
- **End-state: CLEAN** — no defect found. Build 0 warnings / 0 errors; `OncologyAnalyzer_DetectMRD_Tests` = **66 passed, 0 failed** after the 2026-06-24 KDE addition (was 57).
- **2026-06-24 update:** the KDE size-smoothing residual (Assumption #2) is now closed — opt-in `FragmentSizeProfile.FromKernelDensity` (Gaussian KDE, Silverman bandwidth, analytic Gaussian-bin integral) was added with 9 new evidence-based tests (K1–K6, C22–C24), the default discrete profile unchanged. The sole remaining honest residual is the optional, separately-invoked control-derived background re-estimation (`EstimateLocusBackground`), which mirrors INVAR2's own parse/detection separation.
