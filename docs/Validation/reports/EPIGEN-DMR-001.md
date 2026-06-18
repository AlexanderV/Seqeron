# Validation Report: EPIGEN-DMR-001 — Differentially Methylated Region (DMR) detection

- **Validated:** 2026-06-15   **Area:** Epigenetics
- **Canonical method(s):** `EpigeneticsAnalyzer.FindDMRs`, `EpigeneticsAnalyzer.FisherExactProbability`, `EpigeneticsAnalyzer.AnnotateDMRs` (private `FisherExactTwoSided`, `TryBuildRegion`)
- **Stage A verdict:** PASS-WITH-NOTES
- **Stage B verdict:** PASS

## Stage A — Description

### Sources opened this session (retrieved 2026-06-15)
1. **methylKit `tileMethylCounts` man page** — https://github.com/al2na/methylKit/blob/master/man/tileMethylCounts-methods.Rd
   Confirmed defaults `win.size=1000`, `step.size=1000`, `cov.bases=0`, `mc.cores=1`; "summarizes methylated/unmethylated base counts over tilling windows across genome."
2. **methylKit `getMethylDiff` source (diffMeth.R, al2na master)** — https://raw.githubusercontent.com/al2na/methylKit/master/R/diffMeth.R
   Confirmed defaults `difference=25`, `qvalue=0.01`; hyper rule `qvalue<qvalue & meth.diff > difference` (strict `>`); hypo rule `qvalue<qvalue & meth.diff < -1*difference` (strict `<`).
3. **methylKit `getMethylDiff` man page** (RDocumentation/GitHub) — corroborates `difference` is "absolute value of methylation percentage change", default 25; `qvalue` default 0.01; `type` ∈ {hyper, hypo, all}.
4. **Fisher's exact test — Wikipedia** — https://en.wikipedia.org/wiki/Fisher%27s_exact_test
   Confirmed single-table hypergeometric `p = (a+b)!(c+d)!(a+c)!(b+d)! / (a!b!c!d!n!) = C(a+b,a)C(c+d,c)/C(n,a+c)`; worked example cells (1,9,11,3) → **p ≈ 0.001346076**; two-sided rule = sum of probabilities of all same-margin tables with prob ≤ observed.
5. **methylKit `calculateDiffMeth`** (Evidence prior session): one sample per group → Fisher's exact test; replicates → logistic regression. Consistent with the unit's scope.

### Formula check
- Per-site methylation level C/(C+T) and meanDiff = mean(level2 − level1): matches Akalin 2012.
- Strict reporting cutoff `|meanDiff| > minDifference`: matches diffMeth.R `meth.diff > difference` (verbatim).
- Hyper iff diff > 0, hypo iff diff < 0: matches getMethylDiff hyper/hypo rules.
- Hypergeometric single-table probability: matches Wikipedia formula exactly.
- Two-sided p = Σ p(table) over same-margin tables with p ≤ p(observed): matches Wikipedia/standard convention.

### Edge-case semantics
Empty input → no DMRs (no tiles); window with < minCpGCount → not a region; |diff| == cutoff → excluded (strict); degenerate margin (zero row/col total) → Fisher p = 1.0; null sample → ArgumentNullException. All defined and sourced.

### Independent cross-check (numbers, scipy 1.13.1 + math.comb)
- Single-table (1,9,11,3): `comb(10,1)*comb(14,11)/comb(24,12)` = **0.0013460761879122358**; lgamma form = 0.0013460761879122362 → matches Wikipedia 0.001346076. ✓
- Two-sided Fisher via the implementation's own algorithm, replicated in Python and compared to `scipy.stats.fisher_exact` for tables (1,9,11,3), (12,72,108,12), (3,7,8,2), (5,5,5,5), (2,8,7,3), (0,60,60,0), (10,10,3,17): **all match to <1e-9**. ✓
- Symmetric (5,5,5,5) single-table = `252*252/184756` = **0.34371820130334063**. ✓

### Findings / divergences (PASS-WITH-NOTES)
- **Units note (documented, not a defect):** methylKit expresses `meth.diff`/`difference` in percentage points (25); this unit uses fractions in [−1,1] (0.25). The boundary semantics are identical (25% = 0.25), and the description records the equivalence explicitly. No correctness impact.
- **Two-sided convention note:** Wikipedia's worked-example narrative quotes the *one-tailed* sum (0.001379728). The implementation uses the symmetric two-sided rule (sum of all same-margin tables with p ≤ observed), which is the standard Fisher two-sided p and matches scipy's default `fisher_exact` (0.0027594561852… for that table). This is correct, not a divergence from the formula.
- **q-value / multiple-testing intentionally omitted** (raw Fisher p returned); single-sample-per-group only (no logistic regression). Both are documented scope limitations, sourced as the correct regimes.

## Stage B — Implementation

### Code path reviewed
`src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/EpigeneticsAnalyzer.cs`
- `FindDMRs` L604–615 (null guards + deferred iterator); `FindDMRsIterator` L617–657 (dictionary index, sorted union, tiling: new window when `pos - start >= windowSize`); `TryBuildRegion` L659–720 (minCpGCount gate, meanDiff, count reconstruction via `round(level*coverage)`, strict cutoff `Math.Abs(meanDiff) <= minDifference`, sign-based hyper/hypo); `AnnotateDMRs`/Iterator L730–762; `FisherExactProbability` L770–786; `FisherExactTwoSided` L793–829; `LogFactorial` L831–843.

### Formula realised correctly?
Yes. Tiling, strict `>` cutoff, sign-based annotation, hypergeometric single-table prob, and the two-sided enumeration all match the validated description. The two-sided enumeration was independently re-implemented in Python and shown to match the C# logic and scipy across 7 tables.

### Cross-verification table recomputed vs code (via the suite)
| Case | Input | Expected (external) | Result |
|------|-------|---------------------|--------|
| M6 single-table | FisherExactProbability(1,9,11,3) | 0.001346076 (Wikipedia) | PASS (Within 1e-7) |
| M1 PValue | pooled [[0,60],[60,0]] | 2.070073888186964e-35 (scipy) | PASS (added) |
| M8 PValue | pooled [[12,108],[108,12]] | 2.475428262210228e-39 (scipy) | PASS (added) |
| symmetric | FisherExactProbability(5,5,5,5) | 0.34371820130334063 (comb) | PASS (added) |
| zero total | FisherExactProbability(0,0,0,0) | 1.0 | PASS (added) |
| negative cell | FisherExactProbability(-1,9,11,3) | ArgumentOutOfRangeException | PASS (added) |

### Variant/delegate consistency
`FisherExactTwoSided` (private, drives `PValue`) and the public `FisherExactProbability` are consistent and both verified against scipy. `AnnotateDMRs` overlap logic (half-open feature interval) tested for both hit and miss.

### Test quality audit (HARD gate)
- **Defect found & fixed:** M8 (`FindDMRs_ReportedRegion_PValueWithinUnitInterval`) only asserted `InRange(0,1)` — a bounds check a deliberately-wrong implementation would also pass (code-echo / weak assertion). Rewrote to lock the exact externally-computed two-sided p (scipy 2.475428262210228e-39) in addition to the bounds. M1 strengthened to also lock its PValue (scipy 2.07e-35).
- **Coverage gaps closed:** added `FisherExactProbability` error branch (negative cell → throws), degenerate `n==0` → 1.0, and a non-extreme symmetric single-table value; added `AnnotateDMRs` null-input rejection. These exercise previously-untested public branches/edge cases.
- No assertions weakened, no tolerances widened, no tests skipped/ignored. All expected values trace to scipy / `math.comb` / Wikipedia retrieved this session, not to code output.
- **Honest green:** full unfiltered suite `dotnet test` = **Failed: 0, Passed: 6482** (was 6478; +4 new tests). `dotnet build` = 0 errors; the DMR test file compiles warning-free (the 4 build warnings are pre-existing NUnit2007 warnings in the unrelated `ApproximateMatcher_EditDistance_Tests.cs`).

### Findings / defects
- One test-quality defect (weak M8 bounds-only assertion) — **fixed this session**. No implementation defect: the algorithm matches the validated description and an independent reference (scipy) exactly.

## Verdict & follow-ups
- **Stage A:** PASS-WITH-NOTES (units fraction-vs-percentage; two-sided convention note; documented scope limits). Description is biologically/mathematically correct.
- **Stage B:** PASS — code faithfully realises the description; all cross-checked values match scipy/Wikipedia.
- **End-state:** ✅ CLEAN — the test-quality defect was completely fixed; build + full suite green.
- **Test-quality gate:** PASS (after fix).
