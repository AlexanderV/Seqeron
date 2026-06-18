# Validation Report: ONCO-HETERO-001 — Tumor Heterogeneity Analysis

- **Validated:** 2026-06-16   **Area:** Oncology
- **Canonical method(s):** `OncologyAnalyzer.CalculateITH(IReadOnlyList<double>)`, `OncologyAnalyzer.InferSubclones(CcfClustering)`, `OncologyAnalyzer.AnalyzeHeterogeneity(vafs, ccfValues, clusterCount)` (`src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.cs:6432–6607`)
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS

## Stage A — Description

### Sources opened & what they confirm (this session)

1. **maftools `mathScore.R`** (raw GitHub master, WebFetched) — verbatim reference code:
   `abs.med.dev = abs(vaf - median(vaf))`; `pat.mad = median(abs.med.dev) * 100`;
   `pat.math = pat.mad * 1.4826 / median(vaf)`. VAF > 1 ⇒ divide by 100 (percent input).
   Algebraically: MATH = 100 · 1.4826 · median(|VAF − median(VAF)|) / median(VAF).
2. **Mroz & Rocco (2013) / Mroz et al. (2015)** (WebSearch, PMC summaries) — primary formula
   **MATH = 100 × MAD/median** over mutant-allele fractions; **MAD = 1.4826 × median(|MAF − median(MAF)|)**,
   the 1.4826 factor chosen "so that the expected MAD of a sample from a normal distribution equals
   the standard deviation."
3. **Wikipedia *Diversity index*** (WebFetched) — Shannon index **H′ = −Σ pᵢ ln(pᵢ)**; natural log
   conventional; single class (p=1) ⇒ H = 0; k equally-abundant classes ⇒ H = ln(k).
4. Repository constants reused and re-checked: `MadConsistencyConstant = 1.4826`, `MathPercentScale = 100.0`,
   `ClonalCcfThreshold = 0.95` (Landau et al. 2013, validated in ONCO-CLONAL-001).

### Formula check

- **MATH** code (`CalculateITH`, lines 6476–6493): `rawMad = Median(|vᵢ − median|)`; `scaledMad = 1.4826·rawMad`;
  returns `100·scaledMad/median`. Matches the primary formula and maftools cell-for-cell. The maftools
  ordering (`median(absdev)*100` then `*1.4826/median`) is algebraically identical (multiplication/division
  commute) — no operation-order precision concern at these magnitudes.
- **Shannon** code (`AnalyzeHeterogeneity`, lines 6574–6579): `H = −Σ (size/n)·ln(size/n)` over occupied
  CCF clusters, natural log (`Math.Log`). Matches H′ = −Σ pᵢ ln pᵢ. Clone fractions pᵢ = per-cluster
  mutation proportions (documented assumption; standard operationalisation per Liu & Zhang 2017 richness).
- **Subclone count** (`InferSubclones`): number of distinct occupied cluster labels = richness.
- **Subclonal fraction** (lines 6582–6591): #(CCF < 0.95) / n, strict `<` (Landau 2013).

### Edge-case semantics check

- Median = 0 ⇒ MATH divides by zero ⇒ `ArgumentException` (guarded, line 6477). Sourced as undefined.
- All-identical / single VAF ⇒ MAD = 0 ⇒ MATH = 0 (valid minimum).
- Single clone ⇒ H = −1·ln 1 = 0; k equal clones ⇒ H = ln k. Both sourced.
- CCF exactly 0.95 ⇒ clonal (strict `<`), per Landau threshold.

### Independent cross-check (numbers, hand-derived from sourced formulas — NOT from the code)

| Case | Input | Hand computation | Expected |
|------|-------|------------------|----------|
| MATH odd | VAFs {0.1,0.2,0.3,0.4,0.5} | median 0.3; absdev {0.2,0.1,0,0.1,0.2}; rawMAD 0.1; 100·1.4826·0.1/0.3 | **49.42** |
| MATH even | VAFs {0.2,0.4,0.6,0.8} | median (0.4+0.6)/2=0.5; absdev sorted {0.1,0.1,0.3,0.3}; rawMAD (0.1+0.3)/2=0.2; 100·1.4826·0.2/0.5 | **59.304** |
| MATH all-same | {0.3,0.3,0.3} | MAD 0 | **0.0** |
| Shannon 2 equal | p={0.5,0.5} | −ln 0.5 | **0.6931471805599453** |
| Shannon 4 equal | p={0.25×4} | ln 4 | **1.3862943611198906** |
| Shannon 1 clone | p={1} | −ln 1 | **0.0** |
| Subclonal frac | CCFs {0.40,0.50,0.98,1.0} | 2 of 4 below 0.95 | **0.5** |
| Subclonal boundary | CCFs {0.94,0.95,0.96,0.97} | only 0.94 < 0.95 | **0.25** |

All sourced values reproduced by hand from the external formulas; they match the Evidence/TestSpec
tables and the implementation's behaviour.

### Findings / divergences

None. The description matches the primary literature and an independent reference implementation
(maftools) exactly. Two documented modelling assumptions (clone fractions = cluster mutation
proportions; even-count median = mean of central order statistics, per R/maftools) are standard and
disclosed.

## Stage B — Implementation

### Code path reviewed

- `CalculateITH` — `OncologyAnalyzer.cs:6453–6494`
- `InferSubclones` — `OncologyAnalyzer.cs:6507–6522`
- `AnalyzeHeterogeneity` — `OncologyAnalyzer.cs:6540–6593`
- `Median` helper — `OncologyAnalyzer.cs:6600–6607` (even = mean of two central order statistics; matches R)
- Constants: `MadConsistencyConstant=1.4826`, `MathPercentScale=100.0`, `ClonalCcfThreshold=0.95`

### Formula realised correctly?

Yes. MATH and Shannon both computed exactly as the sourced formulas (natural log; 1.4826; ×100;
divide by median; strict-`<` subclonal threshold). Input validation: VAF/CCF finite and in [0,1],
non-empty, aligned lengths, zero-median guard, clusterCount ∈ [1, n] (delegated to `ClusterCcfValues`).

### Cross-verification table recomputed vs code

The full unfiltered suite (6681 tests) passes; the ONCO-HETERO-001 fixture asserts every hand-derived
value above with `.Within(1e-10)`. MATH 49.42 / 59.304, Shannon ln2 / ln4 / 0, subclone counts 1/2/3/4,
subclonal fractions 0.5 / 0.25 / 1.0 all match.

### Variant/delegate consistency

`AnalyzeHeterogeneity` delegates MATH to `CalculateITH` (M10 locks the aggregate MATH = 49.42 equal to
the standalone call) and subclone count to `InferSubclones`; consistent.

### Test quality audit (HARD gate)

- **Sourced, not code-echoes:** expected values trace to maftools / Mroz / Shannon / Landau (verified
  this session by WebFetch/WebSearch + hand computation), not to implementation output.
- **Discriminating:** M2 (even-count 59.304) fails any lower-/upper-middle median bug (→0.4 median).
  M5 (ln2=0.6931) fails a log2 implementation (→1.0). M9b (0.25) fails a `<=` threshold bug
  (CCF=0.95 would then count as subclonal → 0.5).
- **No green-washing:** all exact equalities `.Within(1e-10)`; the only inequality (`C1` MATH ≥ 0) is a
  legitimate COULD-tier INV-1 invariant property test over varied inputs, not a substitute for a known value.
- **Coverage:** all three public methods exercised; Stage-A branches covered — odd/even median, all-same,
  single, zero-median throw, null/empty/out-of-range throws, mismatched lengths, single/2/4 clones,
  subclonal fraction + boundary, empty clustering throw.
- **Gap fixed this session:** the strict-`<` subclonal boundary (CCF = 0.95 ⇒ clonal) was not directly
  asserted by any test (a `<`→`<=` flip would have passed the suite). Added
  `AnalyzeHeterogeneity_SubclonalThresholdBoundary_ExcludesExactly0Point95` (CCFs {0.94,0.95,0.96,0.97}
  ⇒ fraction 0.25). No assertion weakened, no skip, no tolerance widened.
- **Honest green:** full unfiltered `dotnet test` = **Failed: 0, Passed: 6681**; `dotnet build` 0 errors
  (4 pre-existing NUnit2007 warnings only in unrelated `ApproximateMatcher_EditDistance_Tests.cs`; the
  ONCO-HETERO-001 files are warning-free).

**Test-quality gate: PASS** (one coverage gap found and fixed in-session).

### Findings / defects

No code or description defect. One test-coverage gap (subclonal threshold boundary) closed in-session.

## Verdict & follow-ups

- **Stage A: PASS. Stage B: PASS. End-state: CLEAN.**
- Fixture grew 15 → 16 tests; full suite green (6681).
- No production code changed; no follow-ups required.
