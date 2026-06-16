# Validation Report: ONCO-MSI-001 — Microsatellite Instability (MSI) Detection

- **Validated:** 2026-06-16   **Area:** Oncology
- **Canonical method(s):** `OncologyAnalyzer.CalculateMSIScore(int, int)`, `ClassifyMSIStatus(double)`, `ClassifyBethesdaPanel(int, int)`, `DetectMSI(IEnumerable<bool>)` (plus public constants `MsiHighScoreThreshold`, `BethesdaMsiHighMarkerCount`, `BethesdaMsiLowMarkerCount`, enum `MsiStatus`, record `MsiResult`)
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS

## Scope

The unit implements the **scoring-and-classification layer** of MSI detection:
- Continuous MSI score = fraction of unstable microsatellite loci among all valid loci (MSIsensor / MSIsensor2).
- Continuous status classification at the MSIsensor2 ≥20% cutoff (binary MSI-H vs MSS).
- Categorical NCI/Bethesda classification from a marker count (MSS / MSI-L / MSI-H).
- End-to-end `DetectMSI` over per-locus boolean stability flags.

The upstream per-locus instability call (chi-square comparison of tumor-vs-normal repeat-length
distributions, FDR 0.05) is explicitly out of scope and accepted as boolean input — this is a
documented, reasonable simplification.

## Stage A — Description

### Sources opened this session (retrieved 2026-06-16)

1. **niu-lab/msisensor2 README** — https://github.com/niu-lab/msisensor2/blob/master/README.md
   (WebFetch). Verbatim: MSI score = "number of msi sites / all valid sites"; **"The recommended
   msi score cutoff value is 20% (msi high: msi score >= 20%)."** Boundary is inclusive (≥).
2. **Niu et al. (2014) MSIsensor**, Bioinformatics 30(7):1015–1016 —
   https://academic.oup.com/bioinformatics/article/30/7/1015/236553 (WebFetch). MSI score =
   "the percentage of somatic sites"; valid site = ≥20 spanning reads in both normal and tumor.
   Reported cohort separation: **"Among 71 MSI samples, 70 have an MSI score >3.5. In addition,
   165 of 168 MSS samples have a score <3.5"** → a ~3.5(%) decision boundary on the original
   paired-sample cohort.
3. **Boland et al. (1998) NCI Workshop**, Cancer Res 58(22):5248–5257 —
   https://pubmed.ncbi.nlm.nih.gov/9823339/ (PubMed page itself was CAPTCHA-gated this session;
   the criteria were re-confirmed via independent secondary sources citing it verbatim, incl.
   Frontiers Oncol 2013;3:272 and the Promega MSI review). Criteria: **MSI-H ≥2 of 5 markers
   unstable; MSI-L exactly 1 of 5; MSS 0 of 5.** Panel = BAT-25, BAT-26, D2S123, D5S346, D17S250.
4. **Independent cross-check on thresholds** — search of the MANTIS / performance-evaluation
   literature confirmed the **multiplicity of source-backed cutoffs**: MSIsensor (Niu 2014) ≈3.5
   on its paired cohort; mSINGS / MSI-PCR ≈20%; MSIsensor2 (tumor-only) ≥20%. The 20% cutoff used
   by the code is a legitimate, attributed value (MSIsensor2 / mSINGS), not an invented one.

### Formula check
- MSI score = u / n (unstable loci ÷ valid loci), reported as fraction in [0,1] — matches MSIsensor2
  README verbatim and Niu 2014 ("percentage of somatic sites").
- Continuous MSI-H: score ≥ 0.20 (inclusive) — matches MSIsensor2 README verbatim.
- Bethesda: ≥2→MSI-H, 1→MSI-L, 0→MSS — matches Boland 1998 (and revised-Bethesda fraction form).

### Edge-case semantics (sourced)
- n = 0 → score undefined (division by zero) → must throw. ✓
- 0 ≤ u ≤ n required; u > n or u < 0 invalid → must throw. ✓
- score ∈ [0,1] finite; outside-range / NaN / ±∞ invalid → must throw. ✓
- empty / null flags → throw. ✓
- MSS-vs-MSI-L on a 5-marker panel is intrinsically low-confidence (Boland) — documented as a
  limitation, not a defect.

### Independent cross-check (numbers)
- 5/25 = 0.20 (= cutoff, MSI-H, inclusive); 4/25 = 0.16 (MSS); 3/12 = 0.25; 6/20 = 0.30 (MSI-H);
  2/20 = 0.10 (MSS) — all hand-computed and source-consistent.
- Bethesda 0/5→MSS, 1/5→MSI-L, 2/5→MSI-H, 5/5→MSI-H — Boland 1998.

### Findings / divergences (Stage A)
- **Threshold choice is defensible.** The session prompt mentioned a "~30% / ~3.5" threshold. No
  source supports 30% as an MSI cutoff; the real source-backed cutoffs are MSIsensor ≈3.5 (paired,
  Niu 2014) and MSIsensor2/mSINGS = 20%. The implementation uses 20% and attributes it correctly to
  MSIsensor2. The description/spec/code are internally consistent and externally sourced.
  Minor note: the spec/registry could state more explicitly that 20% is the MSIsensor2 (tumor-only)
  cutoff and that the original MSIsensor paired-sample cutoff is ~3.5 — but the Evidence doc and the
  registry §6.2 already record both. No correction required.

**Stage A verdict: PASS** — every formula, threshold, and edge convention traces verbatim to an
authoritative source retrieved this session.

## Stage B — Implementation

### Code path reviewed
`src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.cs:1569–1750`.

- `CalculateMSIScore` (1634): guards `totalLoci ≤ 0`, `unstable < 0 || unstable > total`, returns
  `(double)u/n`. ✓ matches formula exactly.
- `ClassifyMSIStatus` (1664): guards NaN/±∞/out-of-[0,1]; returns `score >= 0.20 ? MSI_High : MSS`.
  ✓ inclusive boundary correct.
- `ClassifyBethesdaPanel` (1688): guards; `>=2 → MSI_High`, `==1 → MSI_Low`, else `MSS`. ✓.
- `DetectMSI` (1725): null guard, single O(n) pass counting unstable/total, empty → throw, composes
  `CalculateMSIScore` + `ClassifyMSIStatus`. ✓ consistent with the canonical methods.
- Constants (1577/1584/1590) = 0.20 / 2 / 1, each with verbatim source citation in XML doc. ✓.

### Cross-verification table (recomputed vs code, via passing tests)
| Input | Expected (source) | Code |
|-------|-------------------|------|
| 5/25 | 0.20 | 0.20 ✓ |
| 3/12 | 0.25 | 0.25 ✓ |
| 0/25 | 0.0 | 0.0 ✓ |
| 25/25 | 1.0 | 1.0 ✓ |
| score 0.20 | MSI-H (≥20% inclusive) | MSI-H ✓ |
| score 0.19999 | MSS | MSS ✓ |
| score 0.16 | MSS | MSS ✓ |
| Bethesda 0/5 | MSS | MSS ✓ |
| Bethesda 1/5 | MSI-L | MSI-L ✓ |
| Bethesda 2/5 | MSI-H | MSI-H ✓ |
| DetectMSI 6/20 | 0.30, MSI-H | 0.30, MSI-H ✓ |
| DetectMSI 2/20 | 0.10, MSS | 0.10, MSS ✓ |

### Variant/delegate consistency
`DetectMSI` reuses `CalculateMSIScore`/`ClassifyMSIStatus` (no duplicated logic) → INV-05 holds by
construction and is exercised by M13/M14.

### Numerical robustness
Integer division avoided via `(double)`; no overflow on realistic loci counts; div-by-zero guarded;
non-finite scores rejected before comparison.

### Test quality audit (test-quality gate)
File: `tests/Seqeron/Seqeron.Genomics.Tests/OncologyAnalyzer_DetectMSI_Tests.cs`.

- **Sourced, not code-echoes:** every expected value is a verbatim/derived source value
  (5/25=0.20, 3/12=0.25, ≥20% inclusive, Bethesda counts). A wrong implementation (e.g. `>` instead
  of `>=`, or a different cutoff) would fail.
- **No green-washing:** exact `Is.EqualTo` (with 1e-10 tolerance only for double equality), no
  Greater/AtLeast/ranges, no skipped/ignored tests, no widened tolerances.
- **Coverage:** all four public methods + all public constants; every Stage-A branch and error case.
  - *Gaps found and fixed this session* (added, all green):
    - `ClassifyMSIStatus(0.19999)` → MSS — locks the `>=` (not `>`) boundary direction.
    - `ClassifyMSIStatus(1.0)` → MSI-H — exercises the valid upper-bound input.
    - `ClassifyMSIStatus(±Infinity)` → throws — the `IsInfinity` guard branch was previously
      untested (only NaN/1.5/-0.1 were).
    - `Constants_MatchSources` — locks 0.20 / 2 / 1 against silent constant drift.
- **Honest green:** FULL unfiltered suite = **6632 passed, 0 failed, 0 skipped** (was 6629 before
  the +3 additions); `dotnet build` 0 errors, no new warnings in changed files.

### Findings / defects (Stage B)
None. The implementation faithfully realises the validated description; tests were strengthened to
cover previously-untested guard/boundary branches but no code behaviour needed changing.

## Verdict & follow-ups
- **Stage A: PASS. Stage B: PASS. End-state: CLEAN.**
- Test-quality gate: **PASS** (sourced values, no green-washing, full coverage of branches/edge
  cases after adding 3 tests, honest full-suite green).
- No defect logged. Optional future enhancement (not required): expose the alternative MSIsensor
  paired-cohort ~3.5% cutoff as a named constant for callers using paired data — currently the unit
  intentionally implements only the MSIsensor2 20% and Bethesda rules (documented in registry §6.2).
