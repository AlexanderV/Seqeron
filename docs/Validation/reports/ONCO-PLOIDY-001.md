# Validation Report: ONCO-PLOIDY-001 — Tumor Ploidy Estimation + Whole-Genome-Doubling Detection

- **Validated:** 2026-06-16   **Area:** Oncology
- **Canonical method(s):** `OncologyAnalyzer.EstimatePloidy(IEnumerable<AlleleSpecificSegment>)`,
  `OncologyAnalyzer.DetectWholeGenomeDoubling(IEnumerable<AlleleSpecificSegment>)`
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS-WITH-NOTES

## Stage A — Description

### Sources opened & what they confirm (all retrieved this session)

1. **Patchwork (Genome Biology, PMC4053982)** — WebFetch of
   https://pmc.ncbi.nlm.nih.gov/articles/PMC4053982/ returned verbatim:
   *"The average ploidy, PloidyTum, is the average total copy number of all genomic segments
   weighted by segment length."* → confirms ψ = Σ(CN_i·L_i)/Σ(L_i) over per-segment **total** CN.

2. **ASCAT reference implementation (Van Loo lab GitHub)** — WebFetch of
   `ASCAT/R/ascat.runAscat.R` returned the exact computation:
   `ploidy = sum((nA+nB) * s[, "length"]) / sum(s[, "length"])`.
   This is an *independent reference implementation* of the same length-weighted-mean formula
   (CN = nA+nB = Major+Minor; weight = segment length). The two sources agree exactly.

3. **Van Loo et al., PNAS 2010 (PMID 20837533)** — Europe PMC core record abstract retrieved:
   *"we observe aneuploidy (>2.7n) in 45% of the cases"* → ploidy is on the n-scale (2n = diploid),
   and >2.7n marks aneuploidy / near-triploid genomes. Title/authors/journal/year/DOI
   (10.1073/pnas.1009843107) all confirmed.

4. **facets-suite `copy-number-scores.R` (MSKCC, PMID 30013179 / Bielski 2018)** — WebFetch of the
   raw source returned the `is_genome_doubled` function verbatim:
   - default `treshold = 0.5`
   - `frac_elevated_mcn = sum(segs$length[segs$mcn >= 2 & segs$chrom %in% 1:22]) / autosomal_genome`
   - `mcn = tcn - lcn` (major-allele CN = total − minor)
   - `wgd = frac_elevated_mcn > treshold` (strict `>`, autosomes 1–22 only).

### Formula check
- Ploidy ψ = Σ((Major+Minor)·Length)/Σ(Length) — matches Patchwork prose AND ASCAT code exactly.
- WGD = (Σ length where major CN ≥ 2) / (Σ length) **> 0.5** strictly — matches facets-suite code.

### Edge-case semantics (all sourced)
- Empty set → Σ(L)=0 → division undefined → reject (Patchwork weighted mean).
- Length ≤ 0 / negative CN → invalid input.
- WGD boundary exactly 0.5 → false (strict `>`, facets-suite).
- WGD keys on **major** CN ≥ 2, not total ≥ 2: a balanced 1:1 (total 2) genome is NOT doubled.

### Independent cross-check (hand-computed, traced to the sourced formula)
| Case | Computation | Expected |
|------|-------------|----------|
| M1 | (2·100 + 4·100 + 3·50)/250 = 750/250 | **3.0** |
| M3 | (2·300 + 4·10)/310 = 640/310 | **2.0645…** |
| S2 | (0·40 + 4·40)/80 = 160/80 | **2.0** |
| C1 | (3·90 + 3·90 + 4·10)/190 = 580/190 | **3.0526…** (>2.7n aneuploid) |
| M8 | 60/100 = 0.60 > 0.5 | **true** |
| M9 | 50/100 = 0.50, not > 0.5 | **false** |
| M10 | 40/100 = 0.40 ≤ 0.5 | **false** |
The independent ASCAT code (`sum((nA+nB)*length)/sum(length)`) reproduces M1=3.0 etc.

### Findings / divergences
- **Note (documented assumption, not a defect):** facets-suite divides by the *autosomal genome*
  length from a chromosome-size table; this unit divides by the **supplied segments' total length**.
  Semantics are identical when the caller supplies autosomal segments covering the genome, and the
  ≥2/>0.5 rule and operator are unchanged. The TestSpec/Evidence record this assumption explicitly.

Stage A is a faithful, source-exact description. **Verdict: PASS.**

## Stage B — Implementation

### Code path reviewed
- `EstimatePloidy` — `OncologyAnalyzer.cs:4487-4511`
- `DetectWholeGenomeDoubling` — `OncologyAnalyzer.cs:4533-4561`
- `ValidateSegment` — `OncologyAnalyzer.cs:2093-2109`
- `AlleleSpecificSegment` record (Length = End − Start) — `OncologyAnalyzer.cs:1886-1895`
- constant `WholeGenomeDoublingFractionThreshold = 0.5` (4466), `WholeGenomeDoublingMajorCopyNumber = 2` (4459)

### Formula realised correctly?
- `EstimatePloidy`: accumulates `weightedCopyNumberSum += (Major+Minor)*length`, `totalLength += length`,
  returns `weightedCopyNumberSum / totalLength`. Rejects empty (totalLength==0). **Matches** ASCAT/Patchwork.
- `DetectWholeGenomeDoubling`: accumulates `elevatedLength` where `MajorCopyNumber >= 2`, returns
  `(double)elevatedLength/totalLength > 0.5` (strict). **Matches** facets-suite exactly (operator, threshold, mcn≥2).
- Validation (null → ArgumentNullException; End ≤ Start or negative CN → ArgumentException) is shared and correct.

### Cross-verification table recomputed vs code
All seven rows above were re-run via the unit's tests (20/20 pass) and match the externally-derived values.

### Variant/delegate consistency
Two independent public static methods; no `*Fast`/instance variants. Both reuse `ValidateSegment` — consistent
error contract with ONCO-LOH-001/HRD-001.

### Test quality audit (HARD gate)
- **Sourced, not code-echoes:** ploidy assertions are exact `Is.EqualTo(x).Within(1e-10)` with x hand-derived
  from the Patchwork formula (3.0; 640/310; 2.0; 580/190), not echoes of code output. WGD assertions are exact
  booleans from the facets-suite strict-`>` rule.
- **Discriminating test present:** M1 (lengths 100/100/50) is non-discriminating — an unweighted mean also gives
  3.0 — but **M3 (640/310 ≠ 3.0)** explicitly distinguishes length-weighting from a plain mean and would fail a
  buggy unweighted implementation. The plain-mean pitfall is therefore covered.
- **Strict boundary covered:** M9 locks 0.5 → false (would catch a `>=` bug).
- **Major-vs-total covered:** M11 (all 1:1, total 2) → false; S1 (2:0 LOH) → true (would catch a total-CN bug).
- **No green-washing:** no Greater/AtLeast/Contains/ranges where an exact value is known; no widened tolerances;
  no skipped/ignored/commented tests.
- **All branches/edges:** empty, non-positive length, negative CN, null — all asserted for **both** methods;
  zero-CN (S2) and single-segment (M4) covered.
- **Honest green:** full unfiltered suite `Failed: 0, Passed: 6661`; unit subset 20/20; build 0 errors.
  (4 pre-existing NUnit2007 warnings live in an unrelated file `ApproximateMatcher_EditDistance_Tests.cs`;
  no file touched this session.)

### Findings / defects
- **No defect.** Code faithfully realises the Stage-A-validated formulas.
- **Note 1 (input contract, not a defect):** WGD uses `MajorCopyNumber` directly as the major allele CN. The
  record documents `MajorCopyNumber` as "the larger of the two" (line 1884) — a contract shared across the
  oncology segment units — so `MajorCopyNumber == max(Major,Minor) == facets-suite mcn`. Ordering is irrelevant
  for ploidy (Major+Minor). Acceptable, documented assumption.
- **Note 2:** autosomal-genome denominator simplified to supplied-segment length (see Stage A note); does not
  change any threshold/operator.

**Verdict: PASS-WITH-NOTES** (two documented, source-consistent assumptions; no behavioural defect).

## Verdict & follow-ups
- Stage A: **PASS**. Stage B: **PASS-WITH-NOTES**.
- Test-quality gate: **PASS** (exact sourced values, discriminating + boundary + branch coverage, honest green).
- End-state: **CLEAN** — no defect found; no code/test changes required.
- Full suite: 6661 passed / 0 failed. Build: 0 errors. Unit: 20/20.
- No findings logged to FINDINGS_REGISTER (no defect).
