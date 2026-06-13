# Test Specification: EPIGEN-DMR-001

**Test Unit ID:** EPIGEN-DMR-001
**Area:** Epigenetics
**Algorithm:** Differentially Methylated Region (DMR) detection — tiling window + Fisher's exact test (methylKit model)
**Status:** ☐ In Progress
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-13

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Akalin et al. (2012) methylKit, Genome Biology 13:R87 | 1/3 | https://doi.org/10.1186/gb-2012-13-10-r87 (PMC3491415) | 2026-06-13 |
| 2 | methylKit `tileMethylCounts` man page | 3 | https://github.com/al2na/methylKit/blob/master/man/tileMethylCounts-methods.Rd | 2026-06-13 |
| 3 | methylKit `getMethylDiff` source (diffMeth.R) | 3 | https://rdrr.io/github/vd4mmind/methylkit/src/R/diffMeth.R | 2026-06-13 |
| 4 | methylKit `calculateDiffMeth` man page | 3 | https://github.com/al2na/methylKit/blob/master/man/calculateDiffMeth-methods.Rd | 2026-06-13 |
| 5 | Fisher's exact test (hypergeometric formula, worked example) | 4 | https://en.wikipedia.org/wiki/Fisher's_exact_test | 2026-06-13 |

### 1.2 Key Evidence Points

1. A DMR is a genomic region of adjacent CpG sites that are differentially methylated; per-site methylation level = C/(C+T) — Source 1.
2. Default cutoffs: q-value < 0.01 and %methylation difference > 25% (strict `>`); `difference=25`, `qvalue=0.01` — Sources 1, 3.
3. Hyper-methylated = higher methylation than control (`meth.diff > difference`); hypo-methylated = lower (`meth.diff < -difference`) — Sources 1, 3.
4. Tiling defaults: `win.size=1000`, `step.size=1000`; tiles summarize methylated/unmethylated counts over the window — Source 2.
5. With one sample per group (no replicates) the Fisher's exact test is applied for differential methylation — Source 4.
6. Fisher's exact single-table probability `p = (a+b)!(c+d)!(a+c)!(b+d)! / (a!b!c!d!n!)`; two-sided p = sum of probabilities of all same-margin tables with probability ≤ observed. Worked example a=1,b=9,c=11,d=3 → p ≈ 0.001346076 — Source 5.

### 1.3 Documented Corner Cases

- Empty input → no tiles, no DMRs (Source 2).
- Window with too few covered cytosines → not a region (Sources 1, 2, `cov.bases`).
- Degenerate group (zero coverage in a group within a window) → Fisher exact p = 1.0 (Source 5, fixed-margin degenerate table).

### 1.4 Known Failure Modes / Pitfalls

1. Using a non-strict `>=` cutoff instead of strict `>` over-calls boundary windows — Source 3 (`meth.diff > difference`).
2. Approximating the p-value with a t-test / normal CDF instead of the exact Fisher test mis-states significance for small counts — Sources 4, 5 (exact test mandated for one sample per group).
3. Treating the difference threshold as a fraction vs percentage points inconsistently — Sources 1, 3 (25% = 0.25 fraction).

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `FindDMRs(sample1, sample2, windowSize, minDifference, minCpGCount)` | EpigeneticsAnalyzer | **Canonical** | Tiling + Fisher's exact; deep evidence-based testing |
| `FisherExactTwoSided(a, b, c, d)` (private, exercised via FindDMRs and one direct-equivalent check) | EpigeneticsAnalyzer | **Internal** | Hypergeometric p-value; verified through worked example via a public helper test path |
| `AnnotateDMRs(dmrs, annotations)` | EpigeneticsAnalyzer | **Delegate** | Maps gene/feature annotations onto called DMRs; smoke verification |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | Every reported DMR has `|MeanDifference| > minDifference` (strict) | Yes | Source 3 (`meth.diff > difference`) |
| INV-2 | Annotation is "Hypermethylated" iff MeanDifference > 0, "Hypomethylated" iff < 0 | Yes | Sources 1, 3 (sign of meth.diff) |
| INV-3 | `0 ≤ PValue ≤ 1` for every reported DMR | Yes | Source 5 (probability) |
| INV-4 | Sites separated by ≥ windowSize are placed in different windows | Yes | Source 2 (win.size tiling) |
| INV-5 | A reported DMR contains ≥ minCpGCount covered sites | Yes | Sources 1, 2 (region of adjacent CpGs) |
| INV-6 | MeanDifference = mean(level2 − level1) over the window, in the same [−1,1] fraction units as MethylationLevel | Yes | Source 1 (meth.diff definition) |
| INV-7 | Output is deterministic and ordered by Start ascending | Yes | Source 2 (deterministic tiling) |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | Hyper DMR | Window: g1 all level 0 cov 20, g2 all level 1 cov 20, 3 sites | 1 DMR, MeanDifference = +1.0, Annotation "Hypermethylated", PValue ≈ 0 | Sources 1,3,5 |
| M2 | Hypo DMR | Window: g1 all level 1, g2 all level 0, 3 sites | 1 DMR, MeanDifference = −1.0, Annotation "Hypomethylated" | Sources 1,3 |
| M3 | Below cutoff | Window mean diff exactly 0.25 with minDifference 0.25 | NOT reported (strict `>`) | Source 3 |
| M4 | Above cutoff | Window mean diff 0.30, minDifference 0.25 | reported | Source 3 |
| M5 | Tiling boundary | Two site-clusters > windowSize apart | reported as ≥2 separate windows | Source 2 |
| M6 | Fisher worked example | 2×2 table a=1,b=9,c=11,d=3 single-table probability | ≈ 0.001346076 (Within 1e-7) | Source 5 |
| M7 | Empty input | both samples empty | no DMRs | Source 2 |
| M8 | PValue bounds | any reported DMR | 0 ≤ PValue ≤ 1 | Source 5 |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | Too few sites | Window with < minCpGCount covered sites | not reported | INV-5; region needs adjacent CpGs |
| S2 | Degenerate group | A window where one group has zero coverage | Fisher p = 1.0; not reported (also fails diff cutoff) | Source 5 degenerate margin |
| S3 | One side missing position | Position present in g1 only → g2 level treated as 0 | diff computed against 0 | methylKit per-base; repo fallback |
| S4 | Null input | null sample → ArgumentNullException | throws | repo input-validation contract |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | Determinism / ordering | repeated calls; DMRs ordered by Start | identical, ascending Start | INV-7 |
| C2 | AnnotateDMRs delegation | DMR overlapping a gene annotation | DMR annotated with the gene label | Delegate smoke |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- No `EpigeneticsAnalyzer_*DMR*` test file exists. Sibling tests: `EpigeneticsAnalyzer_CpGDetection_Tests.cs`, `EpigeneticsAnalyzer_Methylation_Tests.cs`.
- Existing production `FindDMRs` used an ad-hoc t-test approximated through `StatisticsHelper.NormalCDF` (no authoritative basis) and a non-strict `>=` cutoff — both are defects per Evidence; corrected in Phase 5. No `AnnotateDMRs` existed; added in Phase 5.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1 Hyper DMR | ❌ Missing | new unit |
| M2 Hypo DMR | ❌ Missing | new unit |
| M3 Below cutoff (strict) | ❌ Missing | new unit |
| M4 Above cutoff | ❌ Missing | new unit |
| M5 Tiling boundary | ❌ Missing | new unit |
| M6 Fisher worked example | ❌ Missing | new unit |
| M7 Empty input | ❌ Missing | new unit |
| M8 PValue bounds | ❌ Missing | new unit |
| S1 Too few sites | ❌ Missing | new unit |
| S2 Degenerate group | ❌ Missing | new unit |
| S3 One side missing | ❌ Missing | new unit |
| S4 Null input | ❌ Missing | new unit |
| C1 Determinism/order | ❌ Missing | new unit |
| C2 AnnotateDMRs delegation | ❌ Missing | new unit |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/EpigeneticsAnalyzer_DMR_Tests.cs` — all DMR cases here.
- **Remove:** nothing (no prior DMR tests).

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| EpigeneticsAnalyzer_DMR_Tests.cs | Canonical DMR fixture | 14 |

### 5.5 Phase 7 Work Queue

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| 1 | M1 | ❌ Missing | Implemented | ✅ Done |
| 2 | M2 | ❌ Missing | Implemented | ✅ Done |
| 3 | M3 | ❌ Missing | Implemented | ✅ Done |
| 4 | M4 | ❌ Missing | Implemented | ✅ Done |
| 5 | M5 | ❌ Missing | Implemented | ✅ Done |
| 6 | M6 | ❌ Missing | Implemented (public `FisherExactProbability`) | ✅ Done |
| 7 | M7 | ❌ Missing | Implemented | ✅ Done |
| 8 | M8 | ❌ Missing | Implemented | ✅ Done |
| 9 | S1 | ❌ Missing | Implemented | ✅ Done |
| 10 | S2 | ❌ Missing | Implemented | ✅ Done |
| 11 | S3 | ❌ Missing | Implemented | ✅ Done |
| 12 | S4 | ❌ Missing | Implemented | ✅ Done |
| 13 | C1 | ❌ Missing | Implemented | ✅ Done |
| 14 | C2 | ❌ Missing | Implemented | ✅ Done |

**Total items:** 14
**✅ Done:** 14 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ Covered | `FindDMRs_HyperMethylatedWindow_ReportedAsHypermethylated` |
| M2 | ✅ Covered | `FindDMRs_HypoMethylatedWindow_ReportedAsHypomethylated` |
| M3 | ✅ Covered | `FindDMRs_DifferenceEqualToCutoff_NotReported` |
| M4 | ✅ Covered | `FindDMRs_DifferenceAboveCutoff_Reported` |
| M5 | ✅ Covered | `FindDMRs_SitesBeyondWindowSize_SplitIntoSeparateWindows` |
| M6 | ✅ Covered | `FisherExactProbability_WikipediaWorkedExample_MatchesPublishedValue` |
| M7 | ✅ Covered | `FindDMRs_EmptyInput_ReturnsNoRegions` |
| M8 | ✅ Covered | `FindDMRs_ReportedRegion_PValueWithinUnitInterval` |
| S1 | ✅ Covered | `FindDMRs_FewerSitesThanMinCpGCount_NotReported` |
| S2 | ✅ Covered | `FisherExactProbability_DegenerateMargin_ReturnsOne` |
| S3 | ✅ Covered | `FindDMRs_PositionMissingInSecondSample_TreatedAsZero` |
| S4 | ✅ Covered | `FindDMRs_NullSample_Throws` |
| C1 | ✅ Covered | `FindDMRs_RepeatedCalls_DeterministicAndOrderedByStart` |
| C2 | ✅ Covered | `AnnotateDMRs_OverlappingAnnotation_LabelsRegion` |

---

## 6. Assumption Register

**Total assumptions:** 2

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | Per-window pooling of single-sample sites into one 2×2 table (mirrors tileMethylCounts + Fisher) | FindDMRs p-value |
| 2 | numC/numT reconstructed from `MethylationLevel × Coverage` (rounded) since `MethylationSite` stores a fraction, not raw counts | FindDMRs p-value |

---

## 7. Open Questions / Decisions

1. The repo API expresses `minDifference` and `MeanDifference` as fractions in [−1,1] (consistent with `MethylationLevel`), whereas methylKit uses percentage points (25). Decision: keep fraction units; documented equivalence 0.25 fraction = 25 percentage points = methylKit `difference=25`. No correctness impact (same boundary).
2. `AnnotateDMRs` has no single canonical specification; it is a thin mapping of region→feature labels and is tested only for delegation behavior.
