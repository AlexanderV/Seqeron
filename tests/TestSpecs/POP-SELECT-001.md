# Test Specification: POP-SELECT-001

**Test Unit ID:** POP-SELECT-001
**Area:** PopGen
**Algorithm:** Selection Signature Detection (integrated Haplotype Score, iHS)
**Status:** ☐ In Progress
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-13

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Voight, Kudaravalli, Wen & Pritchard (2006), PLoS Biology 4(3):e72 | 1 | https://doi.org/10.1371/journal.pbio.0040072 (PDF: https://web.stanford.edu/group/pritchardlab/publications/VoightEtAl06.pdf) | 2026-06-13 |
| 2 | Sabeti et al. (2002), Nature 419:832–837 | 1 | https://pubmed.ncbi.nlm.nih.gov/12397357/ | 2026-06-13 |
| 3 | Szpiech & Hernandez (2014), selscan, MBE 31(10):2824 | 3 | https://arxiv.org/pdf/1403.6854 | 2026-06-13 |
| 4 | rehh package vignette (Gautier et al.), CRAN | 3 | https://cran.r-project.org/web/packages/rehh/vignettes/rehh.html | 2026-06-13 |

### 1.2 Key Evidence Points

1. EHH_c(x_i) = Σ_h C(n_h,2) / C(n_c,2) over distinct extended haplotypes carrying core allele c — selscan Eq. 3; rehh (equivalent form).
2. iHH = area under the EHH curve, trapezoidal rule, summed over both directions away from the core SNP — Voight (2006) §M&M; selscan Eq. 4.
3. Integration is truncated where EHH first drops below 0.05 — Voight (2006) §M&M; rehh `limehh = 0.05`.
4. Unstandardized iHS = ln(iHH_A / iHH_D) (ancestral numerator, derived denominator) — Voight (2006). selscan uses the reciprocal (opposite sign) and explicitly notes the swap.
5. Standardized iHS = (ln(iHH_A/iHH_D) − E_p[·]) / SD_p[·] within derived-allele-frequency bins; ~ standard normal — Voight (2006).
6. Genome-wide signal = proportion of SNPs with |iHS| > 2 in 50-SNP windows — Voight (2006) §M&M.

### 1.3 Documented Corner Cases

- Balanced EHH decay ⇒ iHH_A/iHH_D ≈ 1 ⇒ unstandardized iHS ≈ 0 (Voight 2006).
- iHS reported only for polymorphic SNPs with ancestral state and MAF > 5% (Voight 2006).
- Long derived haplotype (slow EHH decay on derived) ⇒ large negative iHS (Voight 2006).

### 1.4 Known Failure Modes / Pitfalls

1. Sign-convention confusion between Voight (ln iHH_A/iHH_D) and selscan (ln iHH_1/iHH_0) — selscan (2014) explicit note.
2. Forgetting the 0.05 truncation inflates iHH for slowly-decaying alleles — Voight (2006) §M&M.

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `CalculateEhh(IReadOnlyList<string>)` | PopulationGeneticsAnalyzer | Canonical | EHH = Σ C(n_h,2)/C(n_c,2) |
| `CalculateIHS(IReadOnlyList<string>, IReadOnlyList<int>, int)` | PopulationGeneticsAnalyzer | Canonical | Full iHS pipeline → unstandardized ln(iHH_A/iHH_D) |
| `StandardizeIHS(IReadOnlyList<(double,double)>, int)` | PopulationGeneticsAnalyzer | Canonical | Frequency-binned z-standardization |
| `ScanForSelection(IReadOnlyList<double>, int)` | PopulationGeneticsAnalyzer | Canonical | Proportion of |iHS|>2 per window |
| `IntegrateEhh` / `IntegrateDirection` / `Choose2` | PopulationGeneticsAnalyzer | Internal | Tested via canonical methods |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | EHH ∈ [0,1]; EHH=1 for a single chromosome, 0 for all-distinct | Yes | selscan Eq. 3 |
| INV-2 | iHH ≥ 0 (sum of non-negative trapezoids) | Yes | Voight §M&M |
| INV-3 | Balanced EHH decay ⇒ unstandardized iHS = 0 | Yes | Voight (2006) |
| INV-4 | ln(iHH_A/iHH_D) = −ln(iHH_D/iHH_A) (Voight vs selscan sign symmetry) | Yes | selscan note |
| INV-5 | Standardized scores within a bin have mean 0 (and sd 1 when >1 element) | Yes | Voight standardization |
| INV-6 | ScanForSelection proportion = ExtremeCount / SnpCount ∈ [0,1] | Yes | Voight §M&M |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | EHH worked values | `CalculateEhh` on `{11,11,11,10}` and `{00,00,01,01}` | 0.5 ; 0.333333… | selscan Eq. 3 |
| M2 | EHH single / all-distinct | single chrom; three distinct | 1.0 ; 0.0 | selscan Eq. 3 |
| M3 | iHS constructed panel | 3× `AA1GG` derived (identical), `TC0TC`/`GA0AG`/`CT0CA` ancestral; pos 0,10,20,30,40; core=2 | iHH_A=10, iHH_D=40, iHS=ln(0.25)=−1.386294361, derivedFreq=0.5 | Voight iHH+trapezoid+0.05; ln(A/D) |
| M4 | Balanced decay | symmetric panel, both alleles identical flanks | unstandardized iHS = 0 | Voight (2006) |
| M5 | Sign convention | long derived haplotype | iHS < 0 (Voight: negative ⇒ derived sweep) | Voight (2006) |
| M6 | Standardization | bin {−1.0, 0.0, +1.0} (same freq bin) | mean 0, sd 1 ⇒ standardized {−1, 0, +1} | Voight standardization |
| M7 | Standardization single bin | one score in a bin | standardized = 0 (sd undefined ⇒ 0) | Assumption (SD estimator) |
| M8 | ScanForSelection proportion | scores with 2 of 4 having |iHS|>2, window 4 | ProportionExtreme = 0.5, ExtremeCount=2 | Voight |iHS|>2 |
| M9 | ScanForSelection windows | 5 scores, window 2 | 3 windows of sizes 2,2,1 | Voight 50-SNP windows (generalized) |
| M10 | iHS null haplotypes | `CalculateIHS(null, pos, 0)` | `ArgumentNullException` | input-validation |
| M11 | iHS monomorphic core | all core alleles '1' | `ArgumentException` | Voight (polymorphic SNPs only) |
| M12 | iHS inconsistent length | haplotype shorter than positions | `ArgumentException` | input-validation |
| M13 | iHS invalid allele | core allele not '0'/'1' | `ArgumentException` | EHH defined on alleles 0/1 |
| M14 | iHS coreIndex OOR | coreIndex = positions.Count | `ArgumentOutOfRangeException` | input-validation |
| M15 | EHH null | `CalculateEhh(null)` | `ArgumentNullException` | input-validation |
| M16 | EHH empty | `CalculateEhh([])` | 0.0 | boundary of formula |
| M17 | ScanForSelection null/bad window | null scores; windowSize 0 | `ArgumentNullException`; `ArgumentOutOfRangeException` | input-validation |
| M18 | Standardize null/bad bin | null scores; binCount 0 | `ArgumentNullException`; `ArgumentOutOfRangeException` | input-validation |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | rehh ratio | unstandardized iHS from IHH_A=284429.9, IHH_D=2057107.4 (component-level via ln) | −1.978569274 | reference-implementation worked value |
| S2 | Two-bin standardization | scores split across two freq bins standardized independently | each bin centered at 0 | Voight per-bin |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | Sign-symmetry property | ln(iHH_A/iHH_D) = −ln(iHH_D/iHH_A) over random valid panels | holds within 1e-10 | property-based (INV-4) |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- `tests/Seqeron/Seqeron.Genomics.Tests/PopulationGeneticsAnalyzerTests.cs` contains legacy tests for the pre-existing `CalculateIHS(ehh0, ehh1, positions)` overload and the region-based `ScanForSelection(...)`. These exercise *different* (non-canonical, EHH-curve-as-input) overloads and use permissive assertions (`Is.GreaterThan(0)`); they do not cover the canonical haplotype-based pipeline of this unit.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1 EHH worked values | ❌ Missing | |
| M2 EHH single/all-distinct | ❌ Missing | |
| M3 iHS constructed panel | ❌ Missing | |
| M4 Balanced decay | ❌ Missing | legacy `CalculateIHS_BalancedEHH` tests a different overload |
| M5 Sign convention | ❌ Missing | |
| M6 Standardization | ❌ Missing | |
| M7 Standardization single bin | ❌ Missing | |
| M8 Scan proportion | ❌ Missing | |
| M9 Scan windows | ❌ Missing | |
| M10–M18 edge cases | ❌ Missing | |
| S1 rehh ratio | ❌ Missing | |
| S2 Two-bin | ❌ Missing | |
| C1 Sign-symmetry property | ❌ Missing | |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/PopulationGeneticsAnalyzer_SelectionSignature_Tests.cs` — all canonical iHS-pipeline tests for POP-SELECT-001.
- **Remove:** nothing. The legacy `PopulationGeneticsAnalyzerTests.cs` covers different overloads still consumed by the MCP layer; leaving it untouched keeps scope isolated.

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| `PopulationGeneticsAnalyzer_SelectionSignature_Tests.cs` | Canonical POP-SELECT-001 fixture | 23 |

### 5.5 Phase 7 Work Queue

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| 1 | M1 | ❌ Missing | Implemented | ✅ Done |
| 2 | M2 | ❌ Missing | Implemented | ✅ Done |
| 3 | M3 | ❌ Missing | Implemented | ✅ Done |
| 4 | M4 | ❌ Missing | Implemented | ✅ Done |
| 5 | M5 | ❌ Missing | Implemented | ✅ Done |
| 6 | M6 | ❌ Missing | Implemented | ✅ Done |
| 7 | M7 | ❌ Missing | Implemented | ✅ Done |
| 8 | M8 | ❌ Missing | Implemented | ✅ Done |
| 9 | M9 | ❌ Missing | Implemented | ✅ Done |
| 10 | M10 | ❌ Missing | Implemented | ✅ Done |
| 11 | M11 | ❌ Missing | Implemented | ✅ Done |
| 12 | M12 | ❌ Missing | Implemented | ✅ Done |
| 13 | M13 | ❌ Missing | Implemented | ✅ Done |
| 14 | M14 | ❌ Missing | Implemented | ✅ Done |
| 15 | M15 | ❌ Missing | Implemented | ✅ Done |
| 16 | M16 | ❌ Missing | Implemented | ✅ Done |
| 17 | M17 | ❌ Missing | Implemented | ✅ Done |
| 18 | M18 | ❌ Missing | Implemented | ✅ Done |
| 19 | S1 | ❌ Missing | Implemented | ✅ Done |
| 20 | S2 | ❌ Missing | Implemented | ✅ Done |
| 21 | C1 | ❌ Missing | Implemented (property) | ✅ Done |

**Total items:** 21
**✅ Done:** 21 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ Covered | CalculateEhh_WorkedValues_MatchSelscanFormula |
| M2 | ✅ Covered | CalculateEhh_SingleAndDistinct_ReturnExtremes |
| M3 | ✅ Covered | CalculateIHS_ConstructedPanel_MatchesDerivedValues |
| M4 | ✅ Covered | CalculateIHS_BalancedDecay_ReturnsZero |
| M5 | ✅ Covered | CalculateIHS_LongDerivedHaplotype_ReturnsNegative |
| M6 | ✅ Covered | StandardizeIHS_SingleBin_CentersToZeroUnitSd |
| M7 | ✅ Covered | StandardizeIHS_SingletonBin_ReturnsZero |
| M8 | ✅ Covered | ScanForSelection_HalfExtreme_ReturnsProportionHalf |
| M9 | ✅ Covered | ScanForSelection_FiveScoresWindowTwo_ReturnsThreeWindows |
| M10 | ✅ Covered | CalculateIHS_NullHaplotypes_Throws |
| M11 | ✅ Covered | CalculateIHS_MonomorphicCore_Throws |
| M12 | ✅ Covered | CalculateIHS_InconsistentLength_Throws |
| M13 | ✅ Covered | CalculateIHS_InvalidAllele_Throws |
| M14 | ✅ Covered | CalculateIHS_CoreIndexOutOfRange_Throws |
| M15 | ✅ Covered | CalculateEhh_Null_Throws |
| M16 | ✅ Covered | CalculateEhh_Empty_ReturnsZero |
| M17 | ✅ Covered | ScanForSelection_NullOrBadWindow_Throws |
| M18 | ✅ Covered | StandardizeIHS_NullOrBadBin_Throws |
| S1 | ✅ Covered | CalculateIHS_RehhRatio_MatchesReference |
| S2 | ✅ Covered | StandardizeIHS_TwoBins_StandardizedIndependently |
| C1 | ✅ Covered | CalculateIHS_SignSymmetry_Property |

**In-scope cases:** 21 | **✅:** 21

---

## 6. Assumption Register

**Total assumptions:** 2

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | SD estimator = sample (N−1) standard deviation within a frequency bin (Voight does not specify N vs N−1) | StandardizeIHS, M6, M7 |
| 2 | Frequency-bin width = 0.05 (20 bins) by default, matching rehh `freqbin` | StandardizeIHS, S2 |

Neither assumption affects the canonical, fully evidence-backed value (unstandardized iHS); both affect only the magnitude scaling of standardized scores and are exposed as explicit parameters. They are not correctness-affecting for the core score.

---

## 7. Open Questions / Decisions

1. **Decision:** Followed the Voight et al. (2006) sign convention ln(iHH_A/iHH_D); selscan's reciprocal is documented in Evidence. The conflict is a sign flip only and is resolved in favor of the primary peer-reviewed source.
2. **Decision:** Existing non-canonical `CalculateIHS(ehh0, ehh1, positions)` overload (consumed by the MCP layer) is retained unchanged; the canonical haplotype-based overload is added alongside it.
