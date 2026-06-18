# Test Specification: POP-ROH-001

**Test Unit ID:** POP-ROH-001
**Area:** PopGen
**Algorithm:** Runs of Homozygosity (ROH) detection and genomic inbreeding coefficient F_ROH
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-13

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | McQuillan et al. (2008), Am J Hum Genet 83(3):359-372 | 1 | https://pmc.ncbi.nlm.nih.gov/articles/PMC2556426/ | 2026-06-13 |
| 2 | Chang et al. (2015) — PLINK 1.9 --homozyg docs | 3 | https://www.cog-genomics.org/plink/1.9/ibd | 2026-06-13 |
| 3 | Marras et al. (2015), Anim Genet 46(2):110-121 — detectRUNS consecutive method | 3 | https://cran.r-project.org/web/packages/detectRUNS/vignettes/detectRUNS.vignette.html | 2026-06-13 |

### 1.2 Key Evidence Points

1. F_ROH = ΣL_roh / L_auto, the total ROH length divided by SNP-covered autosomal genome length — McQuillan et al. (2008), verbatim "Froh = ∑Lroh/Lauto".
2. McQuillan used L_auto = 2,673,768 kb — McQuillan et al. (2008).
3. Consecutive-runs method scans SNP by SNP (window-free) and breaks a run when opposite genotypes exceed maxOppRun, missing exceed maxMissRun, or an inter-SNP gap exceeds maxGap — Marras et al. (2015) / detectRUNS vignette.
4. PLINK defaults: --homozyg-snp 100, --homozyg-kb 1000 (1,000,000 bp), --homozyg-window-het 1, --homozyg-gap 1000 kb — Chang et al. (2015).
5. Both SNP count and physical length thresholds must be satisfied for a run to be retained — Chang et al. (2015); Marras et al. (2015).

### 1.3 Documented Corner Cases

- A small number (maxOppRun) of heterozygous calls is tolerated inside a run (Marras 2015).
- An inter-SNP gap larger than maxGap breaks an all-homozygous run (PLINK --homozyg-gap; Marras 2015).
- Stretches failing minSNP or minLengthBps are discarded (PLINK; Marras 2015).

### 1.4 Known Failure Modes / Pitfalls

1. Passing only one of the two thresholds (count OR length) must NOT retain a run — Chang et al. (2015).
2. A heterozygous SNP cannot itself seed a homozygous run — Marras et al. (2015) consecutive scan semantics.

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `FindROH(genotypes, minSnps, minLength, maxHeterozygotes, maxGap)` | PopulationGeneticsAnalyzer | Canonical | Consecutive-runs ROH detection (Marras 2015; PLINK defaults) |
| `CalculateInbreedingFromROH(rohSegments, genomeLength)` | PopulationGeneticsAnalyzer | Canonical | F_ROH = ΣL_roh / L_auto (McQuillan 2008) |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | Every reported run has SnpCount ≥ minSnps and End − Start ≥ minLength | Yes | Chang et al. (2015) |
| INV-2 | A reported run contains at most maxHeterozygotes opposite genotypes and no inter-SNP gap > maxGap | Yes | Marras et al. (2015) |
| INV-3 | Runs are emitted in ascending Start order, Start ≤ End | Yes | consecutive ascending scan |
| INV-4 | 0 ≤ F_ROH ≤ 1 when all ROH lie within [0, genomeLength] | Yes | McQuillan et al. (2008) |
| INV-5 | F_ROH is linear in ΣL_roh: ΣL_roh / L_auto exactly | Yes | McQuillan et al. (2008) |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | Single homozygous run | 100 hom SNPs @20kb | 1 run, Start=0, End=1,980,000, SnpCount=100 | Marras 2015; PLINK defaults |
| M2 | One tolerated het | 101 SNPs, het at idx 50, maxHet=1 | 1 run, Start=0, End=2,000,000, SnpCount=101 | Marras 2015 maxOppRun |
| M3 | Het beyond tolerance splits | 201 SNPs, hets at 50 & 100, maxHet=1 | 2 runs [0..99] and [101..200], 100 SNPs each | Marras 2015 termination |
| M4 | Gap > maxGap breaks | two 60-SNP blocks across 2 Mb jump | 2 runs of 60 SNPs | PLINK --homozyg-gap; Marras 2015 |
| M5 | < minSnps discarded | 20 hom SNPs, default minSnps=100 | empty | PLINK --homozyg-snp |
| M6 | < minLength discarded | 120 SNPs @1kb (119 kb) | empty | PLINK --homozyg-kb |
| M7 | F_ROH two segments | ΣL_roh=20M, L_auto=100M | 0.20 | McQuillan 2008 |
| M8 | F_ROH whole genome | ROH = L_auto | 1.0 | McQuillan 2008 |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | Unsorted input | reversed-position input | same single run, Start=0 End=1,980,000 | internal ordering |
| S2 | Empty input | no genotypes | empty | boundary |
| S3 | Homozygous-alt (2) | all genotype 2 | 1 run, 100 SNPs | encoding |
| S4 | Leading heterozygotes | first two SNPs het | run starts at idx 2 | het cannot seed run |
| S5 | No ROH segments | empty F_ROH input | 0.0 | boundary |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | Invalid FindROH args | null, minSnps=0, neg minLength/maxHet/maxGap | throws | failure modes |
| C2 | Invalid F_ROH input | genomeLength ≤ 0; null segments | 0.0; throws | failure modes |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- Pre-existing weak tests in `tests/Seqeron/Seqeron.Genomics.Tests/PopulationGeneticsAnalyzerTests.cs` (region "Inbreeding Tests"): `CalculateInbreedingFromROH_NoROH_ReturnsZero`, `CalculateInbreedingFromROH_WithROH_CalculatesCorrectly`, `FindROH_LongHomozygousRun_DetectsROH`, `FindROH_TooManyHeterozygotes_NoROH`, `FindROH_ShortRun_NotDetected`.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| FindROH single run (M1) | ⚠ Weak | legacy `FindROH_LongHomozygousRun_DetectsROH` checks Count==1 only, no bounds/SnpCount |
| FindROH tolerated het (M2) | ❌ Missing | not covered |
| FindROH split (M3) | ❌ Missing | not covered |
| FindROH gap break (M4) | ❌ Missing | not covered |
| FindROH < minSnps (M5) | ⚠ Weak | legacy `FindROH_ShortRun_NotDetected` count/existence only |
| FindROH < minLength (M6) | ❌ Missing | not covered |
| F_ROH two segments (M7) | ⚠ Weak | legacy uses `.Within(0.001)` permissive tolerance |
| F_ROH whole genome (M8) | ❌ Missing | not covered |
| Unsorted (S1) | ❌ Missing | not covered |
| Empty (S2) | ❌ Missing | not covered |
| Hom-alt (S3) | ❌ Missing | not covered |
| Leading het (S4) | ❌ Missing | not covered |
| No segments (S5) | ⚠ Weak | legacy `CalculateInbreedingFromROH_NoROH_ReturnsZero` no message |
| Invalid FindROH args (C1) | ❌ Missing | not covered |
| Invalid F_ROH input (C2) | ❌ Missing | not covered |
| Legacy `FindROH_TooManyHeterozygotes_NoROH` | 🔁 Duplicate | superseded by M3 |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/PopulationGeneticsAnalyzer_FindROH_Tests.cs` — all ROH and F_ROH cases with exact evidence-based values.
- **Remove:** the five legacy tests in the "Inbreeding Tests" region of `PopulationGeneticsAnalyzerTests.cs` (weak/duplicate); replaced by a pointer comment.

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| `PopulationGeneticsAnalyzer_FindROH_Tests.cs` | Canonical POP-ROH-001 | 15 |
| `PopulationGeneticsAnalyzerTests.cs` | Legacy (ROH tests removed) | 0 ROH tests |

### 5.5 Phase 7 Work Queue

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| 1 | M1 | ⚠ Weak | rewrote with exact bounds + SnpCount | ✅ Done |
| 2 | M2 | ❌ Missing | implemented | ✅ Done |
| 3 | M3 | ❌ Missing | implemented | ✅ Done |
| 4 | M4 | ❌ Missing | implemented | ✅ Done |
| 5 | M5 | ⚠ Weak | rewrote with message | ✅ Done |
| 6 | M6 | ❌ Missing | implemented | ✅ Done |
| 7 | M7 | ⚠ Weak | rewrote with `.Within(1e-10)` | ✅ Done |
| 8 | M8 | ❌ Missing | implemented | ✅ Done |
| 9 | S1 | ❌ Missing | implemented | ✅ Done |
| 10 | S2 | ❌ Missing | implemented | ✅ Done |
| 11 | S3 | ❌ Missing | implemented | ✅ Done |
| 12 | S4 | ❌ Missing | implemented | ✅ Done |
| 13 | S5 | ⚠ Weak | rewrote with message + tol | ✅ Done |
| 14 | C1 | ❌ Missing | implemented | ✅ Done |
| 15 | C2 | ❌ Missing | implemented | ✅ Done |
| 16 | Legacy het dup | 🔁 Duplicate | removed legacy tests | ✅ Done |

**Total items:** 16
**✅ Done:** 16 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ Covered | `FindROH_SingleHomozygousRun_ReportsExactBoundsAndCount` |
| M2 | ✅ Covered | `FindROH_OneToleratedHeterozygote_KeepsSingleRun` |
| M3 | ✅ Covered | `FindROH_HeterozygoteBeyondTolerance_SplitsIntoTwoRuns` |
| M4 | ✅ Covered | `FindROH_GapExceedsMaxGap_BreaksRun` |
| M5 | ✅ Covered | `FindROH_FewerThanMinSnps_NotReported` |
| M6 | ✅ Covered | `FindROH_ShorterThanMinLength_NotReported` |
| M7 | ✅ Covered | `CalculateInbreedingFromROH_TwoSegments_MatchesFRohFormula` |
| M8 | ✅ Covered | `CalculateInbreedingFromROH_WholeGenomeROH_ReturnsOne` |
| S1 | ✅ Covered | `FindROH_UnsortedInput_OrdersByPosition` |
| S2 | ✅ Covered | `FindROH_EmptyInput_ReturnsEmpty` |
| S3 | ✅ Covered | `FindROH_HomozygousAlternateGenotype_CountsAsHomozygous` |
| S4 | ✅ Covered | `FindROH_LeadingHeterozygotes_RunStartsAtFirstHomozygous` |
| S5 | ✅ Covered | `CalculateInbreedingFromROH_NoSegments_ReturnsZero` |
| C1 | ✅ Covered | `FindROH_InvalidArguments_Throw` |
| C2 | ✅ Covered | `CalculateInbreedingFromROH_InvalidInput_HandledPerContract` |

**Total in-scope cases:** 15 | **✅:** 15

---

## 6. Assumption Register

**Total assumptions:** 2

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | Genotype encoding 0=hom-ref, 1=het, 2=hom-alt (API convention; not correctness-affecting — rule is hom vs opposite) | FindROH input semantics |
| 2 | Missing-genotype handling out of scope (no missing sentinel in input) | FindROH limitation |

---

## 7. Open Questions / Decisions

1. PLINK's sliding-window two-phase scan (--homozyg-window-*) is intentionally not reproduced; the simpler window-free consecutive-runs method (Marras 2015) is implemented, which is itself an authoritative reference method and matches the registry O(n) complexity. Documented in the algorithm doc §5.3.
