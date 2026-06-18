# Test Specification: EPIGEN-BISULF-001

**Test Unit ID:** EPIGEN-BISULF-001
**Area:** Epigenetics
**Algorithm:** Bisulfite Sequencing Analysis (conversion simulation, methylation calling, profile aggregation)
**Status:** ☐ In Progress
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-13

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Frommer et al. (1992) PNAS 89:1827–1831 | 1 | https://doi.org/10.1073/pnas.89.5.1827 | 2026-06-13 |
| 2 | Krueger & Andrews (2011) Bioinformatics 27:1571–1572 | 1/3 | https://doi.org/10.1093/bioinformatics/btr167 | 2026-06-13 |
| 3 | Bismark User Guide v0.15.0 | 3 | https://www.bioinformatics.babraham.ac.uk/projects/bismark/Bismark_User_Guide_v0.15.0.pdf | 2026-06-13 |
| 4 | Schultz et al. (2012) Trends Genet. 28:583–585 | 1 | https://doi.org/10.1016/j.tig.2012.10.012 | 2026-06-13 |

### 1.2 Key Evidence Points

1. Bisulfite converts unmethylated cytosine → uracil (reads as thymine), but 5-methylcytosine remains nonreactive (stays C) — Frommer et al. (1992).
2. Methylation call: a C remaining at a reference-C position = methylated; a T at that position = unmethylated — Krueger & Andrews (2011).
3. Methylation percentage = 100 × methylated / (methylated + unmethylated); as a fraction = meth/(meth+unmeth) — Bismark User Guide v0.15.0.
4. Weighted methylation level = Σ(methylated reads) / Σ(methylated+unmethylated reads) over a context = Σ(level·coverage)/Σ(coverage) — Schultz et al. (2012).
5. Non-cytosine bases (A,G,T) are unaffected by bisulfite — Frommer et al. (1992).

### 1.3 Documented Corner Cases

- Zero-coverage cytosine: percentage undefined (denominator 0) → excluded from output (Bismark).
- Read bases past the reference end / outside reference are ignored; last reference base cannot start a CpG (Bismark/Krueger).
- A read base that is neither C nor T at a reference-C position is not a valid bisulfite call → excluded (Krueger & Andrews 2011).

### 1.4 Known Failure Modes / Pitfalls

1. Treating the percentage as an unweighted mean of per-site fractions under unequal coverage — Schultz et al. (2012) require the weighted (read-pooled) level.
2. Converting protected (methylated) cytosines — Frommer et al. (1992): 5mC must remain C.

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `SimulateBisulfiteConversion(sequence, methylatedPositions)` | EpigeneticsAnalyzer | Canonical | C→T conversion; protected positions stay C; non-C unchanged. |
| `CalculateMethylationFromBisulfite(referenceSequence, bisulfiteReads)` | EpigeneticsAnalyzer | Canonical | Per-CpG level = meth/(meth+unmeth); coverage = valid calls. |
| `GenerateMethylationProfile(sites)` | EpigeneticsAnalyzer | Canonical | Per-context weighted methylation (Schultz). |

<!-- Registry lists CalculateMethylationFromBisulfite(bsSeq, refSeq). The implemented
signature is (referenceSequence, bisulfiteReads): a single converted string cannot yield
per-site coverage, so per-read input is required to compute meth/(meth+unmeth). The
implemented signature is retained as the correct contract; conflict noted in §7. -->

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | Conversion output length equals input length | Yes | Frommer 1992 (per-base substitution) |
| INV-2 | Every protected (methylated) cytosine remains `C`/`c`; every non-protected `C`→`T`, `c`→`t` | Yes | Frommer 1992 |
| INV-3 | Non-cytosine bases are returned unchanged | Yes | Frommer 1992 |
| INV-4 | Every reported methylation level ∈ [0,1]; coverage ≥ 1 | Yes | Bismark (meth/(meth+unmeth)) |
| INV-5 | Profile per-context level = Σ(level·coverage)/Σ(coverage) over that context | Yes | Schultz 2012 |
| INV-6 | Sites with zero coverage are excluded from `CalculateMethylationFromBisulfite` output | Yes | Bismark (undefined percentage) |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | Convert unmethylated, no methylation | `ACGTCGAA`, no protected positions | `ATGTTGAA` (both C → T, rest unchanged) | Frommer 1992 |
| M2 | Convert with one protected C | `ACGTCGAA`, methylated `{1}` | `ACGTTGAA` (C@1 stays C, C@4→T) | Frommer 1992 |
| M3 | Non-cytosine unchanged | `AGTAGT` (no C) | `AGTAGT` (identical) | Frommer 1992 |
| M4 | Lowercase unmethylated c → t | `acgt`, no protected | `atgt` (c@1→t) | Frommer 1992 (case-preserving) |
| M5 | Methylation level half | ref `ACGTACGT`; reads `C`@1 and `T`@1 | site@1 level 0.5, coverage 2 | Bismark % formula |
| M6 | Methylation level zero | ref `ACGTACGT`; read `T`@5 | site@5 level 0.0, coverage 1 | Bismark % formula |
| M7 | Methylation level one | ref `ACGT`; read `C`@1 | site@1 level 1.0, coverage 1 | Bismark % formula |
| M8 | Weighted profile ≠ unweighted | CpG sites (1.0,cov10) and (0.0,cov30) | CpGMethylation = 0.25 | Schultz 2012 |
| M9 | Per-context separation | one CpG (1.0) + one CHG (0.0), equal cov | CpGMethylation 1.0, CHGMethylation 0.0 | Schultz 2012 / Bismark contexts |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | Empty conversion input | `""` | `""` | Input validation |
| S2 | Null conversion input | `null` | `""` | Input validation |
| S3 | Zero-coverage site excluded | ref with CpG but no reads covering it | site absent from output | INV-6 / Bismark |
| S4 | Empty profile | no sites | all-zero `MethylationProfile` | Validation |
| S5 | Non-C/T read base ignored | read has `A` at CpG C position | that call not counted | Krueger 2011 |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | Read extends past reference | read longer than remaining ref | extra bases ignored | Bismark boundary |
| C2 | Conversion output length invariant | random-shaped input | length unchanged (INV-1) | Property test |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- `tests/Seqeron/Seqeron.Genomics.Tests/Snapshots/EpigeneticsSnapshotTests.cs` has a `BisulfiteConversion_KnownSequence_MatchesSnapshot` snapshot test (shape-only, snapshot-based; not evidence-asserted).
- No dedicated `EpigeneticsAnalyzer_Bisulfite_Tests.cs` exists. The three methods in scope have no exact-value evidence-based unit tests.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1 | ❌ Missing | new |
| M2 | ❌ Missing | new |
| M3 | ❌ Missing | new |
| M4 | ❌ Missing | new |
| M5 | ❌ Missing | new |
| M6 | ❌ Missing | new |
| M7 | ❌ Missing | new |
| M8 | ❌ Missing | new |
| M9 | ❌ Missing | new |
| S1 | ❌ Missing | new |
| S2 | ❌ Missing | new |
| S3 | ❌ Missing | new |
| S4 | ❌ Missing | new |
| S5 | ❌ Missing | new |
| C1 | ❌ Missing | new |
| C2 | ❌ Missing | new |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/EpigeneticsAnalyzer_Bisulfite_Tests.cs` — all evidence-based tests for the three methods.
- **Remove:** none. The existing snapshot test in `EpigeneticsSnapshotTests.cs` is a separate regression-snapshot harness for the whole class and is left untouched.

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| `EpigeneticsAnalyzer_Bisulfite_Tests.cs` | Canonical evidence-based unit tests | 16 |
| `EpigeneticsSnapshotTests.cs` | Pre-existing snapshot regression (untouched) | unchanged |

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
| 10 | S1 | ❌ Missing | Implemented | ✅ Done |
| 11 | S2 | ❌ Missing | Implemented | ✅ Done |
| 12 | S3 | ❌ Missing | Implemented | ✅ Done |
| 13 | S4 | ❌ Missing | Implemented | ✅ Done |
| 14 | S5 | ❌ Missing | Implemented | ✅ Done |
| 15 | C1 | ❌ Missing | Implemented | ✅ Done |
| 16 | C2 | ❌ Missing | Implemented | ✅ Done |

**Total items:** 16
**✅ Done:** 16 | **⛔ Blocked:** 0 | **Remaining:** must be 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ Covered | `Convert_UnmethylatedSequence_AllCytosinesBecomeT` |
| M2 | ✅ Covered | `Convert_MethylatedPositionProtected_StaysCytosine` |
| M3 | ✅ Covered | `Convert_NoCytosines_ReturnsUnchanged` |
| M4 | ✅ Covered | `Convert_LowercaseUnmethylated_BecomesLowercaseT` |
| M5 | ✅ Covered | `CalculateMethylation_HalfMethylatedCpG_ReturnsHalf` |
| M6 | ✅ Covered | `CalculateMethylation_AllUnmethylated_ReturnsZero` |
| M7 | ✅ Covered | `CalculateMethylation_FullyMethylated_ReturnsOne` |
| M8 | ✅ Covered | `GenerateProfile_UnequalCoverage_UsesWeightedLevel` |
| M9 | ✅ Covered | `GenerateProfile_SeparatesContexts_ReportsPerContextLevels` |
| S1 | ✅ Covered | `Convert_EmptyInput_ReturnsEmpty` |
| S2 | ✅ Covered | `Convert_NullInput_ReturnsEmpty` |
| S3 | ✅ Covered | `CalculateMethylation_NoReadsCoverSite_SiteExcluded` |
| S4 | ✅ Covered | `GenerateProfile_NoSites_ReturnsZeroProfile` |
| S5 | ✅ Covered | `CalculateMethylation_NonCTReadBase_NotCounted` |
| C1 | ✅ Covered | `CalculateMethylation_ReadPastReferenceEnd_ExtraBasesIgnored` |
| C2 | ✅ Covered | `Convert_AnyInput_PreservesLength` |

---

## 6. Assumption Register

**Total assumptions:** 2

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | Positions are 0-based offsets into the supplied sequence (API convention) | All methods |
| 2 | `SimulateBisulfiteConversion` converts only the supplied strand (no complementary-strand merge) | M1–M4 |

---

## 7. Open Questions / Decisions

1. **Registry signature conflict:** the Registry lists `CalculateMethylationFromBisulfite(bsSeq, refSeq)`. The implemented and correct contract is `(referenceSequence, bisulfiteReads)`, because per-site coverage and methylation fraction require per-read input (a single converted string carries no read multiplicity). Decision: retain the implemented signature; record the conflict here.
2. No correctness-affecting assumptions remain; both assumptions are API-shape conventions.
