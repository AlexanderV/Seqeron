# Test Specification: SEQ-DINUC-001

**Test Unit ID:** SEQ-DINUC-001
**Area:** Statistics
**Algorithm:** Dinucleotide Analysis (frequencies, observed/expected relative abundance, codon frequencies)
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-13

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Karlin S. — Pervasive properties of the genomic signature (PMC126251) | 1 | https://pmc.ncbi.nlm.nih.gov/articles/PMC126251/ | 2026-06-13 |
| 2 | Karlin & Burge (1995) criterion via MBE 19(6):964 | 1 | https://academic.oup.com/mbe/article/19/6/964/1095097 | 2026-06-13 |
| 3 | Gardiner-Garden & Frommer (1987) CpG O/E (restated) | 1/4 | https://link.springer.com/article/10.1007/BF00162972 | 2026-06-13 |
| 4 | Kazusa Codon Usage Database (CUTG) readme | 5/2 | https://www.kazusa.or.jp/codon/readme_codon.html | 2026-06-13 |

### 1.2 Key Evidence Points

1. ρ_XY = f_XY / (f_X · f_Y), with f normalized frequencies; ρ = 1 means no bias — Karlin (PMC126251).
2. Classification (interpretive only): under-represented if ρ ≤ 0.78, over-represented if ρ ≥ 1.23 — Karlin & Burge 1995 (MBE 19(6):964).
3. Dinucleotide frequency f_XY is a normalized frequency = count / (number of dinucleotide positions = N−1) — Karlin (PMC126251).
4. Codon frequency = count of a codon over total codons, reading consecutive **non-overlapping** triplets from a frame; ambiguous/non-ACGT triplets excluded — Kazusa CUTG.

### 1.3 Documented Corner Cases

- ρ is uninformative when a constituent base is absent (expected = 0 ⇒ division by zero); repository returns 0 for that dinucleotide — Karlin (independence baseline).
- Codon counting ignores trailing 1–2 leftover bases and excludes non-ACGT triplets — Kazusa CUTG.

### 1.4 Known Failure Modes / Pitfalls

1. Mixing the (N−1) (Karlin) and N (Gardiner-Garden) normalizations changes the numeric ratio by N/(N−1) — Evidence §Assumptions; repository uses Karlin (N−1).
2. Forgetting to exclude non-ACGT triplets/dinucleotides inflates totals — Kazusa CUTG; implementation filters on {A,T,G,C} (and U for dinucleotides).

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `CalculateDinucleotideRatios(string)` | SequenceStatistics | Canonical | ρ_XY = f_XY/(f_X f_Y), Karlin odds ratio |
| `CalculateDinucleotideFrequencies(string)` | SequenceStatistics | Canonical | normalized dinucleotide frequencies, count/(N−1) |
| `CalculateCodonFrequencies(string,int)` | SequenceStatistics | Canonical | non-overlapping triplet frequencies per frame |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | Dinucleotide frequencies over {A,T,G,C,U} sum to 1.0 (when ≥1 valid dinucleotide) | Yes | Karlin (normalized frequency), PMC126251 |
| INV-2 | Codon frequencies sum to 1.0 (when ≥1 valid codon) | Yes | Kazusa CUTG (count/total) |
| INV-3 | ρ_XY = 1.0 ⇔ observed equals product of base frequencies (no bias); ρ ≥ 0 | Yes | Karlin (r=1 no bias), PMC126251 |
| INV-4 | All returned frequencies and ratios are finite and ≥ 0 | Yes | definitions are non-negative ratios |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | Ratios exact `ATGCGCGT` | ρ for GC,CG,AT,TG,GT on hand-derived sequence | ρ_GC=ρ_CG=64/21=3.047619047619048; ρ_AT=32/7=4.571428571428571; ρ_TG=ρ_GT=32/21=1.523809523809524 | Karlin PMC126251 + dataset |
| M2 | Freqs exact `ATGCGCGT` | dinucleotide frequencies, count/(N−1) | GC=CG=2/7=0.2857142857142857; AT=TG=GT=1/7=0.14285714285714285 | Karlin (normalized freq) |
| M3 | Freqs sum to 1 | INV-1 on `ATGCGCGT` | Σ freq = 1.0 | Karlin (normalized freq) |
| M4 | Codon frame 0 | `ATGATGAAA` frame 0 | ATG=2/3, AAA=1/3; codons={ATG,AAA} | Kazusa CUTG |
| M5 | Codon frame 1 | `ATGATGAAA` frame 1 | TGA=1.0; only key TGA | Kazusa CUTG |
| M6 | Codon sum to 1 | INV-2 on `ATGATGAAA` frame 0 | Σ freq = 1.0 | Kazusa CUTG |
| M7 | No-bias baseline ρ=1 | Homopolymer `AAAA`: f_AA=1, f_A=1 ⇒ ρ_AA=1/(1·1) | ρ_AA = 1.0 exactly | Karlin (r=1 no bias), PMC126251 |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | Ratios null/empty/<2 | null, "", "A" | empty dictionary | input guard |
| S2 | Freqs null/empty/<2 | null, "", "A" | empty dictionary | input guard |
| S3 | Codons null/empty/<3 | null, "", "AT" | empty dictionary | input guard |
| S4 | Division-by-zero guard | sequence missing a base so expected=0 | ratio = 0 for that dinucleotide | expected=0 guard |
| S5 | Codon non-ACGT excluded | `ATGNNNAAA` frame 0 | only ATG, AAA counted; NNN excluded; ATG=AAA=0.5 | Kazusa CUTG |
| S6 | Codon trailing bases ignored | `ATGAA` frame 0 | only ATG; ATG=1.0 | Kazusa CUTG |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | Case-insensitive | `atgcgcgt` ratios equal uppercase | identical to M1 | ToUpperInvariant |
| C2 | RNA U handling | `AUGCGC` frequencies include AU | AU present, U counted | RNA support |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- Methods live in `src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/SequenceStatistics.cs` (lines 599, 633, 670).
- No dedicated `SequenceStatistics_CalculateDinucleotide*_Tests.cs` existed. The legacy `SequenceStatisticsTests.cs` does not cover these three methods. So all planned cases start as ❌ Missing.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1 | ❌ Missing | no ratio test existed |
| M2 | ❌ Missing | no frequency test existed |
| M3 | ❌ Missing | |
| M4 | ❌ Missing | no codon test existed |
| M5 | ❌ Missing | |
| M6 | ❌ Missing | |
| M7 | ❌ Missing | |
| S1 | ❌ Missing | |
| S2 | ❌ Missing | |
| S3 | ❌ Missing | |
| S4 | ❌ Missing | |
| S5 | ❌ Missing | |
| S6 | ❌ Missing | |
| C1 | ❌ Missing | |
| C2 | ❌ Missing | |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/SequenceStatistics_CalculateDinucleotide_Tests.cs` — all cases for the three methods, one `#region` per method.
- **Remove:** nothing (no pre-existing duplicate tests for these methods).

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| `SequenceStatistics_CalculateDinucleotide_Tests.cs` | Canonical, all cases | 15 |

### 5.5 Phase 7 Work Queue

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| 1 | M1 | ❌ Missing | Implemented exact ρ test | ✅ Done |
| 2 | M2 | ❌ Missing | Implemented exact freq test | ✅ Done |
| 3 | M3 | ❌ Missing | Implemented sum-to-1 test | ✅ Done |
| 4 | M4 | ❌ Missing | Implemented codon frame 0 test | ✅ Done |
| 5 | M5 | ❌ Missing | Implemented codon frame 1 test | ✅ Done |
| 6 | M6 | ❌ Missing | Implemented codon sum-to-1 test | ✅ Done |
| 7 | M7 | ❌ Missing | Implemented no-bias baseline test | ✅ Done |
| 8 | S1 | ❌ Missing | Implemented ratios guard test | ✅ Done |
| 9 | S2 | ❌ Missing | Implemented freqs guard test | ✅ Done |
| 10 | S3 | ❌ Missing | Implemented codons guard test | ✅ Done |
| 11 | S4 | ❌ Missing | Implemented div-by-zero guard test | ✅ Done |
| 12 | S5 | ❌ Missing | Implemented non-ACGT exclusion test | ✅ Done |
| 13 | S6 | ❌ Missing | Implemented trailing-bases test | ✅ Done |
| 14 | C1 | ❌ Missing | Implemented case-insensitivity test | ✅ Done |
| 15 | C2 | ❌ Missing | Implemented RNA U test | ✅ Done |

**Total items:** 15
**✅ Done:** 15 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ Covered | exact ρ rationals |
| M2 | ✅ Covered | exact freq rationals |
| M3 | ✅ Covered | INV-1 |
| M4 | ✅ Covered | codon frame 0 |
| M5 | ✅ Covered | codon frame 1 |
| M6 | ✅ Covered | INV-2 |
| M7 | ✅ Covered | no-bias baseline |
| S1 | ✅ Covered | ratios guards |
| S2 | ✅ Covered | freqs guards |
| S3 | ✅ Covered | codons guards |
| S4 | ✅ Covered | div-by-zero |
| S5 | ✅ Covered | non-ACGT excluded |
| S6 | ✅ Covered | trailing bases |
| C1 | ✅ Covered | case-insensitive |
| C2 | ✅ Covered | RNA U |

**In-scope cases:** 15 | **✅:** 15

---

## 6. Assumption Register

**Total assumptions:** 2

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | Karlin (N−1) normalization for dinucleotide frequency (vs Gardiner-Garden N) | M1, M2 — both conventions authoritative; repository uses Karlin |
| 2 | U treated as a fifth base in `CalculateDinucleotideRatios` base-frequency denominator | C2 — RNA support |

---

## 7. Open Questions / Decisions

1. Decision: methods return raw frequencies/ratios; the 0.78/1.23 Karlin & Burge classification thresholds are documentation-only (not applied in code), so no threshold constant is needed in the implementation.
