# Test Specification: CODON-RSCU-001

**Test Unit ID:** CODON-RSCU-001
**Area:** Codon
**Algorithm:** Relative Synonymous Codon Usage (RSCU) and codon counting
**Status:** ☐ In Progress
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-13

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Sharp, Tuohy & Mosurski (1986), NAR 14(13):5125-5143 — introduced RSCU | 1 | https://doi.org/10.1093/nar/14.13.5125 (PMC: https://pmc.ncbi.nlm.nih.gov/articles/PMC311530/) | 2026-06-13 |
| 2 | LIRMM "RSCU RS" — explicit formula RSCU_i = X_i/((1/N_i)ΣX_j) | 3 | https://www.lirmm.fr/~rivals/rscu/ | 2026-06-13 |
| 3 | GenomicSig (CRAN) RSCU — RSCU_{i,j} = (n_i·x_{i,j})/Σx; bounds [0,n_i] | 3 | https://rdrr.io/cran/GenomicSig/man/RSCU.html | 2026-06-13 |
| 4 | seqinr `uco` — definition + no-bias value 1.00 | 3 | https://search.r-project.org/CRAN/refmans/seqinr/html/uco.html | 2026-06-13 |
| 5 | cubar `est_rscu` — pseudocount / zero-count handling | 3 | https://rdrr.io/cran/cubar/man/est_rscu.html | 2026-06-13 |
| 6 | Begomovirus codon usage (PMC2528880) — definition restatement | 1 | https://pmc.ncbi.nlm.nih.gov/articles/PMC2528880/ | 2026-06-13 |

### 1.2 Key Evidence Points

1. RSCU_{i,j} = (n_i × x_{i,j}) / Σ_{k=1}^{n_i} x_{i,k}, where n_i = number of synonymous codons (degeneracy) and x = observed count — source 3 (verbatim), equivalent to source 2's X_i/((1/N_i)ΣX_j).
2. No bias ⇒ RSCU = 1.00; >1 overused, <1 underused — sources 2, 4.
3. RSCU ∈ [0, n_i]; max n_i when only one synonymous codon used — source 3.
4. Definition: observed frequency / expected frequency under equal synonymous usage — sources 4, 6; introduced by source 1.
5. Single-codon families (Met=ATG, Trp=TGG): n_i=1 ⇒ RSCU=1 when present — source 3 (bounds), source 5.
6. Absent family (0/0) is implementation-defined; cubar uses a pseudocount default 1; repository returns 0 — source 5.

### 1.3 Documented Corner Cases

- Absent synonymous family → denominator 0 (0/0), undefined canonically; repository returns 0 (source 5 context).
- Single-codon families always give RSCU = 1 (source 3 bounds).
- CountCodons: non-overlapping triplets from offset 0; trailing 1–2 bases ignored; non-ACGT triplets excluded (repository contract; Kazusa codon-counting convention).

### 1.4 Known Failure Modes / Pitfalls

1. Using `(n−1)`-style or per-thousand normalization instead of the synonymous-family ratio would change values — formula must be the family ratio (sources 2, 3).
2. Treating each codon as its own family (ignoring synonymy) yields RSCU = number-of-codons, not the bias measure — sources 2, 3.

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `CalculateRscu(DnaSequence)` | CodonUsageAnalyzer | Canonical | core RSCU; deep evidence-based tests |
| `CalculateRscu(string)` | CodonUsageAnalyzer | Delegate | string overload (uppercases, empty→empty); smoke |
| `CountCodons(DnaSequence)` | CodonUsageAnalyzer | Canonical | non-overlapping triplet counting |
| `CountCodons(string)` | CodonUsageAnalyzer | Delegate | string overload (uppercases, excludes non-ACGT); smoke |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | For a present family, RSCU = (n_i·x_{i,j})/Σx | Yes | Source 3 (verbatim formula) |
| INV-2 | No bias (equal usage within a family) ⇒ RSCU = 1 | Yes | Sources 2, 4 |
| INV-3 | 0 ≤ RSCU ≤ n_i (family degeneracy) | Yes | Source 3 (bounds) |
| INV-4 | RSCU values within a present family sum to n_i | Yes | Derivation from source 3 formula |
| INV-5 | Single-codon family (Met, Trp) ⇒ RSCU = 1 when present | Yes | Source 3 (bounds), source 5 |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | RSCU Leu 6-fold | `CTGCTGCTGCTA` (CTG×3, CTA×1) | RSCU(CTG)=4.5, RSCU(CTA)=1.5, RSCU(TTA/TTG/CTT/CTC)=0.0 | Source 3 formula; Evidence dataset |
| M2 | RSCU Phe 2-fold biased | `TTTTTTTTC` (TTT×2, TTC×1) | RSCU(TTT)=4/3=1.3333333333333333, RSCU(TTC)=2/3=0.6666666666666666 | Source 3 formula |
| M3 | RSCU unbiased | `TTTTTC` (equal Phe) | RSCU(TTT)=1.0, RSCU(TTC)=1.0 | Sources 2,4 (no-bias value 1) |
| M4 | RSCU single-codon | `ATGATG` (Met) | RSCU(ATG)=1.0 | Source 3 bounds; source 5 |
| M5 | INV-4 family sum | `CTGCTGCTGCTA` | Σ RSCU over 6 Leu codons = 6.0 | Derivation from source 3 |
| M6 | CountCodons basic | `ATGAAATGA` | ATG=1, AAA=1, TGA=1; total 3 | Repository contract / Kazusa convention |
| M7 | CountCodons repeated | `ATGATGATG` | ATG=3 | Triplet counting |
| M8 | CountCodons trailing ignored | `ATGAA` | ATG=1; count=1 (`AA` ignored) | Non-overlapping triplets |
| M9 | CountCodons non-ACGT excluded | `ATGNNNAAA` (string) | ATG=1, AAA=1; `NNN` not counted | IsValidCodon contract |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | RSCU null DnaSequence | `CalculateRscu((DnaSequence)null!)` | throws ArgumentNullException | input guard |
| S2 | RSCU empty string | `CalculateRscu("")` | empty dictionary | input guard |
| S3 | CountCodons null DnaSequence | `CountCodons((DnaSequence)null!)` | throws ArgumentNullException | input guard |
| S4 | CountCodons empty string | `CountCodons("")` | empty dictionary | input guard |
| S5 | RSCU string overload delegation | `CalculateRscu("CTGCTGCTGCTA")` | matches M1 values | Delegate smoke |
| S6 | CountCodons string overload delegation | `CountCodons("atgaaatga")` (lowercase) | ATG=1, AAA=1, TGA=1 | Delegate smoke + case-insensitivity |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | INV-3 bounds check | `CTGCTGCTGCTA` | every RSCU in [0, 6] | bounds invariant |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- `tests/Seqeron/Seqeron.Genomics.Tests/CodonUsageAnalyzerTests.cs` — broad legacy fixture covering CountCodons, RSCU, CAI, ENC, Statistics, reference tables. RSCU and CountCodons cases here overlap this unit; CAI/ENC/Statistics belong to CODON-STATS-001 / CODON-ENC-001 and stay.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1 RSCU Leu 6-fold | ❌ Missing | no 6-fold exact-value test existed |
| M2 RSCU Phe 2-fold biased | ⚠ Weak | legacy `CalculateRscu_BiasedUsage` used `TTTTTT` only, `.Within(0.01)`; no 4/3,2/3 case |
| M3 RSCU unbiased | ⚠ Weak | legacy `CalculateRscu_UnbiasedUsage` uses `.Within(0.01)`, no message |
| M4 RSCU single-codon | ❌ Missing | not covered |
| M5 INV-4 family sum | ❌ Missing | not covered |
| M6 CountCodons basic | ⚠ Weak | legacy `CountCodons_SimpleCodingSequence` lacks messages |
| M7 CountCodons repeated | ⚠ Weak | legacy present, no message |
| M8 CountCodons trailing | ⚠ Weak | legacy `CountCodons_IncompleteCodon` present, no message |
| M9 CountCodons non-ACGT | ❌ Missing | not covered (no N exclusion test) |
| S1 RSCU null | ✅ Covered | legacy `CalculateRscu_NullSequence` (moved to canonical file) |
| S2 RSCU empty | ⚠ Weak | legacy present, no message |
| S3 CountCodons null | ✅ Covered | legacy present (moved) |
| S4 CountCodons empty | ⚠ Weak | legacy present, no message |
| S5 RSCU string delegation | ⚠ Weak | legacy `CalculateRscu_StringOverload` only checks key exists |
| S6 CountCodons string delegation | ⚠ Weak | legacy `CountCodons_StringOverload` no lowercase |
| C1 bounds | ❌ Missing | not covered |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/CodonUsageAnalyzer_CalculateRscu_Tests.cs` — all RSCU + CountCodons cases for this unit, evidence-based, exact values.
- **Remove:** the `#region Codon Counting Tests` and `#region RSCU Tests` blocks and the RSCU/CountCodons null-tests from the `#region Edge Cases` block in `CodonUsageAnalyzerTests.cs` (migrated here). CAI/ENC/Statistics/Reference-table tests remain in the legacy file (other units).

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| `CodonUsageAnalyzer_CalculateRscu_Tests.cs` | Canonical for CODON-RSCU-001 | 16 |
| `CodonUsageAnalyzerTests.cs` | Legacy (CAI/ENC/Stats/Ref only) | reduced |

### 5.5 Phase 7 Work Queue

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| 1 | M1 | ❌ Missing | implemented exact 6-fold test | ✅ Done |
| 2 | M2 | ⚠ Weak | rewritten with 4/3, 2/3 exact values | ✅ Done |
| 3 | M3 | ⚠ Weak | rewritten with messages, Within(1e-10) | ✅ Done |
| 4 | M4 | ❌ Missing | implemented single-codon test | ✅ Done |
| 5 | M5 | ❌ Missing | implemented family-sum invariant test | ✅ Done |
| 6 | M6 | ⚠ Weak | rewritten with Assert.Multiple + messages | ✅ Done |
| 7 | M7 | ⚠ Weak | rewritten with message | ✅ Done |
| 8 | M8 | ⚠ Weak | rewritten with count + message | ✅ Done |
| 9 | M9 | ❌ Missing | implemented non-ACGT exclusion test | ✅ Done |
| 10 | S1 | ✅ Covered | migrated | ✅ Done |
| 11 | S2 | ⚠ Weak | rewritten with message | ✅ Done |
| 12 | S3 | ✅ Covered | migrated | ✅ Done |
| 13 | S4 | ⚠ Weak | rewritten with message | ✅ Done |
| 14 | S5 | ⚠ Weak | rewritten to assert exact delegated values | ✅ Done |
| 15 | S6 | ⚠ Weak | rewritten with lowercase input | ✅ Done |
| 16 | C1 | ❌ Missing | implemented bounds test | ✅ Done |

**Total items:** 16
**✅ Done:** 16 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ | exact 6-fold values |
| M2 | ✅ | 4/3, 2/3 exact |
| M3 | ✅ | unbiased = 1.0 |
| M4 | ✅ | single-codon = 1.0 |
| M5 | ✅ | family sum = 6 |
| M6 | ✅ | basic counting |
| M7 | ✅ | repeated |
| M8 | ✅ | trailing ignored |
| M9 | ✅ | non-ACGT excluded |
| S1 | ✅ | null throws |
| S2 | ✅ | empty → empty |
| S3 | ✅ | null throws |
| S4 | ✅ | empty → empty |
| S5 | ✅ | string delegation exact |
| S6 | ✅ | lowercase delegation |
| C1 | ✅ | bounds [0,6] |

Total in-scope cases: 16. ✅ count: 16.

---

## 6. Assumption Register

**Total assumptions:** 2

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | Absent synonymous family (0/0) returns 0 (no pseudocount) | not exercised by MUST tests (only present families tested); documented |
| 2 | Stop codons grouped as a 3-fold family | does not affect any amino-acid RSCU; documented |

---

## 7. Open Questions / Decisions

1. Pseudocount smoothing (cubar default 1) is intentionally not applied; only affects absent families, which the repository returns as 0. Decision: keep the classic Sharp et al. ratio for present families; document the absent-family convention. No open correctness question for present families.
