# Test Specification: CODON-ENC-001

**Test Unit ID:** CODON-ENC-001
**Area:** Codon
**Algorithm:** Effective Number of Codons (ENC / Nc), Wright 1990
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-13

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Wright F. (1990). The 'effective number of codons' used in a gene. Gene 87(1):23–29. | 1 | https://doi.org/10.1016/0378-1119(90)90491-9 | 2026-06-13 (via refs 2,3) |
| 2 | Fuglsang A. (2004). The 'effective number of codons' revisited. BBRC 317:957–964. | 1 | https://doi.org/10.1016/j.bbrc.2004.03.138 | 2026-06-13 |
| 3 | Fuglsang A. (2006). Estimating the 'effective number of codons'… Genetics 172(2):1301–1307. | 1 | https://academic.oup.com/genetics/article/172/2/1301/5923091 | 2026-06-13 |

### 1.2 Key Evidence Points

1. Codon homozygosity `F̂ = (n·Σp_i² − 1)/(n − 1)`, `p_i = n_i/n`, k synonymous codons — Fuglsang 2004 Eq. 1.
2. Per-amino-acid effective codons `N̂c(aa) = 1/F̂` — Fuglsang 2004 Eq. 2.
3. `N̂c = 2 + 9/F̂₂ + 1/F̂₃ + 5/F̂₄ + 3/F̂₆`; constant 2 = Met + Trp (singlets) — Fuglsang 2004 Eq. 3.
4. Missing amino acid in a class → use the average F of estimable members of that class — Fuglsang 2004 Eq. 4.
5. If N̂c > 61, re-adjust down to 61 — Fuglsang 2004.
6. Isoleucine fallback `F̂₃ = (F̂₂ + F̂₄)/2` when isoleucine unestimable — Fuglsang 2004 Eq. 5a.
7. Range 20 (extreme bias) ≤ Nc ≤ 61 (no bias) — Fuglsang 2004; NCBI standard code degeneracy partition (9/1/5/3 + 2 singlets).

### 1.3 Documented Corner Cases

- Amino acid with n ≤ 1: F undefined (Fuglsang 2004 — "at least two codons for each amino acid").
- Empty degeneracy class: class average undefined; isoleucine has explicit fallback Eq. 5a.
- Overshoot N̂c > 61 → cap at 61.

### 1.4 Known Failure Modes / Pitfalls

1. Using raw counts n_i² instead of frequencies p_i² inside Eq. 1 over-scales F by a factor of n — defect against Fuglsang 2004 Eq. 1.
2. Treating an absent class as contributing its full codon count (e.g. `enc += 9`) instead of the Eq. 4 within-class average / Eq. 5a fallback — defect against Fuglsang 2004 Eqs. 3–5a.

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `CalculateEnc(string)` | CodonUsageAnalyzer | Canonical | Core Wright 1990 computation on raw string. |
| `CalculateEnc(DnaSequence)` | CodonUsageAnalyzer | Delegate | Thin wrapper → `.Sequence`; smoke + null check only. |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | 20 ≤ Nc ≤ 61 for any non-empty coding sequence | Yes | Fuglsang 2004 range + Eq. 3 cap |
| INV-2 | Maximally biased gene (one codon per amino acid) → Nc = 20 | Yes | Fuglsang 2004 range; Eqs. 1–3 (F=1 ⇒ Nc(aa)=1) |
| INV-3 | Near-uniform gene → Nc re-adjusted to exactly 61 | Yes | Fuglsang 2004 cap rule |
| INV-4 | Single two-fold amino acid with counts (3,1) → F=0.5, Nc(aa)=2 | Yes | Fuglsang 2004 Eq. 1, Eq. 2 (hand derivation) |
| INV-5 | Empty / null input handled (0 / ArgumentNullException) | Yes | Contract |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | MaxBias_OneCodonPerAa | Each amino acid uses exactly one codon, ≥2 times | Nc = 20 (Within 1e-9) | Fuglsang 2004 range; F=1 ⇒ Nc(aa)=1, sum=20 |
| M2 | NearUniform_CapsAt61 | All codons equal counts (c=2 per codon) → raw Nc > 61 | Nc = 61 exactly | Fuglsang 2004 cap rule |
| M3 | TwoFold_ExactF | Only Phe: TTT×3, TTC×1 | Nc(Phe)=2 via F=0.5; full Nc derived in test | Fuglsang 2004 Eq. 1, Eq. 2, Eq. 4 |
| M4 | Invariant_Range | Several arbitrary deterministic sequences | 20 ≤ Nc ≤ 61 | Fuglsang 2004 range (INV-1, property test) |
| M5 | IsoleucineAbsent_UsesFallback | Gene with no Ile but Phe(2-fold) and Ala(4-fold) present | F₃ = (F₂+F₄)/2 used; matches hand derivation | Fuglsang 2004 Eq. 5a |
| M6 | Null_Throws | `CalculateEnc((DnaSequence)null!)` | ArgumentNullException | Contract |
| M7 | Empty_ReturnsZero | empty string | 0 | Contract |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | Lowercase_Normalized | lowercase input equals uppercase result | equal | Case normalization |
| S2 | InvalidCodons_Skipped | sequence with an N-containing codon | identical to sequence without it | Non-ACGT codons ignored |
| S3 | DnaSequenceOverload_Delegates | DnaSequence overload equals string overload | equal | Delegate smoke |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | Fuglsang40_5 | Construct gene matching the Nc=40.5 simulation profile | Nc close to 40.5 | Asymptotic; documented as illustrative, not exact-equality |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- Searched `tests/Seqeron/Seqeron.Genomics.Tests/` for `CalculateEnc` / `CodonUsageAnalyzer`. No prior canonical `CodonUsageAnalyzer_CalculateEnc_Tests.cs` existed for CODON-ENC-001.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1 | ❌ Missing | new unit |
| M2 | ❌ Missing | new unit |
| M3 | ❌ Missing | new unit |
| M4 | ❌ Missing | new unit |
| M5 | ❌ Missing | new unit |
| M6 | ❌ Missing | new unit |
| M7 | ❌ Missing | new unit |
| S1 | ❌ Missing | new unit |
| S2 | ❌ Missing | new unit |
| S3 | ❌ Missing | new unit |
| C1 | ❌ Missing | new unit |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/CodonUsageAnalyzer_CalculateEnc_Tests.cs` — all CODON-ENC-001 cases.
- **Remove:** none (no prior dedicated tests).

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| CodonUsageAnalyzer_CalculateEnc_Tests.cs | Canonical CODON-ENC-001 | 11 |

### 5.5 Phase 7 Work Queue

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| 1 | M1 | ❌ Missing | Implemented | ✅ Done |
| 2 | M2 | ❌ Missing | Implemented | ✅ Done |
| 3 | M3 | ❌ Missing | Implemented | ✅ Done |
| 4 | M4 | ❌ Missing | Implemented (property) | ✅ Done |
| 5 | M5 | ❌ Missing | Implemented | ✅ Done |
| 6 | M6 | ❌ Missing | Implemented | ✅ Done |
| 7 | M7 | ❌ Missing | Implemented | ✅ Done |
| 8 | S1 | ❌ Missing | Implemented | ✅ Done |
| 9 | S2 | ❌ Missing | Implemented | ✅ Done |
| 10 | S3 | ❌ Missing | Implemented | ✅ Done |
| 11 | C1 | ❌ Missing | Implemented | ✅ Done |

**Total items:** 11
**✅ Done:** 11 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ Covered | Nc = 20 exact |
| M2 | ✅ Covered | cap to 61 |
| M3 | ✅ Covered | exact F derivation |
| M4 | ✅ Covered | property test, multiple inputs |
| M5 | ✅ Covered | isoleucine fallback |
| M6 | ✅ Covered | ArgumentNullException |
| M7 | ✅ Covered | empty → 0 |
| S1 | ✅ Covered | case normalization |
| S2 | ✅ Covered | invalid codons skipped |
| S3 | ✅ Covered | delegate equality |
| C1 | ✅ Covered | 40.5 profile within tolerance |

---

## 6. Assumption Register

**Total assumptions:** 1

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | Lower clamp at 20 is a defensive structural bound (Wright only prescribes the upper re-adjustment to 61). | INV-1 |

---

## 7. Open Questions / Decisions

1. None. The Eq. 5a isoleucine fallback and Eq. 4 within-class averaging are implemented per Fuglsang 2004; the lower clamp at 20 is documented as a defensive bound consistent with the published range.
