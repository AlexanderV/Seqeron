# Test Specification: CODON-ENC-001

**Test Unit ID:** CODON-ENC-001
**Area:** Codon
**Algorithm:** Effective Number of Codons (ENC / Nc), Wright 1990
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-09-28

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Wright F. (1990). The 'effective number of codons' used in a gene. Gene 87(1):23–29. | 1 | https://doi.org/10.1016/0378-1119(90)90491-9 | 2026-06-13 (via refs 2,3) |
| 2 | Fuglsang A. (2004). The 'effective number of codons' revisited. BBRC 317:957–964. | 1 | https://doi.org/10.1016/j.bbrc.2004.03.138 | 2026-06-13 |
| 3 | Fuglsang A. (2006). Estimating the 'effective number of codons'… Genetics 172(2):1301–1307. | 1 | https://academic.oup.com/genetics/article/172/2/1301/5923091 | 2026-06-13 |
| 4 | Peden J.F. (1999) CodonW 1.4.4 — `codon_us.c` `enc_out`, `README_indices.txt` (source compiled and run) | 1 | https://codonw.sourceforge.net/ | 2026-09-28 |
| 5 | codonbias 0.5.0 (PyPI) `EffectiveNumberOfCodons` (Wright mode: robust=False, pseudocount=0, unweighted) | 2 | https://pypi.org/project/codon-bias/ | 2026-09-28 |

### 1.2 Key Evidence Points

1. Codon homozygosity `F̂ = (n·Σp_i² − 1)/(n − 1)`, `p_i = n_i/n`, k synonymous codons — Fuglsang 2004 Eq. 1.
2. Per-amino-acid effective codons `N̂c(aa) = 1/F̂` — Fuglsang 2004 Eq. 2.
3. `N̂c = 2 + 9/F̂₂ + 1/F̂₃ + 5/F̂₄ + 3/F̂₆`; constant 2 = Met + Trp (singlets) — Fuglsang 2004 Eq. 3.
4. Missing amino acid in a class → use the average F of estimable members of that class — Fuglsang 2004 Eq. 4.
5. If N̂c > 61, re-adjust down to 61 — Fuglsang 2004.
6. Isoleucine fallback `F̂₃ = (F̂₂ + F̂₄)/2` when isoleucine unestimable — Fuglsang 2004 Eq. 5a.
7. Range 20 (extreme bias) ≤ Nc ≤ 61 (no bias) — Fuglsang 2004; NCBI standard code degeneracy partition (9/1/5/3 + 2 singlets).
8. Empty synonymous class (other than a lone Ile class) ⇒ "Nc is not calculated" — CodonW README (citing Wright 1990), `enc_out` prints `*****`. Library returns 0.
9. F̂ = 0 (every observed codon once) is not an estimate: `enc_out` keeps an amino acid only `if (bb > 0.0000001)`.
10. Classes derive from the genetic code (`-enc -code`); for codes with ≠ 61 sense codons the cap is the sense-codon count (codonbias; CodonW hard-codes 61).

### 1.3 Documented Corner Cases

- Amino acid with n ≤ 1: F undefined (Fuglsang 2004 — "at least two codons for each amino acid").
- Empty degeneracy class: class average undefined; isoleucine has explicit fallback Eq. 5a; any other empty class ⇒ Nc not calculated (return 0) — CodonW.
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
| `CalculateEnc(string, GeneticCode)` | CodonUsageAnalyzer | Canonical (code-aware) | Classes from `GeneticCode.CodonTable`; the string overload uses `GeneticCode.Standard`. |
| `CalculateEnc(DnaSequence, GeneticCode)` | CodonUsageAnalyzer | Delegate | Same core. |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | Nc = 0 (not calculable) or 20 ≤ Nc ≤ 61 (table 1; ≤ sense-codon count for other codes) | Yes | Fuglsang 2004 range + Eq. 3 cap; CodonW enc_out |
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
| M3 | FullyPopulatedBiasedGene | All classes estimable (no fallback) | Nc = 41.288461538461526 (independent reference) | Fuglsang 2004 Eq. 1, Eq. 3, Eq. 4 |
| M4 | Invariant_Range | Deterministic genes with every class estimable | 20 ≤ Nc ≤ 61 | Fuglsang 2004 range (INV-1) |
| M5 | IsoleucineAbsent_UsesFallback | All classes present except Ile; Eq. 5a genuinely fires | Nc = 39.47394540942927 (independent reference) | Fuglsang 2004 Eq. 5a; Peden codonW thesis |
| M5b | EmptySynonymousClass_NotCalculated | Phe-only; Ile+4-fold absent; Lys-only; short gene; Met/Trp only; no codon | 0 (CodonW `*****`) | CodonW README + enc_out, binary run (was 29.0 library convention — corrected 2026-09, F15) |
| M8 | ZeroHomozygosityAminoAcid_Excluded | M3 + His CAT,CAC (F̂ = 0) | 41.288461538461526 (CodonW 41.29) | enc_out `bb > 0.0000001` (F16) |
| M9 | AllSenseCodonsOnce_NotCalculated | every sense codon once | 0 (CodonW `*****`; was 20) | enc_out (F16) |
| M10 | DeterministicGene_MatchesCodonW | codon_i × ((3i mod 7)+1) | 57.5614857446809 (CodonW 57.56) | CodonW binary, enc_out port, codonbias |
| M11 | AlternativeGeneticCode_MatchesCodonW | same gene, tables 1/2/3/9; M3 under 2/3/9 | 57.5615/55.3819/57.4784/58.1526; 45.0/43.1853/43.8333 | CodonW `-code 0/1/2/7` (F17) |
| M12 | Overshoot_CappedAtSenseCodonCount | codon_i × ((7i mod 5)+1) | 61/60/62/62 (tables 1/2/3/9) | codonbias `min(len(P),ENC)`; CodonW 61 for table 1 |
| M13 | GeneticCode overload contracts | null code / null DnaSequence; DnaSequence+code equals string+code | ArgumentNullException; equal | Contract |
| M6 | Null_Throws | `CalculateEnc((DnaSequence)null!)` | ArgumentNullException | Contract |
| M7 | Empty_ReturnsZero | empty string | 0 | Contract |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | Lowercase_Normalized | lowercase input equals uppercase result | equal | Case normalization |
| S4 | RnaAndLowerCase_EqualDna | U-spelled and lower-case gene | equal to DNA | CodonW ident_codon reads U as T |
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

1. The Eq. 5a isoleucine fallback and Eq. 4 within-class averaging are implemented per Fuglsang 2004; the lower clamp at 20 is documented as a defensive bound consistent with the published range.
2. **Resolved 2026-09 (review B02, F15–F17).** Whole-class-absent handling now follows CodonW `enc_out` / Wright 1990 ("Nc is not calculated" ⇒ return 0) instead of the former full-count convention (M5b 29.0 and several fuzz/property/combinatorial tests asserting 20 or [20,61] for such genes were corrected); F̂ = 0 amino acids are excluded from class averages (CodonW); classes derive from `GeneticCode` with new `CalculateEnc(..., GeneticCode)` overloads. Cross-checked against the CodonW 1.4.4 binary (6456 cases), a Python port of `enc_out` (27 tables, 21789 cases) and codonbias 0.5.0 (27 tables, 4050 genes).
