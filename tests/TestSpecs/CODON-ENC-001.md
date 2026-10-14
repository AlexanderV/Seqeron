# Test Specification: CODON-ENC-001

**Test Unit ID:** CODON-ENC-001
**Area:** Codon
**Algorithm:** Effective Number of Codons (ENC / Nc), Wright 1990
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-10-14

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Wright F. (1990). The 'effective number of codons' used in a gene. Gene 87(1):23–29. | 1 | https://doi.org/10.1016/0378-1119(90)90491-9 | 2026-06-13 (via refs 2,3) |
| 2 | Fuglsang A. (2004). The 'effective number of codons' revisited. BBRC 317:957–964. | 1 | https://doi.org/10.1016/j.bbrc.2004.03.138 | 2026-06-13 |
| 3 | Fuglsang A. (2006). Estimating the 'effective number of codons'… Genetics 172(2):1301–1307. | 1 | https://academic.oup.com/genetics/article/172/2/1301/5923091 | 2026-06-13 |
| 4 | Peden J.F. (1999) CodonW 1.4.4 — `codon_us.c` `enc_out`, `README_indices.txt` (source compiled and run) | 1 | https://codonw.sourceforge.net/ | 2026-09-28 |
| 5 | codonbias 0.5.0 (PyPI) `EffectiveNumberOfCodons` (Wright mode: robust=False, pseudocount=0, unweighted; reference for the Novembre / Sun–Yang–Xia estimators) | 2 | https://pypi.org/project/codon-bias/ | 2026-10-14 |
| 6 | Novembre J.A. (2002) MBE 19:1390 (Nc′); ENCprime `src/ENCprime.c` (consulted) | 1 | https://github.com/jnovembre/ENCprime | 2026-10-14 (program; paper via refs 5/6 — journal blocked) |
| 7 | Sun X., Yang Q., Xia X. (2013) MBE 30:191 (improved Nc) | 1 | https://doi.org/10.1093/molbev/mss201 | via ref 5 (journal blocked) |

### 1.2 Key Evidence Points

1. Codon homozygosity `F̂ = (n·Σp_i² − 1)/(n − 1)`, `p_i = n_i/n`, k synonymous codons — Fuglsang 2004 Eq. 1.
2. Per-amino-acid effective codons `N̂c(aa) = 1/F̂` — Fuglsang 2004 Eq. 2.
3. `N̂c = 2 + 9/F̂₂ + 1/F̂₃ + 5/F̂₄ + 3/F̂₆`; constant 2 = Met + Trp (singlets) — Fuglsang 2004 Eq. 3.
4. Missing amino acid in a class → use the average F of estimable members of that class — Fuglsang 2004 Eq. 4.
5. If N̂c > 61, re-adjust down to 61 — Fuglsang 2004.
6. Isoleucine fallback `F̂₃ = (F̂₂ + F̂₄)/2` when isoleucine unestimable — Fuglsang 2004 Eq. 5a.
7. Range 20 (extreme bias) ≤ Nc ≤ 61 (no bias) — Fuglsang 2004; NCBI standard code degeneracy partition (9/1/5/3 + 2 singlets).
8. Empty synonymous class (other than a lone Ile class) ⇒ "Nc is not calculated" — CodonW README (citing Wright 1990), `enc_out` prints `*****`. `TryCalculateEnc` returns `false` (explicit not-calculable result; review 2026-09 F36); `CalculateEnc` keeps returning 0 for backward compatibility.
9. F̂ = 0 (every observed codon once) is not an estimate: `enc_out` keeps an amino acid only `if (bb > 0.0000001)`.
10. Classes derive from the genetic code (`-enc -code`); for codes with ≠ 61 sense codons the cap is the sense-codon count (codonbias; CodonW hard-codes 61 — opt-in `EncOptions.CodonWCap`, F51).
11. Opt-in estimators (F50): `EncMethod.Novembre2002` — χ² against background-composition codon frequencies, F = (N·χ² + N − k)/(k(N − 1)), unweighted, no pseudocount; `EncMethod.SunYangXia2013` — pseudocount 1, F = (χ² + 1)/k, N-weighted; empty class ⇒ F̄ = 1/z; no lower clamp — codonbias 0.5.0 `scores.py` `_calc_F`/`_calc_score`.

### 1.3 Documented Corner Cases

- Amino acid with n ≤ 1: F undefined (Fuglsang 2004 — "at least two codons for each amino acid").
- Empty degeneracy class: class average undefined; isoleucine has explicit fallback Eq. 5a; any other empty class ⇒ Nc not calculated (`TryCalculateEnc` → false; `CalculateEnc` → 0) — CodonW.
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
| `TryCalculateEnc(string|DnaSequence[, GeneticCode], out double)` | CodonUsageAnalyzer | Canonical (explicit not-calculable) | `false` + 0 exactly where CodonW prints `*****` (and for input without valid codons); the `CalculateEnc` overloads delegate to it (F36). |
| `CodonUsageStatistics.IsEncCalculable` | CodonUsageStatistics | Derived | `Enc > 0` (a calculated Nc is ≥ the number of amino acids). |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | `TryCalculateEnc` false ⇔ Nc not calculable (then `CalculateEnc` = 0); otherwise 20 ≤ Nc ≤ 61 (table 1; ≤ sense-codon count for other codes) | Yes | Fuglsang 2004 range + Eq. 3 cap; CodonW enc_out |
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
| M30 | TryCalculateEnc_CodonWNotCalculated | Phe-only, short gene (± terminal stop), Met/Trp only, random gene r155 (tables 1, 2) | `false`, enc 0; `CalculateEnc` 0; `IsEncCalculable` false | CodonW 1.4.4 binary `*****` (F36) |
| M31 | TryCalculateEnc_Calculable | M3 gene tables 1/2; random gene r53 (N triplet, in-frame stops) tables 1/2 | `true`; 41.288461538461526 / 45.0; 42.85340746733685 / 40.37735849056604 | CodonW `-enc` 41.29 / 45.00 / 42.85 / 40.38 |
| M32 | TryCalculateEnc_NoCodons | null, "", "A", "NNNRYK" | `false`, 0 | Contract |
| M33 | TryCalculateEnc overload contracts | default-code overloads = table 1; DnaSequence overloads; null code / null DnaSequence | equal; ArgumentNullException | Contract |
| M34 | MCP `effective_number_of_codons.calculable`, `codon_usage_statistics.isEncCalculable` | `*****` gene / M3 gene | false / true | CodonW (F36) |
| M40 | Novembre2002_MatchesCodonbias (+ explicit background, pseudocount 0.5 weighted) | M3, r53, r203, M10, Phe, r15; tables 1/2 | e.g. M3 45.49951134064658 / 46.946215313358586; background M3 38.480446976398305 | codonbias 0.5.0 (F50) |
| M41 | Novembre2002 without background = Wright | M3 | 41.288461538461526 | algebra (χ² = kΣp² − 1) |
| M42 | SunYangXia2013_MatchesCodonbias (+ bg-corrected, pseudocount 0 unweighted) | M3, M10, M12, r53, r203, r15; tables 1/2 | e.g. M3 50.957809330628805 / 50.70834399908201 | codonbias 0.5.0 (F50) |
| M43 | SunYangXia2013 Phe-only hand derivation | TTT×3 TTC | 1028/17 | hand derivation = codonbias 60.470588235294116 |
| M44 | Default `EncOptions` bit-identical to Wright; no sense codon ⇒ false; RNA/lower-case/DnaSequence overloads; invalid option sets / nulls throw | — | — | Contract |
| M45 | CodonWCap on M12 under NCBI 1/2/3/4/5/6/9/10 | default 61/60/62/62/62/63/62/62 | 61 for all | CodonW 1.4.4 `-enc -code 0..7` 61.00 (F51) |
| M46 | CodonWCap r203 table 2 / r15 table 3 | 60 → 60.71285921278657; 61.89427123695478 → 61 | CodonW 60.71 / 61.00 | CodonW 1.4.4 (F51) |
| M47 | MCP `method` / `pseudocount` / `background_correction` / `background` / `codonw_cap` forwarding, defaults, rejections | M3, M12 | as M40–M46 | `EncMethodCodonWCapToolTests.cs` |

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
| M30–M34 | ✅ Covered | `CodonUsageAnalyzer_NotCalculableEnc_CodonWGc_Tests.cs`, MCP `CodonEncCalculableGcConventionTests.cs` (F36) |
| M40–M47 | ✅ Covered | `CodonUsageAnalyzer_EncVariants_Tests.cs`, MCP `EncMethodCodonWCapToolTests.cs` (F50/F51) |

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
3. **Resolved 2026-09 (review B02, F36, audit A1-11).** The 0 sentinel is no longer the only way to see "not calculated": additive `TryCalculateEnc(…, out double)` overloads (repo `Try*` convention, e.g. `OncologyAnalyzer.TryFitPurityPloidy`) return `false` where CodonW prints `*****`; `CalculateEnc` delegates and still returns 0; `CodonUsageStatistics.IsEncCalculable`; MCP `effective_number_of_codons` adds `calculable`, `codon_usage_statistics` adds `isEncCalculable`. Re-verified on the CodonW 1.4.4 binary (codes 0–7 × 304 genes with in-frame stops, N/R/Y, lower-case, RNA, trailing bases): all 273 `*****` cases → `false`; all 2159 calculated cases `true` and equal at 2 dp, except the documented sense-codon cap (CodonW 61 vs 60/62/63).
4. **Resolved 2026-10-14 (review B02, F50/F51, audit A2-7/A2-8).** Novembre (2002) Nc′ and Sun, Yang & Xia (2013) Nc are opt-in via `EncOptions.Method` (additive `CalculateEnc`/`TryCalculateEnc(…, GeneticCode, EncOptions[, out double])`), following codonbias 0.5.0 exactly (57351 gene × table × option cases, ≤ 2.2e-15 relative); `EncOptions.CodonWCap` reproduces CodonW's fixed 61 cap (2448 CodonW 1.4.4 cases). Wright default unchanged (bit-identical).
