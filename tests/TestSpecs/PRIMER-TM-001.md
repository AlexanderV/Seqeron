# Test Specification: PRIMER-TM-001

## Test Unit: Melting Temperature Calculation

**Area:** MolTools
**Status:** Active
**Created:** 2026-01-22
**Last Verified:** 2026-03-04
**Evidence Sources:** Marmur & Doty (1962), Thein & Wallace (1986), Sigma-Aldrich/Merck Technical Docs, Owczarzy et al. (2004)

---

## 1. Overview

Melting temperature (Tm) is the temperature at which 50% of the DNA duplex is dissociated into single strands. Accurate Tm calculation is essential for PCR primer design, probe design, and hybridization experiments.

### Canonical Methods

| Method | Class | Type | Complexity |
|--------|-------|------|------------|
| `CalculateMeltingTemperature(string)` | PrimerDesigner | Canonical | O(n) |
| `CalculateMeltingTemperatureWithSalt(string, double)` | PrimerDesigner | Salt-corrected | O(n) |

### Helper Constants (ThermoConstants)

| Constant/Method | Value/Formula | Description |
|-----------------|---------------|-------------|
| `WallaceMaxLength` | 14 | Threshold: < 14 uses Wallace rule |
| `WallaceAtContribution` | 2 | Tm contribution per A/T base |
| `WallaceGcContribution` | 4 | Tm contribution per G/C base |
| `CalculateWallaceTm(at, gc)` | 2×AT + 4×GC | Wallace rule formula |
| `MarmurDotyBase` | 64.9 | Marmur-Doty base temperature |
| `MarmurDotyGcCoefficient` | 41.0 | GC coefficient |
| `MarmurDotyGcOffset` | 16.4 | GC offset correction |
| `CalculateMarmurDotyTm(gc, len)` | 64.9 + 41×(GC-16.4)/len | Marmur-Doty formula |
| `CalculateSaltCorrection(Na_mM)` | 16.6 × log10(Na/1000) | Salt correction |

---

## 2. Evidence Summary

### 2.1 Wallace Rule (Short Oligonucleotides)

**Source:** Thein & Wallace (1986), widely cited in literature

**Formula:** Tm = 2×(A+T) + 4×(G+C)

**Applicability:** Primers with < 14 valid DNA bases

**Rationale:**
- A-T pairs form 2 hydrogen bonds (weaker)
- G-C pairs form 3 hydrogen bonds (stronger)
- Simple approximation for very short oligos

**External Verification:**
Sigma-Aldrich/Merck uses a modified variant: Tm = 2(A+T) + 4(G+C) − 7 for ≤14 bases.
The −7 correction accounts for "in solution" conditions vs membrane hybridization (Marmur & Doty 1962).
Our implementation uses the original Wallace rule without the −7 correction;
this is a deliberate design choice consistent with published textbook formulations.

**Test Values (Evidence-based):**
| Sequence | Length | A+T | G+C | Expected Tm |
|----------|--------|-----|-----|-------------|
| ATATATAT | 8 | 8 | 0 | 16°C |
| GCGCGCGC | 8 | 0 | 8 | 32°C |
| ACGT | 4 | 2 | 2 | 12°C |
| AAAA | 4 | 4 | 0 | 8°C |
| GGGG | 4 | 0 | 4 | 16°C |
| ACGTACGT | 8 | 4 | 4 | 24°C |

### 2.2 Marmur-Doty Formula (Longer Primers)

**Source:** Marmur & Doty (1962), "Determination of the base composition of deoxyribonucleic acid from its thermal denaturation temperature", J Mol Biol 5:109-118

**Formula:** Tm = 64.9 + 41 × (GC - 16.4) / N

Where:
- GC = number of G and C bases
- N = total valid base count (ACGT only)

**Applicability:** Primers ≥ 14 valid DNA bases

**External Verification:**
Sigma-Aldrich/Merck uses nearest-neighbor (SantaLucia 1998) for ≥15 bases as their primary method.
The Marmur-Doty variant used here is a simplified empirical formula suitable for
basic primer analysis. This is consistent with widely-published bioinformatics references.

**Test Values (Calculated):**
| Sequence (20bp) | GC count | GC% | Expected Tm |
|-----------------|----------|-----|-------------|
| All A/T | 0 | 0% | 64.9 + 41×(-16.4)/20 = 31.28°C |
| 50% GC (10 each) | 10 | 50% | 64.9 + 41×(-6.4)/20 = 51.78°C |
| All G/C | 20 | 100% | 64.9 + 41×(3.6)/20 = 72.28°C |

### 2.3 Salt Correction

**Source:** Owczarzy et al. (2004), general PCR literature

**Formula:** Salt correction = 16.6 × log10([Na⁺]/1000)

Where [Na⁺] is in mM (typical PCR: 50 mM)

**Example:** At 50 mM Na⁺: 16.6 × log10(0.05) = 16.6 × (-1.301) ≈ -21.6°C

**Defined Behavior:** Salt correction returns 0 for empty/null input (no duplex to correct).

---

## 3. Defined Behaviors

### 3.1 Input Alphabet

Only standard DNA bases (A, C, G, T) are recognized. All other characters — including
IUPAC ambiguity codes (N, R, Y, etc.) and RNA bases (U) — are ignored.
Both threshold determination and formula computation use the count of valid ACGT bases only.

| Input | Valid Bases | Behavior |
|-------|------------|----------|
| `"ACNGT"` | A,C,G,T (4 valid) | N ignored; Wallace: 2×2 + 4×2 = 12 |
| `"ACGTNNNNACGT"` | 8 valid | N's ignored; Wallace: 2×4 + 4×4 = 24 |
| `"NNNNN"` | 0 valid | Returns 0 |
| `"ACGUACGU"` | A,C,G,A,C,G (6 valid, U ignored) | Wallace: 2×2 + 4×4 = 20 |

### 3.2 Case Insensitivity

Input is converted to uppercase before processing. `"acgt"` and `"ACGT"` produce identical results.

### 3.3 Empty/Null Handling

Both `CalculateMeltingTemperature` and `CalculateMeltingTemperatureWithSalt`
return 0.0 for empty or null input. Salt correction is not applied to empty/null primers.

---

## 4. Test Cases

### 4.1 Must Tests (Evidence-Based)

#### M1: Empty Input Handling
```
Input: ""
Expected: 0.0
```

#### M2: Null Input Handling
```
Input: null
Expected: 0.0
```

#### M3: Wallace Rule - All A/T (Short)
**Evidence:** Wallace (1986) formula: Tm = 2×AT + 4×GC
```
Input: "ATATATAT" (8 bp, 8 A/T, 0 G/C)
Expected: 2×8 + 4×0 = 16.0
```

#### M4: Wallace Rule - All G/C (Short)
**Evidence:** Wallace (1986) formula
```
Input: "GCGCGCGC" (8 bp, 0 A/T, 8 G/C)
Expected: 2×0 + 4×8 = 32.0
```

#### M5: Wallace Rule - Mixed (Short)
**Evidence:** Wallace (1986) formula
```
Input: "ACGTACGT" (8 bp, 4 A/T, 4 G/C)
Expected: 2×4 + 4×4 = 24.0
```

#### M6: Wallace Rule - Boundary (13 bp, still uses Wallace)
**Evidence:** Implementation uses < 14 valid bases as threshold
```
Input: "ACGTACGTACGTA" (13 bp, 7 A/T, 6 G/C)
Expected: 2×7 + 4×6 = 38.0
```

#### M7: Marmur-Doty - Boundary (14 bp, switches to MD)
**Evidence:** Implementation threshold at 14 valid bases
```
Input: "ACGTACGTACGTAC" (14 bp, 7 A/T, 7 G/C)
Expected: 64.9 + 41×(7-16.4)/14 ≈ 37.36
```

#### M8: Marmur-Doty - Typical Primer (20 bp, 50% GC)
**Evidence:** Marmur & Doty (1962)
```
Input: "ACGTACGTACGTACGTACGT" (20 bp, 10 G/C)
Expected: 64.9 + 41×(10-16.4)/20 ≈ 51.78
```

#### M9: Marmur-Doty - Low GC (20 bp)
**Evidence:** Marmur & Doty (1962)
```
Input: "ATATATATATATATATATAT" (20 bp, 0 G/C)
Expected: 64.9 + 41×(0-16.4)/20 ≈ 31.28
```

#### M10: Marmur-Doty - High GC (20 bp)
**Evidence:** Marmur & Doty (1962)
```
Input: "GCGCGCGCGCGCGCGCGCGC" (20 bp, 20 G/C)
Expected: 64.9 + 41×(20-16.4)/20 ≈ 72.28
```

#### M11: Case Insensitivity
```
Input: "atatatat"
Expected: Same as "ATATATAT" = 16.0 (exact value asserted)
```

#### M12: Salt Correction - Standard 50mM
**Evidence:** Standard PCR salt concentration
```
Input: primer="ACGTACGTACGTACGTACGT", Na=50mM
Base Tm: 51.78
Salt correction: 16.6 × log10(50/1000) ≈ -21.6
Expected: ≈30.2
```

#### M13: Salt Correction - Low Salt (10mM)
**Evidence:** Owczarzy et al. (2004). Lower salt destabilizes duplex.
```
Input: primer="ACGTACGTACGTACGTACGT", Na=10mM
Base Tm: ≈51.78
Salt correction: 16.6 × log10(10/1000) = 16.6 × (-2) = -33.2
Expected: ≈18.58
```

#### M14: Salt Correction - High Salt (200mM)
**Evidence:** Owczarzy et al. (2004). Higher salt stabilizes duplex.
```
Input: primer="ACGTACGTACGTACGTACGT", Na=200mM
Base Tm: ≈51.78
Salt correction: 16.6 × log10(200/1000) ≈ -11.6
Expected: ≈40.18
```

#### M15: Non-ACGT Characters Ignored
**Evidence:** Defined behavior (Section 3.1)
```
Input: "ACNGT" (4 valid bases: A,C,G,T)
Expected: Wallace 2×2 + 4×2 = 12.0
```

#### M16: RNA Base (U) Not Recognized
**Evidence:** Defined behavior — DNA-only tool (Section 3.1)
```
Input: "ACGUACGU" (6 valid bases: A,C,G,A,C,G)
Expected: Wallace 2×2 + 4×4 = 20.0
```

#### M17: All Non-Standard Returns 0
**Evidence:** Defined behavior (Section 3.1)
```
Input: "NNNNN" (0 valid bases)
Expected: 0.0
```

#### M18: Salt Correction - Empty/Null Returns 0
**Evidence:** Defined behavior (Section 3.3)
```
Input: primer="", Na=50mM
Expected: 0.0
```

#### M19: Marmur-Doty - All Same Base (16 bp)
**Evidence:** Edge case — all-same-base above threshold
```
Input: "AAAAAAAAAAAAAAAA" (16 bp, 16 A/T, 0 G/C)
Expected: 64.9 + 41×(0-16.4)/16 = 22.875
```

---

## 5. Invariants

| ID | Invariant | Validation |
|----|-----------|------------|
| INV-1 | Result ≥ 0 for any valid input | Assert.That(tm, Is.GreaterThanOrEqualTo(0)) |
| INV-2 | Higher GC content → Higher Tm | Compare equal-length sequences |
| INV-3 | Salt correction is additive to base Tm | Tm_salt = Tm_base + correction |
| INV-4 | Case insensitivity | toupper(input) == input produces same Tm |

---

## 6. Edge Cases

| Case | Input | Expected Behavior |
|------|-------|-------------------|
| Empty string | "" | Returns 0.0 |
| Null | null | Returns 0.0 |
| Single base | "A" | Wallace: 2.0 |
| Boundary 13bp | 13-char ACGT string | Uses Wallace |
| Boundary 14bp | 14-char ACGT string | Uses Marmur-Doty |
| All same base | "AAAAAAAAAAAAAAAA" (16bp) | Marmur-Doty: 22.875°C (M19) |
| Lowercase | "acgt" | Case-insensitive |
| Non-ACGT (N) | "ACNGT" | Only ACGT counted |
| Only non-ACGT | "NNNNN" | Returns 0.0 |
| RNA (U) | "ACGUACGU" | U ignored; only ACGT counted |

---

## 7. External Source Verification Summary

| Item | Our Implementation | Sigma-Aldrich/Merck | Status |
|------|-------------------|---------------------|--------|
| Short oligo formula | Tm = 2(A+T) + 4(G+C) | Tm = 2(A+T) + 4(G+C) − 7 | **Variant** — −7 correction omitted by design |
| Short oligo threshold | < 14 valid bases | ≤ 14 bases | Aligned (both use 14 as boundary) |
| Long primer formula | Marmur-Doty: 64.9 + 41(GC−16.4)/N | Nearest-neighbor (SantaLucia 1998) | **Simplified** — Marmur-Doty is a simpler, well-published alternative |
| Salt correction | 16.6 × log₁₀(Na_mM/1000) | Integrated into NN formula | Consistent with Owczarzy (2004) |
| Non-ACGT handling | Ignored (only ACGT counted) | Not documented (clean input expected) | Defined behavior |
| RNA (U) | Not supported (ignored) | Not applicable (DNA tool) | Defined behavior |

### Known Variant: −7 Correction Factor

Sigma-Aldrich's "Basic" method includes a −7°C correction factor for
oligonucleotides used in solution (as opposed to membrane hybridization).
Our implementation uses the original Wallace rule without this correction.
This is a **deliberate design choice** — the Wallace rule as published by
Thein & Wallace (1986) does not include the correction.
If in-solution calibration is needed, use `CalculateMeltingTemperatureWithSalt`.

---

## 8. Coverage Classification

All spec test cases verified against `PrimerDesigner_MeltingTemperature_Tests.cs` (34 tests).

| ID | Test Case | Status | Test Method |
|----|-----------|--------|-------------|
| M1 | Empty Input | ✅ Covered | `_EmptyPrimer_Returns0` |
| M2 | Null Input | ✅ Covered | `_NullPrimer_Returns0` |
| M3 | Wallace All A/T | ✅ Covered | `_Wallace_AllAT_Returns16` |
| M4 | Wallace All G/C | ✅ Covered | `_Wallace_AllGC_Returns32` |
| M5 | Wallace Mixed | ✅ Covered | `_Wallace_Mixed_Returns24` |
| M6 | Wallace Boundary 13bp | ✅ Covered | `_Wallace_Boundary13bp_Returns38` |
| M7 | Marmur-Doty Boundary 14bp | ✅ Covered | `_MarmurDoty_Boundary14bp_UsesFormula` |
| M8 | Marmur-Doty 20bp 50%GC | ✅ Covered | `_MarmurDoty_20bp_50GC_ReturnsExpected` |
| M9 | Marmur-Doty 20bp Low GC | ✅ Covered | `_MarmurDoty_20bp_0GC_ReturnsExpected` |
| M10 | Marmur-Doty 20bp High GC | ✅ Covered | `_MarmurDoty_20bp_100GC_ReturnsExpected` |
| M11 | Case Insensitivity | ✅ Covered | `_LowercaseInput_MatchesUppercase` (exact 16.0) |
| M12 | Salt 50mM | ✅ Covered | `_50mM_AppliesCorrection` |
| M13 | Salt 10mM | ✅ Covered | `_10mM_AppliesCorrection` (exact value) |
| M14 | Salt 200mM | ✅ Covered | `_200mM_AppliesCorrection` (exact value) |
| M15 | Non-ACGT Ignored | ✅ Covered | `_NonAcgtIgnored_OnlyValidBasesCounted` |
| M16 | RNA U Ignored | ✅ Covered | `_RnaUracil_NotCountedAsDnaBase` |
| M17 | All Non-Standard → 0 | ✅ Covered | `_AllNonAcgt_Returns0` |
| M18 | Salt Empty/Null → 0 | ✅ Covered | `_EmptyPrimer_Returns0` + `_NullPrimer_Returns0` |
| M19 | All Same Base 16bp | ✅ Covered | `_MarmurDoty_AllSameBase16bp_ReturnsExpected` |
| INV-1 | Non-negative | ✅ Covered | `_AlwaysNonNegative` |
| INV-2 | Higher GC → Higher Tm | ✅ Covered | `_HigherGC_ProducesHigherTm` |
| INV-3 | Salt additive | ✅ Covered | `_50mM_AppliesCorrection` |
| INV-4 | Case insensitive | ✅ Covered | `_LowercaseInput…` + `_MixedCaseInput…` |

**Additional tests** (not spec-required, but valuable):
- `_Wallace_SingleA_Returns2`, `_Wallace_SingleG_Returns4` — single-base edge cases
- `_Wallace_ACGT_Returns12` — 4bp mixed
- `_MarmurDoty_25bp_ReturnsValidRange` — 25bp primer
- `_MixedCaseInput_MatchesUppercase` — mixed case (exact 24.0)
- `_ManyNonAcgt_UsesOnlyValidBases` — many N characters
- 5× `ThermoConstants_*` — direct constant/helper verification

---

## 9. Sign-off

- [x] TestSpec reviewed
- [x] Evidence documented
- [x] External sources verified (Sigma-Aldrich/Merck, Wikipedia, Marmur & Doty 1962)
- [x] All assumptions eliminated (see Section 3: Defined Behaviors)
- [x] Tests implemented (34 tests in PrimerDesigner_MeltingTemperature_Tests.cs)
- [x] All tests pass
