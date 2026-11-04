# Test Specification: PRIMER-TM-001

## Test Unit: Melting Temperature Calculation

**Area:** MolTools
**Status:** ☑ Complete
**Created:** 2026-01-22
**Last Verified:** 2026-10-08 (review-2026-09 B07: salt-adjusted Tm corrected; U read as T per F16; status set complete in audit round 4, A4-4)
**Evidence Sources:** Thein & Wallace (1986), Marmur & Doty (1962), OligoCalc (Kibbe 2007, NAR 35:W43), Schildkraut & Lifson (1965), Biopython `Bio.SeqUtils.MeltingTemp` 1.88, Sigma-Aldrich/Merck Technical Docs

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
| `CalculateSaltCorrection(Na_mM)` | 16.6 × log10(Na/1000) | Absolute Schildkraut–Lifson term (Biopython method 1); NOT used by the primer Tm |
| `CalculateOligoCalcSaltAdjustedTm(at, gc, Na_M)` | N<14: 2AT+4GC+16.6·log10(Na/0.050); N≥14: 100.5+41·GC/N−820/N+16.6·log10(Na) | OligoCalc salt-adjusted Tm |

Guards (B07 audit round 3, A3-16; Biopython `Tm_GC`/`salt_correction` raise `ValueError` for [Na+] ≤ 0): the salt
helpers throw `ArgumentOutOfRangeException` for [Na+] ≤ 0 / NaN / ∞, the count/length helpers for negative
counts/lengths, `CalculateSaltAdjustedTm` also for GC fraction ∉ [0, 1] — `ThermoConstants_SaltHelpers_NonPositiveSodium_Throw`,
`ThermoConstants_CountHelpers_NegativeCountOrLength_Throw`, `ThermoConstants_CalculateSaltAdjustedTm_GcFractionOutOfRange_Throws`,
`ThermoConstants_SaltHelpers_PositiveInputs_Unchanged`.

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

### 2.3 Salt-Adjusted Tm (OligoCalc)

**Source:** OligoCalc (Kibbe 2007, Nucleic Acids Res 35:W43–W46), "Salt Adjusted" Tm; salt slope 16.6·log10 from
Schildkraut & Lifson (1965). The basic formulas of §2.1/§2.2 are OligoCalc's "Basic" Tm and **assume 50 nM primer,
50 mM Na+, pH 7.0**, so an absolute 16.6·log10([Na+]) term must not be added on top of them.

**Formula ([Na+] in mol/L; the API takes mM):**
- N < 14: Tm = 2(A+T) + 4(G+C) − 16.6·log10(0.050) + 16.6·log10([Na+])
- N ≥ 14: Tm = 100.5 + 41·(G+C)/N − 820/N + 16.6·log10([Na+])

Result rounded to 1 decimal. [Na+] ≤ 0 / NaN / ∞ → `ArgumentOutOfRangeException`.

**Reference cross-check:** OligoCalc output for `GAGCAGGATCCCTATAGAGTGACAAAAGGATCTTGGTCC` @ 50 mM: basic 67.6 °C,
salt-adjusted 78 °C (ours 67.633 / 77.852). 20-mer `ACGTACGTACGTACGTACGT` @ 50 mM: ours 58.4; Biopython `Tm_GC(valueset=7)` 50.40;
primer3 `calc_tm` (NN, 50 mM, no Mg) 53.99. **Fixed defect (2026-09):** the previous implementation added
16.6·log10([Na+]/1000) to the basic Tm (double-counting salt): 20-mer 50 %GC @ 50 mM → 30.2 °C.

**Defined Behavior:** returns 0 for empty/null input or no A/C/G/T.

---

## 3. Defined Behaviors

### 3.1 Input Alphabet

Standard DNA bases (A, C, G, T) are recognized; RNA uracil (U) is read as T, as Biopython
`MeltingTemp._check` back-transcribes RNA (review-2026-09 B07 F16 / R24). All other characters — including
IUPAC ambiguity codes (N, R, Y, etc.) — are ignored.
Both threshold determination and formula computation use the count of valid A/C/G/T/U bases only.

| Input | Valid Bases | Behavior |
|-------|------------|----------|
| `"ACNGT"` | A,C,G,T (4 valid) | N ignored; Wallace: 2×2 + 4×2 = 12 |
| `"ACGTNNNNACGT"` | 8 valid | N's ignored; Wallace: 2×4 + 4×4 = 24 |
| `"NNNNN"` | 0 valid | Returns 0 |
| `"ACGUACGU"` | 8 valid (U read as T) | Wallace: 2×4 + 4×4 = 24 (= Biopython `Tm_Wallace`) |

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

#### M12: Salt-Adjusted - 50mM (≥14 nt)
**Evidence:** OligoCalc salt-adjusted formula
```
Input: primer="ACGTACGTACGTACGTACGT", Na=50mM
Expected: 100.5 + 41×10/20 − 820/20 + 16.6×log10(0.05) = 58.403 → 58.4
```

#### M13: Salt-Adjusted - 10mM (≥14 nt)
```
Input: primer="ACGTACGTACGTACGTACGT", Na=10mM
Expected: 100.5 + 20.5 − 41 − 33.2 = 46.8
```

#### M14: Salt-Adjusted - 200mM (≥14 nt)
```
Input: primer="ACGTACGTACGTACGTACGT", Na=200mM
Expected: 68.397 → 68.4
```

#### M14b: Salt-Adjusted short oligo / OligoCalc worked example
```
"ACGTACGT" @ 50 mM → 24.0 (= Wallace); @ 10 mM → 12.4; @ 200 mM → 34.0; "ACGT" @ 1000 mM → 33.6
OligoCalc 39-mer GAGCAGGATCCCTATAGAGTGACAAAAGGATCTTGGTCC @ 50 mM → 77.9 (OligoCalc: 78)
Na ≤ 0, NaN, ∞ → ArgumentOutOfRangeException
```

#### M15: Non-ACGT Characters Ignored
**Evidence:** Defined behavior (Section 3.1)
```
Input: "ACNGT" (4 valid bases: A,C,G,T)
Expected: Wallace 2×2 + 4×2 = 12.0
```

#### M16: RNA Base (U) Read as T
**Evidence:** Biopython 1.88 `Tm_Wallace("ACGUACGU")` = 24.0 (`MeltingTemp._check` back-transcribes U → T; B07 F16)
```
Input: "ACGUACGU" (8 valid bases, U counted as T)
Expected: Wallace 2×4 + 4×4 = 24.0
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
| INV-3 | Salt-adjusted Tm follows the OligoCalc formula; for N<14 at 50 mM it equals the basic Tm; +16.6 °C per decade of [Na+] | property tests |
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
| RNA (U) | "ACGUACGU" | U read as T (24.0, Biopython parity) |

---

## 7. External Source Verification Summary

| Item | Our Implementation | Sigma-Aldrich/Merck | Status |
|------|-------------------|---------------------|--------|
| Short oligo formula | Tm = 2(A+T) + 4(G+C) | Tm = 2(A+T) + 4(G+C) − 7 | **Variant** — −7 correction omitted by design |
| Short oligo threshold | < 14 valid bases | ≤ 14 bases | Aligned (both use 14 as boundary) |
| Long primer formula | Marmur-Doty: 64.9 + 41(GC−16.4)/N | Nearest-neighbor (SantaLucia 1998) | **Simplified** — Marmur-Doty is a simpler, well-published alternative |
| Salt adjustment | OligoCalc salt-adjusted Tm (relative to the 50 mM basis for N<14) | Integrated into NN formula | Matches OligoCalc output |
| Non-ACGT handling | Ignored (only ACGT counted) | Not documented (clean input expected) | Defined behavior |
| RNA (U) | Read as T (Biopython `MeltingTemp._check`) | Not applicable (DNA tool) | Matches Biopython |

### Known Variant: −7 Correction Factor

Sigma-Aldrich's "Basic" method includes a −7°C correction factor for
oligonucleotides used in solution (as opposed to membrane hybridization).
Our implementation uses the original Wallace rule without this correction.
This is a **deliberate design choice** — the Wallace rule as published by
Thein & Wallace (1986) does not include the correction.
For a different [Na+] use `CalculateMeltingTemperatureWithSalt` (OligoCalc salt-adjusted Tm); for accuracy use `CalculateMeltingTemperatureNN`.

---

## 8. Coverage Classification

All spec test cases verified against `PrimerDesigner_MeltingTemperature_Tests.cs` (44 test methods).

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
| M12–M14 | Salt-adjusted ≥14 nt | ✅ Covered | `CalculateMeltingTemperatureWithSalt_Long_OligoCalcSaltAdjusted` (50/10/200/1000 mM) |
| M14b | Salt-adjusted <14 nt, OligoCalc example, invalid Na | ✅ Covered | `_Short_WallaceRelativeCorrection`, `_OligoCalcWorkedExample_39mer`, `_Short_At50mM_EqualsBasicWallace`, `_IncreasesBy16_6PerDecade`, `_InvalidSodium_Throws` |
| M15 | Non-ACGT Ignored | ✅ Covered | `_NonAcgtIgnored_OnlyValidBasesCounted` |
| M16 | RNA U read as T | ✅ Covered | `_RnaUracil_ReadAsThymine_MatchesBiopythonTmWallace` |
| M17 | All Non-Standard → 0 | ✅ Covered | `_AllNonAcgt_Returns0` |
| M18 | Salt Empty/Null → 0 | ✅ Covered | `_EmptyPrimer_Returns0` + `_NullPrimer_Returns0` |
| M19 | All Same Base 16bp | ✅ Covered | `_MarmurDoty_AllSameBase16bp_ReturnsExpected` |
| INV-1 | Non-negative | ✅ Covered | `_AlwaysNonNegative` |
| INV-2 | Higher GC → Higher Tm | ✅ Covered | `_HigherGC_ProducesHigherTm` |
| INV-3 | Salt-adjusted formula | ✅ Covered | `PrimerProbeProperties.MeltingTemperatureWithSalt_MatchesOligoCalcSaltAdjusted`, `_ShortOligoAt50mM_EqualsBaseTm`, `_50mMTo1M_Adds16_6LogTwenty` |
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
- [x] Tests implemented (44 test methods in PrimerDesigner_MeltingTemperature_Tests.cs)
- [x] All tests pass
