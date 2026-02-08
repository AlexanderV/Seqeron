# Test Specification: PRIMER-TM-001

## Test Unit: Melting Temperature Calculation

**Area:** MolTools
**Status:** Draft
**Created:** 2026-01-22
**Evidence Sources:** Wikipedia (Nucleic acid thermodynamics, DNA melting), SantaLucia (1998), Marmur & Doty (1962)

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

**Applicability:** Primers < 14 bp

**Rationale:**
- A-T pairs form 2 hydrogen bonds (weaker)
- G-C pairs form 3 hydrogen bonds (stronger)
- Simple approximation for very short oligos

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

**Source:** Marmur & Doty (1962), "Determination of the base composition of deoxyribonucleic acid from its thermal denaturation temperature"

**Formula:** Tm = 64.9 + 41 × (GC - 16.4) / N

Where:
- GC = number of G and C bases
- N = total sequence length

**Applicability:** Primers ≥ 14 bp

**Test Values (Calculated):**
| Sequence (20bp) | GC count | GC% | Expected Tm |
|-----------------|----------|-----|-------------|
| All A/T | 0 | 0% | 64.9 + 41×(-16.4)/20 = 31.3°C |
| 50% GC (10 each) | 10 | 50% | 64.9 + 41×(-6.4)/20 = 51.8°C |
| All G/C | 20 | 100% | 64.9 + 41×(3.6)/20 = 72.3°C |

### 2.3 Salt Correction

**Source:** Owczarzy et al. (2004), general PCR literature

**Formula:** Salt correction = 16.6 × log10([Na+]/1000)

Where [Na+] is in mM (typical PCR: 50 mM)

**Example:** At 50 mM Na+: 16.6 × log10(0.05) = 16.6 × (-1.301) ≈ -21.6°C

---

## 3. Test Cases

### 3.1 Must Tests (Evidence-Based)

#### M1: Empty Input Handling
**Rationale:** Defensive programming - empty primers should return 0
```
Input: ""
Expected: 0.0
```

#### M2: Null Input Handling
**Rationale:** Defensive programming - null primers should return 0
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
**Evidence:** Implementation uses < 14 as threshold
```
Input: "ACGTACGTACGTA" (13 bp, 7 A/T, 6 G/C)
Expected: 2×7 + 4×6 = 38.0
```

#### M7: Marmur-Doty - Boundary (14 bp, switches to MD)
**Evidence:** Implementation threshold at 14
```
Input: "ACGTACGTACGTAC" (14 bp, 7 A/T, 7 G/C)
Expected: 64.9 + 41×(7-16.4)/14 = 37.4
```

#### M8: Marmur-Doty - Typical Primer (20 bp, 50% GC)
**Evidence:** Marmur & Doty (1962)
```
Input: "ACGTACGTACGTACGTACGT" (20 bp, 10 G/C)
Expected: 64.9 + 41×(10-16.4)/20 = 51.78
```

#### M9: Marmur-Doty - Low GC (20 bp)
**Evidence:** Marmur & Doty (1962), clamped to 0 minimum
```
Input: "ATATATATATATATATATATAT" (20 bp A/T only)
Calculated: 64.9 + 41×(0-16.4)/20 = 31.3 (if > 0)
Expected: Max(0, 31.3) = 31.3
```

#### M10: Marmur-Doty - High GC (20 bp)
**Evidence:** Marmur & Doty (1962)
```
Input: "GCGCGCGCGCGCGCGCGCGC" (20 bp, 20 G/C)
Expected: 64.9 + 41×(20-16.4)/20 = 72.28
```

#### M11: Case Insensitivity
**Rationale:** DNA sequence may be lowercase
```
Input: "atatatat"
Expected: Same as "ATATATAT" = 16.0
```

#### M12: Salt Correction - Standard 50mM
**Evidence:** Standard PCR salt concentration
```
Input: primer="ACGTACGTACGTACGTACGT", Na=50mM
Base Tm: 51.78
Salt correction: 16.6 × log10(50/1000) = -21.6
Expected: ~30.2 (rounded to 1 decimal)
```

#### M13: Salt Correction - Low Salt (10mM)
**Evidence:** Salt affects Tm
```
Salt correction: 16.6 × log10(10/1000) = -33.2
Expected: Base Tm - 33.2
```

#### M14: Salt Correction - High Salt (200mM)
**Evidence:** Higher salt stabilizes duplex
```
Salt correction: 16.6 × log10(200/1000) = -11.6
Expected: Base Tm - 11.6 (higher than 50mM)
```

### 3.2 Should Tests

#### S1: Non-ACGT Characters Ignored
**ASSUMPTION:** Implementation may count only valid bases
```
Input: "ACNGTACNGT" (with N)
Verify: Returns reasonable value (only counts A/C/G/T)
```

#### S2: Single Nucleotide
```
Input: "A"
Expected: 2.0 (Wallace: 2×1 + 4×0)
```

#### S3: Very Long Sequence (Performance)
**ASSUMPTION:** O(n) complexity maintained
```
Input: 1000bp sequence
Verify: Completes in reasonable time, returns valid Tm
```

### 3.3 Could Tests

#### C1: RNA Sequences (U instead of T)
**ASSUMPTION:** May or may not be supported
```
Input: "ACGUACGU"
Verify: Behavior is consistent (either treats U as T or excludes)
```

---

## 4. Invariants

| ID | Invariant | Validation |
|----|-----------|------------|
| INV-1 | Result ≥ 0 for any valid input | Assert.That(tm, Is.GreaterThanOrEqualTo(0)) |
| INV-2 | Higher GC content → Higher Tm | Compare equal-length sequences |
| INV-3 | Longer sequence (same GC%) → Tm approaches asymptote | Compare 20bp vs 100bp at 50% GC |
| INV-4 | Salt correction is additive | Tm_salt = Tm_base + correction |

---

## 5. Edge Cases

| Case | Input | Expected Behavior |
|------|-------|-------------------|
| Empty string | "" | Returns 0.0 |
| Null | null | Returns 0.0 |
| Single base | "A" | Wallace: 2.0 |
| Boundary 13bp | 13-char string | Uses Wallace |
| Boundary 14bp | 14-char string | Uses Marmur-Doty |
| All same base | "AAAAAAAAAAAAAAAA" | Valid calculation |
| Very long | 1000+ bp | O(n) performance |
| Lowercase | "acgt" | Case-insensitive |

---

## 6. Audit of Existing Tests

### Current Coverage (PrimerDesignerTests.cs)

| Test | Coverage | Assessment |
|------|----------|------------|
| CalculateMeltingTemperature_ShortPrimer_UsesWallaceRule | Wallace 8bp A/T | Covered ✓ |
| CalculateMeltingTemperature_ShortAllGC_HighTm | Wallace 8bp G/C | Covered ✓ |
| CalculateMeltingTemperature_LongPrimer_UsesNearestNeighbor | MD 20bp range check | Weak (no exact value) |
| CalculateMeltingTemperature_EmptyPrimer_Returns0 | Empty input | Covered ✓ |
| CalculateMeltingTemperatureWithSalt_AppliesSaltCorrection | Salt ≠ base | Weak (no exact value) |

### Gaps Identified (All Closed)

1. ~~**Missing:** Null input test~~ ✅ Covered
2. ~~**Missing:** Boundary tests (13bp vs 14bp)~~ ✅ Covered
3. ~~**Missing:** Exact Marmur-Doty value verification~~ ✅ Covered
4. ~~**Missing:** Case insensitivity test~~ ✅ Covered
5. ~~**Missing:** Mixed Wallace test~~ ✅ Covered
6. ~~**Missing:** Salt correction exact values~~ ✅ Covered
7. ~~**Weak:** Long primer test only checks range~~ ✅ Strengthened

---

## 7. Consolidation Plan

### Action: Create Dedicated Test File

**File:** `PrimerDesigner_MeltingTemperature_Tests.cs`

**Approach:**
1. Move 5 existing Tm tests from PrimerDesignerTests.cs to new dedicated file
2. Add missing Must tests (exact formula verification)
3. Add boundary tests (13bp/14bp threshold)
4. Add salt correction verification with exact values
5. Remove redundant/weak tests that don't verify exact values

### Tests to Move
- CalculateMeltingTemperature_ShortPrimer_UsesWallaceRule
- CalculateMeltingTemperature_ShortAllGC_HighTm
- CalculateMeltingTemperature_LongPrimer_UsesNearestNeighbor (refactor)
- CalculateMeltingTemperature_EmptyPrimer_Returns0
- CalculateMeltingTemperatureWithSalt_AppliesSaltCorrection (refactor)

### Tests to Add (Must)
- Null input handling
- Mixed Wallace (ACGTACGT = 24)
- Boundary 13bp (Wallace)
- Boundary 14bp (Marmur-Doty)
- Exact Marmur-Doty values
- Case insensitivity
- Exact salt correction values

---

## 8. Open Questions

None - algorithm behavior is well-defined in literature.

---

## 9. Sign-off

- [ ] TestSpec reviewed
- [ ] Evidence documented
- [ ] Tests implemented
- [ ] All tests pass
- [ ] Checklist updated
