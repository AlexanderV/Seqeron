# RNA-ENERGY-001 Test Specification

**Test Unit:** RNA-ENERGY-001  
**Title:** Free Energy Calculation  
**Area:** RnaStructure  
**Created:** 2026-02-05  
**Status:** Complete

---

## 1. Test Unit Definition

### 1.1 Scope

Testing RNA free energy calculation algorithms including stacking energy, hairpin loop energy, minimum free energy (MFE), and related thermodynamic computations based on the Turner 2004 nearest-neighbor parameters.

### 1.2 Methods Under Test

| Method | Class | Type | Priority |
|--------|-------|------|----------|
| `CalculateStemEnergy(sequence, basePairs)` | RnaSecondaryStructure | Canonical | Must |
| `CalculateHairpinLoopEnergy(loop, closing5, closing3)` | RnaSecondaryStructure | Canonical | Must |
| `CalculateMinimumFreeEnergy(rnaSequence, minLoopSize)` | RnaSecondaryStructure | Canonical | Must |
| `CalculateStackingEnergy(bp1, bp2)` | RnaSecondaryStructure | Helper | Should |

---

## 2. Test Categories

### 2.1 Must Tests (Required for Complete)

#### Stem Energy Calculation (Evidence: Turner 2004, NNDB)

| Test ID | Description | Expected |
|---------|-------------|----------|
| SE-001 | Stem with multiple WC pairs returns negative energy | Energy < 0 |
| SE-002 | Single base pair returns zero (no stacking) | Energy = 0 |
| SE-003 | GC-rich stem more stable than AU-rich | GC_energy < AU_energy |
| SE-004 | Empty base pairs list returns zero | Energy = 0 |
| SE-005 | Wobble pairs contribute energy | Energy < 0 (less than WC) |

#### Hairpin Loop Energy (Evidence: Turner 2004, NNDB)

| Test ID | Description | Expected |
|---------|-------------|----------|
| HL-001 | GNRA tetraloop (GAAA) has bonus | Energy(GAAA) < Energy(AAAA) |
| HL-002 | All-C loop has penalty | Energy(CCCC) > Energy(AAAA) |
| HL-003 | Loop size 3 is valid minimum | Returns positive energy |
| HL-004 | Loop size 4 applies tetraloop check | Applies bonus if matched |
| HL-005 | Different loop sizes all positive | All energies > 0 (NOTE: NOT monotonic per Turner 2004) |

#### Minimum Free Energy (Evidence: Zuker 1981, Turner 2004)

| Test ID | Description | Expected |
|---------|-------------|----------|
| MFE-001 | Simple hairpin has negative MFE | MFE < 0 |
| MFE-002 | No structure possible (poly-A) returns 0 | MFE = 0 |
| MFE-003 | Empty/null sequence returns 0 | MFE = 0 |
| MFE-004 | Longer stem more stable than shorter | MFE(long) ≤ MFE(short) |
| MFE-005 | GC-rich hairpin more stable | GC_MFE < AU_MFE |

### 2.2 Should Tests (Recommended)

#### Stacking Energy Individual (Evidence: Turner 2004)

| Test ID | Description | Expected |
|---------|-------------|----------|
| ST-001 | GC/CG stack returns -3.42 kcal/mol (approx) | Near reference value |
| ST-002 | AU/UA stack returns -1.10 kcal/mol (approx) | Near reference value |
| ST-003 | All stacking energies are negative | All < 0 |

#### Energy Parameter Consistency

| Test ID | Description | Expected |
|---------|-------------|----------|
| EC-001 | Multiple identical stacks sum correctly | Linear scaling |
| EC-002 | Mixed stacks combine properly | Sum of individual stacks |
| EC-003 | Temperature at 37°C (310.15 K) | Standard conditions |

### 2.3 Could Tests (Optional)

| Test ID | Description | Expected |
|---------|-------------|----------|
| CT-001 | Very long loop extrapolation | Uses Jacobson-Stockmayer |
| CT-002 | G-U wobble stacking | Reduced stability vs WC |
| CT-003 | Structure probability calculation | Valid range [0,1] |

---

## 3. Edge Cases and Invariants

### 3.1 Edge Cases

| Case | Test | Expected Result |
|------|------|-----------------|
| Empty sequence | MFE("") | 0 |
| Null sequence | MFE(null) | 0 |
| Single base pair | StemEnergy([1 pair]) | 0 |
| Minimum loop size | LoopEnergy("AAA", ...) | Positive value |
| Maximum practical loop | LoopEnergy(30+ nt) | Extrapolated value |

### 3.2 Invariants

1. **I-001**: Stacking energies for Watson-Crick pairs are always negative
2. **I-002**: Hairpin loop initiation energies are always positive
3. **I-003**: MFE for structured RNA is negative; for unstructured is 0
4. **I-004**: Longer stems have more negative (more stable) energy
5. **I-005**: GNRA tetraloops are more stable than non-special tetraloops

---

## 4. Test Data

### 4.1 Reference Sequences

| ID | Sequence | Structure | Notes |
|----|----------|-----------|-------|
| RS-001 | GGGAAAACCC | Hairpin 3bp+4loop | Simple test case |
| RS-002 | GCGCAAAAGCGC | Hairpin 4bp+4loop | GC-rich stem |
| RS-003 | AAAAAAAAA | None | No structure possible |
| RS-004 | GGGCGAAAGCCC | Hairpin 4bp+GNRA | GNRA tetraloop (GAAA loop) |
| RS-005 | GCGCCCCCGCGC | Hairpin 4bp+C-loop | All-C loop penalty |

### 4.2 Reference Energy Values (Turner 2004)

| Stack | ΔG°37 (kcal/mol) | Source |
|-------|------------------|--------|
| GC/CG | -3.42 | NNDB |
| CG/GC | -2.36 | NNDB |
| GG/CC | -3.26 | NNDB |
| AU/UA | -1.10 | NNDB |
| AA/UU | -0.93 | NNDB |

---

## 5. Coverage Classification

### 5.1 Location
- **File:** `RnaSecondaryStructureTests.cs`
- **Regions:** "Energy Calculation Tests", "Turner 2004 NNDB Parameter Validation"

### 5.2 Stem Energy — 12 test methods

| # | Test Method | Spec ID | Status |
|---|-------------|---------|--------|
| 1 | `CalculateStemEnergy_BasePairs_ReturnsNegative` | SE-001 | ✅ |
| 2 | `CalculateStemEnergy_SinglePair_ReturnsZero` | SE-002 | ✅ |
| 3 | `CalculateStemEnergy_GCRichVsAURich_GCMoreStable` | SE-003 | ✅ |
| 4 | `CalculateStemEnergy_EmptyBasePairs_ReturnsZero` | SE-004 | ✅ |
| 5 | `CalculateStemEnergy_WobblePair_ContributesNegativeEnergy` | SE-005 | ✅ |
| 6 | `CalculateStemEnergy_WatsonCrickStacking_MatchesNNDB` (×8) | ST-001,ST-002,ST-003 | ✅ |
| 7 | `CalculateStemEnergy_GUWobbleStacking_IncludesDestabilizing` | CT-002 | ✅ |
| 8 | `CalculateStemEnergy_TerminalAUPenalty_MatchesNNDB` | I-001 | ✅ |
| 9 | `CalculateStemEnergy_GGUC_CUGG_3Stack_MatchesNNDB` | CT-002 | ✅ |
| 10 | `CalculateStemEnergy_GGUC_CUGG_InLongerStem_MatchesNNDB` | EC-002 | ✅ |
| 11 | `CalculateStemEnergy_IdenticalStacks_SumLinearly` | EC-001 | ✅ |

### 5.3 Hairpin Loop Energy — 13 test methods

| # | Test Method | Spec ID | Status |
|---|-------------|---------|--------|
| 12 | `CalculateHairpinLoopEnergy_Tetraloop_HasBonus` | HL-001 | ✅ |
| 13 | `CalculateHairpinLoopEnergy_AllC_HasPenalty` | HL-002 | ✅ |
| 14 | `CalculateHairpinLoopEnergy_MinimumLoop_ReturnsPositive` | HL-003 | ✅ |
| 15 | `CalculateHairpinLoopEnergy_DifferentSizes_AllPositive` | HL-005 | ✅ |
| 16 | `CalculateHairpinLoopEnergy_UNCGTetraloop_HasBonus` | HL-004 | ✅ |
| 17 | `CalculateHairpinLoopEnergy_Initiation_MatchesNNDB` (×5) | HL-003,HL-005,I-002 | ✅ |
| 18 | `CalculateHairpinLoopEnergy_SpecialTetraloop_MatchesNNDB` (×4) | HL-004,I-005 | ✅ |
| 19 | `CalculateHairpinLoopEnergy_NonSpecialTetraloop_UsesStandardModel` | HL-004 | ✅ |
| 20 | `CalculateHairpinLoopEnergy_AllCPenalty_MatchesNNDB` | HL-002 | ✅ |
| 21 | `CalculateHairpinLoopEnergy_MismatchBonuses_MatchesNNDB` | I-005 | ✅ |
| 22 | `CalculateHairpinLoopEnergy_LargeLoop_UsesCorrectExtrapolation` | CT-001 | ✅ |
| 23 | `CalculateHairpinLoopEnergy_TerminalMismatch_IsAdditiveWithBonuses` | HL-004 | ✅ |
| 24 | `CalculateHairpinLoopEnergy_SpecialGUClosure_AppliesBonus` | CT-002 | ✅ |

### 5.4 Minimum Free Energy — 10 test methods

| # | Test Method | Spec ID | Status |
|---|-------------|---------|--------|
| 25 | `CalculateMinimumFreeEnergy_SimpleHairpin_ReturnsNegative` | MFE-001 | ✅ |
| 26 | `CalculateMinimumFreeEnergy_NoStructure_ReturnsZero` | MFE-002 | ✅ |
| 27 | `CalculateMinimumFreeEnergy_EmptySequence_ReturnsZero` | MFE-003 | ✅ |
| 28 | `CalculateMinimumFreeEnergy_LongerStem_MoreStable` | MFE-004,I-004 | ✅ |
| 29 | `CalculateMinimumFreeEnergy_GCRichHairpin_MoreStable` | MFE-005 | ✅ |
| 30 | `CalculateMinimumFreeEnergy_SimpleHairpin_MatchesTurnerManualCalc` | MFE-001 | ✅ |
| 31 | `CalculateMinimumFreeEnergy_FourPairGC_MatchesTurner` | MFE-001 | ✅ |
| 32 | `CalculateMinimumFreeEnergy_AUStem_IncludesTerminalPenalty` | MFE-005 | ✅ |
| 33 | `CalculateMinimumFreeEnergy_InnerAUPenalty_MatchesNNDB` | I-003 | ✅ |
| 34 | `CalculateMinimumFreeEnergy_InnerAUPenalty_NNDBExample_MatchesReference` | I-003 | ✅ |

### 5.5 Standard Conditions — 1 test method

| # | Test Method | Spec ID | Status |
|---|-------------|---------|--------|
| 35 | `EnergyCalculation_UsesStandardTemperature_37C` | EC-003 | ✅ |

### 5.6 Classification Summary

- ✅ Covered: 35 test methods (all spec IDs)
- ❌ Missing: 0
- ⚠ Weak: 0 (EC-001 strengthened with exact linear scaling test; SE-005 added with exact NNDB value)
- 🔁 Duplicate: 0

---

## 6. Verification Status

All parameters hand-verified against NNDB Turner 2004. No assumptions remain.
See Evidence document for full verification matrix.

---

## 7. Test File Location

**Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/RnaSecondaryStructureTests.cs`  
**Region:** Energy Calculation Tests

---

## 8. Definition of Done

- [x] Evidence document created
- [x] Test specification created
- [x] All Must tests implemented
- [x] All Should tests implemented
- [x] All Could tests implemented
- [x] Edge cases covered
- [x] Invariants covered
- [x] Coverage classification complete: 0 missing, 0 weak, 0 duplicate
- [x] Tests pass with zero warnings
- [x] All parameters verified against NNDB Turner 2004
- [x] No assumptions remain — all values sourced from NNDB
- [x] Inner AU/GU terminal penalty in MFE validated
- [x] GGUC/CUGG 3-stack context implemented and tested
- [x] Special GU closure detected in MFE from sequence context
