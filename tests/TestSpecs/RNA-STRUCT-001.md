# RNA-STRUCT-001 Test Specification

**Test Unit:** RNA-STRUCT-001  
**Title:** Secondary Structure Prediction  
**Area:** RnaStructure  
**Created:** 2026-02-05  
**Status:** Complete

---

## 1. Test Unit Definition

### 1.1 Scope

Testing RNA secondary structure prediction algorithms including the Nussinov/Zuker-style dynamic programming approach, stem-loop detection, energy calculations, and notation handling.

### 1.2 Methods Under Test

| Method | Class | Type | Priority |
|--------|-------|------|----------|
| `PredictStructure(sequence, ...)` | RnaSecondaryStructure | Canonical | Must |
| `FindStemLoops(sequence, ...)` | RnaSecondaryStructure | Canonical | Must |
| `CalculateMinimumFreeEnergy(sequence)` | RnaSecondaryStructure | Canonical | Must |
| `ToDotBracket` / `FromDotBracket` | RnaSecondaryStructure | Notation | Must |
| `CanPair(base1, base2)` | RnaSecondaryStructure | Helper | Must |
| `GetBasePairType(base1, base2)` | RnaSecondaryStructure | Helper | Must |
| `CalculateStemEnergy(...)` | RnaSecondaryStructure | Energy | Should |
| `CalculateHairpinLoopEnergy(...)` | RnaSecondaryStructure | Energy | Should |
| `DetectPseudoknots(basePairs)` | RnaSecondaryStructure | Structural | Should |
| `ParseDotBracket(notation)` | RnaSecondaryStructure | Parse | Should |
| `ValidateDotBracket(notation)` | RnaSecondaryStructure | Validation | Should |

---

## 2. Test Categories

### 2.1 Must Tests (Required for Complete)

#### Base Pairing (Evidence: Wikipedia, IUPAC)

| Test ID | Description | Expected |
|---------|-------------|----------|
| BP-001 | Watson-Crick A-U pair | CanPair returns true |
| BP-002 | Watson-Crick U-A pair | CanPair returns true |
| BP-003 | Watson-Crick G-C pair | CanPair returns true |
| BP-004 | Watson-Crick C-G pair | CanPair returns true |
| BP-005 | Wobble G-U pair | CanPair returns true |
| BP-006 | Wobble U-G pair | CanPair returns true |
| BP-007 | Non-pairing A-A | CanPair returns false |
| BP-008 | Non-pairing A-G | CanPair returns false |
| BP-009 | BasePairType for WC pairs | Returns WatsonCrick |
| BP-010 | BasePairType for wobble | Returns Wobble |
| BP-011 | BasePairType for non-pair | Returns null |

#### Stem-Loop Detection (Evidence: Wikipedia, Nussinov 1980)

| Test ID | Description | Expected |
|---------|-------------|----------|
| SL-001 | Simple hairpin GGGAAAACCC | Finds stem-loop |
| SL-002 | No complement poly-A | Returns empty |
| SL-003 | Too short sequence | Returns empty |
| SL-004 | With wobble pairs | Includes G-U pairs |
| SL-005 | Without wobble pairs | Excludes G-U pairs |
| SL-006 | Multiple potential hairpins | Finds multiple |
| SL-007 | Tetraloop detection | Finds 4nt loop |
| SL-008 | Dot-bracket is generated | Contains ( and ) |

#### MFE Calculation (Evidence: Turner 2004, Zuker 1981)

| Test ID | Description | Expected |
|---------|-------------|----------|
| MFE-001 | Simple hairpin | MFE < 0 |
| MFE-002 | No structure (poly-A) | MFE = 0 |
| MFE-003 | Empty sequence | MFE = 0 |
| MFE-004 | Longer stem more stable | MFE(long) ≤ MFE(short) |
| MFE-005 | Stem energy negative | < 0 for valid stem |
| MFE-006 | Single pair energy | = 0 (no stacking) |

#### Dot-Bracket Notation (Evidence: Wikipedia)

| Test ID | Description | Expected |
|---------|-------------|----------|
| DB-001 | Parse simple (((...))) | 3 pairs extracted |
| DB-002 | Parse empty ..... | 0 pairs |
| DB-003 | Parse multiple brackets | All types parsed |
| DB-004 | Validate balanced | Returns true |
| DB-005 | Validate unbalanced | Returns false |
| DB-006 | Validate extra closing | Returns false |

#### Structure Prediction (Evidence: Nussinov 1980, Zuker 1981)

| Test ID | Description | Expected |
|---------|-------------|----------|
| SP-001 | Simple hairpin prediction | Returns valid structure |
| SP-002 | Dot-bracket is valid | Passes validation |
| SP-003 | Has base pairs for structured | Non-empty pairs |
| SP-004 | Empty sequence | Returns empty structure |
| SP-005 | Non-overlapping structures | No overlap in selected |

### 2.2 Should Tests (Recommended)

#### Energy Calculations

| Test ID | Description | Expected |
|---------|-------------|----------|
| EN-001 | Tetraloop GAAA has bonus | Lower energy than AAAA |
| EN-002 | All-C loop has penalty | Higher energy than AAAA |
| EN-003 | Stem stacking is negative | Stabilizing |

#### Pseudoknot Detection

| Test ID | Description | Expected |
|---------|-------------|----------|
| PK-001 | Non-crossing pairs | No pseudoknot |
| PK-002 | Crossing pairs | Pseudoknot detected |

#### Utility Functions

| Test ID | Description | Expected |
|---------|-------------|----------|
| UT-001 | Probability in [0,1] | Valid range |
| UT-002 | MFE structure high probability | > 0.5 |
| UT-003 | Random RNA correct length | Length matches |
| UT-004 | Random RNA valid bases | Only ACGU |
| UT-005 | Random RNA approximate GC | Within tolerance |

### 2.3 Could Tests (Optional)

| Test ID | Description |
|---------|-------------|
| CT-001 | tRNA-like structure analysis |
| CT-002 | Case insensitivity |
| CT-003 | Performance on long sequences |

---

## 3. Test Invariants

### 3.1 Always True

1. `ValidateDotBracket(structure.DotBracket)` = true for any prediction
2. `structure.DotBracket.Length` = `sequence.Length`
3. MFE ≤ 0 for any structurable sequence
4. Selected stem-loops never overlap
5. All base pairs are valid (WC or Wobble)

### 3.2 Boundary Conditions

1. Empty sequence → empty structure, MFE = 0
2. Sequence shorter than minStemLength * 2 + minLoopSize → no stem-loops
3. No complementary bases → no structure

---

## 4. Edge Cases

| Category | Case | Input Example | Expected |
|----------|------|---------------|----------|
| Empty | Empty string | `""` | Empty structure |
| Empty | Null | `null` | Empty structure |
| Length | Too short | `"GC"` | No stem-loops |
| Content | No complement | `"AAAAAAA"` | No structure |
| Content | All same | `"GGGGGGG"` | No structure |
| Case | Lowercase | `"gggaaaaccc"` | Same as uppercase |
| Loop | Minimum loop | 3 nt loop | Should work |
| Loop | All-C loop | `"CCCC"` | Penalty applied |

---

## 5. Audit of Existing Tests

### 5.1 Current Coverage

| Category | Tests | Status | Notes |
|----------|-------|--------|-------|
| Base Pairing | 11 tests | Covered | Complete |
| Stem-Loop Finding | 8 tests | Covered | Good coverage |
| Energy Calculation | 8 tests | Covered | Includes tetraloop |
| Structure Prediction | 5 tests | Covered | Basic coverage |
| Pseudoknot Detection | 2 tests | Covered | Minimal |
| Dot-Bracket | 5 tests | Covered | Validation included |
| Inverted Repeats | 1 test | Smoke | Delegates to RepeatFinder |
| Utility | 4 tests | Covered | Random generation |
| Integration | 3 tests | Covered | tRNA-like, case handling |

### 5.2 Assessment

- **Coverage:** Strong overall
- **Evidence Basis:** Tests align with documented algorithm behavior
- **Missing:** Additional edge cases for energy calculation
- **Weak:** None identified
- **Duplicates:** None

---

## 6. Consolidation Plan

### 6.1 Test File Structure

- **Canonical File:** `RnaSecondaryStructureTests.cs`
- **Location:** `tests/Seqeron/Seqeron.Genomics.Tests/`

### 6.2 Actions

1. Keep existing test structure (well-organized by region)
2. Add edge case tests for empty/null inputs to ensure coverage
3. Add explicit reference data validation tests
4. Ensure all tests have evidence-based comments

---

## 7. Open Questions

None - implementation follows standard Nussinov/Zuker approach.

---

## 8. Decisions

| Decision | Rationale |
|----------|-----------|
| Use Turner 2004 parameters | Standard reference for RNA thermodynamics |
| Minimum loop size = 3 | Biological constraint, widely accepted |
| Allow wobble by default | G-U pairs are valid in RNA |
| Greedy non-overlapping | Practical for simple prediction |

---

## 9. Sign-off

- [x] Evidence document created
- [x] Test specification complete
- [x] Existing tests audited
- [x] Consolidation plan defined
- [x] All Must tests verified
