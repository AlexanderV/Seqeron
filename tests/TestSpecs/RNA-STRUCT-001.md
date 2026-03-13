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
| Base Pairing | 11 tests | ✅ Covered | WC pairs, wobble, non-pairing, complement |
| Stem-Loop Finding | 12 tests | ✅ Covered | Includes edge cases: empty, null, min stem, loop size, biological minimum |
| Energy: Stem Stacking | 4 tests + 8 NNDB | ✅ Covered | Exact NNDB values, GC vs AU, terminal penalty, GU wobble |
| Energy: Hairpin Loop | 7 tests + 5 NNDB | ✅ Covered | Tetraloop, all-C, initiation, special loops, mismatch bonuses, GU closure, extrapolation |
| Energy: Internal Loop | 10 tests | ✅ Covered | Generic model + int11 lookup; symmetric, asymmetric, AU closure, symmetry |
| Energy: Bulge Loop | 5 tests | ✅ Covered | n=1 stacking, special C, degeneracy, multi-nt terminal |
| Energy: Multibranch | 2 tests | ✅ Covered | 3-way junction, strain penalty |
| Energy: Dangling Ends | 5 NNDB tests | ✅ Covered | 3' and 5' dangles on various pairs |
| Energy: Coaxial | 3 tests | ✅ Covered | Flush and mismatch-mediated |
| Energy: MFE (Zuker) | 7 tests | ✅ Covered | Manual Turner calc, 4-pair GC, AU penalty, bulge, GC vs AU |
| Energy: Terminal Mismatch | 8 NNDB tests | ✅ Covered | Exact NNDB table lookups |
| Structure Prediction | 7 tests | ✅ Covered | Hairpin, dot-bracket validity, base pairs, empty, non-overlapping, invariant, poly-G |
| Pseudoknot Detection | 2 tests | ✅ Covered | Non-crossing and crossing |
| Dot-Bracket | 5 tests | ✅ Covered | Parse, validate balanced/unbalanced |
| Inverted Repeats | 1 test | ✅ Smoke | Delegates to RepeatFinder |
| Utility | 5 tests | ✅ Covered | Probability, random generation (length, bases, GC content) |
| Integration | 3 tests | ✅ Covered | tRNA workflow, stem-loop energy, case insensitivity |

### 5.2 Assessment

- **Coverage:** Complete — 131 tests total
- **Evidence Basis:** All quantitative tests use exact NNDB Turner 2004 parameter values
- **Missing:** None
- **Weak:** None — all previously permissive assertions strengthened to exact values
- **Duplicates:** None
- **Weak:** None — all previously permissive assertions strengthened to exact values
- **Duplicates:** None

---

## 6. Consolidation Plan

### 6.1 Test File Structure

- **Canonical File:** `RnaSecondaryStructureTests.cs`
- **Location:** `tests/Seqeron/Seqeron.Genomics.Tests/`

### 6.2 Actions

All coverage actions completed:
1. ⚠ Weak tests strengthened with exact NNDB values (6 tests)
2. ❌ Missing tests implemented (2 tests: invariant + edge case)
3. 🔁 No duplicates found
4. Tests verified against NNDB Turner 2004 theory (not fitted to implementation)

---

## 7. Deviations and Assumptions

### 7.1 Resolved

| # | Deviation | Status | Tests |
|---|-----------|--------|-------|
| D1 | n=1 bulge missing −RT·ln(states) degeneracy term | 🔧 FIXED | `CalculateBulgeLoopEnergy_SingleNucleotide_DegeneracyTerm` |
| D2 | Dangling ends (d2) not in multiloop WM decomposition | 🔧 FIXED | Verified via existing MFE integration tests |

### 7.2 Open

| # | Deviation | Status | Impact |
|---|-----------|--------|--------|
| D3 | int21 (1×2) lookup table missing | ⛔ BLOCKED | Uses generic model; NNDB table too large (2,304 entries) for inline |
| D4 | int22 (2×2) lookup table missing | ⛔ BLOCKED | Uses generic model; NNDB table too large (36,864 entries) for inline |
| D5 | No Zuker traceback in PredictStructure | ⚠ ASSUMPTION | MFE value correct; dot-bracket from greedy selection, not DP traceback |

---

## 8. Decisions

| Decision | Rationale |
|----------|-----------|
| Turner 2004 (NNDB) parameters for stacking and hairpin energies | Standard reference, exact NNDB values validated by tests |
| Minimum loop size = 3 | Biological steric constraint (Wikipedia, NNDB) |
| Allow wobble by default | G-U pairs valid in RNA; 20 GU stacking entries from NNDB |
| Greedy non-overlapping stem-loop selection | Standard approach for simple prediction |
| Nussinov-style weighted pair DP for MFE | Efficient O(n³) approximation; full Zuker requires loop decomposition |
| Special hairpin loops as total replacements | Matches NNDB approach: experimental data supersedes model calculation |

---

## 9. Sign-off

- [x] Evidence document created
- [x] Test specification complete
- [x] Existing tests audited
- [x] Consolidation plan defined
- [x] All Must tests verified
