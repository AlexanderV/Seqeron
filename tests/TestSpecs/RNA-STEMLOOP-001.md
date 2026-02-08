# RNA-STEMLOOP-001 Test Specification

**Test Unit:** RNA-STEMLOOP-001  
**Title:** Stem-Loop Detection  
**Area:** RnaStructure  
**Created:** 2026-02-05  
**Status:** Complete

---

## 1. Test Unit Definition

### 1.1 Scope

Testing stem-loop (hairpin) detection algorithms for RNA sequences. This includes finding stem-loop structures, handling various parameters, and pseudoknot detection from base pair lists.

### 1.2 Methods Under Test

| Method | Class | Type | Priority |
|--------|-------|------|----------|
| `FindStemLoops(sequence, minStem, minLoop, maxLoop, allowWobble)` | RnaSecondaryStructure | Canonical | Must |
| `FindHairpins(sequence, params)` | RnaSecondaryStructure | Variant | Should |
| `FindPseudoknots(sequence)` | RnaSecondaryStructure | Structural | Should |

**Note:** `FindStemLoops` is the canonical method. `FindHairpins` is semantically equivalent (hairpin = stem-loop). `FindPseudoknots(sequence)` extends existing `DetectPseudoknots(basePairs)`.

---

## 2. Test Categories

### 2.1 Must Tests (Required for Complete)

#### Basic Stem-Loop Detection (Evidence: Wikipedia Stem-loop)

| Test ID | Description | Input | Expected |
|---------|-------------|-------|----------|
| SL-001 | Simple hairpin detection | `GGGAAAACCC` | Finds ≥1 stem-loop |
| SL-002 | Stem length validation | `GGGAAAACCC` | Stem.Length ≥ 3 |
| SL-003 | Loop type identification | Any hairpin | Loop.Type = Hairpin |
| SL-004 | No complement returns empty | `AAAAAAAAA` | Empty result |
| SL-005 | Too short returns empty | `GCAUC` | Empty result |

#### Parameter Handling (Evidence: Wikipedia, Biology)

| Test ID | Description | Input/Params | Expected |
|---------|-------------|--------------|----------|
| PH-001 | Wobble pairs included | allowWobble=true | G-U pairs in stem |
| PH-002 | Wobble pairs excluded | allowWobble=false | No G-U pairs |
| PH-003 | Minimum stem length | minStem=4 | Only ≥4bp stems |
| PH-004 | Loop size range | minLoop=4, maxLoop=4 | Only 4nt loops |

#### Edge Cases (Evidence: Implementation)

| Test ID | Description | Input | Expected |
|---------|-------------|-------|----------|
| EC-001 | Empty string | `""` | Empty result |
| EC-002 | Null sequence | `null` | Empty result |
| EC-003 | Lowercase input | `gggaaaaccc` | Same as uppercase |
| EC-004 | Minimum loop constraint | <3 nt loop | No such loop found |

#### Dot-Bracket Notation (Evidence: Wikipedia)

| Test ID | Description | Input | Expected |
|---------|-------------|-------|----------|
| DB-001 | DotBracket is generated | Any hairpin | Non-empty notation |
| DB-002 | Contains opening brackets | Any hairpin | Has '(' |
| DB-003 | Contains closing brackets | Any hairpin | Has ')' |

#### Tetraloop Detection (Evidence: Wikipedia Tetraloop, Woese 1990)

| Test ID | Description | Input | Expected |
|---------|-------------|-------|----------|
| TL-001 | Finds 4nt loop | `GGGGCGAACCCC` | Loop.Size = 4 |
| TL-002 | GNRA pattern recognized | GAAA loop | Structure found |

### 2.2 Should Tests (Recommended)

#### Multiple Structures (Evidence: RNA biology)

| Test ID | Description | Input | Expected |
|---------|-------------|-------|----------|
| MS-001 | Multiple hairpins | Sequence with 2+ | Finds ≥1 |

#### Pseudoknot Detection (Evidence: Wikipedia Pseudoknot, Rivas 1999)

| Test ID | Description | Input | Expected |
|---------|-------------|-------|----------|
| PK-001 | Non-crossing pairs | Nested pairs | No pseudoknot |
| PK-002 | Crossing pairs detected | i<k<j<l pattern | Pseudoknot found |

#### Inverted Repeat Support (Evidence: Wikipedia)

| Test ID | Description | Input | Expected |
|---------|-------------|-------|----------|
| IR-001 | RNA hairpin stems | Complementary regions | Finds repeat |

### 2.3 Could Tests (Optional)

| Test ID | Description |
|---------|-------------|
| CT-001 | tRNA-like structure analysis |
| CT-002 | Energy comparison for different loop sizes |
| CT-003 | Performance on long sequences |

---

## 3. Test Invariants

### 3.1 Always True

1. `stemLoop.Loop.Type == LoopType.Hairpin` for all found hairpins
2. `stemLoop.Stem.Length >= minStemLength` parameter
3. `stemLoop.Loop.Size >= minLoopSize` and `<= maxLoopSize`
4. `stemLoop.DotBracketNotation.Length > 0` for valid structures
5. All base pairs in stem are valid (WatsonCrick or Wobble if allowed)

### 3.2 Boundary Conditions

1. Empty sequence → empty result (no exception)
2. Sequence shorter than minStem*2 + minLoop → empty result
3. No complementary bases → empty result

---

## 4. Audit of Existing Tests

### 4.1 Location

**File:** `tests/Seqeron/Seqeron.Genomics.Tests/RnaSecondaryStructureTests.cs`

### 4.2 Existing Coverage for FindStemLoops

| Test | Method | Coverage | Status |
|------|--------|----------|--------|
| `FindStemLoops_SimpleHairpin_FindsStructure` | SL-001, SL-002, SL-003 | ✓ Covered |
| `FindStemLoops_NoComplement_ReturnsEmpty` | SL-004 | ✓ Covered |
| `FindStemLoops_TooShort_ReturnsEmpty` | SL-005 | ✓ Covered |
| `FindStemLoops_WithWobblePairs_IncludesWobble` | PH-001 | ✓ Covered |
| `FindStemLoops_WithoutWobble_ExcludesWobble` | PH-002 | ✓ Covered |
| `FindStemLoops_MultipleStemLoops_FindsAll` | MS-001 | ✓ Covered |
| `FindStemLoops_Tetraloop_FindsSpecialLoop` | TL-001, TL-002 | ✓ Covered |
| `FindStemLoops_DotBracket_IsGenerated` | DB-001, DB-002, DB-003 | ✓ Covered |
| `DetectPseudoknots_NoCrossing_ReturnsEmpty` | PK-001 | ✓ Covered |
| `DetectPseudoknots_CrossingPairs_DetectsKnot` | PK-002 | ✓ Covered |
| `FindInvertedRepeats_RnaHairpin_SmokeTest` | IR-001 | ✓ Covered |

### 4.3 Missing Tests (All Closed)

| Test ID | Description | Status |
|---------|-------------|--------|
| EC-001 | Empty string handling | ✅ Covered |
| EC-002 | Null handling | ✅ Covered |
| EC-003 | Lowercase input | ✅ Covered |
| PH-003 | MinStem parameter effect | ✅ Covered |
| PH-004 | Loop size range | ✅ Covered |

### 4.4 Assessment

- **Existing Coverage:** Strong for core functionality
- **Gaps:** Parameter validation edge cases
- **Quality:** Tests are well-structured with evidence comments

---

## 5. Consolidation Plan

### 5.1 Strategy

The existing tests in `RnaSecondaryStructureTests.cs` already cover RNA-STEMLOOP-001 methods as part of the RNA-STRUCT-001 Test Unit. The tests are organized in the `#region Stem-Loop Finding Tests` section.

**Decision:** Add missing edge case tests to existing file rather than create new file, as the tests are cohesive and well-organized.

### 5.2 Tests to Add (All Added)

1. ~~`FindStemLoops_EmptyString_ReturnsEmpty` - EC-001~~ ✅ Done
2. ~~`FindStemLoops_NullString_ReturnsEmpty` - EC-002~~ ✅ Done
3. ~~`FindStemLoops_LowercaseInput_HandledCorrectly` - EC-003~~ ✅ Done
4. ~~`FindStemLoops_MinStemParameter_RespectsMinimum` - PH-003~~ ✅ Done
5. ~~`FindStemLoops_LoopSizeRange_RespectsLimits` - PH-004~~ ✅ Done

---

## 6. Open Questions

None - the algorithm is well-defined by standard RNA biology.

---

## 7. Decisions

| Decision | Rationale |
|----------|-----------|
| Use existing test file | Tests for FindStemLoops already exist and are well-organized |
| Add edge case tests only | Core functionality already well-tested |
| No separate FindHairpins | Semantically identical to FindStemLoops |
| DetectPseudoknots takes base pairs | Aligns with existing implementation |

---

## 8. Sign-off

- [x] Test specification created
- [x] Existing tests audited
- [x] Consolidation plan defined
- [x] Missing tests identified
- [x] Evidence linked
