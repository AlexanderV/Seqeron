# RNA-STEMLOOP-001 Test Specification

**Test Unit:** RNA-STEMLOOP-001  
**Title:** Stem-Loop Detection  
**Area:** RnaStructure  
**Created:** 2026-02-05  
**Status:** Complete (re-validated 2026-11-25, B12)

---

## 1. Test Unit Definition

### 1.1 Scope

Testing stem-loop (hairpin) detection algorithms for RNA sequences. This includes finding stem-loop structures, handling various parameters, and pseudoknot detection from base pair lists.

### 1.2 Methods Under Test

| Method | Class | Type | Priority |
|--------|-------|------|----------|
| `FindStemLoops(sequence, minStem, minLoop, maxLoop, allowWobble, maxInternalLoopSize = 0)` | RnaSecondaryStructure | Canonical | Must |
| `CalculateStemEnergy(sequence, basePairs)` | RnaSecondaryStructure | Energy helper (Stem.FreeEnergy) | Must |
| `CalculateHairpinLoopEnergy(loop, c5, c3, specialGU)` | RnaSecondaryStructure | Energy helper (shared with MFE) | Must |
| `DetectPseudoknots(basePairs)` | RnaSecondaryStructure | Structural | Should |

**Note (B12 review 2026-11):** `FindHairpins(sequence, params)` and `FindPseudoknots(sequence)` listed earlier do **not exist** — `FindStemLoops` is the hairpin finder and pseudoknots are detected from base pairs (`DetectPseudoknots`); sequence-level pseudoknot prediction is `PredictStructurePseudoknot` (separate unit). MCP: `find_stem_loops`, `stem_energy`, `hairpin_loop_energy` delegate.

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
| TL-001 | GNRA tetraloop found | `GGGCGAAAGCCC` | Loop.Size = 4, Loop.Sequence = GAAA |
| TL-002 | GNRA pattern verified | GAAA loop | Matches G-N-R-A pattern (Wikipedia Tetraloop) |

#### Turner 2004 energies vs NNDB / ViennaRNA (Evidence §11; tests `RnaSecondaryStructure_StemLoop_Tests.cs`)

| Test ID | Description | Input | Expected (source) |
|---------|-------------|-------|----------|
| SL-E1 | NNDB worked examples (hairpin page §5.3) | `CACAAAAAAAUGUG`, `CACAGAAAGUGUG`, `CACCCGAGGGUG`, `CACACCCCCCUGUG`, `CGGGGAAAGUCCG` | −1.41, −1.91, −4.11, +2.09, −4.80 (NNDB −1.4/−1.9/−4.1/+2.1/−4.8) |
| SL-E2 | A-U-closed special hexaloop: no extra AU end penalty | `GGGACAGUGUUCCC` `((((......))))`; `CGACAGUACUCG` `(((......)))` | −7.07; −1.91 (ViennaRNA eval_structure) |
| SL-E3 | N first mismatch = ViennaRNA N (max over A/C/G/U); DNA T / lower case | `GGGAAANCCC`, `GGGNAAANCCC`, `GGGCNNNNGCCC`; `("AAAN",G,C)`, `("ttcg",c,g)` | −2.02, −1.32, −5.14; 4.50, 3.7 |
| SL-E4 | `CalculateStemEnergy` normalises bases, AU/GU from bases | g·c,g·c; A·T,A·T; G·U labelled WatsonCrick | −3.26; −0.03; = Wobble-labelled value |
| SL-U1 | Opt-in bulge stems = MFE unbranched chain | `GGGACAGGAAACUGCCC`, M = 1 | `(((.(((....))))))` −7.23 (= RNA.fold MFE); M = 0 → `(((....)))` −1.09 |
| SL-U2 | Unbranched optimum per hairpin = ViennaRNA-loop-function DP | `GCAUGCUACGUAGCAUCGAUCAGCUAGCUAGC`, M = 8 | 5 locked rows (e.g. loop 14–21 → −8.86) |
| SL-U3 | Backward compatibility | M = 0 vs legacy call | identical lists |
| SL-U4 | Invalid option | M = −1 | `ArgumentOutOfRangeException` |

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
6. `TotalFreeEnergy == EvaluateStructureEnergy(seq[Start..End], DotBracketNotation)` (ViennaRNA `eval_structure` counterpart of the isolated motif)
7. With `maxInternalLoopSize = M > 0`: the reported stem is the minimum-energy unbranched chain closing that hairpin loop

### 3.2 Boundary Conditions

1. Empty sequence → empty result (no exception)
2. Sequence shorter than minStem*2 + minLoop → empty result
3. No complementary bases → empty result

---

## 4. Audit of Existing Tests

### 4.1 Location

**File:** `tests/Seqeron/Seqeron.Genomics.Tests/RnaSecondaryStructureTests.cs`

### 4.2 Existing Coverage for FindStemLoops

| Test | Spec IDs | Quality | Notes |
|------|----------|---------|-------|
| `FindStemLoops_SimpleHairpin_FindsStructure` | SL-001, SL-002, SL-003, DB-001, DB-002, DB-003 | Exact | Position, stem, loop, dot-bracket assertions |
| `FindStemLoops_NoComplement_ReturnsEmpty` | SL-004 | Exact | |
| `FindStemLoops_TooShort_ReturnsEmpty` | SL-005 | Exact | |
| `FindStemLoops_WithWobblePairs_IncludesWobble` | PH-001 | Exact | Count=2, verifies wobble pair U-G |
| `FindStemLoops_WithoutWobble_ExcludesWobble` | PH-002 | Exact | Count=1, stem=2bp, loop=AAAA, all WC |
| `FindStemLoops_MultipleStemLoops_FindsAll` | MS-001 | Strong | Verifies AAAAA + UUUUU loops present |
| `FindStemLoops_Tetraloop_FindsSpecialLoop` | TL-001, TL-002 | Exact | Count=1, stem=4bp, loop=GAAA, dot-bracket |
| `FindStemLoops_EmptyString_ReturnsEmpty` | EC-001 | Exact | |
| `FindStemLoops_NullString_ReturnsEmpty` | EC-002 | Exact | |
| `FindStemLoops_MinLoopSizeBelowThree_ClampedToThree` | EC-004 | Exact | Uses AAGCAAGCAA (proper clamp test) + positive control |
| `FindStemLoops_MinStemParameter_RespectsMinimum` | PH-003 | Exact | Positive + negative cases |
| `FindStemLoops_LoopSizeRange_RespectsLimits` | PH-004 | Exact | Positive + negative cases |
| `LowerCaseInput_HandlesCorrectly` | EC-003 | Exact | Compares stem, loop, dot-bracket, energy |
| `DetectPseudoknots_NoCrossing_ReturnsEmpty` | PK-001 | Exact | |
| `DetectPseudoknots_CrossingPairs_DetectsKnot` | PK-002 | Exact | |
| `FindInvertedRepeats_RnaHairpin_SmokeTest` | IR-001 | Smoke | |

### 4.3 Assessment

- **Coverage:** Complete — all Must, Should, and edge case tests implemented
- **Quality:** All assertions use exact values derived from RNA biology theory (not fitted to implementation)
- **Consolidation:** DB-001/DB-002/DB-003 merged into SimpleHairpin; MinimumLoopSize_BiologicalConstraint and MinLoopSizeBelowThree merged into one test with proper clamp verification

---

## 5. Consolidation Plan

### 5.1 Strategy

The existing tests in `RnaSecondaryStructureTests.cs` already cover RNA-STEMLOOP-001 methods as part of the RNA-STRUCT-001 Test Unit. The tests are organized in the `#region Stem-Loop Finding Tests` section.

**Decision:** Add missing edge case tests to existing file rather than create new file, as the tests are cohesive and well-organized.

### 5.2 Edge Case Tests

All edge case tests are implemented in the existing file:
- `FindStemLoops_EmptyString_ReturnsEmpty` (EC-001)
- `FindStemLoops_NullString_ReturnsEmpty` (EC-002)
- `LowerCaseInput_HandlesCorrectly` (EC-003) — compares stem, loop, dot-bracket, energy
- `FindStemLoops_MinLoopSizeBelowThree_ClampedToThree` (EC-004) — uses AAGCAAGCAA with positive control
- `FindStemLoops_MinStemParameter_RespectsMinimum` (PH-003)
- `FindStemLoops_LoopSizeRange_RespectsLimits` (PH-004)

### 5.3 Consolidation Actions Taken

| Action | Details |
|--------|---------|
| 🔁 Merged | `FindStemLoops_DotBracket_IsGenerated` → assertions moved to `SimpleHairpin` (exact dot-bracket `(((....))).`) |
| 🔁 Merged | `FindStemLoops_MinimumLoopSize_BiologicalConstraint` + `MinLoopSizeBelowThree_ClampedToThree` → single test with proper sequence `AAGCAAGCAA` |
| ⚠ Strengthened | `SimpleHairpin`: added position, stem coordinates, dot-bracket exact value |
| ⚠ Strengthened | `WithWobblePairs`: count=2, specific wobble U-G pair verification |
| ⚠ Strengthened | `WithoutWobble`: count=1, stem=2bp, loop=AAAA, all-WC assertion |
| ⚠ Strengthened | `MultipleStemLoops`: verifies both AAAAA and UUUUU loops present |
| ⚠ Strengthened | `Tetraloop`: count=1, stem=4bp, dot-bracket `((((....))))` |
| ⚠ Strengthened | `LowerCaseInput`: compares stem length, loop sequence, dot-bracket, energy |

---

## 6. Deviations and Assumptions

- Hairpin loops > 30 nt: NNDB formula 6.4 + 1.75·RT·ln(n/9) unrounded; ViennaRNA truncates to integer dcal/mol (≤ 0.01 kcal/mol).
- NNDB-only terms that ViennaRNA lacks are kept (special GU closure, all-C loops, GGUC/CUGG 3-stack, 1-nt bulge degeneracy / special C).
- Isolated-motif energy: no exterior dangles.

All implementation details are verified against external authoritative sources:

| Aspect | Source | Verification |
|--------|--------|--------------|
| Minimum loop size ≥ 3 nt | Wikipedia Stem-loop: "sterically impossible" | Default `minLoopSize=3` |
| Optimal loop 4–8 nt | Wikipedia Stem-loop: "4-8 bases long" | Default `maxLoopSize=10` is a parameter default, not a biological limit |
| GNRA tetraloop = G-N-R-A | Wikipedia Tetraloop | Test uses GAAA loop, verified pattern |
| GNRA stability via GA mismatch bonus | NNDB Turner 2004 hairpin-special-parameters | No GNRA in special loop table; GA bonus (−0.9) applies |
| Helix = consecutive stacked pairs, AU/GU end 0.45 | NNDB Turner 2004 WC/GU helix pages, hairpin §5.3 | `CalculateStemEnergy`, perfect-stem mode |
| Stem with bulges / internal loops | Zuker & Stiegler 1981 V(i,j) without multiloops; ViennaRNA MAXLOOP | opt-in `maxInternalLoopSize` |
| Energies | ViennaRNA 2.7.2 `eval_structure`, `vrna_E_hairpin` (special loop = total), `mismatchH37` N = max | `EvaluateStructureEnergy`, `HairpinTerm`, `CalculateHairpinLoopEnergy` |
| UNCG most stable tetraloop | Wikipedia Tetraloop, Antao 1991 | NNDB CUUCGG = 3.7 kcal/mol |
| Special loop table (16 tetra + 2 tri + 4 hexa) | NNDB Turner 2004 | All 22 entries verified against source |
| Watson-Crick pairs: A-U, G-C | IUPAC, Wikipedia Base pair | Implemented in `CanPair` |
| Wobble pair: G-U | RNA Biology | Controlled by `allowWobble` parameter |
| Pseudoknot: crossing pairs i < k < j < l | Wikipedia Pseudoknot, Rivas & Eddy 1999 | `DetectPseudoknots` uses exact definition |
| Dot-bracket notation | Wikipedia Nucleic acid secondary structure | Standard `(`, `)`, `.` representation |

## 7. Open Questions

None — the algorithm is well-defined by standard RNA biology.

---

## 8. Decisions

| Decision | Rationale |
|----------|-----------|
| Use existing test file | Tests for FindStemLoops already exist and are well-organized |
| Add edge case tests only | Core functionality already well-tested |
| No separate FindHairpins | Semantically identical to FindStemLoops |
| DetectPseudoknots takes base pairs | Aligns with existing implementation |

---

## 9. Sign-off

- [x] Test specification created
- [x] Existing tests audited
- [x] Consolidation plan defined
- [x] Missing tests identified
- [x] Evidence linked
