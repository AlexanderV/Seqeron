# Test Specification: RNA-PAIR-001

**Test Unit ID:** RNA-PAIR-001
**Area:** RnaStructure
**Algorithm:** RNA Base Pairing (CanPair / GetBasePairType / GetComplement)
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-14

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Crick FHC (1966) Codon–anticodon pairing: the wobble hypothesis, J Mol Biol 19:548–555 | 1 | https://doi.org/10.1016/S0022-2836(66)80022-0 | 2026-06-14 |
| 2 | Wikipedia — Base pair (canonical Watson-Crick A•U, G•C) | 4 | https://en.wikipedia.org/wiki/Base_pair | 2026-06-14 |
| 3 | Wikipedia — Wobble base pair (G–U wobble; cites Crick 1966) | 4 | https://en.wikipedia.org/wiki/Wobble_base_pair | 2026-06-14 |
| 4 | IUPAC-IUB (1970) Biochemistry 9(20):4022–4027, via Nucleic acid notation | 2 | https://en.wikipedia.org/wiki/Nucleic_acid_notation | 2026-06-14 |
| 5 | Biopython — Bio.Seq.complement_rna | 3 | https://biopython.org/docs/latest/api/Bio.Seq.html | 2026-06-14 |

### 1.2 Key Evidence Points

1. Canonical Watson-Crick RNA pairs are A•U (2 H-bonds) and G•C (3 H-bonds) — Source 2.
2. G–U is the standard RNA wobble pair, distinct from Watson-Crick; over the standard alphabet G pairs with C and U, U pairs with A and G — Sources 1, 3.
3. Base pairing is reciprocal/symmetric: A•U ≡ U•A, G•C ≡ C•G — Source 2.
4. RNA complement: A→U, U→A, G→C, C→G, T→A (T treated as U), N→N, with IUPAC degenerate complements — Sources 4, 5.

### 1.3 Documented Corner Cases

- Order independence: `CanPair`/`GetBasePairType` are symmetric in their two arguments (Source 2).
- G–U must be reported as Wobble, never WatsonCrick (Source 3).
- `CanPair`/`GetBasePairType` are defined over the RNA alphabet {A,C,G,U}; T does not pair. For `GetComplement`, T is treated as U → A (Sources 4, 5).

### 1.4 Known Failure Modes / Pitfalls

1. Misclassifying G–U as Watson-Crick — Source 3 (wobble does not follow WC rules).
2. Returning a complement of A as T instead of U in RNA context — Source 5 (complement_rna).
3. Treating the DNA base T as a pairing partner in RNA pairing — not defined by Sources 1–3; `CanPair` returns false for T.

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `CanPair(char, char)` | RnaSecondaryStructure | Canonical | Boolean pairing test; O(1) table lookup |
| `GetBasePairType(char, char)` | RnaSecondaryStructure | Canonical | Returns WatsonCrick / Wobble / null |
| `GetComplement(char)` | RnaSecondaryStructure | Delegate | Delegates to `SequenceExtensions.GetRnaComplementBase` |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | `CanPair(x,y) == CanPair(y,x)` (symmetry) | Yes | Source 2 (reciprocal pairing) |
| INV-2 | `GetBasePairType(x,y) == GetBasePairType(y,x)` (symmetry) | Yes | Source 2 |
| INV-3 | `CanPair(x,y) == (GetBasePairType(x,y) != null)` (consistency) | Yes | Definition: a pair exists iff it has a type |
| INV-4 | `GetBasePairType('G','U')` and `('U','G')` == Wobble (never WatsonCrick) | Yes | Sources 1, 3 |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | CanPair Watson-Crick | A-U, U-A, G-C, C-G | all true | Source 2 |
| M2 | CanPair wobble | G-U, U-G | all true | Sources 1, 3 |
| M3 | CanPair non-pairs | A-A, A-G, A-C, C-U, G-G, C-C | all false | Source 1 (A only U, C only G) |
| M4 | GetBasePairType Watson-Crick | A-U, U-A, G-C, C-G | WatsonCrick | Source 2 |
| M5 | GetBasePairType wobble | G-U, U-G | Wobble | Sources 1, 3 |
| M6 | GetBasePairType non-pairs | A-A, A-G, C-U | null | Source 1 |
| M7 | GetComplement standard bases | A,U,G,C,T | U,A,C,G,A | Sources 4, 5 |
| M8 | Symmetry (INV-1, INV-2) | swap arguments for all pairs | identical results | Source 2 |
| M9 | Consistency (INV-3) | CanPair true ⇔ type non-null | holds for all 16 combinations | definition |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | Case-insensitivity | lowercase a,u,g,c,t | same as uppercase | normalization contract |
| S2 | DNA T in CanPair | T-A, A-T, G-T | false (T not an RNA base) | Sources 1, 2 |
| S3 | GetComplement IUPAC degenerate | N→N, R→Y, Y→R | per IUPAC | Source 4 |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | Out-of-range char | non-ASCII / digit input | CanPair false, type null, no exception | robustness |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- `tests/Seqeron/Seqeron.Genomics.Tests/RnaSecondaryStructureTests.cs` contains pre-existing tests for `CanPair`, `GetBasePairType`, and `GetComplement` (lines ~48–87). These are kept as part of the broader RnaSecondaryStructure fixture but are not the canonical RNA-PAIR-001 file.
- No `RnaSecondaryStructure_CanPair_Tests.cs` existed before this unit.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1 CanPair WC | ❌ Missing | New canonical file does not exist yet |
| M2 CanPair wobble | ❌ Missing | |
| M3 CanPair non-pairs | ❌ Missing | |
| M4 GetBasePairType WC | ❌ Missing | |
| M5 GetBasePairType wobble | ❌ Missing | |
| M6 GetBasePairType non-pairs | ❌ Missing | |
| M7 GetComplement | ❌ Missing | |
| M8 Symmetry | ❌ Missing | not covered by old fixture |
| M9 Consistency | ❌ Missing | not covered by old fixture |
| S1 Case-insensitivity | ❌ Missing | |
| S2 DNA T | ❌ Missing | |
| S3 IUPAC degenerate | ❌ Missing | |
| C1 Out-of-range | ❌ Missing | |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/RnaSecondaryStructure_CanPair_Tests.cs` — all RNA-PAIR-001 cases (M1–M9, S1–S3, C1).
- **Remove:** nothing. Pre-existing `RnaSecondaryStructureTests.cs` covers many other methods (energy, MFE, dot-bracket) and remains the home for those; its few CanPair/GetComplement tests are superseded by the new canonical file but left in place to avoid scope creep into unrelated tests. No duplicate assertions are introduced that weaken evidence.

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| `RnaSecondaryStructure_CanPair_Tests.cs` | Canonical RNA-PAIR-001 | 13 |
| `RnaSecondaryStructureTests.cs` | Pre-existing broad fixture (other methods) | unchanged |

### 5.5 Phase 7 Work Queue

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| 1 | M1 | ❌ Missing | Implemented | ✅ Done |
| 2 | M2 | ❌ Missing | Implemented | ✅ Done |
| 3 | M3 | ❌ Missing | Implemented | ✅ Done |
| 4 | M4 | ❌ Missing | Implemented | ✅ Done |
| 5 | M5 | ❌ Missing | Implemented | ✅ Done |
| 6 | M6 | ❌ Missing | Implemented | ✅ Done |
| 7 | M7 | ❌ Missing | Implemented | ✅ Done |
| 8 | M8 | ❌ Missing | Implemented (property test) | ✅ Done |
| 9 | M9 | ❌ Missing | Implemented (property test) | ✅ Done |
| 10 | S1 | ❌ Missing | Implemented | ✅ Done |
| 11 | S2 | ❌ Missing | Implemented | ✅ Done |
| 12 | S3 | ❌ Missing | Implemented | ✅ Done |
| 13 | C1 | ❌ Missing | Implemented | ✅ Done |

**Total items:** 13
**✅ Done:** 13 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ Covered | CanPair_WatsonCrickPairs_ReturnsTrue |
| M2 | ✅ Covered | CanPair_WobblePairs_ReturnsTrue |
| M3 | ✅ Covered | CanPair_NonPairs_ReturnsFalse |
| M4 | ✅ Covered | GetBasePairType_WatsonCrick_ReturnsWatsonCrick |
| M5 | ✅ Covered | GetBasePairType_Wobble_ReturnsWobble |
| M6 | ✅ Covered | GetBasePairType_NonPairs_ReturnsNull |
| M7 | ✅ Covered | GetComplement_StandardBases_ReturnsRnaComplement |
| M8 | ✅ Covered | CanPair_And_Type_AreSymmetric (property) |
| M9 | ✅ Covered | CanPair_AgreesWith_GetBasePairType (property) |
| S1 | ✅ Covered | CanPair_LowercaseInput_SameAsUppercase |
| S2 | ✅ Covered | CanPair_DnaT_NotAnRnaBase_ReturnsFalse |
| S3 | ✅ Covered | GetComplement_IupacDegenerate_ReturnsExpected |
| C1 | ✅ Covered | CanPair_OutOfRangeChar_ReturnsFalse |

---

## 6. Assumption Register

**Total assumptions:** 1

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | Case-insensitive input is a non-correctness-affecting normalization (lower/upper case denote the same nucleotide) | S1 |

---

## 7. Open Questions / Decisions

1. None — all methods conform to retrieved authoritative sources; no correctness-affecting assumptions remain.
