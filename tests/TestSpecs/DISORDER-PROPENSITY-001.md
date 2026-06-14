# Test Specification: DISORDER-PROPENSITY-001

**Test Unit ID:** DISORDER-PROPENSITY-001
**Area:** ProteinPred
**Algorithm:** Disorder Propensity (TOP-IDP scale lookup + Dunker order/disorder amino-acid classification)
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-14

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Campen et al. (2008) TOP-IDP-Scale, Protein Pept Lett 15(9):956-963 | 1 | https://pmc.ncbi.nlm.nih.gov/articles/PMC2676888/ | 2026-06-14 |
| 2 | Dunker et al. (2001) Intrinsically disordered protein, J Mol Graph Model 19(1):26-59 | 1 | https://pubmed.ncbi.nlm.nih.gov/11381529/ | 2026-06-14 |
| 3 | Wikipedia — Intrinsically disordered proteins (cites Dunker 2001) | 4 | https://en.wikipedia.org/wiki/Intrinsically_disordered_proteins | 2026-06-14 |

### 1.2 Key Evidence Points

1. TOP-IDP per-residue values (all 20 standard residues) extracted verbatim from Campen et al. (2008) Table 2 — source 1.
2. W = -0.884 is the global minimum; P = +0.987 is the global maximum of the TOP-IDP scale — source 1.
3. Disorder-promoting residues = {A, R, G, Q, S, P, E, K}; order-promoting = {W, C, F, I, Y, V, L, N}; ambiguous = {H, M, T, D} — sources 2, 3.

### 1.3 Documented Corner Cases

- Scale defined only for the 20 standard amino acids (Campen 2008 Table 2); no value for B/J/O/U/X/Z or gaps.
- Implementation contract: unknown residue → 0.0 (`GetValueOrDefault(..., 0)`); input is upper-cased before lookup.

### 1.4 Known Failure Modes / Pitfalls

1. Rendered ranking strings in Campen 2008 / Wikipedia place "...Q, K, S, E, P" while the Table 2 *values* give S < K. Use the numeric values, not the rank-string order — source 1.

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `GetDisorderPropensity(char)` | DisorderPredictor | **Canonical** | Returns TOP-IDP value — Campen et al. (2008) Table 2 |
| `IsDisorderPromoting(char)` | DisorderPredictor | **Canonical** | True for {A,R,G,Q,S,P,E,K} — Dunker (2001) |
| `DisorderPromotingAminoAcids` | DisorderPredictor | **Canonical** | 8 AA {A,E,G,K,P,Q,R,S} — Dunker (2001) |
| `OrderPromotingAminoAcids` | DisorderPredictor | **Canonical** | 8 AA {C,F,I,L,N,V,W,Y} — Dunker (2001) |
| `AmbiguousAminoAcids` | DisorderPredictor | **Canonical** | 4 AA {D,H,M,T} — Dunker (2001); completes the classification trio |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | `GetDisorderPropensity` returns the exact Campen 2008 Table 2 value for each of the 20 standard residues | Yes | Campen et al. (2008) Table 2 |
| INV-2 | For all 20 standard residues, -0.884 ≤ propensity ≤ 0.987, with min at W and max at P | Yes | Campen et al. (2008) Table 2 |
| INV-3 | `IsDisorderPromoting(c)` ⇔ c ∈ `DisorderPromotingAminoAcids` (set membership) | Yes | Dunker et al. (2001) |
| INV-4 | Disorder-promoting (8), order-promoting (8), ambiguous (4) sets are pairwise disjoint and cover all 20 standard residues | Yes | Dunker et al. (2001) |
| INV-5 | `GetDisorderPropensity` is case-insensitive (upper-cases input) | Yes | Implementation contract |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | GetDisorderPropensity_AllTwentyAminoAcids_MatchScale | All 20 residues return exact Table 2 value | exact values per dataset | Campen et al. (2008) Table 2 |
| M2 | GetDisorderPropensity_ExtremeAnchors_MinAtTrpMaxAtPro | W = -0.884 (min), P = 0.987 (max) | -0.884 / 0.987 | Campen et al. (2008) Table 2 |
| M3 | IsDisorderPromoting_DisorderPromotingResidues_ReturnsTrue | A,R,G,Q,S,P,E,K → true | true for each | Dunker et al. (2001) |
| M4 | IsDisorderPromoting_OrderPromotingResidues_ReturnsFalse | W,C,F,I,Y,V,L,N → false | false for each | Dunker et al. (2001) |
| M5 | IsDisorderPromoting_AmbiguousResidues_ReturnsFalse | H,M,T,D → false | false for each | Dunker et al. (2001) |
| M6 | DisorderPromotingAminoAcids_EqualsDunkerSet | property = {A,E,G,K,P,Q,R,S}, count 8 | exact set, count 8 | Dunker et al. (2001) |
| M7 | OrderPromotingAminoAcids_EqualsDunkerSet | property = {C,F,I,L,N,V,W,Y}, count 8 | exact set, count 8 | Dunker et al. (2001) |
| M8 | AmbiguousAminoAcids_EqualsDunkerSet | property = {D,H,M,T}, count 4 | exact set, count 4 | Dunker et al. (2001) |
| M9 | ClassificationSets_AreDisjointAndCoverAll20 | 8+8+4 disjoint, union = 20 standard residues | disjoint, union 20 | Dunker et al. (2001) (INV-4) |
| M10 | IsDisorderPromoting_MatchesProperty_AllResidues | INV-3: membership ⇔ predicate for all 20 | consistent for each | Dunker et al. (2001) |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | GetDisorderPropensity_UnknownResidue_ReturnsZero | X, Z, B, '*' → 0.0 | 0.0 | Implementation contract (assumption) |
| S2 | GetDisorderPropensity_LowercaseInput_SameAsUppercase | 'p'/'w'/'e' equal 'P'/'W'/'E' | equal | INV-5 case-insensitivity |
| S3 | IsDisorderPromoting_LowercaseInput_SameAsUppercase | 'p'→true, 'w'→false | true/false | Case-insensitivity |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | Properties_ReturnSorted_StableOrder | the three property lists are returned ascending-sorted | sorted ascending | API stability (cached OrderBy) |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- `tests/Seqeron/Seqeron.Genomics.Tests/DisorderPredictor_DisorderPrediction_Tests.cs` (unit DISORDER-PRED-001) over-claimed the four propensity/classification methods as Canonical and tested them in regions M8, M9, M10, M10b, M11, M12, S6, S7, C5. The Method Index in `ALGORITHMS_CHECKLIST_V2.md` (lines 4814-4817) assigns `GetDisorderPropensity`, `IsDisorderPromoting`, `DisorderPromotingAminoAcids`, `OrderPromotingAminoAcids` to DISORDER-PROPENSITY-001, so those tests belong here.
- No file existed at the canonical path for this unit.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1 (all 20 values) | 🔁 Duplicate | Existed in PRED-001 file (M8 region); migrate to this unit's canonical file |
| M2 (anchors) | ❌ Missing | Not separately asserted |
| M3 (disorder true) | 🔁 Duplicate | PRED-001 M9; migrate |
| M4 (order false) | 🔁 Duplicate | PRED-001 M10; migrate |
| M5 (ambiguous false) | 🔁 Duplicate | PRED-001 M10b; migrate |
| M6 (disorder set) | 🔁 Duplicate | PRED-001 M11; migrate |
| M7 (order set) | 🔁 Duplicate | PRED-001 M12; migrate |
| M8 (ambiguous set) | ❌ Missing | only counts checked via C5 |
| M9 (disjoint/cover) | 🔁 Duplicate | PRED-001 C5; migrate |
| M10 (membership ⇔ predicate) | ❌ Missing | INV-3 not directly tested |
| S1 (unknown → 0) | 🔁 Duplicate | PRED-001 S6; migrate |
| S2 (lowercase prop) | 🔁 Duplicate | PRED-001 S7; migrate |
| S3 (lowercase predicate) | ❌ Missing | Not tested |
| C1 (sorted) | ❌ Missing | Not tested |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/DisorderPredictor_GetDisorderPropensity_Tests.cs` — owns all tests for the four scope methods + `AmbiguousAminoAcids`.
- **Remove:** from `DisorderPredictor_DisorderPrediction_Tests.cs` the regions that test the four scope methods directly (M8, M9, M10, M10b, M11, M12, S6, S7, C5), to eliminate duplicates. The M8b region (tests `PredictDisorder` normalized scores) stays in PRED-001 because it tests the PRED-001 canonical method.

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| `DisorderPredictor_GetDisorderPropensity_Tests.cs` | Canonical (this unit) | 14 |
| `DisorderPredictor_DisorderPrediction_Tests.cs` | PRED-001 (scope methods removed) | reduced |

### 5.5 Phase 7 Work Queue

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| 1 | M1 | 🔁 Duplicate | Migrated to canonical file (exact 20 values) | ✅ Done |
| 2 | M2 | ❌ Missing | Added anchor min/max test | ✅ Done |
| 3 | M3 | 🔁 Duplicate | Migrated | ✅ Done |
| 4 | M4 | 🔁 Duplicate | Migrated | ✅ Done |
| 5 | M5 | 🔁 Duplicate | Migrated | ✅ Done |
| 6 | M6 | 🔁 Duplicate | Migrated (exact set + count) | ✅ Done |
| 7 | M7 | 🔁 Duplicate | Migrated (exact set + count) | ✅ Done |
| 8 | M8 | ❌ Missing | Added ambiguous-set exact test | ✅ Done |
| 9 | M9 | 🔁 Duplicate | Migrated (disjoint + cover 20) | ✅ Done |
| 10 | M10 | ❌ Missing | Added membership⇔predicate property test | ✅ Done |
| 11 | S1 | 🔁 Duplicate | Migrated | ✅ Done |
| 12 | S2 | 🔁 Duplicate | Migrated | ✅ Done |
| 13 | S3 | ❌ Missing | Added lowercase-predicate test | ✅ Done |
| 14 | C1 | ❌ Missing | Added sorted-order test | ✅ Done |
| 15 | (cleanup) | 🔁 Duplicate | Removed duplicate regions from PRED-001 file | ✅ Done |

**Total items:** 15
**✅ Done:** 15 | **⛔ Blocked:** 0 | **Remaining:** must be 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ Covered | `GetDisorderPropensity_AllTwentyAminoAcids_MatchScale` |
| M2 | ✅ Covered | `GetDisorderPropensity_ExtremeAnchors_MinAtTrpMaxAtPro` |
| M3 | ✅ Covered | `IsDisorderPromoting_DisorderPromotingResidues_ReturnsTrue` |
| M4 | ✅ Covered | `IsDisorderPromoting_OrderPromotingResidues_ReturnsFalse` |
| M5 | ✅ Covered | `IsDisorderPromoting_AmbiguousResidues_ReturnsFalse` |
| M6 | ✅ Covered | `DisorderPromotingAminoAcids_EqualsDunkerSet` |
| M7 | ✅ Covered | `OrderPromotingAminoAcids_EqualsDunkerSet` |
| M8 | ✅ Covered | `AmbiguousAminoAcids_EqualsDunkerSet` |
| M9 | ✅ Covered | `ClassificationSets_AreDisjointAndCoverAll20` |
| M10 | ✅ Covered | `IsDisorderPromoting_MatchesProperty_AllStandardResidues` |
| S1 | ✅ Covered | `GetDisorderPropensity_UnknownResidue_ReturnsZero` |
| S2 | ✅ Covered | `GetDisorderPropensity_LowercaseInput_SameAsUppercase` |
| S3 | ✅ Covered | `IsDisorderPromoting_LowercaseInput_SameAsUppercase` |
| C1 | ✅ Covered | `Properties_ReturnSortedAscending` |

---

## 6. Assumption Register

**Total assumptions:** 2

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | Unknown residue → 0.0 (implementation `GetValueOrDefault`, not source-defined) | S1 |
| 2 | S/K rank-string vs Table 2 value discrepancy; numeric Table 2 values are authoritative | M1, INV-1 |

---

## 7. Open Questions / Decisions

1. `AmbiguousAminoAcids` is not assigned to any unit in the Method Index; included here because it completes the Dunker (2001) classification trio that `IsDisorderPromoting` partitions and has no other home.
2. Decision: duplicate tests for the four scope methods were removed from the DISORDER-PRED-001 fixture (which over-claimed them) so one canonical owner remains, per the duplicate-elimination rule.
