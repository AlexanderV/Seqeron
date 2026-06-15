# Test Specification: RNA-PSEUDOKNOT-001

**Test Unit ID:** RNA-PSEUDOKNOT-001
**Area:** RnaStructure
**Algorithm:** Pseudoknot Detection (crossing base pairs)
**Status:** ☐ In Progress
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-14

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Antczak et al. (2018), *Bioinformatics* 34(8):1304–1312 | 1 | https://academic.oup.com/bioinformatics/article/34/8/1304/4721780 | 2026-06-14 |
| 2 | Smit et al. (2008), *RNA* 14(3):410–416 | 1 | https://rnajournal.cshlp.org/content/14/3/410 | 2026-06-14 |
| 3 | biotite.structure.pseudoknots (reference impl.) | 3 | https://www.biotite-python.org/latest/apidoc/biotite.structure.pseudoknots.html | 2026-06-14 |
| 4 | Pseudoknot — Wikipedia (cites Rivas & Eddy 1999) | 4 | https://en.wikipedia.org/wiki/Pseudoknot | 2026-06-14 |

### 1.2 Key Evidence Points

1. Two base pairs (i,j) and (k,l) **cross** (form a pseudoknot) iff **i < k < j < l** — Antczak (2018) [§ crossing/conflict]; biotite; restated symmetrically.
2. **Nested** pairs (i < k < l < j) do not cross — Antczak (2018); Wikipedia "not well nested".
3. **Disjoint** pairs (j < k) do not cross — derived from overlap requirement.
4. A pseudoknot requires ≥ 2 base pairs that cross — Smit (2008) (removal requires crossing pairs).

### 1.3 Documented Corner Cases

- Nested containment is not a pseudoknot (Antczak 2018).
- Disjoint pairs are not a pseudoknot.
- < 2 pairs cannot cross → no pseudoknot (derived).

### 1.4 Known Failure Modes / Pitfalls

1. Treating pair endpoints as ordered (Position1<Position2) without normalizing — would misclassify pairs stored as (close,open). Must min/max before the crossing test — biotite (positions interleave on open/close).

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `DetectPseudoknots(IReadOnlyList<BasePair>)` | RnaSecondaryStructure | Canonical | Returns one `Pseudoknot` per crossing pair-of-pairs. |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | Every reported pseudoknot's two crossing pairs satisfy i < k < j < l (with i<j, k<l after normalization). | Yes | Antczak (2018) crossing condition |
| INV-2 | Nested pairs (i<k<l<j) produce no pseudoknot. | Yes | Antczak (2018); Wikipedia |
| INV-3 | Disjoint pairs (j<k) produce no pseudoknot. | Yes | overlap requirement |
| INV-4 | Fewer than two base pairs → empty result. | Yes | derived from crossing definition |
| INV-5 | Detection is deterministic and order-independent over the same pair set. | Yes | pure combinatorial scan |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | H-type crossing | Pairs (0,2)+(1,3) (`([)]`) | exactly 1 pseudoknot; its CrossingPairs are those 2 pairs | Antczak (2018) i<k<j<l, `([)]` example |
| M2 | Nested control | Pairs (0,5)+(1,4) | 0 pseudoknots | Antczak nested i<k<l<j; Wikipedia |
| M3 | Disjoint control | Pairs (0,2)+(3,5) | 0 pseudoknots | overlap requirement |
| M4 | Endpoint normalization | Pair stored as (2,0)+(3,1) | exactly 1 pseudoknot (same as M1) | biotite open/close min-max |
| M5 | Reported coordinates | (0,2)+(1,3) | result has Start1=0,End1=2,Start2=1,End2=3 | crossing definition |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | Empty input | empty list | 0 pseudoknots | INV-4 |
| S2 | Single pair | one pair | 0 pseudoknots | INV-4 |
| S3 | Mixed set | nested+crossing+disjoint together | only the crossing relations reported | INV-1..3 |
| S4 | Order independence | shuffle input order | identical detected set | INV-5 |
| S5 | Multi-crossing reporting | `([{)]}` pairs (0,3)+(1,4)+(2,5) all mutually cross | exactly 3 pseudoknots, one per pairwise crossing {(0,3,1,4),(0,3,2,5),(1,4,2,5)} | Antczak (2018) crossing condition applied pairwise; documented separate-reporting contract (Assumption 2) |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | Property — every result crosses | random pair sets | every reported pair-of-pairs satisfies i<k<j<l | INV-1, O(n²) property test |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- `DetectPseudoknots` exists in `src/.../RnaSecondaryStructure.cs` (region "Pseudoknot Detection"). No dedicated test file existed before this unit.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1 | ❌ Missing | new |
| M2 | ❌ Missing | new |
| M3 | ❌ Missing | new |
| M4 | ❌ Missing | new |
| M5 | ❌ Missing | new |
| S1 | ❌ Missing | new |
| S2 | ❌ Missing | new |
| S3 | ❌ Missing | new |
| S4 | ❌ Missing | new |
| C1 | ❌ Missing | new |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/RnaSecondaryStructure_DetectPseudoknots_Tests.cs` — all cases above.
- **Remove:** none (no prior pseudoknot tests).

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| `RnaSecondaryStructure_DetectPseudoknots_Tests.cs` | canonical | 12 |

### 5.5 Phase 7 Work Queue

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| 1 | M1 | ❌ Missing | implemented | ✅ Done |
| 2 | M2 | ❌ Missing | implemented | ✅ Done |
| 3 | M3 | ❌ Missing | implemented | ✅ Done |
| 4 | M4 | ❌ Missing | implemented | ✅ Done |
| 5 | M5 | ❌ Missing | implemented | ✅ Done |
| 6 | S1 | ❌ Missing | implemented | ✅ Done |
| 7 | S2 | ❌ Missing | implemented | ✅ Done |
| 8 | S3 | ❌ Missing | implemented | ✅ Done |
| 9 | S4 | ❌ Missing | implemented | ✅ Done |
| 10 | C1 | ❌ Missing | implemented | ✅ Done |

**Total items:** 10
**✅ Done:** 10 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ | implemented |
| M2 | ✅ | implemented |
| M3 | ✅ | implemented |
| M4 | ✅ | implemented |
| M5 | ✅ | implemented |
| S1 | ✅ | implemented |
| S2 | ✅ | implemented |
| S3 | ✅ | implemented |
| S4 | ✅ | implemented |
| C1 | ✅ | implemented |

---

## 6. Assumption Register

**Total assumptions:** 2

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | Fewer than two pairs → empty result (derived from crossing definition). | INV-4, S1, S2 |
| 2 | Each crossing pair-of-pairs reported as one Pseudoknot (DBL order grouping Not Implemented). | M1, M5, output shape |

---

## 7. Open Questions / Decisions

1. Higher-order pseudoknot grouping / DBL order assignment (Antczak 2018) is intentionally out of scope for this unit; documented as Not Implemented.
