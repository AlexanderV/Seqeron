# Test Specification: RESTR-FILTER-001

**Test Unit ID:** RESTR-FILTER-001
**Area:** MolTools
**Algorithm:** Restriction Enzyme Filtering (by recognition-site length, blunt cutters, sticky cutters)
**Status:** ☐ In Progress
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-13

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Sticky and blunt ends (Wikipedia) | 4 | https://en.wikipedia.org/wiki/Sticky_and_blunt_ends | 2026-06-13 |
| 2 | Restriction enzyme (Wikipedia) | 4 | https://en.wikipedia.org/wiki/Restriction_enzyme | 2026-06-13 |
| 3 | List of restriction enzyme cutting sites (Wikipedia) | 4 | https://en.wikipedia.org/wiki/List_of_restriction_enzyme_cutting_sites | 2026-06-13 |
| 4 | KpnI (NEB R0142 / REBASE) | 3 | https://www.neb.com/en/products/r0142-kpni | 2026-06-13 |
| 5 | EcoRI-HF (NEB R3101 / REBASE) | 3 | https://www.neb.com/en/products/r3101-ecori-hf | 2026-06-13 |
| 6 | SfiI interrupted palindrome (PMC548270 / REBASE) | 1 | https://www.ncbi.nlm.nih.gov/pmc/articles/PMC548270/ | 2026-06-13 |

### 1.2 Key Evidence Points

1. A blunt end means "both strands terminate in a base pair"; a sticky/cohesive end is "a stretch of unpaired nucleotides" (a 5' or 3' overhang). Every end is one or the other. — Source 1.
2. Type II enzymes "can either cleave at the center of both strands to yield a blunt end, or at a staggered position leaving overhangs called sticky ends." — Source 2.
3. Type II recognition sites are "4–8 nucleotides in length." — Source 2.
4. SmaI (CCCGGG) is blunt; EcoRI (GAATTC) is a 5' overhang (sticky); KpnI (GGTACC) and PstI (CTGCAG) are 3' overhangs (sticky). — Sources 2, 4, 5.
5. Length categories: EcoRI/BamHI/HindIII/PstI = 6 bp; AluI/HaeIII/TaqI = 4 bp; NotI = 8 bp. — Source 3.
6. SfiI is an interrupted palindrome `5'-GGCCNNNN^NGGCC-3'` (13-nt recognition string), a divided site outside the 4–8 nt undivided range. — Source 6.

### 1.3 Documented Corner Cases

- A length range excluding all of 4–8 nt returns no enzymes; a range covering 4–8 returns the whole library (Source 2, 3).
- Blunt and sticky sets are complementary and disjoint over the library — the partition is total (Source 1).

### 1.4 Known Failure Modes / Pitfalls

1. Treating a palindromic-but-staggered cutter (KpnI, PstI) as blunt because the overhang string is symmetric — geometry, not the string, defines blunt vs sticky. — Source 2.

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `GetEnzymesByCutLength(int minLength, int maxLength)` | RestrictionAnalyzer | **Canonical** | New range overload; inclusive bounds |
| `GetBluntCutters()` | RestrictionAnalyzer | **Canonical** | Enzymes with blunt ends |
| `GetStickyCutters()` | RestrictionAnalyzer | **Canonical** | Enzymes with 5'/3' overhangs |
| `GetEnzymesByCutLength(int length)` | RestrictionAnalyzer | **Delegate** | Existing single-length overload; smoke + consistency with range overload |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | Blunt-cutter set and sticky-cutter set are disjoint and their union equals the full enzyme library (total partition). | Yes | Source 1 (every end is blunt or an overhang) |
| INV-2 | `GetBluntCutters()` returns exactly enzymes whose forward and reverse cut positions are equal (center cut → blunt). | Yes | Source 2 (center cut → blunt) |
| INV-3 | `GetEnzymesByCutLength(min,max)` returns exactly enzymes with `min ≤ RecognitionLength ≤ max` (inclusive). | Yes | Source 2/3 (recognition lengths) + ASSUMPTION (bound inclusivity) |
| INV-4 | `GetEnzymesByCutLength(L)` equals `GetEnzymesByCutLength(L,L)` for any L (overload consistency). | Yes | Derivation |
| INV-5 | When `min > max`, the range filter returns an empty set. | Yes | Empty-interval logic |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | BluntCutters_IncludesKnownBlunt_ExcludesSticky | Blunt set contains SmaI, EcoRV, AluI, HaeIII; excludes EcoRI, KpnI | SmaI/EcoRV/AluI/HaeIII ∈ result; EcoRI/KpnI ∉ result | Source 1, 2 (SmaI blunt) |
| M2 | StickyCutters_IncludesKnownSticky_ExcludesBlunt | Sticky set contains EcoRI, KpnI, PstI; excludes SmaI, EcoRV | EcoRI/KpnI/PstI ∈ result; SmaI/EcoRV ∉ result | Source 2 (EcoRI 5', KpnI/PstI 3'), Source 4 |
| M3 | BluntAndSticky_PartitionLibrary | Blunt ∪ Sticky = all enzymes; Blunt ∩ Sticky = ∅; counts sum to total | Disjoint, union == full library | Source 1 (total partition) |
| M4 | EnzymesByCutLength_Range6_ReturnsAll6Cutters | Range [6,6] returns all and only 6-bp recognition enzymes (EcoRI, BamHI, PstI ∈; AluI, NotI ∉) | All results length 6; EcoRI/BamHI/PstI present; AluI/NotI absent | Source 3 (6-bp set) |
| M5 | EnzymesByCutLength_Range4to8_ReturnsAllUndividedSites_ExcludesSfiI | Range [4,8] returns the full library except the 13-nt interrupted palindrome SfiI | count == total − 1; SfiI absent; all results 4–8 nt | Source 2 (4–8 nt undivided); Source 6 (SfiI interrupted palindrome) |
| M6 | EnzymesByCutLength_Range9to10_ReturnsEmpty | Range above all recognition lengths returns nothing | empty | Source 2 (max 8 nt) |
| M7 | EnzymesByCutLength_Range4_ReturnsOnly4Cutters | Range [4,4] returns only 4-bp enzymes (AluI, HaeIII, TaqI ∈; EcoRI ∉) | All results length 4 | Source 3 (4-bp set) |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | EnzymesByCutLength_MinGreaterThanMax_ReturnsEmpty | min=8, max=4 (inverted) → empty | empty | INV-5 |
| S2 | EnzymesByCutLength_RangeMatchesSingleLengthOverload | Range [L,L] equals single-length overload for L=4,6,8 | set equality | INV-4 |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | EnzymesByCutLength_ZeroOrNegativeBounds_ReturnsEmpty | min=-1, max=0 → empty | empty | No site has length ≤ 0 |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- `tests/Seqeron/Seqeron.Genomics.Tests/RestrictionAnalyzer_FindSites_Tests.cs` (RESTR-FIND-001) — site finding only.
- `tests/Seqeron/Seqeron.Genomics.Tests/RestrictionAnalyzer_Digest_Tests.cs` (RESTR-DIGEST-001) — digest/map/compatibility.
- No existing tests cover `GetEnzymesByCutLength`, `GetBluntCutters`, or `GetStickyCutters`. No prior range overload exists.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1 | ❌ Missing | New unit |
| M2 | ❌ Missing | New unit |
| M3 | ❌ Missing | New unit |
| M4 | ❌ Missing | New unit |
| M5 | ❌ Missing | New unit |
| M6 | ❌ Missing | New unit |
| M7 | ❌ Missing | New unit |
| S1 | ❌ Missing | New unit |
| S2 | ❌ Missing | New unit |
| C1 | ❌ Missing | New unit |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/RestrictionAnalyzer_Filter_Tests.cs` — all filtering tests for RESTR-FILTER-001.
- **Remove:** none (no pre-existing filter tests).

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| RestrictionAnalyzer_Filter_Tests.cs | Canonical for RESTR-FILTER-001 | 10 |

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
| 8 | S1 | ❌ Missing | Implemented | ✅ Done |
| 9 | S2 | ❌ Missing | Implemented | ✅ Done |
| 10 | C1 | ❌ Missing | Implemented | ✅ Done |

**Total items:** 10
**✅ Done:** 10 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ Covered | Exact membership assertions |
| M2 | ✅ Covered | Exact membership assertions |
| M3 | ✅ Covered | Disjoint + union == full library |
| M4 | ✅ Covered | Exact length + membership |
| M5 | ✅ Covered | count == library size − 1; SfiI (interrupted palindrome) excluded |
| M6 | ✅ Covered | empty |
| M7 | ✅ Covered | exact length + membership |
| S1 | ✅ Covered | empty on inverted range |
| S2 | ✅ Covered | set equality with single-length overload |
| C1 | ✅ Covered | empty on non-positive bounds |

**Total in-scope cases:** 10 | **✅:** 10

---

## 6. Assumption Register

**Total assumptions:** 1

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | Range bounds of `GetEnzymesByCutLength(min,max)` are inclusive on both ends (API-shape convention; recognition-length values themselves are source-backed). | INV-3, M4–M7, S1, S2, C1 |

---

## 7. Open Questions / Decisions

1. None. The single assumption is API-shape only (range bound inclusivity); all biological/recognition-length values are source-backed.
