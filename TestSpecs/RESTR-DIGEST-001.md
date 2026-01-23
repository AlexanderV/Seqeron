# Test Specification: RESTR-DIGEST-001

## Test Unit Information

| Field | Value |
|-------|-------|
| **Test Unit ID** | RESTR-DIGEST-001 |
| **Area** | MolTools |
| **Title** | Restriction Digest Simulation |
| **Canonical Class** | `RestrictionAnalyzer` |
| **Canonical Methods** | `Digest`, `GetDigestSummary`, `CreateMap`, `AreCompatible`, `FindCompatibleEnzymes` |
| **Complexity** | O(n + k log k) where k = number of cut sites |
| **Status** | ☑ Complete |
| **Last Updated** | 2026-01-23 |

---

## Evidence Sources

| Source | Type | Key Information |
|--------|------|-----------------|
| Wikipedia: Restriction digest | Academic | Fragment generation, gel electrophoresis verification |
| Wikipedia: Restriction enzyme | Academic | Type II enzymes, palindromic recognition, overhang types |
| Wikipedia: Restriction map | Academic | Site positions, unique cutters identification |
| Addgene: Restriction Digest Protocol | Protocol | Fragment sum = original length, validation by gel electrophoresis |
| Roberts RJ (1976) | Research | Restriction endonuclease classification |
| REBASE | Database | Enzyme compatibility, overhang sequences |

---

## Invariants

### Digest Invariants

1. **Fragment Count**: With k unique cut positions, produces k+1 fragments (Source: Wikipedia Restriction digest)
2. **Fragment Sum**: Sum of all fragment lengths == original sequence length (Source: Addgene Protocol)
3. **Fragment Order**: Fragments are numbered sequentially from 5' to 3' (Source: Implementation)
4. **Fragment Positions**: Fragment start positions are sorted in ascending order (Source: Implementation)
5. **Boundary Fragments**: First fragment has LeftEnzyme=null, last has RightEnzyme=null (Source: Implementation)
6. **Non-Empty Fragments**: All fragments have Length > 0 (Source: Implementation)
7. **No Cuts Returns Whole**: No cut sites returns original sequence as single fragment (Source: Implementation)

### Summary Invariants

8. **Fragment Sizes Sorted**: FragmentSizes in GetDigestSummary are sorted descending (Source: Implementation)
9. **Size Bounds**: SmallestFragment ≤ AverageFragmentSize ≤ LargestFragment (Source: Mathematical)
10. **Enzyme List**: EnzymesUsed contains all input enzyme names (Source: Implementation)

### Map Invariants

11. **Non-Cutters Identified**: Enzymes with zero sites listed in NonCutters (Source: Implementation)
12. **Unique Cutters Identified**: Enzymes with exactly one forward-strand site listed in UniqueCutters (Source: Implementation)
13. **Site Count Consistency**: TotalSites counts forward-strand sites only (Source: Implementation)

### Compatibility Invariants

14. **Blunt Compatibility**: All blunt-end enzymes are compatible with each other (Source: Wikipedia)
15. **Overhang Compatibility**: Enzymes with identical overhangs are compatible (Source: Wikipedia)
16. **Symmetry**: AreCompatible(A, B) == AreCompatible(B, A) (Source: Mathematical)

---

## Test Cases

### Must (Required - Evidence-Based)

| ID | Test Case | Rationale | Source |
|----|-----------|-----------|--------|
| M1 | Single cut produces two fragments | Invariant #1 | Wikipedia |
| M2 | Fragment lengths sum to original sequence length | Invariant #2 | Addgene |
| M3 | No cut sites returns original sequence as single fragment | Invariant #7 | Implementation |
| M4 | Multiple cuts produce correct number of fragments | Invariant #1 | Wikipedia |
| M5 | First fragment has LeftEnzyme=null | Invariant #5 | Implementation |
| M6 | Last fragment has RightEnzyme=null | Invariant #5 | Implementation |
| M7 | Fragment numbers are sequential (1, 2, 3...) | Invariant #3 | Implementation |
| M8 | Fragment start positions increase monotonically | Invariant #4 | Implementation |
| M9 | Multiple enzymes produce combined cut sites | API contract | Implementation |
| M10 | No enzymes provided throws ArgumentException | API contract | Implementation |
| M11 | Null sequence throws ArgumentNullException | API contract | Implementation |
| M12 | GetDigestSummary fragment sizes sorted descending | Invariant #8 | Implementation |
| M13 | GetDigestSummary LargestFragment ≥ SmallestFragment | Invariant #9 | Implementation |
| M14 | GetDigestSummary EnzymesUsed contains input enzymes | Invariant #10 | Implementation |
| M15 | CreateMap identifies non-cutters | Invariant #11 | Implementation |
| M16 | CreateMap identifies unique cutters | Invariant #12 | Implementation |
| M17 | CreateMap TotalSites counts forward-strand only | Invariant #13 | Implementation |
| M18 | CreateMap null sequence throws ArgumentNullException | API contract | Implementation |
| M19 | Blunt enzymes compatible with each other | Invariant #14 | Wikipedia |
| M20 | Same-overhang enzymes compatible (BamHI/BglII) | Invariant #15 | Wikipedia |
| M21 | Different-overhang enzymes not compatible | Invariant #15 | Wikipedia |
| M22 | Unknown enzyme in AreCompatible returns false | API contract | Implementation |
| M23 | FindCompatibleEnzymes includes known pairs | Invariant #15 | Wikipedia |

### Should (Important)

| ID | Test Case | Rationale | Source |
|----|-----------|-----------|--------|
| S1 | Fragment sequence content is correct substring | Verify fragment extraction | Implementation |
| S2 | Very short sequence (shorter than recognition) returns whole | Edge case | Implementation |
| S3 | Adjacent cut sites handled correctly | Edge case | Implementation |
| S4 | AreCompatible is symmetric (A,B) == (B,A) | Invariant #16 | Mathematical |
| S5 | CreateMap with no enzymes specified searches all | API contract | Implementation |

### Could (Optional)

| ID | Test Case | Rationale | Source |
|----|-----------|-----------|--------|
| C1 | Performance with many cut sites | Complexity verification | Implementation |
| C2 | Fragment properties include correct enzyme names | Detailed verification | Implementation |

---

## Audit Results

### Existing Test Coverage (RestrictionAnalyzerTests.cs)

| Test | Status | Coverage |
|------|--------|----------|
| Digest_SingleCut_ReturnsTwoFragments | ✓ Covered | M1 |
| Digest_NoCuts_ReturnsWholeSequence | ✓ Covered | M3 |
| Digest_MultipleCuts_ReturnsCorrectFragments | ✓ Covered | M4 |
| Digest_FragmentsHaveCorrectProperties | ✓ Weak | M5, M7 (partial) |
| Digest_MultipleEnzymes_CutsWithBoth | ✓ Covered | M9 |
| Digest_NoEnzymes_ThrowsException | ✓ Covered | M10 |
| Digest_NullSequence_ThrowsException | ✓ Covered | M11 |
| GetDigestSummary_ReturnsCorrectSummary | ✓ Covered | M14 (partial) |
| GetDigestSummary_FragmentsSortedDescending | ✓ Covered | M12 |
| CreateMap_ReturnsCorrectMap | ✓ Covered | M16 (partial) |
| CreateMap_IdentifiesUniqueCutters | ✓ Covered | M16 |
| CreateMap_IdentifiesNonCutters | ✓ Covered | M15 |
| CreateMap_NullSequence_ThrowsException | ✓ Covered | M18 |
| AreCompatible_BluntEnzymes_AreCompatible | ✓ Covered | M19 |
| AreCompatible_SameOverhang_AreCompatible | ✓ Covered | M20 |
| AreCompatible_DifferentOverhangs_NotCompatible | ✓ Covered | M21 |
| AreCompatible_UnknownEnzyme_ReturnsFalse | ✓ Covered | M22 |
| FindCompatibleEnzymes_FindsPairs | ✓ Covered | M23 |

### Missing Tests

| ID | Test Case | Priority |
|----|-----------|----------|
| M2 | Fragment sum invariant | Must |
| M6 | Last fragment RightEnzyme=null | Must |
| M8 | Fragment start positions monotonic | Must |
| M13 | LargestFragment ≥ SmallestFragment | Must |
| M17 | TotalSites counts forward-strand only | Must |
| S1 | Fragment sequence content | Should |
| S4 | AreCompatible symmetry | Should |

### Weak Tests Needing Strengthening

| Test | Issue | Fix |
|------|-------|-----|
| Digest_FragmentsHaveCorrectProperties | Only checks first two fragments | Add comprehensive invariant checks |
| GetDigestSummary_ReturnsCorrectSummary | Missing some invariants | Add size bounds verification |
| CreateMap_ReturnsCorrectMap | Partial coverage | Add TotalSites verification |

---

## Consolidation Plan

1. **Canonical File**: Rename existing `RestrictionAnalyzerTests.cs` to `RestrictionAnalyzer_Digest_Tests.cs`
2. **Remove Smoke Test**: Remove FindSites smoke test (already covered in RESTR-FIND-001)
3. **Add Missing Tests**: M2, M6, M8, M13, M17, S1, S4
4. **Strengthen Weak Tests**: Add comprehensive invariant assertions with Assert.Multiple
5. **Group by Feature**: Organize into Digest, Summary, Map, Compatibility regions

---

## Open Questions

None - behavior is well-documented in Wikipedia and Addgene protocols.

---

## Assumptions

| ID | Assumption | Justification |
|----|------------|---------------|
| A1 | Only forward-strand cuts considered for fragment generation | Avoid double-counting palindromic sites |
| A2 | Fragments sorted by start position | Natural ordering for molecular biology |
| A3 | Cut position determines fragment boundary | Standard molecular biology convention |

