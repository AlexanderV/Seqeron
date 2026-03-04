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
| **Last Updated** | 2026-03-04 |

---

## Evidence Sources

| Source | Type | Key Information |
|--------|------|-----------------|
| Wikipedia: Restriction digest | Academic | Fragment generation, gel electrophoresis verification, fragment sum = original length |
| Wikipedia: Restriction enzyme | Academic | Type II enzymes, palindromic recognition, overhang types, enzyme Examples table with cut positions |
| Wikipedia: Restriction map | Academic | Site positions, unique cutters identification |
| Wikipedia: Sticky and blunt ends | Academic | Overhang compatibility rules: same type + complementary sequence required; blunt ends always compatible |
| Addgene: Restriction Digest Protocol | Protocol | Fragment sum = original length, validation by gel electrophoresis |
| Roberts RJ (1976) | Research | Restriction endonuclease classification |
| REBASE | Database | Enzyme data, overhang sequences |

---

## Invariants

### Digest Invariants

1. **Fragment Count**: With k unique forward-strand cut positions, produces k+1 fragments (Source: Wikipedia Restriction digest)
2. **Fragment Sum**: Sum of all fragment lengths == original sequence length (Source: Addgene Protocol — "The sum of the individual fragments should equal the size of the original fragment")
3. **Fragment Order**: Fragments are numbered sequentially from 5' to 3' (Source: Standard molecular biology convention)
4. **Fragment Positions**: Fragment start positions are sorted in ascending order (Source: Standard molecular biology convention)
5. **Boundary Fragments**: First fragment has LeftEnzyme=null, last has RightEnzyme=null (Source: DNA ends at sequence termini have no enzyme cut)
6. **Non-Empty Fragments**: All fragments have Length > 0 (Source: Cut positions are at distinct integer positions)
7. **No Cuts Returns Whole**: No cut sites returns original sequence as single fragment (Source: Standard — no double-strand break = no fragmentation)

### Summary Invariants

8. **Fragment Sizes Sorted**: FragmentSizes in GetDigestSummary are sorted descending (Source: Convention for gel electrophoresis comparison)
9. **Size Bounds**: SmallestFragment ≤ AverageFragmentSize ≤ LargestFragment (Source: Mathematical)
10. **Enzyme List**: EnzymesUsed contains all input enzyme names (Source: API contract)

### Map Invariants

11. **Non-Cutters Identified**: Enzymes with zero sites listed in NonCutters (Source: Wikipedia Restriction map)
12. **Unique Cutters Identified**: Enzymes with exactly one forward-strand site listed in UniqueCutters (Source: Wikipedia Restriction map — forward-strand only to avoid double-counting palindromic sites)
13. **Site Count Consistency**: TotalSites counts forward-strand sites only (Source: Palindromic sites appear on both strands at same position; counting both would double-count)

### Compatibility Invariants

14. **Blunt Compatibility**: All blunt-end enzymes are compatible with each other (Source: Wikipedia Sticky and blunt ends — "blunt ends are always compatible with each other")
15. **Overhang Compatibility**: Enzymes with identical overhang **type** (5' or 3') AND identical overhang **sequence** are compatible (Source: Wikipedia Sticky and blunt ends — "overhangs have to be complementary in order for the ligase to work")
15b. **Cross-Type Incompatibility**: A 5' overhang and a 3' overhang are NOT compatible, even if the overhang sequence string is the same palindrome (Source: Wikipedia — strand geometry prevents base-pairing between same-strand extensions)
16. **Symmetry**: AreCompatible(A, B) == AreCompatible(B, A) (Source: Mathematical)

---

## Test Cases

### Must (Required - Evidence-Based)

| ID | Test Case | Rationale | Source |
|----|-----------|-----------|--------|
| M1 | Single cut produces two fragments | Invariant #1 | Wikipedia |
| M2 | Fragment lengths sum to original sequence length | Invariant #2 | Addgene |
| M3 | No cut sites returns original sequence as single fragment | Invariant #7 | Wikipedia |
| M4 | Multiple cuts produce correct number of fragments | Invariant #1 | Wikipedia |
| M5 | First fragment has LeftEnzyme=null | Invariant #5 | Molecular biology |
| M6 | Last fragment has RightEnzyme=null | Invariant #5 | Molecular biology |
| M7 | Fragment numbers are sequential (1, 2, 3...) | Invariant #3 | Convention |
| M8 | Fragment start positions increase monotonically | Invariant #4 | Convention |
| M9 | Multiple enzymes produce combined cut sites | API contract | Implementation |
| M10 | No enzymes provided throws ArgumentException | API contract | Implementation |
| M11 | Null sequence throws ArgumentNullException | API contract | Implementation |
| M12 | GetDigestSummary fragment sizes sorted descending | Invariant #8 | Convention |
| M13 | GetDigestSummary LargestFragment ≥ SmallestFragment | Invariant #9 | Mathematical |
| M14 | GetDigestSummary EnzymesUsed contains input enzymes | Invariant #10 | API contract |
| M15 | CreateMap identifies non-cutters | Invariant #11 | Wikipedia |
| M16 | CreateMap identifies unique cutters (forward-strand) | Invariant #12 | Wikipedia |
| M17 | CreateMap TotalSites counts forward-strand only | Invariant #13 | Palindromic convention |
| M18 | CreateMap null sequence throws ArgumentNullException | API contract | Implementation |
| M19 | Blunt enzymes compatible with each other | Invariant #14 | Wikipedia |
| M20 | Same-type same-overhang enzymes compatible (BamHI/BglII) | Invariant #15 | Wikipedia |
| M21 | Different-overhang enzymes not compatible | Invariant #15 | Wikipedia |
| M22 | Unknown enzyme in AreCompatible returns false | API contract | Implementation |
| M23 | FindCompatibleEnzymes includes known pairs | Invariant #15 | Wikipedia |
| M24 | Cross-type overhangs NOT compatible (HindIII/SacI) | Invariant #15b | Wikipedia |
| M25 | Cross-type overhangs NOT compatible (SphI/NcoI) | Invariant #15b | Wikipedia |
| M26 | GATC overhang family all compatible (MboI/Sau3AI/BamHI/BglII) | Invariant #15 | Wikipedia |
| M27 | CTAG overhang family all compatible (NheI/XbaI/SpeI/AvrII) | Invariant #15 | Wikipedia |
| M28 | TCGA overhang family compatible (SalI/XhoI) | Invariant #15 | Wikipedia |
| M29 | Enzyme database cut positions match Wikipedia Examples table | Data | Wikipedia |
| M30 | UniqueCutters correctly identifies palindromic single-site enzymes | Invariant #12 | Wikipedia |
| M31 | UniqueCutters excludes multi-site enzymes | Invariant #12 | Wikipedia |
| M32 | FindCompatibleEnzymes excludes cross-type pairs | Invariant #15b | Wikipedia |

### Should (Important)

| ID | Test Case | Rationale | Source |
|----|-----------|-----------|--------|
| S1 | Fragment sequence content is correct substring | Verify fragment extraction | Convention |
| S2 | Very short sequence (shorter than recognition) returns whole | Edge case | Convention |
| S3 | Adjacent cut sites handled correctly | Edge case | Convention |
| S4 | AreCompatible is symmetric (A,B) == (B,A) | Invariant #16 | Mathematical |
| S5 | CreateMap with no enzymes specified searches all | API contract | Implementation |

### Could (Optional)

| ID | Test Case | Rationale | Source |
|----|-----------|-----------|--------|
| C1 | Performance with many cut sites | Complexity verification | Implementation |
| C2 | Fragment properties include correct enzyme names | Detailed verification | Implementation |

---

## Deviations and Assumptions

**None.** All behavior is derived from external sources:
- Fragment generation, counting, and sizing: Wikipedia (Restriction digest), Addgene Protocol
- Enzyme cut positions and recognition sequences: Wikipedia (Restriction enzyme) Examples table, REBASE
- Overhang compatibility rules: Wikipedia (Sticky and blunt ends)
- Unique cutter identification: Wikipedia (Restriction map)
- Forward-strand-only counting for palindromic sites: standard molecular biology convention for Type II palindromic enzymes

---

## Audit Results

### Test Coverage (RestrictionAnalyzer_Digest_Tests.cs)

| Test | Status | Coverage |
|------|--------|----------|
| Digest_SingleCut_ReturnsTwoFragmentsWithCorrectSum | ✓ | M1, M2 |
| Digest_NoCuts_ReturnsWholeSequenceAsSingleFragment | ✓ | M3 |
| Digest_MultipleCuts_ReturnsCorrectFragmentCount | ✓ | M4 |
| Digest_FragmentsHaveCorrectProperties | ✓ | M5, M6, M7, M8 |
| Digest_FragmentSequenceContent_MatchesExpectedSubstring | ✓ | S1 |
| Digest_MultipleEnzymes_CutsWithBoth | ✓ | M9, C2 |
| Digest_NoEnzymes_ThrowsArgumentException | ✓ | M10 |
| Digest_NullSequence_ThrowsArgumentNullException | ✓ | M11 |
| Digest_SequenceShorterThanRecognition_ReturnsWholeSequence | ✓ | S2 |
| Digest_AdjacentCutSites_ProducesSmallFragment | ✓ | S3 |
| GetDigestSummary_ReturnsCorrectSummaryWithInvariants | ✓ | M12, M13, M14 |
| CreateMap_ReturnsCorrectMapWithAllFields | ✓ | M13, M17, M31 |
| CreateMap_IdentifiesUniqueCutters | ✓ | M16, M30 |
| CreateMap_IdentifiesNonCutters | ✓ | M15 |
| CreateMap_TotalSites_CountsForwardStrandOnly | ✓ | M17 |
| CreateMap_NullSequence_ThrowsArgumentNullException | ✓ | M18 |
| CreateMap_NoEnzymesSpecified_SearchesAll | ✓ | S5 |
| AreCompatible_BluntEnzymes_AreCompatible | ✓ | M19 |
| AreCompatible_SameOverhang_AreCompatible_BamHI_BglII | ✓ | M20 |
| AreCompatible_SameOverhang_AreCompatible_SalI_XhoI | ✓ | M28 |
| AreCompatible_GATC_OverhangFamily_AllCompatible (3 cases) | ✓ | M26 |
| AreCompatible_CTAG_OverhangFamily_AllCompatible (3 cases) | ✓ | M27 |
| AreCompatible_CrossTypeOverhang_NotCompatible_HindIII_SacI | ✓ | M24 |
| AreCompatible_CrossTypeOverhang_NotCompatible_SphI_NcoI | ✓ | M25 |
| AreCompatible_DifferentOverhangs_NotCompatible | ✓ | M21 |
| AreCompatible_UnknownEnzyme_ReturnsFalse | ✓ | M22 |
| AreCompatible_IsSymmetric (5 cases) | ✓ | S4 |
| EnzymeDatabase_MatchesWikipediaData (20 cases) | ✓ | M29 |
| FindCompatibleEnzymes_FindsKnownPairs | ✓ | M23, M32 |
| FindCompatibleEnzymes_AllReturnedPairsAreActuallyCompatible | ✓ | M23 |

### Coverage Classification Summary

All Must (M1–M32), Should (S1–S5), and Could (C2) test cases are covered.
No missing, weak, or duplicate tests remain.

---

## Open Questions

None — all behavior is sourced from Wikipedia and Addgene protocols.
