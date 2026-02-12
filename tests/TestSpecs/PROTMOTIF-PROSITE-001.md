# Test Specification: PROTMOTIF-PROSITE-001

**Test Unit ID:** PROTMOTIF-PROSITE-001
**Area:** ProteinMotif
**Algorithm:** PROSITE Pattern Matching
**Status:** ☐ In Progress
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-02-12

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | PROSITE User Manual (PA line spec) | 2 | https://prosite.expasy.org/prosuser.html | 2026-02-12 |
| 2 | ScanProsite Documentation | 2 | https://prosite.expasy.org/scanprosite/scanprosite_doc.html | 2026-02-12 |
| 3 | Hulo et al. (2007) | 1 | https://doi.org/10.1093/nar/gkm977 | 2026-02-12 |
| 4 | De Castro et al. (2006) | 1 | https://doi.org/10.1093/nar/gkl124 | 2026-02-12 |
| 5 | PROSITE PS00001 entry | 2 | https://prosite.expasy.org/PS00001 | 2026-02-12 |
| 6 | PROSITE PS00028 entry | 2 | https://prosite.expasy.org/PS00028 | 2026-02-12 |
| 7 | ScanProsite P02787 vs PS00001 | 2 | ScanProsite JSON output | 2026-02-12 |

### 1.2 Key Evidence Points

1. PROSITE pattern syntax defines 10 elements: letters, x, [ABC], {ABC}, -, (n), (n,m), <, >, period — PROSITE User Manual §IV.E
2. Range (n,m) only valid with x, fixed (n) valid with any element — ScanProsite docs
3. Human Transferrin (P02787) has 2 N-glycosylation sites per ScanProsite — ScanProsite output
4. PS00028 zinc finger uses both variable ranges x(2,4) and x(3,5) — PROSITE PS00028

### 1.3 Documented Corner Cases

1. Empty pattern → empty output
2. Trailing period in data file patterns → should be stripped
3. `[G>]` C-terminus inside brackets (PS00267, PS00539) → converts to `(?:G|$)` alternation

### 1.4 Known Failure Modes / Pitfalls

1. Confusing PROSITE repetition with regex quantifiers — PROSITE `A(3)` → regex `A{3}`, not `A(3)`
2. Incorrect exclusion handling — `{P}` must become `[^P]` not `{P}` in regex

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `ConvertPrositeToRegex(prositePattern)` | ProteinMotifFinder | **Canonical** | Converts PROSITE notation to regex |
| `FindMotifByProsite(sequence, pattern, name)` | ProteinMotifFinder | **Canonical** | End-to-end PROSITE matching |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | ConvertPrositeToRegex("") == "" | Yes | Edge case |
| INV-2 | ConvertPrositeToRegex always produces valid .NET regex | Yes | Implementation requirement |
| INV-3 | FindMotifByProsite("", any) yields no matches | Yes | Edge case |
| INV-4 | FindMotifByProsite(any, "") yields no matches | Yes | Edge case |
| INV-5 | FindMotifByProsite is case-insensitive | Yes | Implementation spec |
| INV-6 | Match positions are 0-based inclusive [Start, End] | Yes | Implementation convention |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | ConvertProsite_SimpleLiteral | Convert `R-G-D` (PS00016) | `RGD` | PROSITE PS00016 |
| M2 | ConvertProsite_AnyAminoAcid | Convert `A-x-G` | `A.G` | PROSITE User Manual |
| M3 | ConvertProsite_ExactRepeat | Convert `x(3)-A` | `.{3}A` | PROSITE User Manual |
| M4 | ConvertProsite_RangeRepeat | Convert `A-x(2,4)-G` | `A.{2,4}G` | PROSITE User Manual |
| M5 | ConvertProsite_CharacterClass | Convert `[ST]-x-[RK]` (PS00005) | `[ST].[RK]` | PROSITE PS00005 |
| M6 | ConvertProsite_ExclusionClass | Convert `N-{P}-[ST]-{P}` (PS00001) | `N[^P][ST][^P]` | PROSITE PS00001 |
| M7 | ConvertProsite_NTerminus | Convert `<M-x-K` | `^M.K` | PROSITE User Manual |
| M8 | ConvertProsite_CTerminus | Convert `A-x-G>` | `A.G$` | PROSITE User Manual |
| M9 | ConvertProsite_ElementRepeat | Convert `[RK](2)-x-[ST]` (PS00004) | `[RK]{2}.[ST]` | PROSITE PS00004 |
| M10 | ConvertProsite_ComplexPS00028 | Convert full zinc finger pattern | `C.{2,4}C.{3}[LIVMFYWC].{8}H.{3,5}H` | PROSITE PS00028 |
| M11 | ConvertProsite_ComplexPS00008 | Convert N-myristoylation | `G[^EDRKHPFYW].{2}[STAGCN][^P]` | PROSITE PS00008 |
| M12 | ConvertProsite_ComplexPS00018 | Convert EF-hand pattern | Full expected regex | PROSITE PS00018 |
| M13 | ConvertProsite_EmptyInput | Convert empty string | Empty string | Edge case |
| M14 | ConvertProsite_TrailingPeriod | Convert `R-G-D.` | `RGD` | PROSITE User Manual |
| M15 | FindProsite_NGlycoMatch | Match PS00001 on `AANASAAANGTAAA` | 2 matches at positions 2 and 8 | Manual derivation from PS00001 |
| M16 | FindProsite_RGDMatch | Match `R-G-D` on `AAARGDAAA` | 1 match at position 3 | PS00016 |
| M17 | FindProsite_NoMatch | Match PS00001 on `AANPSAAA` | 0 matches | P at pos 3 violates {P} |
| M18 | FindProsite_EmptySequence | Match any pattern on empty | 0 matches | Edge case |
| M19 | FindProsite_EmptyPattern | Match empty pattern on any | 0 matches | Edge case |
| M20 | FindProsite_CaseInsensitive | Match PS00016 uppercase=lowercase | Same results | INV-5 |
| M21 | FindProsite_MultipleMatches | Match PS00005 on multi-site sequence | Multiple exact positions | PS00005 |
| M22 | ConvertProsite_CTermBrackets_PS00267 | Convert `F-[IVFY]-G-[LM]-M-[G>].` | `F[IVFY]G[LM]M(?:G\|$)` | PROSITE User Manual §IV.E, PS00267 |
| M22b | ConvertProsite_CTermBrackets_PS00539 | Convert `F-[GSTV]-P-R-L-[G>].` | `F[GSTV]PRL(?:G\|$)` | PROSITE User Manual §IV.E, PS00539 |
| M23 | FindProsite_CTermBrackets_MatchG | Match PS00267 with G at final pos | Match FVGLMG | PS00267 |
| M23b | FindProsite_CTermBrackets_MatchEnd | Match PS00267 at C-terminus | Match FVGLM | PS00267 (> branch) |
| M23c | FindProsite_CTermBrackets_NoMidMatch | PS00267 mid-sequence without G | 0 matches | PS00267 (rejects mid-seq) |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | FindProsite_RealProtein | Human Transferrin vs PS00001 | 2 matches | ScanProsite verified |
| S2 | FindProsite_NTerminalAnchor | `<M-x-K` matches only at start | Match at 0 only | PROSITE User Manual |
| S3 | FindProsite_CTerminalAnchor | `A-x-G>` matches only at end | Match at end only | PROSITE User Manual |
| S4 | FindProsite_MatchedSequence | Verify matched sequence string | Exact subsequence | Implementation contract |
| S5 | ConvertProsite_AminoAcidRepeat | Convert `L-x(6)-L` (PS00029) | `L.{6}L` (partial) | PS00029 |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | FindProsite_ZincFingerComplex | Match PS00028 on synthetic | Match at expected position | Complex pattern |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- 7 PROSITE tests in `ProteinMotifFinderTests.cs` (lines 12–77)
- No tests in `ProteinMotifFinder_MotifSearch_Tests.cs` for PROSITE conversion (that file covers PROTMOTIF-FIND-001)

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Canonical Test Method |
|---------------------|--------|----------------------|
| M1: Simple literal | ✅ Covered | `ConvertPrositeToRegex_SimpleLiteralRGD_ProducesExactRegex` |
| M2: Any amino acid | ✅ Covered | `ConvertPrositeToRegex_AnyAminoAcid_ProducesDot` |
| M3: Exact repeat x(n) | ✅ Covered | `ConvertPrositeToRegex_ExactRepeatXN_ProducesQuantifier` |
| M4: Range repeat x(n,m) | ✅ Covered | `ConvertPrositeToRegex_RangeRepeatXNM_ProducesRangeQuantifier` |
| M5: Character class | ✅ Covered | `ConvertPrositeToRegex_CharacterClass_PreservesSquareBrackets` |
| M6: Exclusion class (full PS00001) | ✅ Covered | `ConvertPrositeToRegex_ExclusionClass_ProducesNegatedCharClass` |
| M7: N-terminus | ✅ Covered | `ConvertPrositeToRegex_NTerminusAnchor_ProducesCaret` |
| M8: C-terminus | ✅ Covered | `ConvertPrositeToRegex_CTerminusAnchor_ProducesDollar` |
| M9: Element repetition A(n) | ✅ Covered | `ConvertPrositeToRegex_ElementRepetition_ProducesQuantifier` |
| M10: Complex PS00028 | ✅ Covered | `ConvertPrositeToRegex_ComplexZincFinger_ProducesFullRegex` |
| M11: Complex PS00008 | ✅ Covered | `ConvertPrositeToRegex_ComplexMyristoylation_ProducesFullRegex` |
| M12: Complex PS00018 | ✅ Covered | `ConvertPrositeToRegex_ComplexEFHand_ProducesFullRegex` |
| M13: Empty input | ✅ Covered | `ConvertPrositeToRegex_EmptyString_ReturnsEmpty`, `ConvertPrositeToRegex_Null_ReturnsEmpty` |
| M14: Trailing period | ✅ Covered | `ConvertPrositeToRegex_TrailingPeriod_IsIgnored` |
| M15: N-glyco match positions | ✅ Covered | `FindMotifByProsite_NGlycosylation_MatchesAtCorrectPositions` |
| M16: RGD match position | ✅ Covered | `FindMotifByProsite_RGDPattern_MatchesAtExactPosition` |
| M17: No match | ✅ Covered | `FindMotifByProsite_ExclusionBlocks_ReturnsNoMatches` |
| M18: Empty sequence | ✅ Covered | `FindMotifByProsite_EmptySequence_ReturnsNoMatches` |
| M19: Empty pattern | ✅ Covered | `FindMotifByProsite_EmptyPattern_ReturnsNoMatches` |
| M20: Case insensitivity | ✅ Covered | `FindMotifByProsite_CaseInsensitive_ProducesSameResults` |
| M21: Multiple matches | ✅ Covered | `FindMotifByProsite_MultipleMatches_FindsAllPositions` |
| M22: `[G>]` conversion PS00267 | ✅ Covered | `ConvertPrositeToRegex_PS00267_CTermInsideBrackets_ConvertsCorrectly` |
| M22b: `[G>]` conversion PS00539 | ✅ Covered | `ConvertPrositeToRegex_PS00539_CTermInsideBrackets_ConvertsCorrectly` |
| M23: `[G>]` match with G | ✅ Covered | `FindMotifByProsite_Tachykinin_MatchesWithG` |
| M23b: `[G>]` match at C-terminus | ✅ Covered | `FindMotifByProsite_Tachykinin_MatchesAtCTerminus` |
| M23c: `[G>]` no mid-seq match | ✅ Covered | `FindMotifByProsite_Tachykinin_NoMatchInMiddleWithoutG` |
| S1: Real protein | ✅ Covered | `FindMotifByProsite_HumanTransferrin_FindsExactlyTwoNGlycoSites` |
| S2: N-terminal anchor | ✅ Covered | `FindMotifByProsite_NTerminalAnchor_OnlyMatchesAtStart` |
| S3: C-terminal anchor | ✅ Covered | `FindMotifByProsite_CTerminalAnchor_OnlyMatchesAtEnd` |
| S4: Matched sequence | ✅ Covered | `FindMotifByProsite_MatchedSequence_IsExactSubstring` |
| S5: Leucine zipper partial | ✅ Covered | `ConvertPrositeToRegex_LeucineZipperLongRepeat_ProducesCorrectRegex` |
| C1: Zinc finger complex | ✅ Covered | `FindMotifByProsite_ZincFingerC2H2_MatchesSyntheticDomain` |
| Bonus: Match properties | ✅ Covered | `FindMotifByProsite_MatchProperties_AllFieldsPopulated` |

**Missing:** 0 &emsp; **Weak:** 0 &emsp; **Duplicate:** 0

### 5.3 Consolidation Plan

- **Canonical file:** `ProteinMotifFinder_PrositePattern_Tests.cs` — all PROSITE conversion and matching tests
- **Remove:** 7 PROSITE tests from `ProteinMotifFinderTests.cs` (lines 12–77, `#region PROSITE Pattern Conversion Tests`)
- **Keep:** All non-PROSITE tests in `ProteinMotifFinderTests.cs` unchanged

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| `ProteinMotifFinder_PrositePattern_Tests.cs` | Canonical PROSITE tests | 34 |
| `ProteinMotifFinderTests.cs` | Residual (non-PROSITE): 12 strengthened, 11 OK | 23 |

### 5.5 Phase 7 Work Queue

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| 1 | M1 | ✅ Covered | Rewrite in canonical file | ✅ Done |
| 2 | M2 | ✅ Covered | Rewrite in canonical file | ✅ Done |
| 3 | M3 | ❌ Missing | Implement | ✅ Done |
| 4 | M4 | ✅ Covered | Rewrite in canonical file | ✅ Done |
| 5 | M5 | ✅ Covered | Rewrite in canonical file | ✅ Done |
| 6 | M6 | ⚠ Weak | Rewrite with full PS00001 pattern | ✅ Done |
| 7 | M7 | ⚠ Weak | Rewrite with exact assertion | ✅ Done |
| 8 | M8 | ❌ Missing | Implement | ✅ Done |
| 9 | M9 | ❌ Missing | Implement | ✅ Done |
| 10 | M10 | ❌ Missing | Implement | ✅ Done |
| 11 | M11 | ❌ Missing | Implement | ✅ Done |
| 12 | M12 | ❌ Missing | Implement | ✅ Done |
| 13 | M13 | ❌ Missing | Implement | ✅ Done |
| 14 | M14 | ❌ Missing | Implement | ✅ Done |
| 15 | M15 | ❌ Missing | Implement | ✅ Done |
| 16 | M16 | ❌ Missing | Implement | ✅ Done |
| 17 | M17 | ❌ Missing | Implement | ✅ Done |
| 18 | M18 | ❌ Missing | Implement | ✅ Done |
| 19 | M19 | ❌ Missing | Implement | ✅ Done |
| 20 | M20 | ❌ Missing | Implement | ✅ Done |
| 21 | M21 | ❌ Missing | Implement | ✅ Done |
| 22 | S1 | ❌ Missing | Implement | ✅ Done |
| 23 | S2 | ❌ Missing | Implement | ✅ Done |
| 24 | S3 | ❌ Missing | Implement | ✅ Done |
| 25 | S4 | ❌ Missing | Implement | ✅ Done |
| 26 | S5 | ❌ Missing | Implement | ✅ Done |
| 27 | C1 | ❌ Missing | Implement | ✅ Done |
| 28 | FindMotifByProsite_UsesConversion | ⚠ Weak | Remove (superseded by M15) | ✅ Done |
| 29 | M22 | ❌ Missing | Implement `[G>]` conversion PS00267 | ✅ Done |
| 30 | M22b | ❌ Missing | Implement `[G>]` conversion PS00539 | ✅ Done |
| 31 | M23 | ❌ Missing | Implement `[G>]` match with G | ✅ Done |
| 32 | M23b | ❌ Missing | Implement `[G>]` match at C-terminus | ✅ Done |
| 33 | M23c | ❌ Missing | Implement `[G>]` no mid-seq match | ✅ Done |

**Total items:** 33
**✅ Done:** 33 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1: Simple literal | ✅ | Canonical test in new file |
| M2: Any amino acid | ✅ | Canonical test in new file |
| M3: Exact repeat x(n) | ✅ | New test implemented |
| M4: Range repeat x(n,m) | ✅ | Canonical test in new file |
| M5: Character class | ✅ | Canonical test in new file |
| M6: Exclusion class (full PS00001) | ✅ | Rewritten with exact assertion |
| M7: N-terminus | ✅ | Rewritten with exact assertion |
| M8: C-terminus | ✅ | New test implemented |
| M9: Element repetition | ✅ | New test implemented |
| M10: Complex PS00028 | ✅ | New test implemented |
| M11: Complex PS00008 | ✅ | New test implemented |
| M12: Complex PS00018 | ✅ | New test implemented |
| M13: Empty input | ✅ | New test implemented |
| M14: Trailing period | ✅ | New test implemented |
| M15: N-glyco match positions | ✅ | New test implemented |
| M16: RGD match position | ✅ | New test implemented |
| M17: No match | ✅ | New test implemented |
| M18: Empty sequence | ✅ | New test implemented |
| M19: Empty pattern | ✅ | New test implemented |
| M20: Case insensitivity | ✅ | New test implemented |
| M21: Multiple matches | ✅ | New test implemented |
| S1: Real protein | ✅ | New test implemented |
| S2: N-terminal anchor | ✅ | New test implemented |
| S3: C-terminal anchor | ✅ | New test implemented |
| S4: Matched sequence | ✅ | New test implemented |
| S5: Leucine zipper partial | ✅ | New test implemented |
| C1: Zinc finger complex | ✅ | New test implemented |
| M22: `[G>]` conversion PS00267 | ✅ | New test implemented |
| M22b: `[G>]` conversion PS00539 | ✅ | New test implemented |
| M23: `[G>]` match with G | ✅ | New test implemented |
| M23b: `[G>]` match at C-terminus | ✅ | New test implemented |
| M23c: `[G>]` no mid-seq match | ✅ | New test implemented |

---

## 6. Assumption Register

**Total assumptions:** 0

None. All PROSITE syntax elements are fully implemented per the PROSITE User Manual §IV.E,
including `[G>]` C-terminus inside brackets (PS00267, PS00539).

---

## 7. Open Questions / Decisions

None. All evidence points resolved.
