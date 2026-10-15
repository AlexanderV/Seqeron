# Test Specification: RESTR-FIND-001

## Test Unit Information

| Field | Value |
|-------|-------|
| **Test Unit ID** | RESTR-FIND-001 |
| **Area** | MolTools |
| **Title** | Restriction Site Detection |
| **Canonical Class** | `RestrictionAnalyzer` |
| **Canonical Methods** | `FindSites`, `FindAllSites`, `FindPatternOccurrences`, `GetEnzyme`, `CreateMap` (site reporting) |
| **Complexity** | O(n × k × m) |
| **Status** | ☑ Complete |
| **Last Updated** | 2026-10-14 |

---

## Evidence Sources

| Source | Type | Key Information |
|--------|------|-----------------|
| Wikipedia: Restriction enzyme | Academic | Type II enzymes, palindromic recognition, cut positions, overhang types |
| Wikipedia: Restriction site | Academic | Recognition sequences 4-8 bp, sticky vs blunt ends |
| Wikipedia: EcoRI | Academic | GAATTC recognition, 5' overhang AATT, cut position G↓AATTC |
| Wikipedia: BamHI | Academic | GGATCC recognition, 5' overhang GATC |
| Wikipedia: HindIII | Academic | AAGCTT recognition, 5' overhang AGCT |
| Wikipedia: PstI | Academic | CTGCAG recognition, 3' cohesive termini |
| Wikipedia: EcoRV | Academic | GATATC recognition, blunt ends |
| Wikipedia: List of restriction enzyme cutting sites | Academic | Cut positions for 21 enzymes cross-verified |
| Roberts RJ (1976) | Research | Restriction endonuclease classification |
| REBASE | Database | Comprehensive enzyme database, recognition sequences |
| **Biopython 1.88 `Bio.Restriction`** | Reference implementation (oracle) | Site positions, cut positions (1-based top strand), linear cut-feasibility filter `AbstractCut._drop`, circular search `FormattedSeq.finditer` (`data + data[1:size]`), `Defined._modify` / `_rev_modify`, `Palindromic._search`, `NonPalindromic._search`, `elucidate()`, `is_blunt()`/`is_5overhang()`/`is_3overhang()` |
| **REBASE `allenz.104` (2021), field `<6>`** | Database snapshot | The cognate-methyltransferase site of every built-in (D7 / F52); vendored in `crisprVerse/crisprBase`, pinned in `docs/Evidence/B08-Source-Provenance.md` §3a |
| **REBASE EMBOSS files version 404 (2024)** | Database snapshot | The enzyme table audited field by field; shipped as `Bio/Restriction/Restriction_Dictionary.py` in Biopython 1.88, whose generated header states the REBASE version (rebase.neb.com is unreachable from the validation environment) |

---

## Invariants

1. **Position Range**: 0 ≤ Position ≤ sequence.Length - recognitionLength (linear); 0 ≤ Position < sequence.Length (circular) (Source: Implementation)
2. **Cut Position**: Position ≤ CutPosition ≤ Position + recognitionLength (linear) (Source: Implementation)
3. **Recognition Match**: RecognizedSequence.Length == Enzyme.RecognitionSequence.Length (Source: Implementation)
4. **Case Insensitivity**: GetEnzyme("EcoRI") == GetEnzyme("ecori") (Source: Implementation)
5. **Empty Returns Empty**: Empty sequence returns no sites (Source: Implementation)
6. **Unknown Throws**: Unknown enzyme name throws ArgumentException (Source: Implementation)
7. **Both Strands Searched**: Forward and reverse complement both searched (Source: Biological requirement)
8. **Coordinate convention**: Position and CutPosition are 0-based forward-strand coordinates for BOTH strands; CutPosition is the number of forward-strand bases 5' of the cut made on the strand named by IsForwardStrand; `biopython search() == forward CutPosition + 1` (Source: Biopython `Bio.Restriction`, `Defined._modify`; REBASE)
9. **Reverse-strand cut**: for a reverse-strand site, CutPosition = Position + recognitionLength − CutPositionForward (= Position + CutPositionReverse for every palindromic enzyme) (Source: Biopython `Defined._rev_modify` / `NonPalindromic._search`)
10. **Reverse-strand window**: RecognizedSequence is the matched window as read on the strand named by IsForwardStrand (Source: strand semantics; consistent with CRISPR-PAM-001 F2)
11. **Linear cut feasibility**: a site is reported only when both of its strand cuts lie in [1, Length − 1] (Source: Biopython `AbstractCut._drop`)
12. **Circular search**: origin-spanning windows are found and cut positions are reduced modulo Length; nothing is dropped (Source: Biopython `FormattedSeq.finditer`, `AbstractCut._drop` circular branch)
13. **Overlapping sites**: every occurrence of a self-overlapping recognition sequence is reported (Source: Biopython's lookahead site regexes, `compsite`)
14. **Short molecule**: no site is reported when Length < recognitionLength (Source: biology; deliberate divergence from the reference's unconditional circular padding)
15. **Pattern validation**: an empty or non-IUPAC recognition pattern raises ArgumentException eagerly (Source: API contract)
16. **CreateMap partition**: SitesByEnzyme.Keys ∪ NonCutters covers the considered enzyme set exactly, both keyed by canonical table name (Source: RESTR-DIGEST-001 invariant #11)
17. **Sequence ambiguity (D3)**: per position, a site `N` matches any of the 15 IUPAC codes in the sequence (an `N` included); every other site code matches only a concrete A/C/G/T base it expands to, so a sequence ambiguity code blocks it; a non-IUPAC sequence character matches nothing (Source: Biopython 1.88 `compsite` — site `N` → regex `.`, other codes → their A/C/G/T class; `SfiI.search(Seq("GGCCNNNNNGGCC")) == [9]`, `EcoRI.search(Seq("AGAANTCA")) == []`)
18. **Cognate methylation site (D7)**: every built-in's `RebaseMethylationField` equals REBASE `allenz.104` field `<6>` verbatim, and `CognateMethylationSites` is its parse per the file's format header (`X(Y)`, negative X = complementary strand, `?` = base undetermined, `(4)/(5)/(6)` = m4C/m5C/m6A); a known top-strand position names an `A` for m6A and a `C` for m4C/m5C. The field is the cognate MTase's protection site — no Dam/Dcm/CpG sensitivity is modelled (Source: REBASE `allenz.104`, sha256 `43b12432…69388`)

---

## Test Cases

### Must (Required - Evidence-Based)

| ID | Test Case | Rationale | Source |
|----|-----------|-----------|--------|
| M1 | EcoRI finds GAATTC at correct position | Known recognition sequence | Wikipedia (EcoRI) |
| M2 | BamHI finds GGATCC at correct position | Known recognition sequence | Wikipedia |
| M3 | FindSites returns multiple sites for repeated pattern | Sequence may contain multiple sites | Wikipedia |
| M4 | FindSites returns empty for sequence without recognition site | No false positives | Implementation |
| M5 | GetEnzyme case-insensitive lookup works | Invariant #4 | Implementation |
| M6 | GetEnzyme returns null for unknown enzyme | API contract | Implementation |
| M7 | FindSites throws ArgumentException for unknown enzyme | Invariant #6 | Implementation |
| M8 | FindSites with empty sequence returns empty | Invariant #5 | Implementation |
| M9 | CutPosition calculated correctly (Position + CutPositionForward) | Invariant #2 | Implementation |
| M10 | RecognizedSequence matches the actual sequence at position | Invariant #3 | Implementation |
| M11 | FindAllSites finds sites from multiple enzymes | API contract | Implementation |
| M12 | FindSites with multiple enzyme names finds all | API contract | Implementation |
| M13 | Custom enzyme support works | API extensibility | Implementation |
| M14 | 5' overhang enzyme (EcoRI) properties correct | Overhang type | Wikipedia (EcoRI) |
| M15 | 3' overhang enzyme (PstI) properties correct | Overhang type | Wikipedia |
| M16 | Blunt cutter enzyme (EcoRV) properties correct | Overhang type | Wikipedia |
| M35 | Every built-in's REBASE methylation-site field equals allenz.104 `<6>` | D7: the field was claimed unobtainable | REBASE allenz.104 |
| M36 | REBASE `<6>` notation parsed (complementary strand, lists, `?`); malformed values rejected | Faithful typed representation | REBASE allenz.104 format header |

### Should (Important)

| ID | Test Case | Rationale | Source |
|----|-----------|-----------|--------|
| S1 | 4-cutter, 6-cutter, 8-cutter all detected | Different recognition lengths | Wikipedia |
| S2 | Enzyme database contains 30+ enzymes | Comprehensive coverage | REBASE |
| S3 | GetEnzymesByCutLength returns correct enzymes | Filtering by length | Implementation |
| S4 | GetBluntCutters returns only blunt enzymes | Filtering by type | Implementation |
| S5 | GetStickyCutters returns only sticky enzymes | Filtering by type | Implementation |

### Could (Optional)

| ID | Test Case | Rationale | Source |
|----|-----------|-----------|--------|
| C1 | IUPAC ambiguous recognition sequences work | Degenerate sequences (HincII: GTYRAC) | REBASE, Implementation |
| C2 | Both strands searched for palindromic sites | Biological correctness | Wikipedia |

---

## Coverage Classification

### Covered (✅)

| ID | Test Case | Covering Tests |
|----|-----------|----------------|
| M1 | EcoRI finds GAATTC at correct position | `GetEnzyme_EcoRI_*`, `FindSites_EcoRI_FindsSiteAtCorrectPosition`, `EnzymeDatabase_CutPositions_MatchWikipedia` |
| M2 | BamHI finds GGATCC at correct position | `GetEnzyme_BamHI_*`, `FindSites_BamHI_FindsSiteAtCorrectPosition`, `EnzymeDatabase_CutPositions_MatchWikipedia` |
| M3 | Multiple sites found | `FindSites_MultipleSites_FindsAllAtCorrectPositions` |
| M4 | No sites returns empty | `FindSites_NoSites_ReturnsEmptyCollection` |
| M5 | Case insensitive lookup | `GetEnzyme_CaseInsensitive_ReturnsCorrectEnzyme` (×4 casings) |
| M6 | Unknown returns null | `GetEnzyme_UnknownOrInvalidName_ReturnsNull` (×4 inputs) |
| M7 | Unknown enzyme throws | `FindSites_UnknownEnzyme_ThrowsArgumentException` |
| M8 | Empty sequence returns empty | `FindSites_EmptyStringSequence_ReturnsEmptyCollection` |
| M9 | CutPosition calculation | `FindSites_EcoRI_CutPositionCalculatedCorrectly` (exact: 4 = 3+1) |
| M10 | RecognizedSequence matches | `FindSites_EcoRI_RecognizedSequenceMatchesEnzyme` |
| M11 | FindAllSites multi-enzyme | `FindAllSites_FindsSitesFromMultipleEnzymes` |
| M12 | FindSites multiple names | `FindSites_MultipleEnzymes_FindsSitesForAllEnzymes` |
| M13 | Custom enzyme support | `FindSites_CustomEnzyme_FindsSiteCorrectly`, `FindSites_CustomEnzymeWithAsymmetricCut_CutPositionCorrect` |
| M14 | 5' overhang (EcoRI) | `GetEnzyme_EcoRI_*`, `FindSites_EcoRI_OverhangSequenceIsAATT` |
| M15 | 3' overhang (PstI) | `RestrictionEnzyme_PstI_HasThreePrimeOverhang`, `FindSites_PstI_ProducesThreePrimeOverhangOfFourBases` |
| M16 | Blunt (EcoRV) | `RestrictionEnzyme_EcoRV_IsBluntCutter`, `FindSites_EcoRV_ProducesBluntEndCuts` |
| S1 | 4/6/8-cutters detected | `GetEnzymesByCutLength_ReturnsEnzymesOfCorrectLength` (×3 lengths) |
| S2 | 30+ enzymes in database | `Enzymes_ContainsAtLeast30CommonEnzymes` |
| S3 | GetEnzymesByCutLength | `GetEnzymesByCutLength_ReturnsEnzymesOfCorrectLength` |
| S4 | GetBluntCutters | `GetBluntCutters_ReturnsOnlyBluntEndEnzymes` |
| S5 | GetStickyCutters | `GetStickyCutters_ReturnsOnlyStickyEndEnzymes` |
| C1 | IUPAC degenerate recognition | `FindSites_HincII_MatchesAllDegenerateCombinations` (×4), `*_RejectsNonMatchingDegenerateBases` (×4), `*_DegenerateEnzymeProducesCorrectCutPosition` |
| C2 | Palindromic both strands | `FindSites_PalindromicSite_FoundOnBothStrandsAtSamePosition` |
| M17 | Table audited against REBASE 404 | `EnzymeDatabase_AllEntries_MatchRebase404` (×39), `EnzymeDatabase_AuditCoversEveryEntry`, `EnzymeDatabase_EveryEntry_HasSymmetricCutGeometry` |
| M18 | Cut positions vs Biopython | `FindSites_EcoRI_SiteFillsWholeSequence_IsStillReported`, `FindSites_Sau3AI_EdgeSites_AreNotCuttable_MatchesBiopython`, `FindSites_HincII_ReverseRecognizedSequenceIsReadOnItsOwnStrand` (plus the 34,906-comparison sweep recorded in the Evidence file) |
| M19 | Linear cut feasibility | `FindSites_Sau3AI_EdgeSites_AreNotCuttable_MatchesBiopython`, `FindSites_GatcEdgeCutter_WholeSequenceIsTheSite_ReportsNothing` (×2), `FindSites_Sau3AI_SiteAtThreePrimeEnd_ReportsNothing` |
| M20 | Circular origin-spanning | `FindSites_EcoRI_Circular_FindsOriginSpanningSite`, `FindSites_BamHI_Circular_FindsOriginSpanningSite`, `FindSites_Sau3AI_Circular_ReportsEveryGatc`, `FindAllSites_Circular_FindsOriginSpanningSite`, `CreateMap_Circular_FindsOriginSpanningSite` |
| M21 | Circular modulo cut | `FindSites_Circular_CutPositionIsReducedModuloLength` |
| M22 | Non-palindromic reverse cut | `FindSites_CustomNonPalindromicEnzyme_CdiI_ReverseStrandCutMatchesReference`, `…_SimI_…` |
| M23 | Reverse-strand window | `FindSites_HincII_ReverseRecognizedSequenceIsReadOnItsOwnStrand` |
| M24 | Overlapping sites | `FindPatternOccurrences_ReportsOverlappingOccurrences` |
| M25 | Pattern validation | `FindSites_CustomEnzyme_InvalidRecognitionSequence_ThrowsArgumentException` (×3), `FindSites_CustomEnzyme_LowerCaseRecognitionSequence_StillMatches` |
| M26 | Eager enzyme-name validation | `FindSites_RawString_UnknownEnzyme_ThrowsAtCallTime` |
| M27-M30 | `FindPatternOccurrences` | `FindPatternOccurrences_ReportsZeroBasedStartIndices`, `…_ReportsOverlappingOccurrences`, `…_DegeneratePattern_MatchesEveryExpansion` (×4), `…_DegeneratePattern_RejectsNonMatchingBases`, `…_AmbiguousSequenceBase_NeverMatches` (×2), `…_AcceptsSequenceThatDnaSequenceRejects`, `…_IsCaseInsensitive`, `…_StartIndex_SkipsEarlierOccurrences` (×5), `…_SequenceTooShort_ReturnsEmpty` (×2), `…_InvalidArguments_Throw`, `…_ValidatesEagerly`, `…_AgreesWithFindSitesForwardPositions` |
| M31-M32 | CreateMap site reporting | `CreateMap_LowerCaseEnzymeNames_DoesNotListACutterAsNonCutter`, `CreateMap_NoEnzymesSpecified_NonCuttersCoverTheRestOfTheTable`, `CreateMap_RepeatedEnzymeName_IsConsideredOnce`, `CreateMap_UnknownEnzymeName_ThrowsArgumentException` |
| INV-14 | Short molecule | `FindSites_Circular_MoleculeShorterThanSite_ReportsNothing` |
| INV-17 / M33 | Sequence ambiguity (D3) | `FindSites_SfiI_SequenceNInsideSpacer_MatchesBiopython` (×2), `FindSites_SfiI_SequenceN_TopStrandCutsMatchBiopython` (×3), `FindSites_SequenceAmbiguityAtConcreteSiteCode_StillBlocks` (×5), `SiteN_AcceptsEveryIupacCode_InFindSitesAndFindPatternOccurrences`, `FindSites_NonIupacCharacterInSpacer_MatchesNothing`, plus `FindPatternOccurrences_AmbiguousSequenceBase_NeverMatches` (×2, unchanged) |
| INV-17 / M34 | ASCII-only fold (F49): 'ſ' (U+017F) never matches and a pattern containing it is rejected | `RawStringEntryPoints_LongS_IsNotFoldedOntoIupacS`, `Pattern_WithLongS_IsRejectedLikeAnyNonIupacCharacter` |
| INV-18 / M35-M36 | Cognate methylation site (D7 / F52): table locked to allenz.104 `<6>` (M35); notation parse + rejection of malformed / `(99)` values (M36) | `EnzymeDatabase_MethylationField_MatchesRebaseAllenz104` (×39), `EnzymeDatabase_MethylationAuditCoversEveryEntry`, `CognateMethylationSites_BuiltIns_ParsedToTypedSites`, `CognateMethylationSites_RebaseNotation_ComplementaryStrandAndLists`, `RestrictionEnzyme_MalformedMethylationField_Throws` (×7), `RestrictionEnzyme_NullMethylationField_Throws`; MCP `GetEnzyme_ReportsRebaseCognateMethylationSite` |

### External Source Verification Tests

| Test | Verified Against | Source |
|------|-----------------|--------|
| FindSites_EcoRI_OverhangSequenceIsAATT | Overhang = AATT (4 bp, 5') | Wikipedia: EcoRI |
| FindSites_BamHI_OverhangSequenceIsGATC | Overhang = GATC (4 bp, 5') | Wikipedia: BamHI |
| FindSites_HindIII_OverhangSequenceIsAGCT | Overhang = AGCT (4 bp, 5') | Wikipedia: HindIII |
| FindSites_PstI_ProducesThreePrimeOverhangOfFourBases | 3' cohesive termini, 4 bp | Wikipedia: PstI |
| FindSites_EcoRV_ProducesBluntEndCuts | Blunt ends (fwd == rev cut) | Wikipedia: EcoRV |
| FindSites_PalindromicSite_FoundOnBothStrandsAtSamePosition | Both strands recognized | Wikipedia: Restriction enzyme |
| EnzymeDatabase_CutPositions_MatchWikipedia (×21) | Recognition sequences and cut positions | Wikipedia: Restriction enzyme (Examples table) |
| EnzymeDatabase_AllEntries_MatchRebase404 (×39) | site, fst5, len(site)+fst3, blunt/5'/3' class for every table entry — 39/39 matched, 0 wrong entries | REBASE 404 via Biopython 1.88 |
| FindSites_Sau3AI_EdgeSites_AreNotCuttable_MatchesBiopython | `Sau3AI.search(Seq("GATCGATCGATC")) == [5]` | Biopython 1.88 |
| FindSites_GatcEdgeCutter_WholeSequenceIsTheSite_ReportsNothing | `Sau3AI.search(Seq("GATC")) == []`, `MboI.search(Seq("GATC")) == []` | Biopython 1.88 |
| FindSites_Sau3AI_SiteAtThreePrimeEnd_ReportsNothing | `Sau3AI.search(Seq("AAAAAAAAGATC")) == []` | Biopython 1.88 |
| FindSites_EcoRI_SiteFillsWholeSequence_IsStillReported | `EcoRI.search(Seq("GAATTC")) == [2]` | Biopython 1.88 |
| FindSites_EcoRI_Circular_FindsOriginSpanningSite | `EcoRI.search(Seq("TTCAAAAAAAGAA"), linear=False) == [12]`, linear `== []` | Biopython 1.88 |
| FindSites_BamHI_Circular_FindsOriginSpanningSite | `BamHI.search(Seq("CCAAAAAAAAGGAT"), linear=False) == [12]`, linear `== []` | Biopython 1.88 |
| FindSites_Sau3AI_Circular_ReportsEveryGatc | `Sau3AI.search(Seq("GATCGATCGATC"), linear=False) == [1, 5, 9]` | Biopython 1.88 |
| FindSites_Circular_CutPositionIsReducedModuloLength | `Sau3AI.search(Seq("AAAAAAAAGATC"), linear=False) == [9]` | Biopython 1.88 |
| FindSites_Circular_MoleculeShorterThanSite_ReportsNothing | `ClaI.search(Seq("GATC"), linear=False) == [4]` — deliberate divergence | Biopython 1.88 |
| FindSites_CustomNonPalindromicEnzyme_CdiI_ReverseStrandCutMatchesReference | CdiI CATCG 4/4, `elucidate() == "CATC^_G"`, `ovhg == 0`, `search(Seq("AAACATCGAAACGATGAAA")) == [8, 13]` | REBASE 404 / Biopython 1.88 |
| FindSites_CustomNonPalindromicEnzyme_SimI_ReverseStrandCutMatchesReference | SimI GGGTC 2/5, `elucidate() == "GG^GTC_N"`, `ovhg == -3`, `search(Seq("AAAGGGTCAAAGACCCAAA")) == [6, 12]` | REBASE 404 / Biopython 1.88 |
| FindSites_HincII_ReverseRecognizedSequenceIsReadOnItsOwnStrand | `HincII.search(Seq("AAGTCAACAA")) == [6]` | Biopython 1.88 |
| FindPatternOccurrences_ReportsOverlappingOccurrences | `PacI.search(Seq("TTAATTAATTAA")) == [6, 10]` (sites at 0 and 4) | Biopython 1.88 |
| FindSites_SfiI_SequenceNInsideSpacer_MatchesBiopython | `SfiI.search(Seq("GGCCNNNNNGGCC")) == [9]`, linear and circular | Biopython 1.88 |
| FindSites_SfiI_SequenceN_TopStrandCutsMatchBiopython | `SfiI.search(Seq("AGGCCNNNNNGGCCA")) == [10]`, `SfiI.search(Seq("ggccnnnnnggcc")) == [9]`, `SfiI.search(Seq("NNNNNGGCCAGGCC"), linear=False) == [5]` | Biopython 1.88 |
| FindSites_SequenceAmbiguityAtConcreteSiteCode_StillBlocks | `EcoRI` on `AGAANTCA`, `HincII` on `AGTYRACA`/`AGTNNACA`, `SfiI` on `GGCNNNNNNGGCC`, `EcoRI` on `N`×16 — all `[]`, both topologies | Biopython 1.88 |
| SiteN_AcceptsEveryIupacCode_InFindSitesAndFindPatternOccurrences | `SfiI.search(Seq("GGCCRYSWKGGCC")) == [9]`; `Hpy188I.compsite` `TC.GA`, `Hpy8I` `GT..AC`, `AccI` `GT[AC][GT]AC` | Biopython 1.88 |
| FindSites_NonIupacCharacterInSpacer_MatchesNothing | `SfiI.search(Seq("GGCCXXXXXGGCC")) == [9]` — deliberate divergence | Biopython 1.88 |
| EnzymeDatabase_MethylationField_MatchesRebaseAllenz104 (×39) | field `<6>` of every table entry, strings extracted programmatically (34 with a site, 5 empty) | REBASE `allenz.104` |

### Actions Taken

| Action | Tests | Reason |
|--------|-------|--------|
| ❌ Implemented | `FindSites_HincII_MatchesAllDegenerateCombinations`, `*_RejectsNonMatchingDegenerateBases`, `*_DegenerateEnzymeProducesCorrectCutPosition` | C1 was missing — IUPAC degenerate recognition untested |
| ⚠ Strengthened | `GetEnzyme_BamHI_ReturnsEnzymeWithCorrectProperties` | Added Name, CutPositionForward/Reverse, IsBluntEnd (was only checking RecognitionSequence + Length + OverhangType) |
| ⚠ Strengthened | `FindSites_CustomEnzyme_FindsSiteCorrectly` | Hardcoded exact position (2) and cut position (4) instead of computing via IndexOf |
| 🔁 Removed | `RestrictionEnzyme_EcoRI_HasFivePrimeOverhang` | Fully covered by `GetEnzyme_EcoRI_ReturnsEnzymeWithCorrectProperties` (OverhangType + IsBluntEnd already asserted) |
| 🔁 Removed | `RestrictionEnzyme_NotI_IsEightCutter` | Fully covered by `EnzymeDatabase_CutPositions_MatchWikipedia("NotI", "GCGGCCGC", 2, 6)` |
| ❌ Implemented (2026-10-14 review) | the whole `Reference-Locked Site Search`, `CreateMap Site Reporting`, `FindPatternOccurrences` and `Enzyme Table Audit against REBASE 404` regions | M17-M32 and invariants #8-#16 were all missing; they lock the fixes F16-F21 |
| ⚠ Corrected (2026-10-14 review) | `CreateMap_NoEnzymesSpecified_SearchesAll` (RESTR-DIGEST-001 test file) | It asserted `NonCutters` is empty when no enzyme names are given, which restated an implementation artefact and contradicted RESTR-DIGEST-001 invariant #11; it now asserts the cutter/non-cutter partition |
| ❌ Implemented (B08 WP9 / D7, F52) | the `D7 — REBASE allenz.104 methylation-site field` region (50 cases) + 1 MCP test | the cognate-MTase site was absent although REBASE `allenz.104` field `<6>` is reachable |
| ❌ Implemented (B08 WP3 / D3, F40) | the `D3 — a site N matches a sequence N` region (12 cases) | a sequence `N` used to block every site, so `SfiI` on `GGCCNNNNNGGCC` reported nothing where Biopython reports `[9]`; the old `FindPatternOccurrences` doc comment saying ambiguity is read in the pattern only was corrected |

---

## Open Questions

None — all behavior verified against external sources.

---

## Deviations and Assumptions

Two deliberate divergences from the oracle, with their reasons:

- **Molecule shorter than the recognition site, circular.** Biopython pads a circular sequence with
  `data[1:size]` unconditionally, so a molecule shorter than the site is wrapped more than once and a
  site is reported anyway (`ClaI.search(Seq("GATC"), linear=False) == [4]` for the 6-bp site `ATCGAT`
  on a 4-bp circle). This library reports nothing, because a 6-bp recognition window cannot be present
  on a 4-bp circle. Locked by `FindSites_Circular_MoleculeShorterThanSite_ReportsNothing`.
- **Non-IUPAC character under a site `N` (D3).** Biopython compiles a site `N` to the regex `.`,
  so it matches any letter its `FormattedSeq` keeps (`A`-`Z`), e.g. `SfiI.search(Seq("GGCCXXXXXGGCC"))
  == [9]`. This library lets a site `N` match only the 15 IUPAC nucleotide codes: the raw-string
  overloads keep their documented contract that a non-IUPAC character never matches. Locked by
  `FindSites_NonIupacCharacterInSpacer_MatchesNothing`. Over the 15-code alphabet the two agree
  exactly (Evidence §2.5).

Two scope limits, not assumptions:

- ~~A `RestrictionEnzyme` carries one cut offset per strand **inside** the recognition window, so Type IIS /
  Type IIG enzymes that cut at a distance cannot be expressed.~~ **Withdrawn (review 2026-09 B08 WP11a,
  F54/F55):** the offsets are signed `int`s and were never restricted to the window; a **single-sided**
  distal cutter is a plain custom record — `CutPositionForward = fst5`, `CutPositionReverse = size + fst3`
  (BsaI GGTCTC(1/5) = `("GGTCTC", 7, 11)`, FokI `("GGATG", 14, 18)`, MlyI `("GAGTC", 10, 10)`). F54 fixed
  the circular reverse-strand cut, which could be negative; F55 made `Digest`/`PartialDigest` accept
  custom enzymes and cut at the reverse-strand sites of a non-palindromic enzyme. Differential vs
  Biopython 1.88 over all 153 single-sided distal REBASE enzymes (+ 576 in-site, as control):
  0 mismatches (Evidence RESTR-FIND-001 §2.6). **Two-sided** Type IIB cutters (BaeI, BcgI — Biopython
  `TwoCuts`, second pair `scd5`/`scd3`) are expressible since WP11b/F56 through the additive
  `RestrictionEnzyme.SecondCut = (scd5, size + scd3)`: `FindSites` reports one record per cut
  (`RestrictionSite.IsSecondCut`, `CutOffsets`), each cut kept or dropped on its own on a linear
  molecule; 0 mismatches vs Biopython over all 25 two-cut enzymes × both topologies (Evidence
  RESTR-FIND-001 §2.7; `TwoCut_*` tests in `RestrictionAnalyzer_DistalCutter_Tests`). The multi-enzyme
  `FindSites(DnaSequence, MoleculeTopology, IEnumerable<RestrictionEnzyme>)` takes custom records. The 39 table entries all cut inside their site.
- Cases: `RestrictionAnalyzer_DistalCutter_Tests` lock Biopython's `search` /
  `catalyse` values for BsaI, BsmBI, BbsI, SapI, FokI, MlyI, AcuI and the TspRI divergence (Biopython
  stores `ovhg = 10` for TscAI/TspRI against its own `fst5`/`fst3`/`ovhgseq` = 9, so it drops a site
  whose bottom-strand cut is the bond after base 0; this library keeps it).
- Sequence and pattern case folding on the raw-string entry points (`FindSites(string, …)`,
  `FindPatternOccurrences`, the pattern validator) is ASCII-only (`SequenceExtensions.ToUpperAscii`,
  review 2026-09 B08 F49, superseding F43's `ToUpperInvariant`): U+017F 'ſ' is NOT folded onto the IUPAC
  code 'S', so it matches nothing in a sequence and is rejected in a pattern, like any other non-IUPAC
  character. The `DnaSequence` overloads are unaffected (that type rejects non-ACGT input outright).

Previous items A1 ("Forward strand position is primary") and A2 ("Reverse complement sites reported
with forward strand coordinates") were reclassified as **standard conventions** — not assumptions.
Both are the universal convention in molecular biology (REBASE, GenBank, NCBI) and are directly
verified by tests `FindSites_PalindromicSite_FoundOnBothStrandsAtSamePosition` and
`EnzymeDatabase_CutPositions_MatchWikipedia`.
