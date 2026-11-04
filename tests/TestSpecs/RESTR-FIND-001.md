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
| **EMBOSS `restrict -methylation -mfile Emethylsites.dat`** (`emboss/data/Emethylsites.dat`, `nucleus/embpat.c`; binary 6.6.0 executed) | Reference tool, source + execution | Methylase sites and the masking / window rule behind the opt-in methylation mode (E2 / F68); Evidence §7 |
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
19. **REBASE set (D4, F59)**: `RebaseEnzymes` holds the 754 Biopython 1.88 REBASE 404 enzymes with a defined cut (729 OneCut + 25 TwoCuts; the 334 NoCut/NotDefined are excluded and unresolvable); the 39 common enzymes are the same instances; `GetEnzyme` (and every name-taking member) resolves against it; every member that enumerates enzymes without a list or an `EnzymeSet` stays on the 39 common enzymes (Source: Biopython 1.88 `Restriction_Dictionary.py`, REBASE allenz.104; API contract)
20. **Readings vs. windows (F60)**: Biopython `search` == top-strand cuts of `FindSites`' forward records plus those reverse records whose window does not also read as the site on the forward strand (Source: Biopython `Palindromic._search`, `NonPalindromic._search` compsite alternation)
21. **Opt-in methylation (E2 / F68)**: with a methylase set (`MethylaseSets` = EMBOSS `Emethylsites.dat`, `m` before each methylated base) every base found at an `m` mark or at a site `N`, on either strand, is masked, and an enzyme window containing a masked base is not a site; entries are visited last-file-entry-first and each sees earlier masks (EMBOSS `patRestrictMethylMod`/`embPatRestrictScan`). `null`/empty set == the overload without it (byte-identical). Window rule, not a sensitivity table (BamHI is blocked by Dam). Circular: methylase sites spanning the origin are found (EMBOSS does not) (Source: EMBOSS `restrict -methylation`, executed; Evidence §7)

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
| M37 | All 754 REBASE enzymes, by name, linear + circular: `FindSites` cuts (invariant 20) = Biopython 1.88 `search`; `Digest` fragments = gaps of the sorted `search` cuts | D4: enzymes outside the 39 were "not in our table" | Biopython 1.88 (embedded oracle `biopython_rebase_search.tsv`) |
| M38 | Census 754 / 25 two-cut / 137 blunt / 388 5' / 229 3'; 39 common enzymes = same instances, vendored rows agree with them; `Dpi3084I` (no defined cut) unresolvable | Invariant 19 | Biopython 1.88 `is_blunt()`/`is_5overhang()`/`is_3overhang()` |
| M39 | Default-set members unchanged (39 / 24 six-cutters / 10 blunt / 29 sticky; no non-common enzyme in `FindAllSites`/`CreateMap`/`FindCompatibleEnzymes`); `EnzymeSet.Rebase` filters and schizomers = Biopython (385 six-cutters, 676 in 4–8, 137 blunt, 617 sticky; BsaI iso/equi = Bso31I, BspTNI, Eco31I; SmaI neo = Cfr9I, TspMI, XmaI) | Invariant 19 | Biopython 1.88 `.size`, `isoschizomers()`, `equischizomers()`, `neoschizomers()` |
| M40 | MCP `enzyme_set` (F61): default/`"common"` output byte-identical to the common-set library path (7 tools); `"rebase"`: 137 blunt / 617 sticky / 385 six-cutters, BsaI equi `Bso31I, BspTNI, Eco31I`, SmaI neo `Cfr9I, TspMI, XmaI`, enzymes with a site in `find_all_restriction_sites` (linear 141, circular 144) and `restriction_map` keys = Biopython `RestrictionBatch(754).search`, `compatible_enzymes("BamHI"/"EcoRI")` = `compatible_end` (59 / 51); unknown value / `"rebase"` without `enzyme_name` → `ArgumentException("enzyme_set")`; `enzymes_compatible("BamHI","BclI")` = (true, `GATC`) (crashed before: partner searched in the common set) | D4 part 2 | Biopython 1.88 (`EnzymeSetToolsTests`) |
| M41 | JSON of a two-cut enzyme carries `SecondCutForward = 23`, `SecondCutReverse = 18` (BaeI) and no `SecondCut`/`cutOffsets` empty object; one-cut enzymes unchanged | F62 (`"secondCut":{}` before) | Biopython 1.88 BaeI `scd5`, `size + scd3` |

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
| S1 | 4/6/8-cutters detected | `GetEnzymesByCutLength_ReturnsExactlyTheEnzymesOfThatLength` (×3 lengths, exact sets; B08 H54) |
| S2 | 30+ enzymes in database | `Enzymes_ContainsExactlyThe39CommonEnzymes` (B08 H54) |
| S3 | GetEnzymesByCutLength | `GetEnzymesByCutLength_ReturnsExactlyTheEnzymesOfThatLength` |
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
| INV-14 | Short molecule | `FindSites_Circular_MoleculeShorterThanSite_FindsTheWrappedSite` (circular: reported, as Biopython; B08 F100) |
| INV-17 / M33 | Sequence ambiguity (D3) | `FindSites_SfiI_SequenceNInsideSpacer_MatchesBiopython` (×2), `FindSites_SfiI_SequenceN_TopStrandCutsMatchBiopython` (×3), `FindSites_SequenceAmbiguityAtConcreteSiteCode_StillBlocks` (×5), `SiteN_AcceptsEveryIupacCode_InFindSitesAndFindPatternOccurrences`, `FindSites_NonIupacCharacterInSpacer_MatchesNothing`, plus `FindPatternOccurrences_AmbiguousSequenceBase_NeverMatches` (×2, unchanged) |
| INV-17 / M34 | ASCII-only fold (F49): 'ſ' (U+017F) never matches and a pattern containing it is rejected | `RawStringEntryPoints_LongS_IsNotFoldedOntoIupacS`, `Pattern_WithLongS_IsRejectedLikeAnyNonIupacCharacter` |
| INV-18 / M35-M36 | Cognate methylation site (D7 / F52): table locked to allenz.104 `<6>` (M35); notation parse + rejection of malformed / `(99)` values (M36) | `EnzymeDatabase_MethylationField_MatchesRebaseAllenz104` (×39), `EnzymeDatabase_MethylationAuditCoversEveryEntry`, `CognateMethylationSites_BuiltIns_ParsedToTypedSites`, `CognateMethylationSites_RebaseNotation_ComplementaryStrandAndLists`, `RestrictionEnzyme_MalformedMethylationField_Throws` (×7), `RestrictionEnzyme_NullMethylationField_Throws`; MCP `GetEnzyme_ReportsRebaseCognateMethylationSite` |
| INV-21 / M37-M42 | Opt-in methylation (E2 / F68): EMBOSS `restrict -methylation` oracle — L1 all Dam/Dcm windows blocked (MboI, Sau3AI, BamHI, BclI, EcoRII), L2 EcoRI 0,18 / SmaI 27 / TaqI 38,44 kept, L3 EcoKI spacer (KpnI/Acc65I 47 only), mask equality via the `NN` probe (6 sequences), circular cases, default-off identity | `RestrictionAnalyzer_Methylation_Tests` (42 cases) |
| M43 (F83, round 3 G3) | `PartialDigest` with a methylase set: F53 constructed set over the cuts that Biopython `search` + EMBOSS `restrict -methylation` leave — L2 Dam `{1,19,30,39,45}` (21 fragments, unmethylated 36), L1 Dam intact molecule, L3 emboss `{22,48,52}`, 500-bp random sequence none/emboss/Dam/Dcm, linear + circular, name + record forms, `MissedCuts = 0` prefix = methylation-aware `Digest`, default-off identity; MCP `partial_digest` `methylation` | `RestrictionAnalyzer_MethylationEntryPoints_Tests` (32 cases) + 4 MCP test methods (6 cases) |
| M44 (F84, round 3 G4) | Methylation mode on `FindAllSites(EnzymeSet)` (one mask per molecule), raw-string `FindSites`, `FindSites(…, enzymes, methylases)`, MCP `find_all_restriction_sites`: EMBOSS `restrict -methylation` over the 39 common enzymes for L1/L2/L3/500-bp x none/dam/dcm/emboss (16 rows, linear = circular), Rebase set = per-enzyme F68 path, raw string = DnaSequence form (case fold, empty, eager unknown enzyme / malformed methylase, `N` sequence), default-off identity | `RestrictionAnalyzer_MethylationEntryPoints_Tests` (G4 region) + 4 MCP test methods (6 cases) |

### External Source Verification Tests

| Test | Verified Against | Source |
|------|-----------------|--------|
| FindSites_EcoRI_OverhangSequenceIsAATT | Overhang = AATT (4 bp, 5') | Wikipedia: EcoRI |
| FindSites_BamHI_OverhangSequenceIsGATC | Overhang = GATC (4 bp, 5') | Wikipedia: BamHI |
| FindSites_HindIII_OverhangSequenceIsAGCT | Overhang = AGCT (4 bp, 5') | Wikipedia: HindIII |
| FindSites_PstI_ProducesThreePrimeOverhangOfFourBases | 3' cohesive termini, 4 bp | Wikipedia: PstI |
| FindSites_EcoRV_ProducesBluntEndCuts | Blunt ends (fwd == rev cut) | Wikipedia: EcoRV |
| FindSites_PalindromicSite_FoundOnBothStrandsAtSamePosition | Both strands recognized | Wikipedia: Restriction enzyme |
| EnzymeDatabase_CutPositions_MatchWikipedia (×21) | Recognition sequences and cut positions | Wikipedia: Restriction enzyme (Examples table) for 19; BglII and ApaI (not on that table) from Biopython 1.88 `Bio.Restriction` (B08 F64) |
| EnzymeDatabase_AllEntries_MatchRebase404 (×39) | site, fst5, len(site)+fst3, blunt/5'/3' class for every table entry — 39/39 matched, 0 wrong entries | REBASE 404 via Biopython 1.88 |
| FindSites_Sau3AI_EdgeSites_AreNotCuttable_MatchesBiopython | `Sau3AI.search(Seq("GATCGATCGATC")) == [5]` | Biopython 1.88 |
| FindSites_GatcEdgeCutter_WholeSequenceIsTheSite_ReportsNothing | `Sau3AI.search(Seq("GATC")) == []`, `MboI.search(Seq("GATC")) == []` | Biopython 1.88 |
| FindSites_Sau3AI_SiteAtThreePrimeEnd_ReportsNothing | `Sau3AI.search(Seq("AAAAAAAAGATC")) == []` | Biopython 1.88 |
| FindSites_EcoRI_SiteFillsWholeSequence_IsStillReported | `EcoRI.search(Seq("GAATTC")) == [2]` | Biopython 1.88 |
| FindSites_EcoRI_Circular_FindsOriginSpanningSite | `EcoRI.search(Seq("TTCAAAAAAAGAA"), linear=False) == [12]`, linear `== []` | Biopython 1.88 |
| FindSites_BamHI_Circular_FindsOriginSpanningSite | `BamHI.search(Seq("CCAAAAAAAAGGAT"), linear=False) == [12]`, linear `== []` | Biopython 1.88 |
| FindSites_Sau3AI_Circular_ReportsEveryGatc | `Sau3AI.search(Seq("GATCGATCGATC"), linear=False) == [1, 5, 9]` | Biopython 1.88 |
| FindSites_Circular_CutPositionIsReducedModuloLength | `Sau3AI.search(Seq("AAAAAAAAGATC"), linear=False) == [9]` | Biopython 1.88 |
| FindSites_Circular_MoleculeShorterThanSite_FindsTheWrappedSite | `ClaI.search(Seq("GATC"), linear=False) == [4]` (this library: forward record at 1, cut 3, 0-based) — the former deliberate divergence was withdrawn in F100 | Biopython 1.88 |
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
| ❌ Implemented (B08 round 2 E2, F68) | the `RestrictionAnalyzer_Methylation_Tests` file (42 cases) + 6 MCP tests | the recorded proof "no source tabulates sensitivity" was defeated: EMBOSS `restrict -methylation -mfile Emethylsites.dat` is driven by methylase sites |

---

## Open Questions

None — all behavior verified against external sources.

---

## Deviations and Assumptions

One deliberate divergence from the oracle remains (the circular short-molecule divergence of earlier rounds was
withdrawn in B08 F100: this library reports the wrapped site exactly as Biopython does, locked by
`FindSites_Circular_MoleculeShorterThanSite_FindsTheWrappedSite`), with its reason:

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

## Round 5 (B08 A5R-2, A5R-3, A5R-5 .. A5R-8, mutation survivors)

Tests in `RestrictionAnalyzer_Round5_Tests.cs` (library), `RestrictionWrapperRound5Tests.cs` and `B07ToolDocsContractTests` (MCP).
Each failed on the audited code (32 of the 45 library tests; the rest kill the auditor's surviving mutants, re-verified by mutation here).

| ID | Test | Pins | Oracle |
|---|---|---|---|
| A5R-2 | `FindPatternOccurrences_HugeStartIndex_IsAnEmptyResult_NotAnOverflow` (3), `..._StartIndexAtTheLastWindow_StillFindsIt` | a start index near `int.MaxValue` = empty, not `IndexOutOfRange`; the last window still found | the method's own documented contract |
| A5R-3 | `EnzymeNameArray_NullIsAMissingArgument_EmptyIsAWrongValue_UnknownAndNullNamesAreRefusedAtTheCall` (13 forms incl. the methylation `FindSites`/`Digest`/`PartialDigest`/`GetDigestSummary`) | null array = `ArgumentNullException`, empty array / empty or unknown name = `ArgumentException`, null element = `ArgumentNullException`, all eager | BCL `ThrowIfNullOrEmpty` split (B08 H19/H24) |
| A5R-5 | `UndefinedTopology_IsRefusedAtTheCall` (18 members; superseded in round 6 by the reflection test over all 41), `RestrictionTools_RefuseAnUndefinedTopology` (8 tools) | `ArgumentOutOfRangeException("topology")` everywhere | enum validation as the sibling wrappers do |
| A5R-6 | `RawStringOverloads_RejectAMalformedMethylase_ForNullAndEmptySequencesToo` | the documented `ArgumentException` whatever the sequence | the sibling overload that always validated |
| A5R-7 | `CustomEnzyme_HugeCutOffsets_OnACircle_ReduceInLongArithmetic` | forward cut 10 and reverse cut 2 for offsets of `int.MaxValue` on a 12-mer | exact arithmetic: (3 + 2^31 - 1) mod 12 = 10; (9 - (2^31 - 1)) mod 12 = 2 |
| A5R-8 | `RawStringOutputs_AreUpperCase_ForASoftMaskedMolecule` | raw-string outputs are upper-case, coordinates unaffected (documented in the tool descriptions) | ASCII fold of `Normalize` |
| T1 | `CreateMap_OneTypeIibSite_IsOneCleavage_AndAUniqueCutter` | one AlfI site: `TotalSites 1`, `UniqueCutters [AlfI]` | F108 rule |
| T2/T3 | `Methylation_BlocksAWindowThatCrossesTheCircularOrigin` | Dam blocks a GATC that spans the origin, in the DnaSequence and the raw-string `FindAllSites` forms | EMBOSS `restrict -methylation` |
| T4 | `RawGetDigestSummary_EnzymesUsed_AreTheNamesAsSupplied` | `["ecori"]`, not `["EcoRI"]` | `DigestSummary.EnzymesUsed` doc |
| T5 | `DigestFragment_Equality_ComparesTheEnzymeListsByContent` | element-wise equality | the record's doc |
| T6/T7 | `RestrictionDigestCustomEnzyme_ForwardsSecondCut_AndMethylation`, `FindSitesCustomEnzyme_ForwardsMethylation_AndTopology`, `CompatibleEnzymes_ReportsBothEnds_ForARevcompOnlyPair_InThePerEnzymeAndTheAllPairsForm` | every argument of the custom tools reaches the algorithm; `end2` in both forms | BaeI geometry (10/33/8 vs 10/41); Biopython compatible_end |
| A5R-1 | `EnzymesCompatible_AcceptsEverySpellingTheLibraryResolves` (5) | trimmed / hyphen / case spellings of the partner | F103 |
| A5R-4 | `CustomEnzymeTools_RefuseANullOrEmptySequence` | same guard as the named wrappers | `find_restriction_sites("")` |
| T8 | `PartialDigest_CountExactlyAtTheLimit_IsAccepted_OneStepOverIsRefused` | limit inclusive: 1 000 cuts x (99 + 1) = 100 000 accepted, x 101 refused | `MaxPartialDigestFragments` |
| docs | `McpJson_Examples_ReproduceTheTool` (25 tools), `McpJson_InputSchema_...`, `McpJson_OutputSchema_...` | every restriction doc example is the tool's real output (prose summaries skipped); schemas equal the wrapper | the tool itself |

## Round 6 (B08 A6-2 .. A6-5, A6-7, U3)

Tests in `RestrictionAnalyzer_Round6_Tests.cs`, `IupacPattern_Tests.cs` (library) and `RestrictionCrisprRound6WrapperTests.cs` (MCP); heavy tier in
`B08RestrictionAuditFuzzTests.cs` / `B08RestrictionAuditProperties.cs` (rotation equivariance of the circular sites, digest and map; partial-digest counts
(k+1)(k+2)/2 and k^2+1; IUPAC symbols in the sequence; methylation only blocks; undefined topology and the name-array contract over every member).
27 of the 59 library tests and 15 of the 18 MCP tests fail on the code the auditor read (the rest are controls, or the reflection tests that kill the
auditor's guard-deletion mutants: 11 of 11 killed, and the `Modulo` int mutant).

| ID | Test | Pins | Oracle |
|---|---|---|---|
| A6-2 | `RawSequence_AUracilInsideASite_IsNotASite_OnEitherStrand` (5), `..._AUracil_IsAlsoNotReadAcrossTheCircularOrigin`, `Methylation_AUracil_DoesNotCreateAMethylaseContextOnTheReverseStrand`, `Overhang_OfASiteWhoseWindowHoldsAUracil_KeepsTheUracil`, `IupacPattern_Tests.ReverseComplement_*` (9), MCP `FindRestrictionSites_AUracilInTheSite_IsNoSite` | a `U` is outside the 15 IUPAC codes: it is part of no site on either strand, no methylase context and no complemented overhang base | the documented rule "a character outside the 15 codes matches nothing"; the forward scan and `FindPatternOccurrences` |
| A6-3 | `CustomEnzyme_StaggerBeyondTheLimit_IsRefusedAtTheCall_ByEveryMember` (8 offset pairs x 11 members), `..._StaggerAtTheLimit_IsAccepted_EvenAtTheEdgesOfTheIntRange`, `EveryRebaseEnzyme_IsWithinTheLimit_AndTheLargestStaggerIsNine`, `SecondCutStagger_IsComparedInLong_NotWrappedInInt`, MCP `RestrictionDigestCustomEnzyme_StaggerBeyondTheLimit_IsRefused` | `|forward - reverse| <= MaxCutStagger` (64), computed in `long`, refused eagerly | REBASE: largest of the 754 enzymes is 9 (TspRI, TscAI); exact arithmetic |
| A6-4 | `EveryPublicMemberThatTakesATopology_RefusesAnUndefinedOne_AtTheCall` (41 members found by reflection), `..._StillWorks_ForTheDefinedOnes`, `EveryParamsEnzymeNameMember_SharesOneContract_AtTheCall` (21 members x 6 shapes) | no member can be added or lose its guard unnoticed | the 11 guards the auditor deleted without a test failing |
| A6-5 | `CustomEnzyme_MinValueOffsets_OnACircle_ReduceInLongArithmetic` | forward cut 7 and reverse cut 5 for `int.MinValue` offsets on a 12-mer | exact arithmetic: (3 + int.MinValue) mod 12 = 7; (9 - int.MinValue) mod 12 = 5 |
| A6-7 | `IupacPattern_Tests.MatchesAt_ReferenceRows` (int.MaxValue rows) | a window start near `int.MaxValue` does not fit: false, no exception | the method's documented contract |
| U3 | MCP `CompatibleEnzymes_PerEnzymeForm_EchoesTheCanonicalName` (3), `..._ResolvesAHyphenSpelling_ToTheTableName` | `enzyme1` carries the table spelling, as `enzymes_compatible` does since A5R-1 | `get_enzyme` |

## Round 7 (B08 A7-2 .. A7-5, O7-1)

Tests in `RestrictionAnalyzer_Round7_Tests.cs` (library) and `RestrictionCrisprRound7WrapperTests.cs` (MCP); heavy tier: `U` added to the junk list of `B08RestrictionAuditProperties.NRegion_AcceptsEveryIupacSymbol_AndOnlyThose`.
11 of the 16 library tests and 2 of the 3 MCP tests fail on the code the auditor read (the new constants stubbed in, because the tests name them); the 3' distal overhang test is a mutation row
(reverting the 3' branch of `DownstreamEnd` to the RNA-aware complement gives `CGAT`).

| ID | Test | Pins | Oracle |
|---|---|---|---|
| A7-2 | `CircularDigest_CustomEnzymeAtTheIntExtremes_CutsAtBothReadings` (3 offset pairs, DnaSequence and raw overload), `CircularMap_AndPartialDigest_..._CountTwoCleavages` (3) | `GAATTC` at 3 of a 12-mer, offsets `(Min+4, Min+2)`, `(Min+3, Min+3)`, `(Max-3, Max-1)`: two cleavages, fragments `3+8 / 11+4`, `2+8 / 10+4`, `3+4 / 7+8`; `TotalSites 2`, no unique cutter, 5 partial fragments | exact arithmetic with 2^31 = 8 (mod 12), checked outside the library; the reverse reading's top cut is `3 + 6 - rev` |
| A7-3 | `CustomEnzyme_RecognitionSequenceBeyondTheLimit_...` (65 and 200 000 `N`), `CustomEnzyme_NameBeyondTheLimit_...`, `CustomEnzyme_RecognitionSequenceAndNameAtTheLimit_AreAccepted`, `EnzymeLimits_ArePinnedLiterally`, MCP `CustomEnzymeTools_*` (3) | `MaxRecognitionLength` = `MaxEnzymeNameLength` = `MaxCutStagger` = 64; refused at the call by every member that takes an enzyme; an `N` pattern of 64 over 100 bases gives 2 x 37 records | REBASE: longest site 15, longest name 9, largest stagger 9 |
| A7-4 | `Digest_ThreePrimeDistalOverhangWindowHoldingAUracil_KeepsTheUracil` | 3' overhang window `AUCG` read as `AUCG` / `CGUT` (reverse complement, the U kept) | hand-derived reverse complement |
| O7-1 | `GetOverhangSequence_AndAreCompatible_Refuse...` (4), `..._ReadALowerCaseRecognitionSequenceAsUpperCase` | an empty or non-IUPAC recognition sequence is an `ArgumentException` here too; `gaattc` gives `AATT` and is compatible with EcoRI | `FindSites` folds and validates the same pattern |

## Round 8 (B08 A8-1, A8-2, A8-3, A8-7, A8-8)

Tests in `RestrictionAnalyzer_Round8_Tests.cs` (library) and `RestrictionCrisprRound7WrapperTests.cs` (MCP, the `FindPatternOccurrences_*` pair). 5 of the 13 library tests and 1 of the 2 MCP tests fail on the code the auditor read
(A8-7, A8-8); the others pin mutants of the round-7 change (A8-1, A8-2, A8-3), re-killed by mutation here.

| ID | Test | Pins | Oracle |
|---|---|---|---|
| A8-1 | `CircularDigestAndMap_TwoSidedCustomEnzymeAtTheIntExtremes_CutAtTheUnionOfBothReadings` | `GAATTC` at 3 of a 12-mer, cuts `(5, 1)` + second pair `(Min+5, Min+1)`: top cuts {0, 4, 8}, fragments `0+4 AAAG / 4+4 AATT / 8+4 CAAA`, `TotalSites 2` | exact arithmetic, 2^31 = 8 (mod 12): forward 8 and 12 = 0; reverse 9 - 1 = 8 and 9 - (Min+1) = 16 = 4 |
| A8-2 | `AreCompatible_RefusesABadDefinition_EvenWhenTheEndsDecideTheAnswerFirst` | a non-IUPAC site is an `ArgumentException` for blunt x blunt, 5' x 3' and blunt x 5' (it returned True / False) | the definition check of `FindSites` |
| A8-3 | `EnzymeLimits_AreCheckedBeforeThePattern_AndBlameTheCallersParameter` | the size limit is the reported reason for an oversized non-IUPAC site; `enzyme` / `enzymes` / `enzyme1` / `enzyme2` are blamed; the name limit precedes the pattern | the documented order |
| A8-7 | `FindPatternOccurrences_PatternBeyondTheLimit_IsRefusedAtTheCall`, `..._PatternAtTheLimit_Works`, MCP `FindPatternOccurrences_PatternBeyondTheLimit_IsRefusedAtOnce` / `..._PatternAtTheLimit_Works` | a pattern over `MaxRecognitionLength` (64) is refused, 64 works (an `N` pattern over 100 bases: 37 hits); the scan costs n x m | measured 0.07 / 0.20 / 0.63 / 2.50 s for n = 10k .. 80k |
| A8-8 | `CircularMap_OffsetsCongruentModuloTheLength_CountOneCleavage_AsTheDigestDoes` (5), `CircularMap_AnEnzymeWhoseReverseReadingCutsOtherBonds_StillCountsTwoCleavages` (2), `LinearMap_OffsetsAreAbsolute_ThereIsNoModulo` | `(1,5)` `(13,17)` `(-11,-7)` `(25,29)` `(1200001,1200005)` on a 12-mer all cut the same two bonds: one cleavage in the map, one cut in the digest; `(1,3)` and `(13,15)` differ and count 2; linear stays absolute | the documented rule "a reverse record counts unless its backwards reading cuts the same bonds" |

## Round 9 (B08 A9-1, A9-2, A9-3)

Tests in `RestrictionAnalyzer_Round9_Tests.cs`. No defect was found in the round-8 code; every test kills a mutant the auditor built (and the order test also kills the round-8 survivor M26).

| ID | Test | Pins | Oracle |
|---|---|---|---|
| A9-1 | `LinearMoleculeShorterThanTheOffsets_TheReverseReadingIsItsOwnCleavage` | `GAATTC` at 0 of `GAATTCAAAAAA`, offsets `(-3,-3)`: the forward cut -3 is outside, the reverse reading cuts at 9; digest `0+9 GAATTCAAA / 9+3 AAA`, `TotalSites 1`, unique cutter | linear cuts are absolute: -3 and 9 are congruent modulo 12 but different places |
| A9-2 | `CircularTwoSidedEnzyme_ReadingsThatCutTheSameBondsModuloTheLength_AreOneCleavage` (3 rows) | `(1,5)+(7,11)`, `(13,17)+(7,11)`, `(1,3)+(15,17)`: one cleavage, fragments `4+6 / 10+6` and `4+2 / 6+10` | forward top cuts 3 + f, reverse top cuts 9 - r, modulo 12 |
| A9-3 | `FindPatternOccurrences_OversizedPattern_IsRefusedBeforeTheAlphabetIsScanned`, `AreCompatible_NameLimitIsReported_BeforeABadPattern` | a 2 001-character pattern: the length is the reason and the message is not an echo (under 300 characters); the name limit precedes the pattern in both argument positions | the documented order |

