# Test Specification: RESTR-FIND-001

## Test Unit Information

| Field | Value |
|-------|-------|
| **Test Unit ID** | RESTR-FIND-001 |
| **Area** | MolTools |
| **Title** | Restriction Site Detection |
| **Canonical Class** | `RestrictionAnalyzer` |
| **Canonical Methods** | `FindSites`, `FindAllSites`, `GetEnzyme` |
| **Complexity** | O(n × k × m) |
| **Status** | ☑ Complete |
| **Last Updated** | 2026-03-04 |

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

---

## Invariants

1. **Position Range**: 0 ≤ Position ≤ sequence.Length - recognitionLength (Source: Implementation)
2. **Cut Position**: Position ≤ CutPosition ≤ Position + recognitionLength (Source: Implementation)
3. **Recognition Match**: RecognizedSequence.Length == Enzyme.RecognitionSequence.Length (Source: Implementation)
4. **Case Insensitivity**: GetEnzyme("EcoRI") == GetEnzyme("ecori") (Source: Implementation)
5. **Empty Returns Empty**: Empty sequence returns no sites (Source: Implementation)
6. **Unknown Throws**: Unknown enzyme name throws ArgumentException (Source: Implementation)
7. **Both Strands Searched**: Forward and reverse complement both searched (Source: Biological requirement)

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

### Actions Taken

| Action | Tests | Reason |
|--------|-------|--------|
| ❌ Implemented | `FindSites_HincII_MatchesAllDegenerateCombinations`, `*_RejectsNonMatchingDegenerateBases`, `*_DegenerateEnzymeProducesCorrectCutPosition` | C1 was missing — IUPAC degenerate recognition untested |
| ⚠ Strengthened | `GetEnzyme_BamHI_ReturnsEnzymeWithCorrectProperties` | Added Name, CutPositionForward/Reverse, IsBluntEnd (was only checking RecognitionSequence + Length + OverhangType) |
| ⚠ Strengthened | `FindSites_CustomEnzyme_FindsSiteCorrectly` | Hardcoded exact position (2) and cut position (4) instead of computing via IndexOf |
| 🔁 Removed | `RestrictionEnzyme_EcoRI_HasFivePrimeOverhang` | Fully covered by `GetEnzyme_EcoRI_ReturnsEnzymeWithCorrectProperties` (OverhangType + IsBluntEnd already asserted) |
| 🔁 Removed | `RestrictionEnzyme_NotI_IsEightCutter` | Fully covered by `EnzymeDatabase_CutPositions_MatchWikipedia("NotI", "GCGGCCGC", 2, 6)` |

---

## Open Questions

None — all behavior verified against external sources.

---

## Deviations and Assumptions

None.

Previous items A1 ("Forward strand position is primary") and A2 ("Reverse complement sites reported
with forward strand coordinates") were reclassified as **standard conventions** — not assumptions.
Both are the universal convention in molecular biology (REBASE, GenBank, NCBI) and are directly
verified by tests `FindSites_PalindromicSite_FoundOnBothStrandsAtSamePosition` and
`EnzymeDatabase_CutPositions_MatchWikipedia`.
