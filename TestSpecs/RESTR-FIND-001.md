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
| **Last Updated** | 2026-01-23 |

---

## Evidence Sources

| Source | Type | Key Information |
|--------|------|-----------------|
| Wikipedia: Restriction enzyme | Academic | Type II enzymes, palindromic recognition, cut positions, overhang types |
| Wikipedia: Restriction site | Academic | Recognition sequences 4-8 bp, sticky vs blunt ends |
| Wikipedia: EcoRI | Academic | GAATTC recognition, 5' overhang AATT, cut position G↓AATTC |
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
| C1 | IUPAC ambiguous recognition sequences work | Degenerate sequences | Implementation |
| C2 | Both strands searched for palindromic sites | Biological correctness | Wikipedia |

---

## Audit Results

### Existing Test Coverage (RestrictionAnalyzerTests.cs)

| Test | Status | Notes |
|------|--------|-------|
| GetEnzyme_EcoRI_ReturnsCorrectEnzyme | Covered | M1 partially |
| GetEnzyme_CaseInsensitive_Works | Covered | M5 |
| GetEnzyme_UnknownEnzyme_ReturnsNull | Covered | M6 |
| FindSites_EcoRI_FindsSite | Covered | M1 partially |
| FindSites_MultipleSites_FindsAll | Covered | M3 |
| FindSites_NoSites_ReturnsEmpty | Covered | M4 |
| FindSites_UnknownEnzyme_ThrowsException | Covered | M7 |
| FindSites_MultipleEnzymes_FindsAllSites | Covered | M12 |
| FindAllSites_FindsMultipleEnzymes | Covered | M11 |

### Missing Tests

| ID | Test Case | Priority |
|----|-----------|----------|
| M8 | Empty sequence handling | Must |
| M9 | CutPosition calculation verification | Must |
| M10 | RecognizedSequence verification | Must |
| M13 | Custom enzyme support | Must |

### Weak Tests

| Test | Issue | Fix |
|------|-------|-----|
| FindSites_EcoRI_FindsSite | Only checks Any() | Verify exact position |
| GetEnzyme_EcoRI_ReturnsCorrectEnzyme | No cut position verification | Add invariant assertions |

---

## Consolidation Plan

1. **Canonical File**: Create `RestrictionAnalyzer_FindSites_Tests.cs` with FindSites/FindAllSites/GetEnzyme tests
2. **Keep in RestrictionAnalyzerTests.cs**: Digest, Map, Compatibility tests (for RESTR-DIGEST-001)
3. **Move to Canonical**: All FindSites and GetEnzyme tests
4. **Add Missing Tests**: M8, M9, M10, M13, S1-S5
5. **Strengthen Tests**: Use Assert.Multiple for invariant grouping

---

## Open Questions

None - behavior is well-documented in Wikipedia and REBASE.

---

## Assumptions

| ID | Assumption | Justification |
|----|------------|---------------|
| A1 | Forward strand position is primary | Convention in molecular biology |
| A2 | Reverse complement sites reported with forward strand coordinates | User convenience |
