# TestSpec: CRISPR-PAM-001

## Test Unit Information

| Field | Value |
|-------|-------|
| **ID** | CRISPR-PAM-001 |
| **Area** | Molecular Tools (MolTools) |
| **Algorithm** | PAM Site Detection |
| **Canonical Method** | `CrisprDesigner.FindPamSites(DnaSequence, CrisprSystemType)` |
| **Complexity** | O(n) |
| **Status** | ☑ Complete |

## Methods Under Test

| Method | Class | Type | Test Depth |
|--------|-------|------|------------|
| `FindPamSites(DnaSequence, CrisprSystemType)` | CrisprDesigner | Canonical | Deep |
| `FindPamSites(string, CrisprSystemType)` | CrisprDesigner | Overload | Deep |
| `GetSystem(CrisprSystemType)` | CrisprDesigner | Helper | Deep |

## Evidence Sources

| Source | Type | Used For |
|--------|------|----------|
| Wikipedia: Protospacer adjacent motif | Reference | PAM definitions, behavior |
| Wikipedia: CRISPR | Reference | System types, PAM sequences |
| Jinek et al. (2012), Science | Academic | SpCas9 NGG PAM |
| Zetsche et al. (2015), Cell | Academic | Cas12a TTTV PAM |
| Anders et al. (2014), Nature | Academic | PAM recognition mechanics |

## Test Categories

### MUST Tests (Required for Completion)

| ID | Test Case | Rationale | Source |
|----|-----------|-----------|--------|
| M1 | GetSystem returns correct PAM for SpCas9 (NGG) | Canonical PAM sequence | Wikipedia PAM |
| M2 | GetSystem returns correct PAM for SaCas9 (NNGRRT) | System-specific PAM | Wikipedia CRISPR |
| M3 | GetSystem returns correct PAM for Cas12a (TTTV) | Cas12a uses T-rich PAM | Wikipedia CRISPR, Zetsche 2015 |
| M4 | GetSystem returns correct guide lengths | System-specific guide lengths | Wikipedia CRISPR |
| M5 | GetSystem indicates PAM position (before/after target) | Cas9=after, Cas12a=before | Wikipedia PAM |
| M6 | FindPamSites finds NGG on forward strand | Core functionality | Wikipedia PAM |
| M7 | FindPamSites matches all NGG variants (AGG, CGG, TGG, GGG) | N = any nucleotide | IUPAC nomenclature |
| M8 | FindPamSites searches reverse strand | Both strands must be searched | Implementation spec |
| M9 | FindPamSites returns empty for sequences without PAM | No false positives | Core functionality |
| M10 | FindPamSites returns empty for empty sequence | Edge case | Core functionality |
| M11 | FindPamSites handles case-insensitive input | Usability requirement | Implementation spec |
| M12 | FindPamSites returns correct target sequence | Target must be extractable | Core functionality |
| M13 | FindPamSites Cas12a detects TTTA, TTTC, TTTG (V variants) | V = A, C, G not T | IUPAC nomenclature |
| M14 | FindPamSites does not detect TTTT for Cas12a | T is excluded from V | IUPAC nomenclature |
| M15 | FindPamSites SaCas9 matches NNGRRT pattern | R = A or G | IUPAC nomenclature |

### SHOULD Tests (Important but not blocking)

| ID | Test Case | Rationale | Source |
|----|-----------|-----------|--------|
| S1 | FindPamSites excludes sites where target would be out of bounds | Boundary safety | Implementation spec |
| S2 | FindPamSites returns correct position for reverse strand | Coordinate conversion | Implementation spec |
| S3 | FindPamSites handles overlapping PAM sites | Multiple valid sites | ASSUMPTION |
| S4 | GetSystem throws for unknown system type | Input validation | Implementation spec |
| S5 | FindPamSites with null DnaSequence throws ArgumentNullException | Input validation | .NET conventions |

### COULD Tests (Nice to have)

| ID | Test Case | Rationale | Source |
|----|-----------|-----------|--------|
| C1 | FindPamSites performance is linear in sequence length | O(n) complexity | Complexity spec |
| C2 | All seven CRISPR system types return valid configurations | Complete coverage | Implementation spec |

## Invariants

1. **PAM pattern matching**: A site is returned if and only if the PAM pattern matches using IUPAC rules
2. **Target extraction**: Target sequence length equals the system's guide length
3. **Strand correctness**: Forward strand sites have PAM at reported position; reverse strand sites are correctly transformed
4. **No self-targeting**: Empty input returns empty output (no crashes, no false positives)

## Audit of Existing Tests

### Current State (CrisprDesignerTests.cs)

| Test | Coverage | Quality | Action |
|------|----------|---------|--------|
| GetSystem_SpCas9_ReturnsCorrectSystem | M1, M4, M5 | Good | Keep |
| GetSystem_SaCas9_ReturnsCorrectSystem | M2, M4 | Good | Keep |
| GetSystem_Cas12a_ReturnsCorrectSystem | M3, M4, M5 | Good | Keep |
| FindPamSites_SpCas9_NGG_FindsSites | M6 | Weak (only AGG) | Enhance |
| FindPamSites_SpCas9_CGG_FindsSite | M7 partial | OK | Keep |
| FindPamSites_SpCas9_TGG_FindsSite | M7 partial | OK | Keep |
| FindPamSites_NoPam_ReturnsEmpty | M9 partial | Weak assertion | Enhance |
| FindPamSites_BothStrands_FindsSites | M8 | Weak | Enhance |
| FindPamSites_StringOverload_Works | M6 | OK | Keep |
| FindPamSites_EmptySequence_ReturnsEmpty | M10 | Good | Keep |
| FindPamSites_ReturnsTargetSequence | M12 partial | OK | Enhance |
| FindPamSites_Cas12a_TTTV_FindsSites | M13 partial | Weak (only TTTA) | Enhance |
| FindPamSites_SaCas9_NNGRRT_FindsSites | M15 partial | Weak | Enhance |

### Gaps Identified (All Closed)

- ~~M7: Need GGG variant test~~ ✅ Covered
- ~~M13: Need TTTC, TTTG variant tests~~ ✅ Covered
- ~~M14: Missing test for TTTT exclusion in Cas12a~~ ✅ Covered
- ~~M15: Need comprehensive NNGRRT pattern tests (R = A, G)~~ ✅ Covered
- ~~S1-S5: Most SHOULD tests missing~~ ✅ Covered

### Consolidation Plan

1. ~~**Create** `CrisprDesigner_PAM_Tests.cs`~~ ✅ Done
2. ~~**Extract** PAM-related tests from `CrisprDesignerTests.cs`~~ ✅ Done
3. ~~**Enhance** weak tests with stronger assertions~~ ✅ Done
4. ~~**Add** missing Must tests (M7, M13-M15)~~ ✅ Done
5. **Keep** remaining tests in `CrisprDesignerTests.cs` (guide RNA, off-target)

## Open Questions

None - all PAM sequences and behaviors are well-documented in literature.

## ASSUMPTIONS

| ID | Assumption | Justification |
|----|------------|---------------|
| A1 | Overlapping PAM sites are all reported | Implementation choice - not specified in sources |
