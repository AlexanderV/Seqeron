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
| Wikipedia: Protospacer adjacent motif | Reference | PAM definitions, canonical NGG, PAM position (3'/5') |
| Wikipedia: CRISPR | Reference | System types overview |
| Wikipedia: Cas9 | Reference | SpCas9 NGG 20 nt spacer, NAG as tolerated secondary PAM |
| Wikipedia: Cas12a | Reference | TTTV PAM (V=A,C,G), 5' PAM before target, Cpf1 naming |
| Jinek et al. (2012), Science 337:816 | Academic | SpCas9 NGG PAM, 20 nt guide RNA |
| Zetsche et al. (2015), Cell 163:759 | Academic | Cas12a/Cpf1 TTTV PAM, ~23 nt spacer, AsCas12a/LbCas12a variants |
| Anders et al. (2014), Nature 513:569 | Academic | Structural basis of PAM-dependent target DNA recognition |
| Ran et al. (2015), Nature 520:186 | Academic | SaCas9 NNGRRT PAM, 21 nt guide (20-24 nt range) |
| Hsu et al. (2013), Nat Biotechnol 31:827 | Academic | NAG as secondary SpCas9 PAM with lower activity |
| Liu et al. (2019), Nature 566:218 | Academic | CasX (Cas12e) TTCN PAM, 20 nt guide |

## Test Categories

### MUST Tests (Required for Completion)

| ID | Test Case | Rationale | Source |
|----|-----------|-----------|--------|
| M1-M3 | GetSystem returns correct Name and PAM for all 7 systems | Canonical PAM sequences | Jinek 2012, Ran 2015, Zetsche 2015, Hsu 2013, Liu 2019 |
| M4 | GetSystem returns correct guide lengths for all 7 systems | System-specific guide lengths | Jinek 2012, Ran 2015, Zetsche 2015, Liu 2019 |
| M5 | GetSystem indicates PAM position (before/after target) for all 7 systems | Cas9=3' after, Cas12a/CasX=5' before | Wikipedia PAM, Zetsche 2015, Liu 2019 |
| M6 | FindPamSites finds NGG on forward strand with correct position and target start | Core detection with position verification | Wikipedia PAM |
| M7 | FindPamSites matches all NGG variants (AGG, CGG, TGG, GGG) with position | N = any nucleotide per IUPAC; position verified | IUPAC nomenclature |
| M8 | FindPamSites searches reverse strand with exact PAM, position, and count | Both strands searched; PAM reverse-complemented; exact coordinates | Biology: guide RNA targets both strands |
| M9 | FindPamSites returns no sites on either strand for PAM-free sequence | No false positives on forward or reverse | Core functionality |
| M10 | FindPamSites returns empty for empty sequence | Edge case | Core functionality |
| M11 | FindPamSites handles case-insensitive input | Usability requirement | Implementation spec |
| M12 | FindPamSites returns target with correct length AND content | Target extraction verified by value, not just length | Core functionality |
| M13 | FindPamSites Cas12a detects TTTA, TTTC, TTTG with position, target start, and content | V = A, C, G not T | IUPAC, Zetsche 2015 |
| M14 | FindPamSites does not detect TTTT for Cas12a | T is excluded from IUPAC V | IUPAC nomenclature |
| M15 | FindPamSites SaCas9 detects NNGRRT with R=A, R=G, and mixed R; verified by PAM content | R = A or G; negative test with C at R position | Ran 2015, IUPAC |
| M16 | SpCas9_NAG detects NAG pattern with all N variants with position verification | Secondary SpCas9 PAM | Hsu 2013, Wikipedia Cas9 |
| M17 | CasX detects TTCN with all N variants (TTCA, TTCC, TTCG, TTCT) with position and target | CasX/Cas12e 5' PAM | Liu et al. 2019 |

### SHOULD Tests (Important but not blocking)

| ID | Test Case | Rationale | Source |
|----|-----------|-----------|--------|
| S1 | FindPamSites excludes sites where target would be out of bounds (Cas9 and Cas12a) | Boundary safety for both PAM-before and PAM-after systems | Implementation spec |
| S2 | FindPamSites returns correct forward-strand coordinates for reverse strand hits | Coordinate conversion: revcomp position → forward position | Implementation spec |
| S3 | FindPamSites handles overlapping PAM sites — exact count and positions verified | Each PAM position defines a unique guide RNA | Wikipedia PAM, Hsu 2013 |
| S4 | GetSystem throws ArgumentException for unknown system type | Input validation | .NET conventions |
| S5 | FindPamSites with null DnaSequence throws ArgumentNullException | Input validation | .NET conventions |
| S6 | AsCas12a FindPamSites produces 23bp target end-to-end | Variant-specific guide length | Zetsche 2015 |
| S7 | LbCas12a FindPamSites produces 24bp target end-to-end | Variant-specific guide length | Zetsche 2015 |

### COULD Tests (Nice to have)

| ID | Test Case | Rationale | Source |
|----|-----------|-----------|--------|
| C1 | FindPamSites performance is linear in sequence length | O(n) complexity | Complexity spec |

## Invariants

1. **PAM pattern matching**: A site is returned if and only if the PAM pattern matches using IUPAC rules
2. **Target extraction**: Target sequence length equals the system's guide length; content matches the genomic region
3. **Strand correctness**: Forward strand sites have PAM at reported position; reverse strand positions are converted to forward-strand coordinates
4. **No self-targeting**: Empty input returns empty output (no crashes, no false positives)

## Coverage Classification

**Total: 58 tests in canonical file (was 56 before classification — 7 duplicates removed, 2 missing tests added, 7 parametric cases added)**

### Summary

| Category | Count |
|----------|-------|
| ❌ Missing → Implemented | 3 (AsCas12a PAM name, LbCas12a PAM name via parametric; AsCas12a FindPamSites S6; LbCas12a FindPamSites S7) |
| ⚠ Weak → Strengthened | 4 (M7, M8, M16 variants, S3) |
| 🔁 Duplicate → Removed | 7 (M1/M2/M3/SpCas9NAG PAM/CasX PAM → parametric; C2 subsumed; PamBeforeTarget → merged into DetectsTTTA) |
| ✅ Covered | All remaining |

### Classification Detail

#### ❌ Missing (Implemented)

| ID | Test | Action |
|----|------|--------|
| M1-M3 | AsCas12a/LbCas12a PAM + Name not verified | Added to parametric `GetSystem_ReturnsCorrectNameAndPam` (7 cases for all systems) |
| S6 | AsCas12a not exercised through FindPamSites | Implemented `FindPamSites_AsCas12a_DetectsTTTV_With23bpTarget` — verifies position, target start, length=23, content |
| S7 | LbCas12a not exercised through FindPamSites | Implemented `FindPamSites_LbCas12a_DetectsTTTV_With24bpTarget` — verifies position, target start, length=24, content |

#### ⚠ Weak (Strengthened)

| ID | Test | Before | After |
|----|------|--------|-------|
| M7 | MatchesAllNGG_Variants (4 cases) | `Is.Not.Null` only | `Is.Not.Null` + `Position == 20` |
| M8 | SearchesReverseStrand | `Is.Not.Empty` + `Position >= 0` | `Count == 1` + `Position == 0` + `PamSequence == "CCA"` |
| M16 | MatchesAllNAG_Variants (4 cases) | `Is.Not.Null` only | `Is.Not.Null` + `Position == 20` |
| S3 | OverlappingPamSites_AllReported | `Does.Contain` + `Count >= 2` | `Count == 2` + `[0].Position == 20` + `[1].Position == 23` |

#### 🔁 Duplicate (Removed)

| Test | Reason |
|------|--------|
| `GetSystem_SpCas9_ReturnsNGG_Pam` | Absorbed into parametric `GetSystem_ReturnsCorrectNameAndPam` SpCas9 case |
| `GetSystem_SaCas9_ReturnsNNGRRT_Pam` | Absorbed into parametric SaCas9 case |
| `GetSystem_Cas12a_ReturnsTTTV_Pam` | Absorbed into parametric Cas12a case |
| `GetSystem_SpCas9NAG_ReturnsNAG_Pam` | Absorbed into parametric SpCas9_NAG case; GuideLength/PamAfterTarget assertions redundant with M4+M5 |
| `GetSystem_CasX_ReturnsTTCN_Pam` | Absorbed into parametric CasX case; GuideLength/PamAfterTarget assertions redundant with M4+M5 |
| `GetSystem_AllSystemTypes_ReturnValidConfigurations` (C2) | Fully subsumed by parametric PAM test (Name+PAM exact) + M4 (GuideLength exact) + M5 (PamAfterTarget exact) |
| `FindPamSites_Cas12a_PamBeforeTarget` | Same input sequence and PAM as `DetectsTTTA`; assertions merged into M13 TTTA test |

### Canonical File (`CrisprDesigner_PAM_Tests.cs`) — 58 tests

| # | Test Method | Spec ID | Status |
|---|-------------|---------|--------|
| 1 | `GetSystem_ReturnsCorrectNameAndPam` (7 cases) | M1-M3 | ✅ |
| 2 | `GetSystem_ReturnsCorrectGuideLength` (7 cases) | M4 | ✅ |
| 3 | `GetSystem_ReturnsCorrectPamPosition` (7 cases) | M5 | ✅ |
| 4 | `GetSystem_UnknownType_ThrowsArgumentException` | S4 | ✅ |
| 5 | `FindPamSites_SpCas9_DetectsNGG_OnForwardStrand` | M6 | ✅ |
| 6 | `FindPamSites_SpCas9_MatchesAllNGG_Variants` (4 cases) | M7 | ✅ |
| 7 | `FindPamSites_SpCas9_SearchesReverseStrand` | M8 | ✅ |
| 8 | `FindPamSites_SpCas9_ReverseStrand_PositionConversion` | S2 | ✅ |
| 9 | `FindPamSites_SpCas9_CaseInsensitive` | M11 | ✅ |
| 10 | `FindPamSites_SpCas9_ReturnsTargetSequence_WithCorrectContent` | M12 | ✅ |
| 11 | `FindPamSites_NoPamPresent_ReturnsEmpty` | M9 | ✅ |
| 12 | `FindPamSites_EmptySequence_ReturnsEmpty` | M10 | ✅ |
| 13 | `FindPamSites_TargetOutOfBounds_Excluded` | S1 | ✅ |
| 14 | `FindPamSites_Cas12a_TargetOutOfBounds_Excluded` | S1 | ✅ |
| 15 | `FindPamSites_NullSequence_ThrowsArgumentNullException` | S5 | ✅ |
| 16 | `FindPamSites_Cas12a_DetectsTTTA` | M13 | ✅ |
| 17 | `FindPamSites_Cas12a_DetectsTTTC` | M13 | ✅ |
| 18 | `FindPamSites_Cas12a_DetectsTTTG` | M13 | ✅ |
| 19 | `FindPamSites_Cas12a_DoesNotDetectTTTT` | M14 | ✅ |
| 20 | `FindPamSites_AsCas12a_DetectsTTTV_With23bpTarget` | S6 | ✅ |
| 21 | `FindPamSites_LbCas12a_DetectsTTTV_With24bpTarget` | S7 | ✅ |
| 22 | `FindPamSites_SaCas9_DetectsNNGAAT` | M15 | ✅ |
| 23 | `FindPamSites_SaCas9_DetectsNNGAGT` | M15 | ✅ |
| 24 | `FindPamSites_SaCas9_DetectsNNGGGT` | M15 | ✅ |
| 25 | `FindPamSites_SaCas9_Returns21bp_TargetSequence` | M15 | ✅ |
| 26 | `FindPamSites_SaCas9_RejectsNonR_AtRPosition` | M15 | ✅ |
| 27 | `FindPamSites_SpCas9NAG_DetectsNAG` | M16 | ✅ |
| 28 | `FindPamSites_SpCas9NAG_MatchesAllNAG_Variants` (4 cases) | M16 | ✅ |
| 29 | `FindPamSites_CasX_DetectsTTCN_Variants` (4 cases) | M17 | ✅ |
| 30 | `FindPamSites_OverlappingPamSites_AllReported` | S3 | ✅ |
| 31 | `FindPamSites_StringOverload_WorksIdentically` | M6 | ✅ |

### Classification Summary

- ✅ Covered: 58 tests (31 methods, 26 parametric cases)
- ❌ Missing: 0
- ⚠ Weak: 0
- 🔁 Duplicate: 0

## Open Questions

None — all PAM sequences, guide lengths, and PAM positions are verified against peer-reviewed literature.

## ASSUMPTIONS

None — all behaviors are grounded in external sources.

## Validation Checklist

- [x] All MUST tests have evidence source
- [x] Invariants specified and tested
- [x] Edge cases documented (empty, null, out-of-bounds, PAM-free)
- [x] API consistency verified (string vs DnaSequence overload)
- [x] All 7 CRISPR systems tested for Name, PAM, GuideLength, PamAfterTarget
- [x] Both strands tested with exact coordinates and PAM content
- [x] IUPAC ambiguity codes tested (N, R, V) with positive and negative cases
- [x] No assumptions — all design decisions backed by external sources
- [x] No duplicates — each test serves a distinct purpose
- [x] Coverage classification complete: 0 missing, 0 weak, 0 duplicate
- [x] Tests passing (58/58)
