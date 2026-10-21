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
| CRISPOR `crispor.py` (maximilianh/crisporWebsite, master, fetched 2026-09-28) | Reference implementation | Coordinate/orientation convention (`findAllPams`, `flankSeqIter`), PAM + guide-length table (`setupPamInfo`), IUPAC matcher (`patMatch`); 136-site numerical cross-check — see `docs/Evidence/CRISPR-PAM-001-Evidence.md` |
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
| Kleinstiver et al. (2015), Nature 523:481 (PMC4540238, opened 2026-10-21) | Academic | VQR SpCas9 variant: strongly depleted `NGAN` PAMs; NGA as a non-canonical SpCas9 PAM (B08 round 3 G8) |
| CRISPOR `crispor.py`, `pamDesc` L188 (`('NGA','20bp-NGA - Cas9 S. Pyogenes mutant VQR')`), `offtargetPams` L243-249 (`"NGA" : ["NGG"]`, `"NNGRRT" : ["NNGRRN"]`), `setupPamInfo` L458-520 (guide length 20 for NGA; 21 for NNGRRT) | Reference implementation | `SpCas9_NGA` / `SaCas9_NNGRRN` PAM motifs, guide lengths, PAM-after-target (B08 round 3 G8) |

## Test Categories

### MUST Tests (Required for Completion)

| ID | Test Case | Rationale | Source |
|----|-----------|-----------|--------|
| M1-M3 | GetSystem returns correct Name and PAM for all 9 systems | Canonical PAM sequences | Jinek 2012, Ran 2015, Zetsche 2015, Hsu 2013, Liu 2019, Kleinstiver 2015, CRISPOR `offtargetPams` |
| M4 | GetSystem returns correct guide lengths for all 9 systems | System-specific guide lengths | Jinek 2012, Ran 2015, Zetsche 2015, Liu 2019, CRISPOR `setupPamInfo` |
| M5 | GetSystem indicates PAM position (before/after target) for all 9 systems | Cas9=3' after, Cas12a/CasX=5' before | Wikipedia PAM, Zetsche 2015, Liu 2019 |
| M6 | FindPamSites finds NGG on forward strand with correct position and target start | Core detection with position verification | Wikipedia PAM |
| M7 | FindPamSites matches all NGG variants (AGG, CGG, TGG, GGG) with position | N = any nucleotide per IUPAC; position verified | IUPAC nomenclature |
| M8 | FindPamSites searches reverse strand with exact PAM, position, target start and count | Both strands searched; PAM read on the protospacer strand; forward-strand coordinates | CRISPOR `flankSeqIter`; biology: guide RNA targets both strands |
| M9 | FindPamSites returns no sites on either strand for PAM-free sequence | No false positives on forward or reverse | Core functionality |
| M10 | FindPamSites returns empty for empty sequence | Edge case | Core functionality |
| M11 | FindPamSites handles case-insensitive input | Usability requirement | Implementation spec |
| M12 | FindPamSites returns target with correct length AND content | Target extraction verified by value, not just length | Core functionality |
| M13 | FindPamSites Cas12a detects TTTA, TTTC, TTTG with position, target start, and content | V = A, C, G not T | IUPAC, Zetsche 2015 |
| M14 | FindPamSites does not detect TTTT for Cas12a | T is excluded from IUPAC V | IUPAC nomenclature |
| M15 | FindPamSites SaCas9 detects NNGRRT with R=A, R=G, and mixed R; verified by PAM content | R = A or G; negative test with C at R position | Ran 2015, IUPAC |
| M16 | SpCas9_NAG detects NAG pattern with all N variants with position verification | Secondary SpCas9 PAM | Hsu 2013, Wikipedia Cas9 |
| M17 | CasX detects TTCN with all N variants (TTCA, TTCC, TTCG, TTCT) with position and target | CasX/Cas12e 5' PAM | Liu et al. 2019 |
| M18 | `SpCas9_NGA` (`NGA`, 20 nt, PAM after target) finds every site on both strands: strand, `Position`, `TargetStart`, PAM text, 20-nt target (3 forward + 6 reverse on the constructed sequence); NGG (`TGG`) is not an NGA site | VQR PAM; both-strand convention | Kleinstiver 2015; CRISPOR `findAllPams`/`flankSeqIter`; independent brute-force oracle (round 3 G8) |
| M19 | `SaCas9_NNGRRN` (`NNGRRN`, 21 nt, PAM after target) finds every site on both strands with 6-nt PAM text and 21-nt target (4 forward + 6 reverse); all 16 R/R × last-base PAMs accepted, non-R / non-G rejected; NNGRRN is a superset of NNGRRT (`CCGAAC` is NNGRRN only) | SaCas9 relaxed PAM | Ran 2015; CRISPOR `offtargetPams` `NNGRRT → NNGRRN`; IUPAC; independent oracle (round 3 G8) |

### SHOULD Tests (Important but not blocking)

| ID | Test Case | Rationale | Source |
|----|-----------|-----------|--------|
| S1 | FindPamSites excludes sites where target would be out of bounds (Cas9 and Cas12a) | Boundary safety for both PAM-before and PAM-after systems | Implementation spec |
| S2 | FindPamSites returns correct forward-strand coordinates for reverse strand hits | Coordinate conversion: revcomp position → forward position | Implementation spec |
| S2b | Reverse-strand SpCas9 sites match CRISPOR field-by-field (pamStart, guideStart, pamSeq, guideSeq) | Locks the reference convention: forward-strand coordinates, protospacer-strand sequences | CRISPOR `findAllPams` + `flankSeqIter` |
| S2c | Reverse-strand Cas12a site matches CRISPOR (`guideStart = Position - guideLength`) | PAM-before systems mirror the guide side on the reverse strand | CRISPOR `flankSeqIter`, `pamIsFirst` branch |
| S2d | For every site (either strand, 6 systems: SpCas9, SpCas9_NGA, SaCas9, SaCas9_NNGRRN, Cas12a, CasX) `TargetSequence` = forward window `[TargetStart, +guideLength)` (revcomp on reverse) and `PamSequence` matches the motif | Invariants INV-04/05/06 hold jointly | CRISPOR convention, IUPAC |
| S3 | FindPamSites handles overlapping PAM sites — exact count and positions verified | Each PAM position defines a unique guide RNA | Wikipedia PAM, Hsu 2013 |
| S4 | GetSystem throws ArgumentException for unknown system type | Input validation | .NET conventions |
| S5 | FindPamSites with null DnaSequence throws ArgumentNullException | Input validation | .NET conventions |
| S6 | AsCas12a FindPamSites produces 23bp target end-to-end | Variant-specific guide length | Zetsche 2015 |
| S7 | LbCas12a FindPamSites produces 24bp target end-to-end | Variant-specific guide length | Zetsche 2015 |
| S8 | `SpCas9_NGA` / `SaCas9_NNGRRN`: a PAM whose target would run off the sequence start (forward) or end (reverse) is excluded; the exact-fit case (target flush with the edge) is kept | Boundary rule on both strands for the new systems | CRISPOR `findPams` `minPosPlus` / `maxPosMinus`; independent oracle (round 3 G8) |
| S9 | `SpCas9_NGA` / `SaCas9_NNGRRN` with `N` in the genome: an `N` inside the target is carried through on both strands (CRISPOR `flankSeqIter`, `doFilterNs=False`); an `N` at a concrete PAM position (`G`, `A`, `R`, leading `T` of the reverse reading) is not a PAM | Ambiguity handling; agrees with CRISPOR `patMatch` | CRISPOR `patMatch` L2347-2394; independent oracle (round 3 G8) |

### COULD Tests (Nice to have)

| ID | Test Case | Rationale | Source |
|----|-----------|-----------|--------|
| C1 | FindPamSites performance is linear in sequence length | O(n) complexity | Complexity spec |

## Review-fix traceability (B08 F1, F2)

| Fix | Contract | Tests | Source |
|-----|----------|-------|--------|
| F1 | Reverse-strand `TargetStart` is a forward-strand 0-based coordinate (`Position + pamLen` for PAM-after, `Position − guideLength` for PAM-before) — CRISPOR minus-strand `guideStart` | M8 `FindPamSites_SpCas9_SearchesReverseStrand`, S2b `FindPamSites_SpCas9_ReverseStrand_MatchesCrisporConvention`, S2c `FindPamSites_Cas12a_ReverseStrand_MatchesCrisporConvention`, S2d `FindPamSites_AllSites_AreConsistentInForwardCoordinates` (4 cases); downstream `GuideRnaCandidate.Position` on a reverse guide: CRISPR-GUIDE-001 S-012 | CRISPOR `crispor.py` `findAllPams`, `flankSeqIter`, `findPams` doctests; 136-site cross-check (`docs/Evidence/CRISPR-PAM-001-Evidence.md`) |
| F2 | `PamSequence` is read on the protospacer strand and satisfies the system's PAM motif (`TGG`, `TTTA`, not `CCA`, `TAAA`) | M8, S2b, S2c, S2d (PAM motif asserted per site) | CRISPOR `pamSeq` (same run) |

## Invariants

1. **PAM pattern matching**: A site is returned if and only if the PAM pattern matches using IUPAC rules
2. **Target extraction**: Target sequence length equals the system's guide length; content matches the genomic region
3. **Strand correctness / coordinate convention** (CRISPOR `findAllPams` + `flankSeqIter`): `Position` and `TargetStart` are forward-strand 0-based coordinates on **both** strands; `PamSequence` and `TargetSequence` are read 5'→3' on the protospacer strand, so `PamSequence` always satisfies the system's PAM motif and, for reverse-strand hits, `TargetSequence == revcomp(seq[TargetStart, TargetStart+guideLength))`
4. **No self-targeting**: Empty input returns empty output (no crashes, no false positives)

## Coverage Classification

**Total: 81 tests in canonical file (58 + 6 added 2026-09-28 by the review campaign: S2b, S2c and S2d with 4 parametric cases; +17 on 2026-10-21, B08 round 3 G8: 9 new test methods, 6 new GetSystem parametric cases, 2 new S2d cases)**

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
| M8 | SearchesReverseStrand | `Is.Not.Empty` + `Position >= 0` | `Count == 1` + `Position == 0` + `PamSequence == "TGG"` + `TargetStart == 3` + `TargetSequence` (CRISPOR row; was `"CCA"` before B08 F2) |
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

### Canonical File (`CrisprDesigner_PAM_Tests.cs`) — 43 test methods (81 cases)

| # | Test Method | Spec ID | Status |
|---|-------------|---------|--------|
| 1 | `GetSystem_ReturnsCorrectNameAndPam` (9 cases) | M1-M3 | ✅ |
| 2 | `GetSystem_ReturnsCorrectGuideLength` (9 cases) | M4 | ✅ |
| 3 | `GetSystem_ReturnsCorrectPamPosition` (9 cases) | M5 | ✅ |
| 4 | `GetSystem_UnknownType_ThrowsArgumentException` | S4 | ✅ |
| 5 | `FindPamSites_SpCas9_DetectsNGG_OnForwardStrand` | M6 | ✅ |
| 6 | `FindPamSites_SpCas9_MatchesAllNGG_Variants` (4 cases) | M7 | ✅ |
| 7 | `FindPamSites_SpCas9_SearchesReverseStrand` | M8 | ✅ |
| 8 | `FindPamSites_SpCas9_ReverseStrand_PositionConversion` | S2 | ✅ |
| 8b | `FindPamSites_SpCas9_ReverseStrand_MatchesCrisporConvention` | S2b | ✅ |
| 8c | `FindPamSites_Cas12a_ReverseStrand_MatchesCrisporConvention` | S2c | ✅ |
| 8d | `FindPamSites_AllSites_AreConsistentInForwardCoordinates` (6 cases) | S2d | ✅ |
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
| 32 | `FindPamSites_SpCas9NGA_BothStrands_MatchesIndependentScanAndCrispor` | M18 | ✅ |
| 33 | `FindPamSites_SpCas9NGA_DoesNotMatchNGG` | M18 | ✅ |
| 34 | `FindPamSites_SpCas9NGA_SequenceEdge_ExcludesOverrunAndKeepsExactFit` | S8 | ✅ |
| 35 | `FindPamSites_SpCas9NGA_AmbiguityN_InTargetKept_AtConcretePamPositionRejected` | S9 | ✅ |
| 36 | `FindPamSites_SaCas9NNGRRN_BothStrands_MatchesIndependentScanAndCrispor` | M19 | ✅ |
| 37 | `FindPamSites_SaCas9NNGRRN_AcceptsAllRrCombinationsAndAnyLastBase_RejectsNonR` | M19 | ✅ |
| 38 | `FindPamSites_SaCas9NNGRRN_IsSupersetOfNNGRRT` | M19 | ✅ |
| 39 | `FindPamSites_SaCas9NNGRRN_SequenceEdge_ExcludesOverrunAndKeepsExactFit` | S8 | ✅ |
| 40 | `FindPamSites_SaCas9NNGRRN_AmbiguityN_InTargetKept_AtConcretePamPositionRejected` | S9 | ✅ |

### Classification Summary

- ✅ Covered: 81 tests (43 methods, 45 parametric cases)
- ❌ Missing: 0
- ⚠ Weak: 0
- 🔁 Duplicate: 0

## Open Questions

- `LbCas12a` is configured with a 24-nt guide; CRISPOR uses 23 nt for every Cpf1 PAM and Zetsche et al. (2015) treat
  23 nt as canonical (20-24 functional). Left unchanged in the 2026-09 review because no opened source supports
  either value over the other for LbCas12a specifically — recorded for the ledger.

- **Genome `N` against a PAM `N`.** CRISPOR's `patMatch` lets a pattern `N` match anything, including a genome `N`;
  `IupacHelper.MatchesIupac('N','N')` accepts only `A/C/G/T`, so `FindPamSites` does not report a PAM whose
  "any base" position is itself an `N` in the genome (e.g. genome `NGG` for the NGG motif, `CCGAAN` for NNGRRN). This is
  a recorded divergence that is NOT pinned by a test (the library and the reference disagree; no source ranks one
  over the other); the S9 tests use only the cases where both agree. CRISPOR's optional `doFilterNs` (drop guides
  containing `N`) is likewise not applied: the target keeps its `N`.

## ASSUMPTIONS

None — all behaviors are grounded in external sources.

## Validation Checklist

- [x] All MUST tests have evidence source
- [x] Invariants specified and tested
- [x] Edge cases documented (empty, null, out-of-bounds, PAM-free)
- [x] API consistency verified (string vs DnaSequence overload)
- [x] All 9 CRISPR systems tested for Name, PAM, GuideLength, PamAfterTarget (`SpCas9_NGA` / `SaCas9_NNGRRN` added in B08 round 3 G8)
- [x] Both strands tested with exact coordinates and PAM content
- [x] IUPAC ambiguity codes tested (N, R, V) with positive and negative cases
- [x] No assumptions — all design decisions backed by external sources
- [x] No duplicates — each test serves a distinct purpose
- [x] Coverage classification complete: 0 missing, 0 weak, 0 duplicate
- [x] Tests passing (64/64, 2026-09-28; 81/81, 2026-10-21 after round 3 G8)
- [x] Reverse-strand coordinate/orientation convention cross-checked numerically against CRISPOR (136 sites, 0 mismatches)
- [x] B08 F1/F2 traced to tests (B08 WP6/D16, 2026-10-15; 64/64 passing)
- [x] `SpCas9_NGA` / `SaCas9_NNGRRN`: both strands, target length, PAM text, coordinates, `N` in genome, sequence edge — expected values from an independent Python brute-force scan cross-checked against CRISPOR's own functions (B08 round 3 G8, 2026-10-21; 81/81 passing)
