# Test Specification: PROBE-DESIGN-001

## Test Unit Information

| Field | Value |
|-------|-------|
| **Test Unit ID** | PROBE-DESIGN-001 |
| **Area** | MolTools |
| **Title** | Hybridization Probe Design |
| **Canonical Class** | `ProbeDesigner` |
| **Canonical Methods** | `DesignProbes`, `DesignTilingProbes`, `ScoreProbe` (via EvaluateProbe), `EvaluateTaqManProbe`, `SelectTaqManStrand` |
| **Complexity** | O(n²) |
| **Status** | ☐ Pending re-validation (TaqMan opt-in rules added) |
| **Last Updated** | 2026-06-24 |

---

## Evidence Sources

| Source | Type | Key Information |
|--------|------|-----------------|
| Wikipedia: Nucleic acid thermodynamics | Academic | Tm calculation, nearest-neighbor method, GC content effects |
| Wikipedia: Hybridization probe | Academic | Probe design principles, applications (15-10000 nt) |
| Wikipedia: FISH | Academic | BAC probes ~100 kb; oligo FISH probes 10-25 nt; smFISH 20-50 oligos per target |
| Wikipedia: DNA microarray | Academic | Microarray probe design: Affymetrix 25-mer, Agilent 60-mer |
| Wikipedia: Molecular beacon | Academic | Loop 18-30 bp, stem 5-7 nt each side, typical total 25 nt |
| SantaLucia (1998) | Research | Unified nearest-neighbor thermodynamics for Tm |
| Breslauer et al. (1986) | Research | Predicting DNA duplex stability |
| PREMIER Biosoft "TaqMan probe design tips" | Vendor doc | TaqMan: length 18-22 nt, GC 30-80%, more Cs than Gs and no G at 5' end, no ≥4-G runs, probe Tm 10 °C above primer Tm |
| Applied Biosystems / Thermo Fisher "Designing a TaqMan Gene Expression Assay" | Manufacturer | No 5' G (interferes with reporter fluorescence); probe Tm ~10 °C above primer; antisense fallback when 5' G unavoidable |
| ScienceDirect "TaqMan — an overview" | Reference work | 5' G adjacent to reporter quenches fluorescence even after cleavage (hard rule) |

---

## Invariants

1. **Score Range**: 0.0 ≤ score ≤ 1.0 (Source: Implementation)
2. **GC Range**: 0.0 ≤ GC content ≤ 1.0 (Source: Mathematical definition)
3. **Tm Positivity**: Tm > 0 for valid probes (Source: Physical law)
4. **Coordinate Validity**: 0 ≤ Start < End < sequence.Length (Source: Implementation)
5. **Probe Substring**: probe.Sequence == input.Substring(probe.Start, probe.End - probe.Start + 1) (Source: Implementation)

---

## Test Cases

### Must (Required - Evidence-Based)

| ID | Test Case | Rationale | Source |
|----|-----------|-----------|--------|
| M1 | Empty sequence returns empty result | Boundary condition | Implementation spec |
| M2 | Null sequence returns empty result | Boundary condition | Implementation spec |
| M3 | Sequence shorter than MinLength returns empty | Length constraint | Implementation spec |
| M4 | Valid sequence produces probes with score in [0,1] | Invariant #1 | Implementation |
| M5 | All probes have GC content in [0,1] | Invariant #2 | Mathematical |
| M6 | All probes have Tm > 0 | Invariant #3 | Physical law |
| M7 | Probe coordinates are valid (Start ≥ 0, End < seq.Length) | Invariant #4 | Implementation |
| M8 | Probe sequence matches substring at coordinates | Invariant #5 | Implementation |
| M9 | Tiling probes cover expected positions | Coverage guarantee | Algorithm spec |
| M10 | Tiling probes all have Type = Tiling | Type consistency | Implementation |
| M11 | Microarray defaults: length 50-60 bp | Application param | Wikipedia (DNA microarray) |
| M12 | FISH defaults: length 200-500 bp | Application param | Standard molecular biology practice |
| M13 | High GC content (100%) results in GcContent ≈ 1.0 | Edge case | Mathematical |
| M14 | Low GC content (all A/T) results in low GcContent | Edge case | Mathematical |
| M15 | maxProbes parameter limits returned count | API contract | Implementation |
| TM1 | TaqMan probe with 5' G is flagged (`NoGuanineAt5Prime == false`) and `PassesAll == false` | 5' G quenches reporter even after cleavage | ABI / ScienceDirect |
| TM2 | Run of ≥4 consecutive Gs flagged (`NoRunOfFourOrMoreG == false`) | No ≥4-G runs | PREMIER Biosoft |
| TM3 | More G than C flagged (`MoreCytosineThanGuanine == false`; C=1, G=9) | More Cs than Gs | PREMIER Biosoft |
| TM4 | Length outside 18-22 flagged (`LengthInRange == false`; 16 nt) | Length 18-22 nt | PREMIER Biosoft |
| TM5 | GC outside 30-80% flagged (`GcContentInRange == false`; GC = 1.0) | GC 30-80% | PREMIER Biosoft |
| TM6 | Probe Tm < primerTm + 10 flagged (Tm 49.35, primer 45) | Probe Tm ≥ primer + 10 °C | PREMIER Biosoft / ABI |
| TM7 | Fully compliant probe accepted (`PassesAll == true`, no violations) | All rules satisfied | PREMIER Biosoft |
| TM8 | Null primerTm skips the Tm gate (reported satisfied) | Optional primer Tm | Implementation contract |
| TM9 | `SelectTaqManStrand` picks antisense when sense has 5' G / more G | Antisense fallback | ABI |
| TM10 | `SelectTaqManStrand` keeps the sense strand when already compliant | No needless RC | ABI |

### Should (Important)

| ID | Test Case | Rationale | Source |
|----|-----------|-----------|--------|
| S1 | Homopolymer runs generate warnings | Quality check | General practice |
| S2 | Case-insensitive input handling | Usability | Implementation |
| S3 | DesignAntisenseProbes returns Antisense type | Type correctness | Implementation |
| S4 | MolecularBeacon has stem sequences | Structure check | Implementation |
| S5 | Tiling probes calculate mean Tm correctly | Statistics | Implementation |
| S6 | Probes are sorted by score descending | Ranking | Implementation |

### Could (Optional)

| ID | Test Case | Rationale | Source |
|----|-----------|-----------|--------|
| C1 | qPCR defaults produce 20-30 bp probes | Application param | Standard practice |
| C2 | Self-complementarity detection works correctly | Quality metric | Implementation |
| C3 | Secondary structure detection identifies hairpins | Quality metric | Implementation |

---

## Coverage Classification

Canonical file (generic designer): `ProbeDesigner_ProbeDesign_Tests.cs` (29 tests).
Canonical file (TaqMan opt-in rules): `ProbeDesigner_TaqMan_Tests.cs` (12 tests).
Supplementary file: `ProbeDesignerTests.cs` (6 tests — smoke/utility, no PROBE-DESIGN-001 scope).

| ID | Test | Status |
|----|------|--------|
| M1 | DesignProbes_EmptySequence_ReturnsEmpty | ✅ Covered |
| M2 | DesignProbes_NullSequence_ReturnsEmpty | ✅ Covered |
| M3 | DesignProbes_ShortSequence_ReturnsEmpty | ✅ Covered |
| M4 | DesignProbes_ValidSequence_ProbesHaveScoreInValidRange | ✅ Covered |
| M5 | DesignProbes_ValidSequence_ProbesHaveGcContentInValidRange | ✅ Covered |
| M6 | DesignProbes_ValidSequence_ProbesHavePositiveTm | ✅ Covered |
| M7 | DesignProbes_ValidSequence_ProbesHaveValidCoordinates | ✅ Covered |
| M8 | DesignProbes_ValidSequence_ProbeSequenceMatchesSubstring | ✅ Covered |
| M9 | DesignTilingProbes_CoversExpectedPositions | ✅ Covered |
| M10 | DesignTilingProbes_AllProbesHaveTilingType | ✅ Covered |
| M11 | DesignProbes_MicroarrayDefaults_ProducesCorrectLengthProbes | ✅ Covered |
| M12 | DesignProbes_FISHDefaults_ProducesCorrectLengthProbes | ✅ Covered |
| M13 | DesignProbes_AllGC_ReturnsProbesWithHighGcContent | ✅ Covered |
| M14 | DesignProbes_AllAT_ReturnsProbesWithLowGcContent | ✅ Covered |
| M15 | DesignProbes_MaxProbesParameter_LimitsResultCount | ✅ Covered |
| S1 | DesignProbes_HomopolymerSequence_GeneratesWarnings | ✅ Covered |
| S2 | DesignProbes_CaseInsensitiveInput_ProducesConsistentResults | ✅ Covered |
| S3 | DesignAntisenseProbes_ReturnsAntisenseType | ✅ Covered |
| S4 | DesignMolecularBeacon_CreatesBeaconWithStem | ✅ Covered |
| S5 | DesignTilingProbes_CalculatesTmStatisticsCorrectly | ✅ Covered |
| S6 | DesignProbes_ProbesAreSortedByScoreDescending | ✅ Covered |
| C1 | DesignProbes_qPCRDefaults_ProducesCorrectLengthProbes | ✅ Covered |
| C2 | ValidateProbe_SelfComplementarity_DetectsCorrectly | ✅ Covered |
| C3 | ValidateProbe_SecondaryStructure_IdentifiesHairpins | ✅ Covered |
| — | DesignMolecularBeacon_ShortSequence_ReturnsNull | ✅ Boundary |
| — | DesignMolecularBeacon_AtRichTarget_ScorePenalizedForGcAndTm | ✅ Mutation |
| — | DesignMolecularBeacon_GcRichTarget_ScorePenalizedForGcAndTm | ✅ Mutation |
| — | DesignProbes_WithSuffixTree_FiltersNonUniqueProbes | ✅ Specificity |
| — | DesignProbes_WithSuffixTree_PerformanceImprovement | ✅ Integration |
| TM1 | EvaluateTaqManProbe_FivePrimeGuanine_FlaggedAndRejected | ✅ Covered |
| TM2 | EvaluateTaqManProbe_RunOfFourGuanines_Rejected | ✅ Covered |
| TM3 | EvaluateTaqManProbe_MoreGuanineThanCytosine_FlagsCgRule | ✅ Covered |
| TM4 | EvaluateTaqManProbe_LengthOutsideRange_FlagsLengthRule | ✅ Covered |
| TM5 | EvaluateTaqManProbe_GcContentOutsideRange_FlagsGcRule | ✅ Covered |
| TM6 | EvaluateTaqManProbe_ProbeTmNotTenAbovePrimer_FlagsTmGate | ✅ Covered |
| TM7 | EvaluateTaqManProbe_AllRulesSatisfied_Accepted | ✅ Covered |
| TM8 | EvaluateTaqManProbe_NoPrimerTm_TmGateReportedSatisfied | ✅ Covered |
| TM9 | SelectTaqManStrand_SenseHas5PrimeGAndMoreG_PicksAntisense | ✅ Covered |
| TM10 | SelectTaqManStrand_SenseAlreadyCompliant_KeepsSense | ✅ Covered |
| — | EvaluateTaqManProbe_NullSequence_Throws | ✅ Edge |
| — | SelectTaqManStrand_NullSequence_Throws | ✅ Edge |

---

## Open Questions

None - behavior is well-documented in implementation and sources.

