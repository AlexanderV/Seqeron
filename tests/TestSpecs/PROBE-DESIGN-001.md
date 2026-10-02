# Test Specification: PROBE-DESIGN-001

## Test Unit Information

| Field | Value |
|-------|-------|
| **Test Unit ID** | PROBE-DESIGN-001 |
| **Area** | MolTools |
| **Title** | Hybridization Probe Design |
| **Canonical Class** | `ProbeDesigner` |
| **Canonical Methods** | `DesignProbes`, `DesignProbesPrimer3`, `DesignTilingProbes`, `DesignMolecularBeacon`, `EvaluateTaqManProbe`, `SelectTaqManStrand`, `AnalyzeOligo`, `CalculateMolecularWeight`, `CalculateExtinctionCoefficient(NearestNeighbor)` |
| **Complexity** | O(n²) |
| **Status** | Reviewed 2026-10 (B07 F17–F21) |
| **Last Updated** | 2026-10-01 |

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
| primer3 `libprimer3.cc` + primer3-py 2.3.1 | Reference implementation | Probe Tm (`seqtm`), internal-oligo defaults (50 nM, 50 mM, 0 Mg, 0 dNTP; size 18/20/27; Tm 57/60/63; GC 20–80; poly-X 5; ntthal limits 47 °C), `pick_hyb_probe_only` selection and ordering |
| Biopython 1.88 `molecular_weight` | Reference implementation | Single-stranded DNA/RNA molecular weight |
| Cantor, Warshaw & Shapiro 1970; Warshaw & Tinoco 1966 | Research | Nearest-neighbour ε260 tables |
| Tyagi & Kramer 1996 / Marras et al. | Research | Molecular beacon: stem 5–7 bp; probe and stem Tm 7–10 °C above the detection temperature |
| Applied Biosystems Primer Express guidelines | Manufacturer | TaqMan probe Tm 68–70 °C, G+C 30–80 % (qPCR preset) |

---

## Invariants

1. **Score Range**: 0.0 ≤ score ≤ 1.0 (Source: Implementation)
2. **GC Range**: 0.0 ≤ GC content ≤ 1.0 (Source: Mathematical definition)
3. **Tm**: Tm = Primer3 seqtm at the ProbeParameters conditions (primer3-py calc_tm); > 0 for the tested ≥ 20-nt probes
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

| P1 | `DesignProbesPrimer3` defaults = primer3-py `pick_hyb_probe_only` (positions, Tm, penalty, SELF_ANY/END/HAIRPIN_TH, order) | Reference parity | primer3 |
| P2 | `DesignProbesPrimer3` non-default settings + PCR buffer = primer3-py | Reference parity | primer3 |
| P3 | ntthal limits reject 41 of 99 windows of a self-complementary template (primer3 explain) | Thermodynamic screen | primer3 |
| P4 | Probe Tm in `DesignProbes` = Primer3 seqtm at stated conditions | Probe Tm | primer3 |
| P5 | Thermodynamic self-dimer screen flags a GC palindrome | Self-structure | primer3 ntthal |
| P6 | Lazy (branch-and-bound) ranking = exhaustive ranking | Algorithm contract | Implementation |
| P7 | Beacon: loop Tm, ntthal stem-loop Tm, 7 °C rules for a detection temperature | Beacon design | Tyagi & Kramer |
| P8 | NN ε260 = Cantor/Warshaw tables (ACGT 40300, ACGU RNA 41300) | Oligo property | Cantor 1970 |
| P9 | MW = Biopython molecular_weight (DNA/RNA, U = UMP) | Oligo property | Biopython |
| TM11 | TaqMan Tm at stated conditions = primer3-py calc_tm; non-ACGT → gate fails | Probe Tm | primer3 |
| P10 | `DesignProbesPrimer3` with PRIMER_INTERNAL_MISHYB_LIBRARY: probe positions, penalties (incl. PRIMER_INTERNAL_WT_LIBRARY_MISHYB term) and PRIMER_INTERNAL_n_LIBRARY_MISHYB score + entry (dpal unanchored LOCAL, IUPAC per PRIMER_LIB_AMBIGUITY_CODES_CONSENSUS, short entries score their length) = primer3-py `mishyb_lib`, both alignment modes; `CalculateLibraryMishyb` = check_primers; weight without library / limit > 32767 in alignment mode throw (`ProbeDesigner_MishybLibrary_Tests`) | Mishyb library | primer3 libprimer3.c |
| P11 | `DesignProbesPrimer3` / MCP `design_probes_primer3` reject numReturn < 1 (Primer3 "PRIMER_NUM_RETURN < 1") | Data control | primer3 _pr_data_control |
| P12 | `DesignProbesPrimer3` with PRIMER_ANNEALING_TEMP: PRIMER_INTERNAL_n_BOUND, rejection outside PRIMER_INTERNAL_MIN/MAX_BOUND, PRIMER_INTERNAL_WT_BOUND_GT/LT terms (ungated: bound −999999.9999 without an annealing temperature) and the opt-bound / annealing-temperature data control = primer3-py pick_hyb_probe_only (`PrimerDesigner_BoundAndPosition_Tests.DesignProbesPrimer3_*`) | Fraction bound | primer3 libprimer3.c, oligotm.c |
| P13 | `DesignProbesPrimer3` with SEQUENCE_QUALITY: PRIMER_INTERNAL_n_MIN_SEQ_QUALITY, rejection below PRIMER_INTERNAL_MIN_QUALITY, PRIMER_INTERNAL_WT_SEQ_QUAL term (WT_END_QUAL inert) and the quality data control = primer3-py pick_hyb_probe_only (`PrimerDesigner_SequenceQuality_Tests.DesignProbesPrimer3_Quality_MatchesPrimer3`, `DataControl_MatchesPrimer3Messages`; MCP `DesignProbesPrimer3_SequenceQuality_MatchesPrimer3`) | Sequence quality | primer3 libprimer3.c |

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
Canonical file (TaqMan opt-in rules): `ProbeDesigner_TaqMan_Tests.cs`.
Canonical file (Primer3 picker, probe Tm/structure, beacon, oligo properties): `ProbeDesigner_Primer3Probe_Tests.cs` (P1–P8); MW/AnalyzeOligo Tm reference values: `Mutation/ProbeDesignerMutationTests.cs` (P9).
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
| P1 | DesignProbesPrimer3_Defaults_MatchPrimer3PickHybProbeOnly | ✅ Covered |
| P2 | DesignProbesPrimer3_TaqManLikeSettingsAndPcrBuffer_MatchPrimer3 | ✅ Covered |
| P3 | DesignProbesPrimer3_ThermodynamicLimits_RejectSelfComplementaryWindows | ✅ Covered |
| — | DesignProbesPrimer3_InvalidArguments_Throw | ✅ Edge |
| P4 | DesignProbes_ProbeTm_IsPrimer3SeqtmAtParameterConditions | ✅ Covered |
| P5 | DesignProbes_ThermodynamicScreen_FlagsSelfDimerByNtthalTm | ✅ Covered |
| P6 | DesignProbes_LazyStructureScreen_EqualsExhaustiveRanking | ✅ Covered |
| P7 | DesignMolecularBeacon_DetectionTemperature_AppliesTyagiKramerRules | ✅ Covered |
| P8 | ExtinctionCoefficientNearestNeighbor_MatchesCantorTable | ✅ Covered |
| P9 | CalculateMolecularWeight_MatchesBiopython, AnalyzeOligo_Tm_MatchesPrimer3CalcTmAtProbeDefaults | ✅ Covered |
| TM11 | EvaluateTaqManProbe_TmAtStatedConditions_MatchesPrimer3CalcTm | ✅ Covered |

---

## Open Questions

None - behavior is well-documented in implementation and sources.

