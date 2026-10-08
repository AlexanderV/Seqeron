# Test Specification: PRIMER-STRUCT-001

## Test Unit

| Field | Value |
|-------|-------|
| **ID** | PRIMER-STRUCT-001 |
| **Title** | Primer Structure Analysis |
| **Area** | Molecular Tools |
| **Status** | ☑ Complete |
| **Last Updated** | 2026-10-01 |
| **Owner** | GitHub Copilot |

## Methods Under Test

| Method | Class | Type | Complexity |
|--------|-------|------|------------|
| `HasHairpinPotential(seq, minStemLength, minLoopLength)` | PrimerDesigner | Canonical | O(n²) <100bp, O(n) ≥100bp* |
| `HasPrimerDimer(primer1, primer2, minComp)` | PrimerDesigner | Canonical | O(n·m) |
| `CalculatePrimerDimerEndComplementarity(p1, p2)` / `CalculatePrimerSelfEndComplementarity(p)` | PrimerDesigner | Canonical (Primer3 alignment compl_end / self_end) | O(n·m) |
| `CalculatePrimer3OligoStructure(p)` / `CalculatePrimer3PairComplementarity(l, r)` | PrimerDesigner | Canonical (Primer3 *_TH, ntthal) | O(n·m·L²) |
| `CalculateDimerThermodynamicsNtthal(a, b, mode, mv, dv, dntp, c)` / `CalculateHairpinThermodynamicsNtthal(s, mv, dv, dntp)` | PrimerDesigner | ntthal END1/END2 + divalent salt | O(n·m·L²) |
| `Calculate3PrimeStability(seq)` | PrimerDesigner | Canonical | O(1) |
| `FindLongestHomopolymer(seq)` | PrimerDesigner | Canonical | O(n) |
| `FindLongestDinucleotideRepeat(seq)` | PrimerDesigner | Canonical | O(n) |

*Uses suffix tree optimization for long sequences (≥100bp)

## Evidence Sources

| Source | Type | URL/Reference |
|--------|------|---------------|
| Wikipedia - Primer (molecular biology) | Encyclopedia | https://en.wikipedia.org/wiki/Primer_(molecular_biology) |
| Wikipedia - Primer dimer | Encyclopedia | https://en.wikipedia.org/wiki/Primer_dimer |
| Wikipedia - Stem-loop | Encyclopedia | https://en.wikipedia.org/wiki/Stem-loop |
| Wikipedia - Nucleic acid thermodynamics | Encyclopedia | https://en.wikipedia.org/wiki/Nucleic_acid_thermodynamics |
| Primer3 Manual | Tool Documentation | https://primer3.org/manual.html |
| SantaLucia (1998) | Primary Literature | PNAS 95:1460-65 |
| Primer3 source (libprimer3.cc, oligotm.c, dpal.c, thal.c) | Reference implementation | https://github.com/primer3-org/primer3 (oligotm.c and dpal.c compiled locally as oracles) |
| primer3-py 2.3.1 | Reference implementation | `design_primers` (check_primers), `calc_hairpin`, `calc_homodimer`, `calc_heterodimer`, `calc_end_stability` |

## Invariants

1. **Homopolymer invariant:** Result ≥ 1 for non-empty sequences, 0 for empty
2. **Dinucleotide invariant:** Result ≥ 0; 0 for sequences < 4 bp (unsourced library screen, no Primer3 counterpart — audit round 3, A3-8)
3. **Hairpin invariant:** Requires minimum length (2×stem + loop) to return true (`HasHairpinPotential` is a sequence-only library screen: default loop 3 = Primer3 `thal.c` `min_hrpn_loop`, default stem 4 unsourced; it can disagree with Primer3's ntthal hairpin, see M10b)
4. **Stability invariant:** GC-rich 3' ends have more negative (stable) ΔG
5. **Primer-dimer invariant:** Returns false for empty primers; the compl_end score is symmetric and ≥ 0
6. **Primer3 structure invariant:** every *_TH value is ≥ 0 (negative Tm / no structure → 0)

## Test Cases

### MUST Tests (Required)

| ID | Category | Test Case | Expected | Source |
|----|----------|-----------|----------|--------|
| M1 | Homopolymer | Null/empty sequence | 0 | Primer3 default behavior |
| M2 | Homopolymer | No run (ACGT) | 1 | Primer3 PRIMER_MAX_POLY_X |
| M3 | Homopolymer | All same (AAAAAA) | 6 | Primer3 PRIMER_MAX_POLY_X |
| M4 | Homopolymer | Mixed case (AaAaAa) | 6 (case insensitive) | Universal DNA convention |
| M5 | Dinucleotide | Null/empty/short (<4 bp) | 0 | Implementation bounds |
| M6 | Dinucleotide | No repeat (ACGT) | 1 | Library screen contract (unsourced; Primer3 has no dinucleotide limit — A3-8) |
| M7 | Dinucleotide | ACACACAC | 4 | Library screen contract (unsourced; Primer3 has no dinucleotide limit — A3-8) |
| M8 | Hairpin | Null/empty/too-short | false | Stem-loop theory (min 2×stem + loop) |
| M9 | Hairpin | Non-self-complementary | false | Wikipedia Stem-loop |
| M10 | Hairpin | Self-complementary | true | Wikipedia Stem-loop |
| M10b | Hairpin | Library screen ≠ Primer3: `CAGTAAAACCCTTTTGCAGC` flagged; `AAAACCCTTTT` flagged | true / true, while Primer3 ntthal hairpin Tm 37.65 °C (< 47) / no structure | primer3-py 2.3.1 `calc_hairpin` (50 mM / 1.5 mM / 0.6 mM / 50 nM) |
| M11 | Primer-dimer | Null/empty primer (either side) | false | Null guard |
| M12 | Primer-dimer | Non-complementary 3' ends (A₈ vs AAAACCCC) | false, compl_end 0 | primer3-py check_primers (alignment mode) |
| M13 | Primer-dimer | Complementary 3' ends (A₈ vs T₈); identical poly-A (A₈ vs A₈) is NOT a dimer; GGCC/GGCC 3' ends ARE | true / false / true (8 / 0 / 4) | Primer3 dpal.c (compiled), primer3-py check_primers |
| M13b | Primer-dimer | compl_end and self_end match Primer3 | PRIMER_PAIR_0_COMPL_END 1/8/4/0; SELF_END 6/4/0/0 | primer3-py check_primers (alignment mode) |
| M14 | 3' Stability | Null/empty → 0; < 5 bp scored whole (ACGT −2.56, GC +0.15, A +2.06, AT +1.61); N (ACGTN −3.62, NNNNN −0.36, AAANA −1.40); invalid char in window → NaN | as listed | Primer3 end_oligodg (oligotm.c compiled) |
| M18 | Thermodynamic screen | SELF_ANY_TH / SELF_END_TH / HAIRPIN_TH / COMPL_ANY_TH / COMPL_END_TH for 4 primer pairs | primer3-py values to 1e-9 | primer3-py design_primers check_primers |
| M19 | ntthal modes | END1 / END2 / ANY at dv 1.5, dntp 0.6; hairpin at dv 1.5 and 0 | calc_end_stability / calc_heterodimer / calc_hairpin to 1e-9 | primer3-py |
| M20 | Homopolymer N | NNG 3, ANA 3, CNN 3, TNNG 3, GNGNG 5, ANGNG 4 + 5 primers | Primer3 _pr_violates_poly_x | libprimer3.cc comment + primer3-py check_primers |
| M21 | EvaluatePrimer | Thermodynamic screen populates SelfAnyTh/SelfEndTh/HairpinTh, rejects > 47 °C; Heuristic screen → null; MaxStructureTm 100 → no structure issue | ACGTACGTACGTACGTACGT 60.49/60.49/71.76 | primer3-py check_primers |
| M15 | 3' Stability | GC-rich vs AT-rich (exact values) | GCGCG = -6.86, TATAT = -0.86 | SantaLucia (1998) + Primer3 Manual |
| M16 | 3' Stability | GCGCG (most stable 5mer) | -6.86 kcal/mol | Primer3 Manual + SantaLucia (1998) |
| M17 | 3' Stability | TATAT (least stable 5mer) | -0.86 kcal/mol | Primer3 Manual + SantaLucia (1998) |

### SHOULD Tests (Recommended)

| ID | Category | Test Case | Expected | Source |
|----|----------|-----------|----------|--------|
| S1 | Homopolymer | Internal run (ACAAAAGT) | 4 | Common pattern |
| S2 | Dinucleotide | ATATATAT pattern | 4 | Common microsatellite |
| S3 | Hairpin | Custom minStemLength | Respects parameter | API contract |
| S4 | Hairpin | Custom minLoopLength | Respects parameter | API contract / Wikipedia Stem-loop |
| S5 | Primer-dimer | Custom minComplementarity (ACGTACGT self, score 8: min 8 true, min 9 false) | Respects parameter | dpal.c (compiled) |
| S6 | 3' Stability | Exact 5-base input (TACGT) | -3.57 kcal/mol | SantaLucia (1998) |
| S7 | 3' Stability | Case insensitive (GCGCG vs gcgcg) | Both = -6.86 | Universal DNA convention |

### COULD Tests (Optional)

| ID | Category | Test Case | Expected | Source |
|----|----------|-----------|----------|--------|
| C1 | Homopolymer | Run at end / multiple runs | Detected correctly | Edge case |
| C2 | Dinucleotide | Multiple repeat types | Returns longest | Logic verification |
| C3 | Hairpin | Long sequence (>100bp) suffix tree path | Correct detection | Performance optimization |
| C6 | DesignPrimers | Default thermodynamic screen = primer3-py design_primers (random template: right AGGAACGGATCGAGGACTGC, pair penalty 2.425289002682007); hairpin 48.215 °C right primer → no valid primers; Heuristic / MaxStructureTm 100 → pair 1.869300477208128 | Primer3 values | primer3-py design_primers |
| C7 | Alignment mode (A1, audit round 2) | dpal LOCAL `CalculatePrimerSelfAnyComplementarity` / `CalculatePrimerDimerAnyComplementarity` (+ self_end) = compiled dpal.c + `align()` (len < 3 rule, N −0.25, 136-nt input); `EvaluatePrimer` with `Primer3Alignment` limits 8/3 | dpal.c values (`PrimerDesigner_AlignmentMode_Tests`) | dpal.c, libprimer3.cc |
| C4 | Integration | Well-designed primer exact metrics | Homopolymer=1, dinuc=1, ΔG=-3.57 | Combined verification |
| C5 | Integration | Problematic primer exact metrics | Homopolymer=20, ΔG=-5.40 | Primer3 failure modes |

## Coverage Classification

**Total: 31 test runs (was 35 before classification — 4 duplicates merged, 5 weak strengthened, 4 missing added)**

| Test Method | Runs | Classification | Action |
|-------------|------|---------------|--------|
| `FindLongestHomopolymer_EmptySequence_ReturnsZero` | 1 | ✅ Covered | M1 |
| `FindLongestHomopolymer_NullSequence_ReturnsZero` | 1 | ✅ Covered | M1 |
| `FindLongestHomopolymer_NoRun_ReturnsOne` | 1 | ✅ Covered | M2 |
| `FindLongestHomopolymer_InternalRun_ReturnsRunLength` | 1 | ✅ Covered | S1 |
| `FindLongestHomopolymer_AllSame_ReturnsFullLength` | 1 | ✅ Covered | M3 |
| `FindLongestHomopolymer_MixedCase_IsCaseInsensitive` | 1 | ✅ Covered | M4 |
| `FindLongestHomopolymer_RunAtEnd_ReturnsRunLength` | 1 | ✅ Covered | C1 |
| `FindLongestHomopolymer_MultipleRuns_ReturnsLongest` | 1 | ✅ Covered | C1 |
| `FindLongestHomopolymer_NIsWorstCaseWildcard_MatchesPrimer3PolyX` | 11 | ✅ Covered | M20 |
| `FindLongestDinucleotideRepeat_InvalidInput_ReturnsZero` | 3 | ✅ Covered | M5 (null+empty+short merged) |
| `FindLongestDinucleotideRepeat_NoRepeat_ReturnsOne` | 1 | ✅ Covered | M6 (was ⚠ Weak: `≤1` → exact `1`) |
| `FindLongestDinucleotideRepeat_AcRepeat_ReturnsCount` | 1 | ✅ Covered | M7 |
| `FindLongestDinucleotideRepeat_AtRepeat_ReturnsCount` | 1 | ✅ Covered | S2 |
| `FindLongestDinucleotideRepeat_MultipleRepeats_ReturnsLongest` | 1 | ✅ Covered | C2 |
| `HasHairpinPotential_InvalidOrTooShort_ReturnsFalse` | 4 | ✅ Covered | M8 (null+empty+short+borderline merged) |
| `HasHairpinPotential_NonSelfComplementary_ReturnsFalse` | 1 | ✅ Covered | M9 |
| `HasHairpinPotential_SelfComplementary_ReturnsTrue` | 1 | ✅ Covered | M10 |
| `HasHairpinPotential_CustomMinStem_RespectsParameter` | 1 | ✅ Covered | S3 |
| `HasHairpinPotential_CustomMinLoopLength_RespectsParameter` | 1 | ✅ Covered | S4 (was ❌ Missing) |
| `HasHairpinPotential_LibraryScreen_DiffersFromPrimer3NtthalHairpin` | 2 | ✅ Covered | M10b |
| `HasHairpinPotential_LongSequence_UsesSuffixTreeOptimization` | 1 | ✅ Covered | C3 |
| `HasHairpinPotential_LongSequenceNoHairpin_ReturnsFalse` | 1 | ✅ Covered | C3 |
| `HasHairpinPotential_NonAcgtBlocks_NeverPair_OnBothPaths` | 2 | ✅ Covered | C3: only A·T/G·C pair on both paths (97 nt O(n²), 113 nt suffix tree; the suffix path used the IUPAC N→N complement before audit round 5, A5-2) |
| `HasHairpinPotential_UracilAndLowercase_SameOnBothPaths` | 5 | ✅ Covered | C3: U never pairs with A and lowercase pairs like uppercase on both paths (14 nt O(n²), 134 nt suffix tree; the suffix path mapped U→A before audit round 6, A6-1) |
| `HasHairpinPotential_ShortScanAndSuffixTree_AgreeOnRandomSequences` | 1 | ✅ Covered | C3: 2000 seeded random < 100-nt sequences (ACGTUNS + lowercase) give the same result as their 100-N-padded suffix-tree extension (A6-1) |
| `HasPrimerDimer_NullOrEmptyPrimer_ReturnsFalse` | 4 | ✅ Covered | M11 (4 cases: null/empty × both sides) |
| `HasPrimerDimer_NonComplementary3Ends_ReturnsFalse` | 1 | ✅ Covered | M12 (fixture replaced 2026-10-01: the old AAAACCCCCCCC/GGGGGGGGTTTT pair is fully reverse-complementary, compl_end 12) |
| `HasPrimerDimer_Complementary3Ends_ReturnsTrue` | 1 | ✅ Covered | M13 (A₈/T₈) |
| `HasPrimerDimer_IdenticalPolyA_IsNotADimer` | 1 | ✅ Covered | M13 regression (A₈/A₈ was flagged) |
| `HasPrimerDimer_SelfComplementaryGgccEnds_Detected` | 1 | ✅ Covered | M13 regression (GGCC ends) |
| `CalculatePrimerDimerEndComplementarity_MatchesPrimer3ComplEnd` | 4 | ✅ Covered | M13b |
| `CalculatePrimerSelfEndComplementarity_MatchesPrimer3SelfEnd` | 4 | ✅ Covered | M13b |
| `HasPrimerDimer_CustomMinComplementarity_RespectsParameter` | 1 | ✅ Covered | S5 (old assertion min 8 → false was wrong: ACGTACGT self compl_end = 8) |
| `Calculate3PrimeStability_InvalidInput_ReturnsZero` | 2 | ✅ Covered | M14 (null+empty) |
| `Calculate3PrimeStability_MatchesPrimer3EndOligoDg` | 8 | ✅ Covered | M14 (short, N) |
| `Calculate3PrimeStability_InvalidCharacterInWindow_ReturnsNaN` | 1 | ✅ Covered | M14 |
| `Calculate3PrimeStability_Exact5Bases_ProducesCorrectDeltaG` | 1 | ✅ Covered | S6 (was ❌ Missing) |
| `Calculate3PrimeStability_GcRich_MoreNegativeThanAtRich` | 1 | ✅ Covered | M15 (was ⚠ Weak: now exact -6.86/-0.86) |
| `Calculate3PrimeStability_MixedCase_ReturnsSameExactValue` | 1 | ✅ Covered | S7 (was ⚠ Weak: now checks -6.86) |
| `Calculate3PrimeStability_MostStable5mer_MatchesPrimer3` | 1 | ✅ Covered | M16 |
| `Calculate3PrimeStability_LeastStable5mer_MatchesPrimer3` | 1 | ✅ Covered | M17 |
| `PrimerStructureAnalysis_WellDesignedPrimer_ExactMetrics` | 1 | ✅ Covered | C4 (was ⚠ Weak: `True.Or.False` → exact values) |
| `PrimerStructureAnalysis_ProblematicPrimer_ExactMetrics` | 1 | ✅ Covered | C5 (was ⚠ Weak: `<-5.0` → exact -5.40) |
| `Primer3ThermodynamicStructure_MatchesPrimer3CheckPrimers` | 4 | ✅ Covered | M18 |
| `CalculateDimerThermodynamicsNtthal_AlignmentModes_MatchPrimer3Py` | 4 | ✅ Covered | M19 |
| `CalculateHairpinThermodynamicsNtthal_Divalent_MatchesPrimer3Py` | 1 | ✅ Covered | M19 |
| `Primer3ThermodynamicStructure_NoStructureAndInvalidInput` | 1 | ✅ Covered | invariant 6 |
| `EvaluatePrimer_ThermodynamicScreen_ReportsPrimer3Values` | 1 | ✅ Covered | M21 |
| `DesignPrimers_RandomTemplate_MatchesPrimer3DesignPrimers`, `DesignPrimers_HeuristicScreen_KeepsSequenceOnlyChecks`, `DesignPrimers_RightPrimerHairpinAbovePrimer3Limit_NoValidPrimers` (PrimerDesigner_PrimerDesign_Tests) | 3 | ✅ Covered | C6 |

### Classification Summary

| Status | Count | Details |
|--------|-------|---------|
| ❌ Missing → Implemented | 4 | null hairpin, null dinuc, custom minLoopLength, exact 5-base stability |
| ⚠ Weak → Strengthened | 5 | dinuc `≤1`→`=1`, GC/AT ordering→exact values, mixed case→exact value, integration→exact metrics, problematic→exact ΔG |
| 🔁 Duplicate → Merged | 4 | dinuc empty+short→TestCase, hairpin empty+short+borderline→TestCase, dimer empty+null→TestCase, dimer `ComplementaryEnds_FormsDimer` removed (identical to `Complementary3Ends_ReturnsTrue`) |
| ✅ Covered | 22 | Already had correct exact assertions |

**Result: 0 missing, 0 weak, 0 duplicate**

## External Source Verification

All implementation details have been verified against external sources. No residual assumptions remain.

| Item | Verification | Source |
|------|-------------|--------|
| NN ΔG°37 values (16 dinucleotides) | Exact match with SantaLucia (1998) Table 1 unified parameters | SantaLucia (1998) PNAS 95:1460-65, Table 1 |
| Initiation parameters (+0.98 G·C, +1.03 A·T) | Included in Calculate3PrimeStability | SantaLucia (1998) Table 1 |
| GCGCG = -6.86 kcal/mol | Exact match with Primer3 reference value | Primer3 Manual PRIMER_MAX_END_STABILITY |
| 3' stability formula | Ours = −end_oligodg for 2000/2000 random 1–12-mers incl. N | oligotm.c compiled (2026-10-01) |
| Primer-dimer compl_end | 3000/3000 vs dpal.c GLOBAL_END; 300/300 PRIMER_PAIR_0_COMPL_END; 600/600 SELF_END | dpal.c compiled; primer3-py check_primers |
| Thermodynamic screen | Formulas = libprimer3.cc characterize_pair / oligo_compl_thermod / oligo_hairpin / align_thermod; values match design_primers where the ntthal engines match primer3-py | primer3-py 2.3.1 |
| TATAT = -0.86 kcal/mol | Exact match with Primer3 reference value | Primer3 Manual PRIMER_MAX_END_STABILITY |
| Minimum hairpin loop = 3 nt | "loops fewer than three bases long are sterically impossible" | Wikipedia Stem-loop |
| Case-insensitive matching | Universal convention across all DNA tools | Standard bioinformatics practice |
| 3' end complementarity for primer-dimers | primer2's 3'-terminal bases must be the REVERSE COMPLEMENT of primer1's (the former window comparison tested parallel complementarity and flagged A₈/A₈) | Primer3 dpal GLOBAL_END / ntthal END1 |

### Design Parameters (Configurable, Not Assumptions)

| Parameter | Default | Rationale |
|-----------|---------|-----------|
| `minStemLength` | 4 bp | Configurable; stems < 4 bp are generally unstable at PCR temperatures |
| `minLoopLength` | 3 nt | Sterically required minimum (Wikipedia Stem-loop) |
| `minComplementarity` | 4 bp | Configurable; controls primer-dimer detection sensitivity |
| `MaxStructureTm` | 47 °C | Primer3 PRIMER_MAX_*_TH / PRIMER_PAIR_MAX_COMPL_*_TH defaults |

### Cross-Spec Note

The `EvaluatePrimer` threshold (`stability3Prime < -9`) in PRIMER-DESIGN-001 is unreachable
(most stable 5-mer GCGCG = -6.86), consistent with Primer3's default PRIMER_MAX_END_STABILITY = 100.
The structure screen of `EvaluatePrimer`/`DesignPrimers` is Primer3's thermodynamic one by default (2026-10-01).

## Test File Location

- **Canonical:** `Seqeron.Genomics.Tests/PrimerDesigner_PrimerStructure_Tests.cs`
- **Smoke tests:** Remain in `PrimerDesignerTests.cs` for integration testing
