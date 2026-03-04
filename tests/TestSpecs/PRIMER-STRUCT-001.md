# Test Specification: PRIMER-STRUCT-001

## Test Unit

| Field | Value |
|-------|-------|
| **ID** | PRIMER-STRUCT-001 |
| **Title** | Primer Structure Analysis |
| **Area** | Molecular Tools |
| **Status** | ☑ Complete |
| **Last Updated** | 2026-03-04 |
| **Owner** | GitHub Copilot |

## Methods Under Test

| Method | Class | Type | Complexity |
|--------|-------|------|------------|
| `HasHairpinPotential(seq, minStemLength, minLoopLength)` | PrimerDesigner | Canonical | O(n²) <100bp, O(n) ≥100bp* |
| `HasPrimerDimer(primer1, primer2, minComp)` | PrimerDesigner | Canonical | O(n) |
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

## Invariants

1. **Homopolymer invariant:** Result ≥ 1 for non-empty sequences, 0 for empty
2. **Dinucleotide invariant:** Result ≥ 0; 0 for sequences < 4 bp
3. **Hairpin invariant:** Requires minimum length (2×stem + loop) to return true
4. **Stability invariant:** GC-rich 3' ends have more negative (stable) ΔG
5. **Primer-dimer invariant:** Returns false for empty primers

## Test Cases

### MUST Tests (Required)

| ID | Category | Test Case | Expected | Source |
|----|----------|-----------|----------|--------|
| M1 | Homopolymer | Null/empty sequence | 0 | Primer3 default behavior |
| M2 | Homopolymer | No run (ACGT) | 1 | Primer3 PRIMER_MAX_POLY_X |
| M3 | Homopolymer | All same (AAAAAA) | 6 | Primer3 PRIMER_MAX_POLY_X |
| M4 | Homopolymer | Mixed case (AaAaAa) | 6 (case insensitive) | Universal DNA convention |
| M5 | Dinucleotide | Null/empty/short (<4 bp) | 0 | Implementation bounds |
| M6 | Dinucleotide | No repeat (ACGT) | 1 | Primer3 behavior |
| M7 | Dinucleotide | ACACACAC | 4 | Primer3 behavior |
| M8 | Hairpin | Null/empty/too-short | false | Stem-loop theory (min 2×stem + loop) |
| M9 | Hairpin | Non-self-complementary | false | Wikipedia Stem-loop |
| M10 | Hairpin | Self-complementary | true | Wikipedia Stem-loop |
| M11 | Primer-dimer | Null/empty primer (either side) | false | Null guard |
| M12 | Primer-dimer | Non-complementary 3' ends | false | Wikipedia Primer-dimer |
| M13 | Primer-dimer | Complementary 3' ends (A₈ vs A₈) | true | Wikipedia Primer-dimer |
| M14 | 3' Stability | Null/empty/short (<5 bp) | 0 | Primer3 5-mer standard |
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
| S5 | Primer-dimer | Custom minComplementarity | Respects parameter | API contract |
| S6 | 3' Stability | Exact 5-base input (TACGT) | -3.57 kcal/mol | SantaLucia (1998) |
| S7 | 3' Stability | Case insensitive (GCGCG vs gcgcg) | Both = -6.86 | Universal DNA convention |

### COULD Tests (Optional)

| ID | Category | Test Case | Expected | Source |
|----|----------|-----------|----------|--------|
| C1 | Homopolymer | Run at end / multiple runs | Detected correctly | Edge case |
| C2 | Dinucleotide | Multiple repeat types | Returns longest | Logic verification |
| C3 | Hairpin | Long sequence (>100bp) suffix tree path | Correct detection | Performance optimization |
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
| `HasHairpinPotential_LongSequence_UsesSuffixTreeOptimization` | 1 | ✅ Covered | C3 |
| `HasHairpinPotential_LongSequenceNoHairpin_ReturnsFalse` | 1 | ✅ Covered | C3 |
| `HasPrimerDimer_NullOrEmptyPrimer_ReturnsFalse` | 4 | ✅ Covered | M11 (4 cases: null/empty × both sides) |
| `HasPrimerDimer_NonComplementary3Ends_ReturnsFalse` | 1 | ✅ Covered | M12 |
| `HasPrimerDimer_Complementary3Ends_ReturnsTrue` | 1 | ✅ Covered | M13 |
| `HasPrimerDimer_CustomMinComplementarity_RespectsParameter` | 1 | ✅ Covered | S5 |
| `Calculate3PrimeStability_InvalidInput_ReturnsZero` | 3 | ✅ Covered | M14 (null+empty+short merged) |
| `Calculate3PrimeStability_Exact5Bases_ProducesCorrectDeltaG` | 1 | ✅ Covered | S6 (was ❌ Missing) |
| `Calculate3PrimeStability_GcRich_MoreNegativeThanAtRich` | 1 | ✅ Covered | M15 (was ⚠ Weak: now exact -6.86/-0.86) |
| `Calculate3PrimeStability_MixedCase_ReturnsSameExactValue` | 1 | ✅ Covered | S7 (was ⚠ Weak: now checks -6.86) |
| `Calculate3PrimeStability_MostStable5mer_MatchesPrimer3` | 1 | ✅ Covered | M16 |
| `Calculate3PrimeStability_LeastStable5mer_MatchesPrimer3` | 1 | ✅ Covered | M17 |
| `PrimerStructureAnalysis_WellDesignedPrimer_ExactMetrics` | 1 | ✅ Covered | C4 (was ⚠ Weak: `True.Or.False` → exact values) |
| `PrimerStructureAnalysis_ProblematicPrimer_ExactMetrics` | 1 | ✅ Covered | C5 (was ⚠ Weak: `<-5.0` → exact -5.40) |

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
| TATAT = -0.86 kcal/mol | Exact match with Primer3 reference value | Primer3 Manual PRIMER_MAX_END_STABILITY |
| Minimum hairpin loop = 3 nt | "loops fewer than three bases long are sterically impossible" | Wikipedia Stem-loop |
| Case-insensitive matching | Universal convention across all DNA tools | Standard bioinformatics practice |
| 3' end complementarity for primer-dimers | "two primers anneal at their respective 3' ends" | Wikipedia Primer-dimer |

### Design Parameters (Configurable, Not Assumptions)

| Parameter | Default | Rationale |
|-----------|---------|-----------|
| `minStemLength` | 4 bp | Configurable; stems < 4 bp are generally unstable at PCR temperatures |
| `minLoopLength` | 3 nt | Sterically required minimum (Wikipedia Stem-loop) |
| `minComplementarity` | 4 bp | Configurable; controls primer-dimer detection sensitivity |
| `checkLength` (primer-dimer) | 8 bp | 3' region window; 6-10 bp needed for stable hybridization at PCR temperatures |

### Cross-Spec Note

The `EvaluatePrimer` threshold (`stability3Prime < -9`) in PRIMER-DESIGN-001 is now unreachable
(most stable 5-mer GCGCG = -6.86 with initiation). Review needed in PRIMER-DESIGN-001.

## Test File Location

- **Canonical:** `Seqeron.Genomics.Tests/PrimerDesigner_PrimerStructure_Tests.cs`
- **Smoke tests:** Remain in `PrimerDesignerTests.cs` for integration testing
