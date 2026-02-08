# Test Specification: PRIMER-STRUCT-001

## Test Unit

| Field | Value |
|-------|-------|
| **ID** | PRIMER-STRUCT-001 |
| **Title** | Primer Structure Analysis |
| **Area** | Molecular Tools |
| **Status** | ☑ Complete |
| **Last Updated** | 2026-01-22 |
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
| M1 | Homopolymer | Empty sequence | 0 | Primer3 default behavior |
| M2 | Homopolymer | No run (ACGT) | 1 | Primer3 PRIMER_MAX_POLY_X |
| M3 | Homopolymer | All same (AAAAAA) | 6 | Primer3 PRIMER_MAX_POLY_X |
| M4 | Homopolymer | Mixed case (AaAa) | 4 (case insensitive) | Implementation |
| M5 | Dinucleotide | Sequence < 4 bp | 0 | Implementation bounds |
| M6 | Dinucleotide | No repeat (ACGT) | 1 | Primer3 behavior |
| M7 | Dinucleotide | ACACACAC | 4 | Primer3 behavior |
| M8 | Hairpin | Short sequence | false | Stem-loop theory (min 3 bp loop) |
| M9 | Hairpin | Non-self-complementary | false | Wikipedia Stem-loop |
| M10 | Hairpin | Self-complementary | true | Wikipedia Stem-loop |
| M11 | Primer-dimer | Empty primer | false | Null guard |
| M12 | Primer-dimer | Non-complementary 3' ends | false | Wikipedia Primer-dimer |
| M13 | Primer-dimer | Complementary 3' ends | true | Wikipedia Primer-dimer |
| M14 | 3' Stability | Short sequence (<5 bp) | 0 | Primer3 5-mer standard |
| M15 | 3' Stability | GC-rich vs AT-rich | GC more negative | SantaLucia (1998) |

### SHOULD Tests (Recommended)

| ID | Category | Test Case | Expected | Source |
|----|----------|-----------|----------|--------|
| S1 | Homopolymer | Internal run (ACAAAAGT) | 4 | Common pattern |
| S2 | Dinucleotide | ATATATAT pattern | 4 | Common pattern |
| S3 | Hairpin | Custom minStemLength | Respects parameter | API contract |
| S4 | Primer-dimer | Self-dimer detection | true for A₈ vs A₈ | Self-complementarity |
| S5 | 3' Stability | Returns negative ΔG | Negative for valid inputs | Thermodynamics |

### COULD Tests (Optional)

| ID | Category | Test Case | Expected | Source |
|----|----------|-----------|----------|--------|
| C1 | Homopolymer | Long run at end | Detected correctly | Edge case |
| C2 | Dinucleotide | Multiple repeat types | Returns longest | Logic verification |
| C3 | Hairpin | Borderline length | Correct boundary behavior | Edge case |

## Audit Notes

### Existing Tests (PrimerDesignerTests.cs)

Located in `Seqeron.Genomics.Tests/PrimerDesignerTests.cs`:

| Test | Status | Action |
|------|--------|--------|
| FindLongestHomopolymer_* (4 tests) | Adequate | Move to canonical file |
| FindLongestDinucleotideRepeat_* (3 tests) | Adequate | Move to canonical file |
| HasHairpinPotential_* (3 tests) | Adequate | Move to canonical file |
| HasPrimerDimer_* (3 tests) | Adequate | Move to canonical file |
| Calculate3PrimeStability_* (2 tests) | Adequate | Move to canonical file |

### Consolidation Plan

1. **Create canonical test file:** `PrimerDesigner_PrimerStructure_Tests.cs`
2. **Move structure-related tests:** From PrimerDesignerTests.cs
3. **Leave smoke tests:** In PrimerDesignerTests.cs for EvaluatePrimer integration
4. ~~**Add missing tests:** Case insensitivity, boundary conditions, additional ΔG tests~~ ✅ All added

## Assumptions

| ID | Assumption | Rationale |
|----|------------|-----------|
| A1 | Case-insensitive matching | Standard DNA sequence handling |
| A2 | Simplified ΔG values used | Implementation uses subset of SantaLucia table |
| A3 | 3' check length = 8 for primer-dimer | Implementation-specific optimization |
| A4 | minStemLength default = 4 | Common bioinformatics threshold |

## Open Questions

*None - all critical behaviors are documented in sources.*

## Test File Location

- **Canonical:** `Seqeron.Genomics.Tests/PrimerDesigner_PrimerStructure_Tests.cs`
- **Smoke tests:** Remain in `PrimerDesignerTests.cs` for integration testing
