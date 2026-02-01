# CHROM-CENT-001 Evidence: Centromere Analysis

## Test Unit
- **ID:** CHROM-CENT-001
- **Area:** Chromosome
- **Canonical Method:** `ChromosomeAnalyzer.AnalyzeCentromere(...)`

## Sources Consulted

### Primary Sources

1. **Wikipedia - Centromere**
   - URL: https://en.wikipedia.org/wiki/Centromere
   - Retrieved: 2026-01-31
   - Key information:
     - Centromere links sister chromatids during cell division
     - Creates short arm (p) and long arm (q)
     - Position classifications: Metacentric, Submetacentric, Acrocentric, Telocentric
     - Regional centromeres contain large arrays of repetitive DNA (alpha-satellite)
     - Human centromeric repeat unit is called α-satellite (alphoid)

2. **Wikipedia - Karyotype**
   - URL: https://en.wikipedia.org/wiki/Karyotype
   - Retrieved: 2026-01-31
   - Key information:
     - Human chromosome groups based on centromere position
     - Group A (1-3): Large, metacentric or submetacentric
     - Group D (13-15): Medium-sized, acrocentric
     - Group G (21-22, Y): Very small, acrocentric

3. **Wikipedia - Chromosome**
   - URL: https://en.wikipedia.org/wiki/Chromosome
   - Retrieved: 2026-01-31
   - Key information:
     - Centromere region has constitutive heterochromatin with repetitive sequences
     - p arm (short) named from French "petit"
     - q arm (long) follows p in Latin alphabet

### Key Definitions (from Sources)

| Centromere Position | Arm Ratio | Classification | Description |
|---------------------|-----------|----------------|-------------|
| Medial (1.0-1.7)    | ~1.0      | Metacentric    | p and q arms approximately equal |
| Submedial (≤3.0)    | ~2.0      | Submetacentric | Arms close in length but not equal |
| Subterminal (3.1-6.9) | ~5.0    | Subtelocentric | Arms unequal, centromere toward end |
| Terminal (≥7.0)     | ~7.0+     | Acrocentric    | One arm much shorter than other |
| Terminal (∞)        | ∞         | Telocentric    | Centromere at very end |

Source: Levan A, Fredga K, Sandberg AA (1964). "Nomenclature for centromeric position on chromosomes". Hereditas.

### Human Chromosome Centromere Positions (from Wikipedia)

| Chromosome | Centromere Position (Mb) | Type |
|------------|-------------------------|------|
| 1          | 125.0                   | Metacentric |
| 2          | 93.3                    | Submetacentric |
| 3          | 91.0                    | Metacentric |
| 13         | 17.9                    | Acrocentric |
| 21         | 13.2                    | Acrocentric |
| Y          | 12.5                    | Acrocentric |

### Alpha-Satellite DNA (from Wikipedia)

- Primary centromeric repeat in humans
- ~171 bp monomer repeat
- Forms higher-order repeat arrays
- Associated with heterochromatin
- High repetitive content with low GC variability

## Corner Cases (Documented)

1. **Empty/null sequence**: No centromere can be detected
2. **Sequence shorter than window size**: Cannot perform analysis
3. **No repetitive regions**: Unknown centromere type
4. **Homogeneous sequence**: May falsely detect as centromeric if repetitive

## Testing Methodology

Based on the literature:

1. **Position-based classification test**: Verify classification thresholds match established nomenclature
2. **Alpha-satellite detection**: Verify recognition of repetitive patterns characteristic of centromeric regions
3. **Boundary conditions**: Empty, short, and edge-case sequences
4. **Invariants**:
   - Start ≤ End when centromere found
   - Length = End - Start
   - Type is one of: Metacentric, Submetacentric, Acrocentric, Telocentric, Unknown

## Implementation Notes

The implementation uses:
- Sliding window approach with k-mer frequency analysis
- GC content variability as discriminating feature (centromeres have low GC variability)
- Repeat content estimation using k-mer counting (k=15)
- Position-based classification matching standard nomenclature

## ASSUMPTIONS

1. **ASSUMPTION**: Window size of 100,000 bp is appropriate for detecting centromeric regions in typical genomic sequences
2. **ASSUMPTION**: Alpha-satellite content threshold of 0.3 is reasonable for identifying centromeric candidates
3. **ASSUMPTION**: The DetermineCentromereType position thresholds (< 0.15 = Acrocentric, 0.35-0.65 = Metacentric, etc.) align with biological conventions
