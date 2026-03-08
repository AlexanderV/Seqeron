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
     - Position classifications: Metacentric, Submetacentric, Subtelocentric, Acrocentric, Telocentric
     - Regional centromeres contain large arrays of repetitive DNA (alpha-satellite)
     - Human centromeric repeat unit is called α-satellite (alphoid), ~171 bp monomer
     - "Telocentric chromosomes... are not present in humans but can form through cellular chromosomal errors"

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

### Levan (1964) Classification — Arm Ratio Thresholds

Source: Levan A, Fredga K, Sandberg AA (1964). "Nomenclature for centromeric position on chromosomes". Hereditas. 52(2):201-220.

As cited on Wikipedia (Centromere article), the classification uses the arm ratio q/p (long arm / short arm):

| Centromere Position    | Arms length ratio (q/p) | Sign | Classification   |
|------------------------|-------------------------|------|------------------|
| Medial sensu stricto   | 1.0 – 1.6              | M    | Metacentric      |
| Medial region          | 1.7                     | m    | Metacentric      |
| Submedial              | 3.0                     | sm   | Submetacentric   |
| Subterminal            | 3.1 – 6.9              | st   | Subtelocentric   |
| Terminal region        | 7.0                     | t    | Acrocentric      |
| Terminal sensu stricto | ∞                       | T    | Telocentric      |

Implementation boundaries (derived from Levan table):
- ratio ≤ 1.7 → Metacentric (M + m)
- ratio (1.7, 3.0] → Submetacentric (sm)
- ratio (3.0, 7.0) → Subtelocentric (st)
- ratio ≥ 7.0 → Acrocentric (t)
- p = 0 → Telocentric (T)

**Boundary values:** 1.7 → Metacentric, 3.0 → Submetacentric, 7.0 → Acrocentric (per Levan table entries).

### Human Chromosome Centromere Positions (from Wikipedia)

| Chromosome | Centromere Position (Mb) | Type |
|------------|-------------------------|------|
| 1          | 125.0                   | Metacentric |
| 2          | 93.3                    | Submetacentric |
| 3          | 91.0                    | Metacentric |
| 13         | 17.9                    | Acrocentric |
| 21         | 13.2                    | Acrocentric |
| Y          | 12.5                    | Acrocentric |

Note: Practical human karyotype classifications are based on cytogenetic (microscopic) observation, not genomic coordinate ratios. The Levan thresholds applied to sequence-derived positions may yield slightly different classifications for borderline chromosomes.

### Alpha-Satellite DNA (from Wikipedia)

- Primary centromeric repeat in humans
- ~171 bp monomer repeat
- Forms higher-order repeat arrays
- Associated with heterochromatin

## Testing Methodology

Based on the literature:

1. **Arm-ratio-based classification test**: Verify classification thresholds match Levan (1964) nomenclature using arm ratio (q/p)
2. **Alpha-satellite detection**: Verify recognition of repetitive patterns characteristic of centromeric regions
3. **Boundary conditions**: Empty, short, and edge-case sequences
4. **Invariants**:
   - Start ≤ End when centromere found
   - Length = End - Start
   - Type is one of: Metacentric, Submetacentric, Subtelocentric, Acrocentric, Telocentric, Unknown

## Implementation Notes

The implementation uses:
- Sliding window approach with k-mer frequency analysis
- GC content variability as discriminating feature (centromeres have low GC variability)
- Repeat content estimation using k-mer counting (k=15)
- Arm ratio classification matching Levan (1964) nomenclature exactly
