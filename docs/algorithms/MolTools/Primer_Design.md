# Primer Pair Design

## Overview

Primer pair design is a computational process for selecting two oligonucleotide primers (forward and reverse) that can effectively amplify a target DNA region via PCR. The design process evaluates candidate primers against multiple quality criteria and selects optimal pairs.

## Algorithm Description

### Core Algorithm: `DesignPrimers`

The `DesignPrimers` method performs the following steps:

1. **Define Search Regions**
   - Forward primer region: up to 200 bp upstream of target start
   - Reverse primer region: up to 200 bp downstream of target end

2. **Generate Candidates**
   - Enumerate all possible primers within length constraints (18-25 bp)
   - Evaluate each candidate using `EvaluatePrimer`
   - Filter to retain only valid candidates

3. **Select Best Pair**
   - Rank candidates by score (higher = better)
   - Select highest-scoring forward primer
   - Select highest-scoring reverse primer
   - Verify pair compatibility:
     - Tm difference ≤ 5°C
     - No primer-dimer formation

4. **Return Result**
   - Return `PrimerPairResult` with selected primers, validity status, and product size

### Complexity

- **Time**: O(n²) where n is template length (iterating over positions × lengths)
- **Space**: O(k) where k is number of valid candidates

## Quality Criteria

### Standard Parameters (from Primer3, Addgene)

| Parameter | Min | Optimal | Max | Source |
|-----------|-----|---------|-----|--------|
| Length (bp) | 18 | 20 | 25 | Primer3, Wikipedia |
| GC Content (%) | 40 | 50 | 60 | Addgene, Implementation |
| Melting Temp (°C) | 55 | 60 | 65 | Implementation |
| Homopolymer Run | - | - | 4 | Implementation (Primer3: 5) |
| Dinucleotide Repeats | - | - | 4 | Implementation |

### Pair Compatibility

- **Tm Difference**: ≤ 5°C between forward and reverse primers
- **Primer-Dimer**: No significant 3' complementarity

## Scoring Algorithm

The `CalculatePrimerScore` method assigns a quality score:

```
score = 100
score -= |length - optimalLength| × 2
score -= |Tm - optimalTm| × 2
score -= |GC% - 50| × 0.5
score -= homopolymerLength × 5
score += 5 if 3' ends with G or C (GC clamp bonus)
score = max(0, score)
```

## Implementation Details

### Forward Primer Selection

Forward primers are extracted directly from the template sequence in the 5'→3' direction.

### Reverse Primer Selection

Reverse primers are:
1. Extracted from downstream region
2. Converted to reverse complement
3. This represents the primer sequence that will bind to the template's reverse strand

### 3' Stability Calculation

Uses nearest-neighbor thermodynamic parameters (SantaLucia, 1998) to calculate ΔG of the last 5 bases. More negative values indicate higher stability, which can cause mispriming.

## Related Algorithms

- **PRIMER-TM-001**: Melting temperature calculation (prerequisite)
- **PRIMER-STRUCT-001**: Hairpin and dimer detection (used in evaluation)

## References

1. **Wikipedia**: [Primer (molecular biology)](https://en.wikipedia.org/wiki/Primer_(molecular_biology)) - Standard primer design criteria
2. **Addgene**: [How to Design a Primer](https://www.addgene.org/protocols/primer-design/) - Protocol guidelines
3. **Primer3 Manual**: [primer3.org/manual.html](https://primer3.org/manual.html) - Comprehensive parameter documentation
4. **SantaLucia JR (1998)**: "A unified view of polymer, dumbbell and oligonucleotide DNA nearest-neighbor thermodynamics", PNAS 95:1460-65 - Thermodynamic calculations
