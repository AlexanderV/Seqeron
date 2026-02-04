# CODON-OPT-001: Codon Optimization Evidence

## Test Unit
- **ID:** CODON-OPT-001
- **Area:** Codon
- **Algorithm:** Sequence Optimization (OptimizeSequence)

## Sources

### Primary Sources

1. **Wikipedia - Codon Usage Bias**
   - URL: https://en.wikipedia.org/wiki/Codon_usage_bias
   - Key concepts:
     - Codon optimization adjusts codons to match host tRNA abundances for heterologous expression
     - Strategies include: local mRNA folding, codon pair bias, codon ramp, codon harmonization
     - Rare codons can lead to inefficient use and depletion of ribosomes
     - Protein folding is affected by translation rates (cotranslational folding)

2. **Wikipedia - Codon Adaptation Index (CAI)**
   - URL: https://en.wikipedia.org/wiki/Codon_Adaptation_Index
   - Reference: Sharp, P.M. & Li, W.H. (1987). "The codon adaptation index-a measure of directional synonymous codon usage bias, and its potential applications." Nucleic Acids Research. 15(3): 1281–1295.
   - Key concepts:
     - CAI measures deviation from reference gene set
     - Relative adaptiveness (wi) = fi / max(fj) for synonymous codons
     - CAI = geometric mean of weights: CAI = (∏wi)^(1/L)
     - Range: 0 < CAI ≤ 1

3. **Kazusa Codon Usage Database**
   - URL: https://www.kazusa.or.jp/codon/
   - Provides organism-specific codon usage tables

### Key Academic References

4. **Sharp & Li (1987)**
   - "The codon adaptation index-a measure of directional synonymous codon usage bias"
   - Original CAI definition
   - PMC 340524, PMID 3547335

5. **Plotkin & Kudla (2011)**
   - "Synonymous but not the same: The causes and consequences of codon bias"
   - Nature Reviews Genetics 12(1): 32-42
   - Comprehensive review of codon optimization strategies

6. **Mignon et al. (2018)**
   - "Codon harmonization - going beyond the speed limit for protein expression"
   - FEBS Letters 592(9): 1554-1564
   - Describes HarmonizeExpression strategy

## Documented Behaviors

### Optimization Strategies (from Wikipedia & implementation)

1. **MaximizeCAI**: Use most frequent codons for each amino acid
   - Source: Sharp & Li (1987), CAI definition

2. **BalancedOptimization**: Balance CAI with GC content constraints
   - Source: Wikipedia - mRNA secondary structure affects translation

3. **HarmonizeExpression**: Match host codon usage distribution
   - Source: Mignon et al. (2018)

4. **AvoidRareCodons**: Replace only rare codons (frequency < threshold)
   - Source: Wikipedia - consecutive rare codons inhibit translation

5. **MinimizeSecondary**: Avoid mRNA secondary structures
   - Source: Wikipedia - 5' secondary structure inhibits translation

### Invariants (from theory)

1. **Protein preservation**: Optimization must not change encoded protein
   - Source: Definition of synonymous codon substitution

2. **CAI range**: 0 < CAI ≤ 1
   - Source: Sharp & Li (1987) - geometric mean of values 0 < wi ≤ 1

3. **Methionine/Tryptophan unchanged**: AUG and UGG are unique codons
   - Source: Standard genetic code

## Test Datasets

### Organism-Specific Codon Tables (from Kazusa/implementation)

1. **E. coli K12** preferred codons:
   - Leu: CUG (0.47)
   - Arg: CGU/CGC (0.36 each), AGA/AGG rare (0.07/0.04)
   - Pro: CCG (0.49)
   
2. **S. cerevisiae (Yeast)** preferred codons:
   - Leu: UUA/UUG (0.28/0.29)
   - Arg: AGA (0.48)
   - Pro: CCA (0.42)

3. **H. sapiens (Human)** preferred codons:
   - Leu: CUG (0.41)
   - Arg: AGA/AGG (0.20 each)
   - Pro: CCC (0.33)

### Edge Cases (from theory)

1. **Empty sequence**: Should return empty result
2. **Incomplete codons**: Trim to complete codons (length % 3 == 0)
3. **DNA input (T)**: Convert to RNA (U)
4. **Lowercase input**: Case-insensitive processing
5. **Stop codons**: Preserved, not optimized

## Testing Methodology

1. **Unit tests**: Verify individual method behaviors
2. **Invariant tests**: Verify protein preservation across all strategies
3. **CAI mathematical tests**: Verify CAI calculation formula
4. **Organism-specific tests**: Verify different organisms yield different optimizations
5. **Edge case coverage**: Empty, single codon, all-same codons

## Known Failure Modes

1. **Invalid codon**: Unknown codon should translate to 'X' or error
2. **Non-RNA characters**: Should handle gracefully
3. **Stop codon in middle**: May terminate protein prematurely

## Implementation Notes

- Implementation uses RNA notation (U not T)
- Automatically converts T → U
- Trims to complete codons
- GC content balancing in BalancedOptimization strategy (40-60% target)

## Date
2026-02-04
