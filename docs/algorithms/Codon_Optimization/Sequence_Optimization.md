# Codon Optimization (Sequence Optimization)

## Overview

Codon optimization is a technique used in molecular biology to enhance heterologous gene expression by replacing codons in a coding sequence with synonymous codons that are more frequently used by the target host organism. The method preserves the encoded protein while improving translation efficiency.

## Biological Background

Different organisms exhibit codon usage bias—a preference for certain synonymous codons over others. This bias correlates with tRNA abundance and affects translation rate and accuracy. When expressing a gene from one organism in another (heterologous expression), the original codon usage may be suboptimal for the host's tRNA pool.

## Algorithm Description

### Input
- Coding sequence (DNA or RNA)
- Target organism codon usage table
- Optimization strategy
- Optional: GC content constraints

### Output
- Optimized sequence
- Original and optimized CAI values
- GC content metrics
- List of codon changes

### Optimization Strategies

| Strategy | Description | Use Case |
|----------|-------------|----------|
| MaximizeCAI | Select most frequent codon for each amino acid | Maximum expression level |
| BalancedOptimization | Balance CAI with GC content (40-60%) | Avoid extreme GC |
| HarmonizeExpression | Match host codon usage distribution | Proper protein folding |
| AvoidRareCodons | Replace only codons below threshold | Minimal changes |
| MinimizeSecondary | Reduce mRNA secondary structure potential | Improved translation initiation |

## Mathematical Foundation

### Codon Adaptation Index (CAI)

The CAI measures how well a gene's codon usage matches the preferred codons of highly expressed genes.

**Relative Adaptiveness** for codon $i$:
$$w_i = \frac{f_i}{\max(f_j)}$$

where $f_i$ is the frequency of codon $i$ and $\max(f_j)$ is the maximum frequency among synonymous codons.

**CAI** (geometric mean):
$$CAI = \left(\prod_{i=1}^{L} w_i\right)^{1/L}$$

where $L$ is the number of codons (excluding stop codons in some implementations).

### Properties

- **Range**: $0 < CAI \leq 1$
- **Maximum (1.0)**: All codons are optimal
- **Invariant**: Protein sequence is unchanged

## Implementation Details

```
OptimizeSequence(codingSequence, targetOrganism, strategy):
    1. Convert DNA to RNA (T → U)
    2. Trim to complete codons (length % 3 == 0)
    3. For each codon:
       a. Identify amino acid
       b. If stop codon: preserve
       c. Select replacement based on strategy
    4. Calculate original and optimized CAI
    5. Return optimization result
```

### Complexity

- **Time**: O(n) where n = sequence length
- **Space**: O(n) for storing optimized sequence

## Organism-Specific Codon Preferences

### E. coli K12 (selected codons)

| Amino Acid | Preferred | Frequency | Rare | Frequency |
|------------|-----------|-----------|------|-----------|
| Leucine | CUG | 0.47 | CUA | 0.04 |
| Arginine | CGU | 0.36 | AGG | 0.04 |
| Proline | CCG | 0.49 | CCC | 0.13 |

### S. cerevisiae (Yeast)

| Amino Acid | Preferred | Frequency | Rare | Frequency |
|------------|-----------|-----------|------|-----------|
| Leucine | UUG | 0.29 | CUC | 0.06 |
| Arginine | AGA | 0.48 | CGG | 0.04 |
| Proline | CCA | 0.42 | CCG | 0.12 |

## Edge Cases

1. **Empty sequence**: Returns empty result with CAI = 0
2. **Single amino acid codons**: Met (AUG) and Trp (UGG) cannot be optimized
3. **Stop codons**: Preserved as-is, not optimized
4. **Incomplete codons**: Trimmed to complete triplets

## References

1. Sharp, P.M. & Li, W.H. (1987). "The codon adaptation index-a measure of directional synonymous codon usage bias, and its potential applications." *Nucleic Acids Research* 15(3): 1281–1295.

2. Plotkin, J.B. & Kudla, G. (2011). "Synonymous but not the same: The causes and consequences of codon bias." *Nature Reviews Genetics* 12(1): 32–42.

3. Mignon, C. et al. (2018). "Codon harmonization - going beyond the speed limit for protein expression." *FEBS Letters* 592(9): 1554–1564.

4. Athey, J. et al. (2017). "A new and updated resource for codon usage tables." *BMC Bioinformatics* 18: 391.

## See Also

- [CODON-CAI-001](../CODON-CAI-001.md) - CAI Calculation
- [CODON-RARE-001](../CODON-RARE-001.md) - Rare Codon Detection
- [CODON-USAGE-001](../CODON-USAGE-001.md) - Codon Usage Analysis
