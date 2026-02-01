# Linkage Disequilibrium

## Overview

Linkage disequilibrium (LD) is a measure of the non-random association of alleles at two or more loci in a population. When alleles at different loci occur together more or less frequently than expected by chance, they are said to be in linkage disequilibrium.

**Algorithm Group:** Population Genetics  
**Test Unit:** POP-LD-001  
**Complexity:** O(n) per variant pair, where n = number of genotypes

---

## Mathematical Foundation

### Coefficient D

The coefficient of linkage disequilibrium (D) measures the deviation from random association:

$$D = p_{AB} - p_A \times p_B$$

Where:
- $p_{AB}$ = observed frequency of haplotype AB
- $p_A$ = frequency of allele A at locus 1
- $p_B$ = frequency of allele B at locus 2

**Interpretation:**
- D = 0: Linkage equilibrium (alleles are independent)
- D > 0: Positive LD (alleles co-occur more than expected)
- D < 0: Negative LD (alleles co-occur less than expected)

### D' (Lewontin's Normalized D)

D' normalizes D by its theoretical maximum, allowing comparison across different allele frequencies:

$$D' = \frac{D}{D_{max}}$$

Where:

$$D_{max} = \begin{cases} 
\min(p_A \times q_B, q_A \times p_B) & \text{if } D < 0 \\
\min(p_A \times p_B, q_A \times q_B) & \text{if } D \geq 0
\end{cases}$$

With $q_A = 1 - p_A$ and $q_B = 1 - p_B$.

**Range:** $-1 \leq D' \leq 1$ (typically |D'| is used, range 0 to 1)

**Source:** Lewontin (1964)

### r² (Correlation Coefficient Squared)

r² measures the squared correlation between alleles at two loci:

$$r^2 = \frac{D^2}{p_A \times q_A \times p_B \times q_B}$$

**Range:** $0 \leq r^2 \leq 1$

**Interpretation:**
- r² = 0: No correlation (linkage equilibrium)
- r² = 1: Perfect correlation (complete LD)
- r² ≥ 0.8: Strong LD (commonly used threshold)

**Source:** Hill & Robertson (1968)

---

## Haplotype Block Detection

A haplotype block is a region of the genome with high internal LD and low recombination. Variants within a block tend to be inherited together.

### Gabriel Method (Simplified)

The implementation uses a simplified version of the Gabriel et al. (2002) method:

1. Order variants by genomic position
2. Calculate pairwise LD (r²) between adjacent variants
3. Group consecutive variants with r² ≥ threshold into blocks
4. Blocks must contain ≥ 2 variants

**Default threshold:** r² ≥ 0.7

**Source:** Gabriel et al. (2002)

---

## Implementation

### CalculateLD

```
Location: Seqeron.Genomics/PopulationGeneticsAnalyzer.cs
Method:   PopulationGeneticsAnalyzer.CalculateLD(...)
```

**Signature:**
```csharp
public static LinkageDisequilibrium CalculateLD(
    string variant1Id,
    string variant2Id,
    IEnumerable<(int Geno1, int Geno2)> genotypes,
    int distance)
```

**Parameters:**
- `variant1Id`: Identifier for first variant
- `variant2Id`: Identifier for second variant
- `genotypes`: Pairs of genotype values (0=hom.major, 1=het, 2=hom.minor)
- `distance`: Physical distance between variants (bp)

**Returns:** `LinkageDisequilibrium` record with D', r², and distance

**Algorithm:**
1. Calculate allele frequencies from genotypes
2. Estimate haplotype frequency (assuming random phase for heterozygotes)
3. Calculate D = p11 - p1 × p2
4. Normalize to get D' and r²

### FindHaplotypeBlocks

```
Location: Seqeron.Genomics/PopulationGeneticsAnalyzer.cs
Method:   PopulationGeneticsAnalyzer.FindHaplotypeBlocks(...)
```

**Signature:**
```csharp
public static IEnumerable<HaplotypeBlock> FindHaplotypeBlocks(
    IEnumerable<(string VariantId, int Position, IReadOnlyList<int> Genotypes)> variants,
    double ldThreshold = 0.7)
```

**Parameters:**
- `variants`: Ordered list of variants with positions and genotypes
- `ldThreshold`: Minimum r² to consider variants in same block (default: 0.7)

**Returns:** Enumerable of `HaplotypeBlock` records

---

## Edge Cases

| Case | Handling |
|------|----------|
| Empty genotype list | Returns r² = 0, D' = 0 |
| Monomorphic locus | Division by zero protected, returns r² = 0 |
| Single variant | No blocks returned |
| All variants in perfect LD | Single block spanning all variants |
| Distance preservation | Distance value passed through unchanged |

---

## Invariants

### CalculateLD
1. $0 \leq r^2 \leq 1$
2. $0 \leq |D'| \leq 1$
3. Variant IDs preserved in output
4. Distance preserved in output

### FindHaplotypeBlocks
1. Each block contains ≥ 2 variants
2. Block.Start ≤ Block.End
3. Blocks are non-overlapping
4. Blocks ordered by position

---

## References

1. Lewontin, R.C. (1964). "The interaction of selection and linkage." *Genetics* 49(1): 49–67.
2. Hill, W.G. & Robertson, A. (1968). "Linkage disequilibrium in finite populations." *Theor. Appl. Genet.* 38(6): 226–231.
3. Gabriel, S.B. et al. (2002). "The Structure of Haplotype Blocks in the Human Genome." *Science* 296(5576): 2225–2229.
4. Wikipedia: Linkage disequilibrium - https://en.wikipedia.org/wiki/Linkage_disequilibrium
5. Wikipedia: Haplotype block - https://en.wikipedia.org/wiki/Haplotype_block

---

## Implementation Notes

1. **Phasing assumption:** The implementation estimates haplotype frequencies from unphased genotypes using a probabilistic approach for heterozygotes.

2. **Genotype coding:** Uses standard coding where 0 = homozygous major allele, 1 = heterozygous, 2 = homozygous minor allele.

3. **Adjacent-pair LD:** The haplotype block finder uses only adjacent variant pairs for efficiency. This is a simplification of the full four-gamete test.
