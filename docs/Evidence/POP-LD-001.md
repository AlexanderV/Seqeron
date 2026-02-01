# Evidence: POP-LD-001 - Linkage Disequilibrium

## Test Unit
- **ID:** POP-LD-001
- **Area:** Population Genetics (PopGen)
- **Methods:** CalculateLD, FindHaplotypeBlocks

---

## Sources

### Primary Sources

1. **Wikipedia: Linkage disequilibrium**
   - URL: https://en.wikipedia.org/wiki/Linkage_disequilibrium
   - Accessed: 2026-02-01
   - Key information:
     - Formal definition: D = p_AB - p_A × p_B
     - D' (Lewontin's D prime) normalization formula
     - r² (correlation coefficient squared) formula
     - Decay of LD with recombination

2. **Wikipedia: Haplotype block**
   - URL: https://en.wikipedia.org/wiki/Haplotype_block
   - Accessed: 2026-02-01
   - Key information:
     - Definition based on high LD threshold
     - Gabriel et al. (2002) block definition method
     - Patil et al. (2001) percentage-based definition

3. **Lewontin, R.C. (1964)**
   - "The interaction of selection and linkage"
   - Genetics 49(1): 49–67
   - Introduced D' normalization

4. **Hill, W.G. & Robertson, A. (1968)**
   - "Linkage disequilibrium in finite populations"
   - Theoretical and Applied Genetics 38(6): 226–231
   - Introduced r² correlation measure

5. **Gabriel, S.B. et al. (2002)**
   - "The Structure of Haplotype Blocks in the Human Genome"
   - Science 296(5576): 2225–2229
   - Defined haplotype block detection using LD thresholds

---

## Mathematical Definitions

### Coefficient D (Linkage Disequilibrium)
From Wikipedia (Linkage disequilibrium):
```
D = p_AB - p_A × p_B
```
Where:
- p_AB = observed frequency of haplotype AB
- p_A = frequency of allele A
- p_B = frequency of allele B

### D' (Normalized D)
From Wikipedia (Lewontin 1964):
```
D' = D / D_max

D_max = {
  min(p_A × q_B, q_A × p_B)  when D < 0
  min(p_A × p_B, q_A × q_B)  when D ≥ 0
}
```
Where q_A = 1 - p_A and q_B = 1 - p_B

Range: -1 ≤ D' ≤ 1 (|D'| typically used, range 0 to 1)

### r² (Correlation Coefficient Squared)
From Wikipedia (Hill & Robertson 1968):
```
r² = D² / (p_A × q_A × p_B × q_B)
```
Range: 0 ≤ r² ≤ 1

---

## Test Datasets from Sources

### Perfect LD (Complete Association)
When alleles always co-occur:
- Genotypes: (0,0), (0,0), (1,1), (1,1), (2,2), (2,2)
- Expected: High r² (close to 1), high D'

### No LD (Random Association)
When alleles are independent:
- Genotypes: (0,2), (2,0), (1,1), (0,1), (2,1), (1,0)
- Expected: r² ≈ 0, D' ≈ 0

### Edge Cases from Theory
1. **Empty genotype list** → D = 0, r² = 0 (no data)
2. **Single genotype** → Cannot estimate LD meaningfully
3. **All homozygous major (0,0)** → No variation, r² = 0
4. **All homozygous minor (2,2)** → No variation, r² = 0
5. **Monomorphic locus** → Division by zero prevention (denominator = 0)

---

## Haplotype Block Detection

### Gabriel et al. (2002) Method
- Consecutive variants with r² ≥ threshold form a block
- Common threshold: 0.7–0.8 for strong LD
- Minimum 2 variants required for a block

### Edge Cases
1. **Single variant** → No blocks possible
2. **All variants in strong LD** → Single block spanning all
3. **No variants above threshold** → No blocks
4. **Non-contiguous blocks** → Multiple separate blocks

---

## Invariants

### CalculateLD Invariants
1. 0 ≤ r² ≤ 1
2. 0 ≤ |D'| ≤ 1
3. Distance is preserved from input
4. Variant IDs are preserved from input
5. Empty input → r² = 0, D' = 0

### FindHaplotypeBlocks Invariants
1. Block.Start ≤ Block.End
2. Block contains ≥ 2 variants
3. Blocks are non-overlapping
4. Blocks are ordered by position
5. All variants in a block have pairwise r² ≥ threshold with adjacent variants

---

## Implementation Notes

The current implementation:
1. Estimates haplotype frequency from genotypes (phasing assumption)
2. Uses genotype coding: 0 = homozygous major, 1 = heterozygous, 2 = homozygous minor
3. Assumes biallelic variants
4. FindHaplotypeBlocks uses adjacent-pair LD (simplified Gabriel method)

---

## Corner Cases Requiring Tests

| Case | Source | Expected Behavior |
|------|--------|-------------------|
| Empty genotypes | Mathematical definition | r² = 0, D' = 0 |
| Single genotype pair | Statistical requirement | Valid but unreliable |
| All identical genotypes | No variation | r² = 0 (no polymorphism) |
| Perfect correlation | Mathematical limit | r² = 1, D' = 1 |
| Perfect anti-correlation | Mathematical limit | r² = 1, D' = -1 or 1 |
| Monomorphic locus 1 | Division by zero | r² = 0 (protected) |
| Monomorphic locus 2 | Division by zero | r² = 0 (protected) |
| Single variant for blocks | Block definition | No blocks |
| Two variants, high LD | Block definition | One block |
| Two variants, low LD | Block definition | No blocks |

---

## Last Updated
2026-02-01
