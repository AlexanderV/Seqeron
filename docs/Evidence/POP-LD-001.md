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

### CalculateLD Method (Hill 1974 Composite LD Estimator)

The implementation uses the **Hill (1974) composite linkage disequilibrium estimator**,
which calculates r² directly from genotype correlation without requiring phase information:

**Reference:** Hill WG (1974) "Estimation of linkage disequilibrium in randomly mating populations" Heredity 33:229-239

**Formula:**
```
r² = Cov(X₁, X₂)² / (Var(X₁) × Var(X₂))
```

Where:
- X₁, X₂ are genotype values (0, 1, 2 = count of minor alleles)
- Cov = Σ(X₁ᵢ - μ₁)(X₂ᵢ - μ₂) / n
- Var = Σ(Xᵢ - μ)² / n

**Advantages of this approach:**
1. Does not require haplotype phase information
2. Mathematically equivalent to squared Pearson correlation
3. For identical genotypes: r² = 1.0 (perfect LD)
4. For independent genotypes: r² = 0.0 (no LD)
5. Handles both positive and negative correlations (r² captures magnitude)

**Genotype coding:**
- 0 = homozygous reference (AA)
- 1 = heterozygous (Aa)
- 2 = homozygous alternate (aa)

**D' calculation:**
D' is estimated using the Lewontin normalization based on allele frequencies
derived from genotype means.

### FindHaplotypeBlocks Method
Uses adjacent-pair LD (simplified Gabriel et al. 2002 method):
- Consecutive variants with r² ≥ threshold form a block
- Default threshold: 0.7

---

## Test Coverage

### Reference Data Tests (Added 2026-02-05)
Tests validated against published literature:

**HapMap Consortium (2005):**
- `CalculateLD_HapMapInterpretation_HighLDThreshold` - validates r² = 1.0 for identical genotypes
- Tests reflect HapMap criterion: r² > 0.8 for tag SNP selection

**Hill & Robertson (1968):**
- `CalculateLD_HillRobertsonFormula_CorrectRange` - validates 0 ≤ r² ≤ 1
- Tests confirm r² formula: D² / (pA × qA × pB × qB)

**Lewontin (1964):**
- `CalculateLD_LewontinDPrime_NormalizedCorrectly` - validates |D'| ≤ 1

**Gabriel et al. (2002):**
- `FindHaplotypeBlocks_GabrielCriteria_HighLDBlocksDetected` - validates block detection algorithm

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
2026-02-05

## Change History
- **2026-02-05**: Added reference data tests from HapMap, Hill & Robertson, Gabriel et al.
- **2026-02-05**: Fixed CalculateLD implementation. Replaced haplotype frequency estimation 
  (which incorrectly gave r²≈0.44 for identical genotypes) with Hill (1974) composite LD estimator
  using genotype correlation. Now correctly returns r²=1.0 for perfect LD.
- **2026-02-01**: Initial documentation.
