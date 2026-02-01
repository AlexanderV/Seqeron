# Evidence Document: POP-FREQ-001 (Allele Frequencies)

**Test Unit ID:** POP-FREQ-001  
**Algorithm:** Allele Frequency Calculation (Major/Minor Frequencies, MAF, Filtering)  
**Date:** 2026-02-01  
**Status:** Evidence Gathered  

---

## 1. Sources Consulted

### Primary Sources

| # | Source | Type | URL | Accessed |
|---|--------|------|-----|----------|
| 1 | Wikipedia: Allele frequency | Encyclopedia | https://en.wikipedia.org/wiki/Allele_frequency | 2026-02-01 |
| 2 | Wikipedia: Minor allele frequency | Encyclopedia | https://en.wikipedia.org/wiki/Minor_allele_frequency | 2026-02-01 |
| 3 | Wikipedia: Genotype frequency | Encyclopedia | https://en.wikipedia.org/wiki/Genotype_frequency | 2026-02-01 |
| 4 | Gillespie (2004) - Population Genetics: A Concise Guide | Textbook | ISBN 978-0-8018-8008-7 | 2026-02-01 |
| 5 | NDSU Population Genetics | Educational | https://www.ndsu.edu/pubweb/~mcclean/plsc431/popgen/popgen2.htm | 2026-02-01 |

---

## 2. Algorithm Definitions

### 2.1 Allele Frequency from Genotype Counts

**Definition (Wikipedia):**
> "Allele frequency, or gene frequency, is the relative frequency of an allele (variant of a gene) at a particular locus in a population, expressed as a fraction or percentage."

**Formula for Diploids (Wikipedia):**
For a locus with two alleles A and B:
- If f(AA), f(AB), f(BB) are genotype counts:
  - p = f(AA) + ½f(AB) = frequency of A
  - q = f(BB) + ½f(AB) = frequency of B
  - **Invariant: p + q = 1**

**Alternative formulation (using allele counting):**
- Total alleles = 2 × (homozygous_major + heterozygous + homozygous_minor)
- Major alleles = 2 × homozygous_major + heterozygous
- Minor alleles = 2 × homozygous_minor + heterozygous

**Example (Wikipedia Genotype Frequency):**
> "Consider a population of 100 four-o-clock plants with genotypes: 49 AA, 42 Aa, 9 aa"
> - Allele frequency of 'a' = (42 + 2×9) / (2×49 + 2×42 + 2×9) = 60/200 = 0.30
> - Allele frequency of 'A' = (2×49 + 42) / 200 = 140/200 = 0.70

### 2.2 Minor Allele Frequency (MAF)

**Definition (Wikipedia):**
> "Minor allele frequency (MAF) is the frequency at which the second most common allele occurs in a given population."

**Formula:**
- MAF = min(alt_frequency, 1 - alt_frequency)
- **Invariant: 0 ≤ MAF ≤ 0.5**

**Clinical/Research Significance (Wikipedia):**
> "SNPs with a minor allele frequency of 0.05 (5%) or greater were targeted by the HapMap project."
> "Rare variants (MAF < 0.05) appeared more frequently in coding regions than common variants (MAF > 0.05)."

### 2.3 Genotype Encoding

**Standard VCF/PLINK Convention:**
- 0 = homozygous reference (AA)
- 1 = heterozygous (Aa)
- 2 = homozygous alternate (aa)

**Allele Counting from Genotypes:**
- Alt alleles per individual: genotype value (0, 1, or 2)
- Total alt alleles = sum of all genotype values
- Total alleles = 2 × number of individuals
- Alt frequency = total_alt / total_alleles
- MAF = min(alt_freq, 1 - alt_freq)

---

## 3. Test Cases from Sources

### 3.1 Allele Frequency Calculation

| # | Test Case | Source | Input | Expected |
|---|-----------|--------|-------|----------|
| AF-1 | Wikipedia flower example | Genotype frequency article | AA=49, Aa=42, aa=9 | p=0.70, q=0.30 |
| AF-2 | All homozygous major | Derived from formula | AA=100, Aa=0, aa=0 | p=1.0, q=0.0 |
| AF-3 | All homozygous minor | Derived from formula | AA=0, Aa=0, aa=100 | p=0.0, q=1.0 |
| AF-4 | All heterozygous | Derived from formula | AA=0, Aa=100, aa=0 | p=0.5, q=0.5 |
| AF-5 | Equal genotypes (25-50-25) | Derived from HWE expectation | AA=25, Aa=50, aa=25 | p=0.5, q=0.5 |
| AF-6 | Zero samples | Edge case | AA=0, Aa=0, aa=0 | Handle gracefully (0,0) |

### 3.2 Minor Allele Frequency

| # | Test Case | Source | Input | Expected |
|---|-----------|--------|-------|----------|
| MAF-1 | Alt is minor (freq < 0.5) | MAF definition | genotypes with alt_freq=0.3 | MAF=0.3 |
| MAF-2 | Ref is minor (freq > 0.5) | MAF definition | genotypes with alt_freq=0.7 | MAF=0.3 |
| MAF-3 | Monomorphic ref (alt=0) | Edge case | all 0 genotypes | MAF=0 |
| MAF-4 | Monomorphic alt (alt=1) | Edge case | all 2 genotypes | MAF=0 |
| MAF-5 | Perfect 50/50 | Boundary | genotypes with alt_freq=0.5 | MAF=0.5 |
| MAF-6 | Empty genotypes | Edge case | empty array | MAF=0 |

### 3.3 MAF Filtering

| # | Test Case | Source | Rationale |
|---|-----------|--------|-----------|
| FLT-1 | Filter rare variants (MAF < 0.01) | HapMap convention | Common threshold |
| FLT-2 | Filter very rare (MAF < 0.05) | HapMap threshold | Common variant definition |
| FLT-3 | Filter common (MAF > 0.4) | Range filtering | Upper bound filtering |
| FLT-4 | Empty input | Edge case | Graceful handling |
| FLT-5 | All filtered | Edge case | Empty result |
| FLT-6 | None filtered | Edge case | All variants pass |

---

## 4. Edge Cases and Corner Cases

### 4.1 Numerical Edge Cases

| # | Case | Expected Behavior |
|---|------|-------------------|
| E-1 | Zero samples | Return (0, 0) or handle gracefully |
| E-2 | Single sample | Valid frequency calculation |
| E-3 | Very large population | No overflow |
| E-4 | Negative counts | ASSUMPTION: Not validated (caller responsibility) |

### 4.2 MAF Boundary Conditions

| # | Case | Expected Behavior |
|---|------|-------------------|
| E-5 | MAF at exact threshold | Include or exclude consistently |
| E-6 | MAF = 0 (monomorphic) | Should work, represents fixed allele |
| E-7 | MAF = 0.5 (balanced) | Maximum MAF value |

---

## 5. Invariants (Mathematical Properties)

### From Wikipedia Allele Frequency:
1. **Sum to one:** p + q = 1 (for biallelic locus)
2. **Range constraint:** 0 ≤ p ≤ 1, 0 ≤ q ≤ 1
3. **MAF range:** 0 ≤ MAF ≤ 0.5
4. **Symmetry:** MAF = min(freq, 1-freq)

### From Genotype Counting:
5. **Total alleles:** total = 2 × (hom_maj + het + hom_min)
6. **Allele accounting:** major_alleles + minor_alleles = total_alleles

---

## 6. Implementation Notes

### Current Implementation Analysis

**CalculateAlleleFrequencies:**
```
Total alleles = 2 × (hom_maj + het + hom_min)
Major alleles = 2 × hom_maj + het
Minor alleles = 2 × hom_min + het
Frequencies = alleles / total
```
✓ Matches Wikipedia formula

**CalculateMAF:**
```
Alt freq = sum(genotypes) / (2 × count)
MAF = min(alt_freq, 1 - alt_freq)
```
✓ Matches MAF definition

**FilterByMAF:**
```
For each variant:
  maf = min(AF, 1-AF)
  Include if minMAF ≤ maf ≤ maxMAF
```
✓ Standard filtering approach

---

## 7. Open Questions

None - algorithm is well-defined by sources.

---

## 8. References

1. Gillespie, J.H. (2004). Population genetics: a concise guide (2nd ed.). Johns Hopkins University Press. ISBN 978-0-8018-8008-7.
2. The International HapMap Consortium (2005). "A haplotype map of the human genome". Nature. 437 (7063): 1299–1320.
3. Sidore, C., et al. (2015). "Genome sequencing elucidates Sardinian genetic architecture". Nature Genetics. 47 (11): 1272–1281.
