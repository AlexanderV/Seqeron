# Allele Frequency

**Algorithm Group:** Population Genetics
**Test Unit:** POP-FREQ-001
**Last Updated:** 2026-02-01

---

## Overview

Allele frequency (also called gene frequency) is a fundamental concept in population genetics that measures the relative frequency of a specific allele at a genetic locus within a population. This measurement is essential for studying genetic variation, evolution, and the genetic structure of populations.

---

## Theory

### Definition

Allele frequency is the proportion of all copies of a gene that are of a particular allele type, expressed as a fraction or percentage (Gillespie, 2004).

### Calculation for Diploid Organisms

For a biallelic locus with alleles A and a in a diploid population:

Given genotype counts:
- $n_{AA}$ = number of homozygous major individuals
- $n_{Aa}$ = number of heterozygous individuals
- $n_{aa}$ = number of homozygous minor individuals

The allele frequencies are:

$$p = \frac{2n_{AA} + n_{Aa}}{2(n_{AA} + n_{Aa} + n_{aa})}$$

$$q = \frac{2n_{aa} + n_{Aa}}{2(n_{AA} + n_{Aa} + n_{aa})}$$

**Invariant:** $p + q = 1$

### Minor Allele Frequency (MAF)

MAF is the frequency of the less common allele at a locus:

$$MAF = \min(p, q) = \min(f, 1-f)$$

**Invariant:** $0 \leq MAF \leq 0.5$

MAF is commonly used to:
- Distinguish common variants (MAF ≥ 0.05) from rare variants (MAF < 0.05)
- Filter variants in genome-wide association studies
- Assess population diversity at specific loci

---

## Implementation

### Methods

#### CalculateAlleleFrequencies

```csharp
public static (double MajorFreq, double MinorFreq) CalculateAlleleFrequencies(
    int homozygousMajor,
    int heterozygous,
    int homozygousMinor)
```

**Parameters:**
- `homozygousMajor`: Count of AA individuals
- `heterozygous`: Count of Aa individuals
- `homozygousMinor`: Count of aa individuals

**Returns:** Tuple of (major allele frequency, minor allele frequency)

**Complexity:** O(1)

**Edge Cases:**
- Zero total samples → returns (0, 0)

#### CalculateMAF

```csharp
public static double CalculateMAF(IEnumerable<int> genotypes)
```

**Parameters:**
- `genotypes`: Collection where 0=hom ref, 1=het, 2=hom alt (VCF/PLINK convention)

**Returns:** Minor allele frequency (0 to 0.5)

**Complexity:** O(n) where n = number of genotypes

**Edge Cases:**
- Empty collection → returns 0
- Monomorphic locus → returns 0

#### FilterByMAF

```csharp
public static IEnumerable<Variant> FilterByMAF(
    IEnumerable<Variant> variants,
    double minMAF = 0.01,
    double maxMAF = 0.5)
```

**Parameters:**
- `variants`: Variants to filter
- `minMAF`: Minimum MAF threshold (default 0.01)
- `maxMAF`: Maximum MAF threshold (default 0.5)

**Returns:** Variants with MAF in range [minMAF, maxMAF]

**Complexity:** O(n) where n = number of variants

---

## Mathematical Invariants

1. **Sum Property:** For biallelic loci, $p + q = 1$
2. **Range Constraints:** $0 \leq p \leq 1$, $0 \leq q \leq 1$
3. **MAF Range:** $0 \leq MAF \leq 0.5$
4. **Allele Accounting:** $\text{major alleles} + \text{minor alleles} = \text{total alleles}$

---

## Example

From Wikipedia (Genotype Frequency):

A population of 100 four-o'-clock plants:
- 49 AA (red flowers)
- 42 Aa (pink flowers)
- 9 aa (white flowers)

Calculations:
- Total alleles = 2 × 100 = 200
- A alleles = 2 × 49 + 42 = 140
- a alleles = 2 × 9 + 42 = 60
- p(A) = 140/200 = 0.70
- q(a) = 60/200 = 0.30
- MAF = 0.30

---

## References

1. Gillespie, J.H. (2004). *Population genetics: a concise guide* (2nd ed.). Johns Hopkins University Press.
2. Wikipedia: Allele frequency. https://en.wikipedia.org/wiki/Allele_frequency
3. Wikipedia: Minor allele frequency. https://en.wikipedia.org/wiki/Minor_allele_frequency
4. Wikipedia: Genotype frequency. https://en.wikipedia.org/wiki/Genotype_frequency
5. The International HapMap Consortium (2005). A haplotype map of the human genome. *Nature*, 437(7063), 1299-1320.
