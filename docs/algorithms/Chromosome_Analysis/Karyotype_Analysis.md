# Karyotype Analysis

**Algorithm Category:** Chromosome Analysis  
**Test Unit:** CHROM-KARYO-001  
**Date:** 2026-02-01  

---

## Overview

Karyotype analysis examines the complete set of chromosomes in an organism to determine chromosome count, identify sex chromosomes, detect aneuploidies (abnormal chromosome numbers), and estimate ploidy level.

---

## Biological Background

### Karyotype Definition
A **karyotype** is the complete set of chromosomes in a cell, characterized by:
- Total chromosome number (e.g., 2n = 46 in humans)
- Chromosome sizes and shapes
- Position of centromeres
- Banding patterns (G-banding, etc.)

*Source: Wikipedia - Karyotype*

### Ploidy Levels
- **Haploid (n):** Single chromosome set (gametes)
- **Diploid (2n):** Two chromosome sets (most somatic cells)
- **Polyploid (>2n):** Multiple chromosome sets (common in plants)

*Source: Wikipedia - Ploidy*

### Aneuploidy
Abnormal chromosome number, typically caused by nondisjunction during meiosis:
- **Monosomy:** Missing one chromosome (2n - 1)
- **Trisomy:** Extra chromosome (2n + 1)

Examples:
- Down syndrome: Trisomy 21 (47 chromosomes)
- Turner syndrome: Monosomy X (45 chromosomes)

*Source: Wikipedia - Aneuploidy*

---

## Algorithms

### 1. AnalyzeKaryotype

**Purpose:** Analyze chromosome data to produce a karyotype description.

**Input:**
- `chromosomes`: List of (Name, Length, IsSexChromosome) tuples
- `expectedPloidyLevel`: Expected ploidy (default: 2 for diploid)

**Output:**
- `Karyotype` record containing:
  - TotalChromosomes
  - AutosomeCount
  - SexChromosomes
  - TotalGenomeSize
  - MeanChromosomeLength
  - PloidyLevel
  - HasAneuploidy
  - Abnormalities list

**Algorithm:**
```
1. Separate sex chromosomes from autosomes
2. Group autosomes by base chromosome name
   - Strip copy suffixes (e.g., "chr1_1" → "chr1")
3. For each autosome group:
   - If count < expectedPloidy: Record "Monosomy {chr}"
   - If count > expectedPloidy: Record "Trisomy {chr}"
4. Calculate statistics:
   - TotalGenomeSize = Σ(lengths)
   - MeanLength = TotalGenomeSize / TotalChromosomes
5. Return Karyotype record
```

**Complexity:** O(n) where n = number of chromosomes

### 2. DetectPloidy

**Purpose:** Estimate ploidy level from read depth data.

**Input:**
- `normalizedDepths`: Sequence of normalized read depth values
- `expectedDiploidDepth`: Expected depth for diploid (default: 1.0)

**Output:**
- `PloidyLevel`: Estimated ploidy (1-8)
- `Confidence`: Confidence score (0-1)

**Algorithm:**
```
1. If empty input: Return (2, 0)  // Default diploid, zero confidence
2. Calculate median depth from sorted values
3. Compute ratio = medianDepth / expectedDiploidDepth
4. Estimate ploidy = round(ratio × 2)
5. Clamp ploidy to [1, 8]
6. Calculate confidence:
   - fractionalPart = |ratio × 2 - ploidy|
   - confidence = 1.0 - fractionalPart × 2
7. Return (ploidy, max(0, confidence))
```

**Complexity:** O(n log n) due to median calculation

---

## Implementation Notes

### Chromosome Naming Convention
The implementation expects chromosome names with optional copy suffixes:
- `chr1_1`, `chr1_2` for diploid copies of chromosome 1
- `chrX`, `chrY` for sex chromosomes (no suffix)

### Ploidy Detection Rationale
Read depth-based ploidy detection is standard in NGS analysis:
- Diploid regions have 1× normalized depth
- Tetraploid regions have 2× normalized depth
- Aneuploid regions deviate from expected ratios

### Confidence Calculation
Confidence reflects how close the observed ratio is to a clean integer ploidy:
- Ratio of exactly 1.0 → 100% confidence for diploid
- Ratio of 1.25 → reduced confidence (between diploid and triploid)

---

## References

1. Wikipedia. "Karyotype." https://en.wikipedia.org/wiki/Karyotype
2. Wikipedia. "Ploidy." https://en.wikipedia.org/wiki/Ploidy
3. Wikipedia. "Aneuploidy." https://en.wikipedia.org/wiki/Aneuploidy
4. Tjio, J.H.; Levan, A. (1956). "The chromosome number of man." Hereditas 42: 1-6.

---

## Deviations from Theory

None identified. Implementation follows standard biological definitions.

---

## See Also

- [Centromere_Analysis.md](Centromere_Analysis.md)
- [Telomere_Analysis.md](Telomere_Analysis.md)
