# Diversity Statistics

**Algorithm Group:** Population Genetics
**Test Unit:** POP-DIV-001
**Last Updated:** 2026-02-01

---

## Overview

Diversity statistics are fundamental measures in population genetics used to quantify genetic variation within populations. The three primary metrics covered in this module are:

1. **Nucleotide Diversity (π)** - average pairwise differences between sequences
2. **Watterson's Theta (θ_W)** - diversity estimate based on segregating sites
3. **Tajima's D** - test statistic comparing the two diversity estimators

These metrics are essential for understanding evolutionary processes, detecting natural selection, and analyzing population history.

---

## Theory

### Nucleotide Diversity (π)

**Definition (Nei & Li, 1979):**
Nucleotide diversity is the average number of nucleotide differences per site between two randomly chosen sequences from a population.

**Formula:**
$$\pi = \frac{\sum_{i<j} d_{ij}}{\binom{n}{2} \cdot L}$$

Where:
- $d_{ij}$ = number of nucleotide differences between sequences $i$ and $j$
- $n$ = number of sequences
- $L$ = sequence length
- $\binom{n}{2} = \frac{n(n-1)}{2}$ = number of pairwise comparisons

**Properties:**
- Range: $0 \leq \pi \leq 1$
- π = 0 for monomorphic populations (all sequences identical)
- π = 1 when all positions differ between any two sequences

### Watterson's Theta (θ_W)

**Definition (Watterson, 1975):**
An estimator of the population mutation rate (4Nμ) based on the number of segregating sites.

**Formula:**
$$\theta_W = \frac{S}{a_n \cdot L}$$

Where:
- $S$ = number of segregating (polymorphic) sites
- $L$ = sequence length
- $a_n = \sum_{i=1}^{n-1} \frac{1}{i}$ = $(n-1)$-th harmonic number

**Properties:**
- Unbiased under neutral model assumptions
- Simpler to compute than nucleotide diversity
- Sensitive to rare variants

### Tajima's D

**Definition (Tajima, 1989):**
A test statistic that compares nucleotide diversity (π) with Watterson's theta (θ_W) to detect departures from neutral evolution.

**Formula:**
$$D = \frac{\pi - \theta_W}{\sqrt{Var(\pi - \theta_W)}}$$

The variance calculation involves:
- $a_1 = \sum_{i=1}^{n-1} \frac{1}{i}$, $a_2 = \sum_{i=1}^{n-1} \frac{1}{i^2}$
- $b_1 = \frac{n+1}{3(n-1)}$
- $b_2 = \frac{2(n^2 + n + 3)}{9n(n-1)}$
- $c_1 = b_1 - \frac{1}{a_1}$
- $c_2 = b_2 - \frac{n+2}{a_1 n} + \frac{a_2}{a_1^2}$
- $e_1 = \frac{c_1}{a_1}$
- $e_2 = \frac{c_2}{a_1^2 + a_2}$
- $Var = e_1 S + e_2 S(S-1)$

**Interpretation:**

| D Value | Interpretation |
|---------|----------------|
| D ≈ 0 | Neutral evolution, mutation-drift equilibrium |
| D < 0 | Excess of rare alleles (selective sweep, population expansion) |
| D > 0 | Deficit of rare alleles (balancing selection, population contraction) |

**Significance threshold:** Values beyond ±2 are generally considered significant.

### Segregating Sites

A segregating site is a position in a sequence alignment where at least one sequence differs from others. The count of segregating sites (S) is used in both θ_W and Tajima's D calculations.

### Heterozygosity

**Expected Heterozygosity:**
$$H_e = 1 - \sum_{i=1}^{k} f_i^2$$

Where $f_i$ is the frequency of the $i$-th allele. For sequence data, this is calculated per position and averaged.

**Observed Heterozygosity:**
For sequence data, this is adapted as the fraction of polymorphic sites.

---

## Implementation

### Location

`Seqeron.Genomics.PopulationGeneticsAnalyzer`

### Methods

| Method | Description | Complexity |
|--------|-------------|------------|
| `CalculateNucleotideDiversity(seqs)` | Computes π (nucleotide diversity) | O(n² × m) |
| `CalculateWattersonTheta(S, n, L)` | Computes θ_W from segregating sites | O(n) |
| `CalculateTajimasD(π, θ, S, n)` | Computes Tajima's D test statistic | O(n) |
| `CalculateDiversityStatistics(seqs)` | Computes all diversity metrics | O(n² × m) |

### Return Type

```csharp
public readonly record struct DiversityStatistics(
    double NucleotideDiversity,    // π
    double WattersonTheta,         // θ_W
    double TajimasD,               // D
    int SegregratingSites,         // S
    int SampleSize,                // n
    double HeterozygosityObserved, // H_o
    double HeterozygosityExpected  // H_e
);
```

### Edge Cases

| Case | Behavior |
|------|----------|
| Empty input (n = 0) | Returns zeros |
| Single sequence (n = 1) | Returns zeros (undefined for comparison) |
| n = 2 for Tajima's D | Returns 0 (requires n ≥ 3) |
| No segregating sites (S = 0) | π = 0, θ = 0, D = 0 |
| Identical sequences | S = 0, all metrics = 0 |

---

## Example

Using the Wikipedia Tajima's D example with 5 sequences of length 20 and 4 segregating sites:

```
Seq Y: 00000 00000 00000 00000  (reference)
Seq A: 00100 00000 00100 00010  (3 differences)
Seq B: 00000 00000 00100 00010  (2 differences)
Seq C: 00000 01000 00000 00010  (2 differences)
Seq D: 00000 01000 00100 00010  (3 differences)
```

**Calculations:**
- S = 4 (positions 3, 7, 13, 19)
- n = 5
- L = 20
- $a_1$ = 1 + 0.5 + 0.333 + 0.25 = 2.083
- Average pairwise differences = 2.0
- π = 2.0 / 20 = 0.1
- θ_W = 4 / (2.083 × 20) = 0.096
- d = π - θ_W = 0.004
- D (normalized) ≈ 0.08

---

## References

1. Nei, M.; Li, W.-H. (1979). "Mathematical Model for Studying Genetic Variation in Terms of Restriction Endonucleases". *PNAS*. 76(10): 5269-73.

2. Watterson, G.A. (1975). "On the number of segregating sites in genetical models without recombination". *Theoretical Population Biology*. 7(2): 256-276.

3. Tajima, F. (1989). "Statistical method for testing the neutral mutation hypothesis by DNA polymorphism". *Genetics*. 123(3): 585-95.

4. Hartl, D.L.; Clark, A.G. (2007). *Principles of Population Genetics* (4th ed.). Sinauer Associates.

5. Wikipedia contributors. "Nucleotide diversity", "Watterson estimator", "Tajima's D". *Wikipedia, The Free Encyclopedia*.
