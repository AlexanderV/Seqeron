# Phylogenetic Distance Matrix Calculation

## Overview

Phylogenetic distance matrices quantify the evolutionary divergence between pairs of aligned sequences. These distances are used as input for distance-based phylogenetic tree construction methods (UPGMA, Neighbor-Joining) and for comparative evolutionary analysis.

## Methods

### 1. Hamming Distance (Raw Count)

The simplest measure: count of positions where sequences differ.

**Formula:**
```
d_H(s1, s2) = Σ(s1[i] ≠ s2[i]) for all comparable positions i
```

**Properties:**
- Returns integer count
- Does not account for multiple substitutions
- Not corrected for evolutionary saturation

### 2. p-distance (Uncorrected Distance)

Proportion of differing sites among comparable sites.

**Formula:**
```
p = (number of differences) / (number of comparable sites)
```

**Properties:**
- Range: [0, 1]
- Underestimates true evolutionary distance for divergent sequences
- Simple and widely used for closely related sequences

### 3. Jukes-Cantor Distance (JC69)

Corrects for multiple substitutions assuming equal base frequencies and equal substitution rates.

**Formula:**
```
d_JC = -3/4 × ln(1 - 4p/3)
```

Where p is the p-distance.

**Properties:**
- Always ≥ p-distance
- Becomes undefined (→ ∞) when p ≥ 0.75 (saturation)
- Assumes: π_A = π_G = π_C = π_T = 0.25
- Reference: Jukes & Cantor (1969)

### 4. Kimura 2-Parameter Distance (K2P/K80)

Distinguishes between transitions (purine↔purine, pyrimidine↔pyrimidine) and transversions (purine↔pyrimidine).

**Formula:**
```
d_K2P = -0.5 × ln((1 - 2S - V) × √(1 - 2V))
```

Where:
- S = proportion of transitions
- V = proportion of transversions

**Transitions:** A↔G (purines), C↔T (pyrimidines)
**Transversions:** A↔C, A↔T, G↔C, G↔T

**Properties:**
- Accounts for transition/transversion bias
- More accurate for sequences with unequal mutation rates
- Reference: Kimura (1980)

## Distance Matrix Properties

A valid phylogenetic distance matrix must satisfy:

| Property | Definition |
|----------|------------|
| **Symmetry** | d(i,j) = d(j,i) |
| **Zero Diagonal** | d(i,i) = 0 |
| **Non-Negativity** | d(i,j) ≥ 0 |

## Gap Handling

Standard practice for aligned sequences with gaps ('-'):

- **Pairwise deletion:** Exclude positions with gaps in either sequence from comparison
- Only positions where both sequences have a nucleotide (A, C, G, T) are compared
- This maximizes usable data for each pairwise comparison

## Implementation (Seqeron.Genomics)

### CalculateDistanceMatrix

```csharp
double[,] CalculateDistanceMatrix(
    IReadOnlyList<string> alignedSequences,
    DistanceMethod method = DistanceMethod.JukesCantor)
```

**Parameters:**
- `alignedSequences`: List of equal-length aligned sequences
- `method`: Distance method (Hamming, PDistance, JukesCantor, Kimura2Parameter)

**Returns:** Symmetric n×n matrix of pairwise distances

**Complexity:** O(n² × m) where n = sequence count, m = sequence length

### CalculatePairwiseDistance

```csharp
double CalculatePairwiseDistance(
    string seq1,
    string seq2,
    DistanceMethod method = DistanceMethod.JukesCantor)
```

**Parameters:**
- `seq1`, `seq2`: Aligned sequences (must be equal length)
- `method`: Distance method

**Returns:** Distance value

**Throws:** `ArgumentException` if sequences have different lengths

## Edge Cases

| Scenario | Behavior |
|----------|----------|
| Identical sequences | Returns 0 |
| All gaps | Returns 0 (no comparable sites) |
| High divergence (p ≥ 0.75) | JC69 returns +∞ |
| Unequal lengths | Throws ArgumentException |
| Case difference | Case-insensitive comparison |

## References

1. Jukes TH, Cantor CR (1969). "Evolution of Protein Molecules". *Mammalian Protein Metabolism*. Academic Press. pp. 21–132.

2. Kimura M (1980). "A simple method for estimating evolutionary rates of base substitutions through comparative studies of nucleotide sequences". *Journal of Molecular Evolution*. 16: 111–120.

3. Felsenstein J (2004). *Inferring Phylogenies*. Sinauer Associates.

4. Wikipedia. "Models of DNA evolution". https://en.wikipedia.org/wiki/Models_of_DNA_evolution

5. Wikipedia. "Distance matrices in phylogeny". https://en.wikipedia.org/wiki/Distance_matrices_in_phylogeny
