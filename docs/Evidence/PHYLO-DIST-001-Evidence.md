# PHYLO-DIST-001: Phylogenetic Distance Matrix - Evidence Document

## Test Unit
**ID:** PHYLO-DIST-001
**Title:** Phylogenetic Distance Matrix Calculation
**Status:** Evidence Collected
**Last Updated:** 2026-02-01

---

## 1. Sources Consulted

### 1.1 Wikipedia - Models of DNA Evolution
**URL:** https://en.wikipedia.org/wiki/Models_of_DNA_evolution
**Accessed:** 2026-02-01

**Key Information Extracted:**
- **p-distance (Hamming proportion):** Proportion of differing sites between two sequences; d = (differences / comparable sites)
- **Jukes-Cantor 1969 (JC69):** d = -3/4 × ln(1 - 4p/3) where p is the p-distance
  - Assumes equal base frequencies (π_A = π_G = π_C = π_T = 0.25)
  - Assumes equal mutation rates between all nucleotides
  - Formula corrects for multiple substitutions at the same site
- **Kimura 2-Parameter (K80/K2P):** d = -0.5 × ln((1 - 2S - V) × √(1 - 2V))
  - S = proportion of transition differences
  - V = proportion of transversion differences
  - Distinguishes transitions (purine↔purine: A↔G, pyrimidine↔pyrimidine: C↔T) from transversions
- Gaps should be ignored in distance calculations

### 1.2 Wikipedia - Substitution Model
**URL:** https://en.wikipedia.org/wiki/Substitution_model
**Accessed:** 2026-02-01

**Key Information Extracted:**
- Distance matrix is symmetric: d(i,j) = d(j,i)
- Diagonal elements are zero: d(i,i) = 0
- Time-reversible models allow evolutionary distance calculation between any two sequences
- JC69 reference: Jukes TH, Cantor CR (1969) "Evolution of Protein Molecules"
- K80 reference: Kimura M (1980) "A simple method for estimating evolutionary rates"

### 1.3 Wikipedia - Distance Matrices in Phylogeny
**URL:** https://en.wikipedia.org/wiki/Distance_matrices_in_phylogeny
**Accessed:** 2026-02-01

**Key Information Extracted:**
- Distance is often defined as fraction of mismatches at aligned positions
- Gaps either ignored or counted as mismatches (implementation choice)
- Raw distance values (Hamming distance) can be calculated by counting pairwise differences
- Distance correction (Jukes-Cantor) accounts for back mutations

---

## 2. Published Test Cases

### 2.1 Identity Test Case
**Source:** Mathematical definition
**Test:** Two identical sequences should have distance = 0 for all methods

### 2.2 Complete Difference Test Case
**Source:** Wikipedia JC69 formula behavior
**Test:** When p ≥ 0.75 (75% differences), JC69 returns +∞ (saturation)
- Formula: arg = 1 - 4p/3 → when p ≥ 0.75, arg ≤ 0, ln undefined

### 2.3 Known Distance Values
**Source:** JC69 formula
**Examples:**
- p = 0 → d = 0
- p = 0.125 (1/8) → d = -0.75 × ln(1 - 4×0.125/3) = -0.75 × ln(0.833) ≈ 0.137
- p = 0.25 → d = -0.75 × ln(1 - 1/3) = -0.75 × ln(0.667) ≈ 0.304

### 2.4 Symmetry Property
**Source:** Time-reversibility (Wikipedia Substitution Model)
**Test:** d(seq1, seq2) = d(seq2, seq1)

### 2.5 Triangle Inequality (expected to hold in most cases)
**Source:** Metric space properties
**Test:** d(A,C) ≤ d(A,B) + d(B,C) (may not always hold due to correction)

---

## 3. Corner Cases from Sources

| Corner Case | Expected Behavior | Source |
|-------------|-------------------|--------|
| Identical sequences | Distance = 0 | Definition |
| Single different base | Small positive distance | Formula calculation |
| All gaps in alignment | Distance = 0 (no comparable sites) | Implementation choice |
| Unequal length sequences | Throw ArgumentException | Pre-condition |
| Empty sequences | Distance = 0 or throw | Implementation choice |
| Case-insensitive | Upper/lowercase treated same | Standard practice |
| High divergence (p ≥ 0.75) | JC69 returns +∞ | JC69 formula saturation |

---

## 4. Algorithm Invariants

### 4.1 Distance Matrix Properties
1. **Symmetry:** Matrix[i,j] = Matrix[j,i] for all i,j
2. **Zero diagonal:** Matrix[i,i] = 0 for all i
3. **Non-negative:** Matrix[i,j] ≥ 0 for all i,j (for corrected distances)
4. **Dimensions:** n×n matrix for n sequences

### 4.2 Distance Method Relationships
1. **JC69 ≥ p-distance:** Corrected distance always larger than or equal to raw proportion
2. **K2P ≥ p-distance:** Kimura correction also increases distance
3. **Hamming = differences count:** Raw integer count, not a proportion

### 4.3 Gap Handling
- Gaps ('-') are excluded from comparison (not counted in comparable sites)
- Position with gap in either sequence is skipped

---

## 5. Test Datasets

### 5.1 Minimal Test Set
```
Seq1: ACGTACGT (8 bp)
Seq2: ACGTACGT (identical)
Expected: d = 0 for all methods
```

### 5.2 Single Difference Test
```
Seq1: ACGTACGT (8 bp)
Seq2: TCGTACGT (1 difference at pos 0)
Expected:
- Hamming: 1
- p-distance: 0.125
- JC69: ≈ 0.137
```

### 5.3 Transition vs Transversion Test
```
Seq1: AAAA
Seq2: AGGG (1 transition A→G at pos 1, 2 transversions A→G at pos 2,3)
Note: A→G is a transition (purine to purine)
```

### 5.4 With Gaps Test
```
Seq1: ACG-ACGT (8 chars, 1 gap)
Seq2: ACGTACGT (8 chars)
Expected: Gap position ignored, 7 comparable sites
```

---

## 6. Implementation Notes

### 6.1 Current Implementation Observations
- `CalculateDistanceMatrix(seqs, method)`: Computes pairwise distances for all sequences
- `CalculatePairwiseDistance(s1, s2, method)`: Single pair calculation
- Supports 4 methods: Hamming, PDistance, JukesCantor, Kimura2Parameter
- Gaps are skipped (not counted in comparable sites)
- Case-insensitive comparison (ToUpperInvariant)
- Throws on unequal length sequences

### 6.2 JC69 Formula Implementation
```
d = -0.75 × ln(1 - 4p/3)
where p = differences / comparableSites
```

### 6.3 K2P Formula Implementation
```
d = -0.5 × ln((1 - 2S - V) × √(1 - 2V))
where S = transitions / comparableSites
      V = transversions / comparableSites
```

---

## 7. References

1. Jukes TH, Cantor CR (1969). "Evolution of Protein Molecules". In Munro HN (ed.). Mammalian Protein Metabolism. Academic Press. pp. 21–132.
2. Kimura M (1980). "A simple method for estimating evolutionary rates of base substitutions through comparative studies of nucleotide sequences". Journal of Molecular Evolution. 16: 111–120.
3. Felsenstein J (2004). Inferring Phylogenies. Sinauer Associates.
4. Wikipedia contributors. "Models of DNA evolution." Wikipedia, The Free Encyclopedia.
5. Wikipedia contributors. "Distance matrices in phylogeny." Wikipedia, The Free Encyclopedia.
