# CODON-CAI-001: Codon Adaptation Index (CAI) Calculation Evidence

## Test Unit
- **ID:** CODON-CAI-001
- **Area:** Codon
- **Algorithm:** CAI Calculation

## Sources

### Primary Sources

1. **Wikipedia - Codon Adaptation Index**
   - URL: https://en.wikipedia.org/wiki/Codon_Adaptation_Index
   - Key concepts:
     - CAI is "the most widespread technique for analyzing codon usage bias"
     - Measures "deviation of a given protein coding gene sequence with respect to a reference set of genes"
     - Used as "a quantitative method of predicting the level of expression of a gene based on its codon sequence"
     - Reference set ideally composed of highly expressed genes

2. **Sharp, P.M. & Li, W.H. (1987)** - Original CAI Paper
   - Title: "The codon adaptation index-a measure of directional synonymous codon usage bias, and its potential applications"
   - Journal: Nucleic Acids Research, 15(3): 1281-1295
   - PMC: 340524, PMID: 3547335
   - DOI: 10.1093/nar/15.3.1281

### Mathematical Definition (from Wikipedia)

**Relative Adaptiveness (w_i):**
```
w_i = f_i / max(f_j)  where i,j ∈ [synonymous codons for amino acid]
```

Where:
- f_i = observed frequency of codon i
- max(f_j) = frequency of the most frequent synonymous codon for that amino acid

**CAI Calculation:**
```
CAI = (∏_{i=1}^{L} w_i)^{1/L}
```

Equivalent to:
```
CAI = exp((1/L) × Σ ln(w_i))
```

Where L = number of codons (excluding stop codons per implementation)

### Key Properties (from Sharp & Li 1987 / Wikipedia)

1. **Range:** 0 < CAI ≤ 1
   - CAI = 1 when all codons are the most frequent for their amino acids
   - CAI < 1 when non-optimal codons are used

2. **Geometric Mean:** CAI uses geometric mean, which is sensitive to low values
   - A single rare codon significantly lowers overall CAI

3. **Single-Codon Amino Acids:** Methionine (AUG) and Tryptophan (UGG) have w=1.0 always
   - Only one codon exists, so it's always the "most frequent"

4. **Stop Codons:** Typically excluded from CAI calculation
   - Source: Implementation standard practice

### Documented Edge Cases

1. **Empty Sequence:**
   - No codons to evaluate → CAI = 0 (by convention)
   - Source: Implementation convention

2. **Sequence with only Met/Trp:**
   - All codons have w=1.0 → CAI = 1.0
   - Source: Mathematical definition

3. **All Optimal Codons:**
   - Every codon is the most frequent for its amino acid → CAI = 1.0
   - Source: CAI definition

4. **All Rare Codons:**
   - Low w values for all codons → CAI approaches 0
   - Source: Geometric mean properties

5. **Different Organisms:**
   - Same sequence has different CAI values for different organisms
   - Source: Codon usage varies by organism

### Implementation Notes (from source code analysis)

The implementation:
1. Converts T→U and handles case-insensitively
2. Splits sequence into codons
3. For each non-stop codon:
   - Finds the amino acid
   - Calculates relative adaptiveness: codon_freq / max_synonym_freq
   - Accumulates log(w) sum
4. Returns exp(sum / count)

## Test Datasets

### Reference Codon Tables (from implementation)

**E. coli K12 Leucine Codons:**
| Codon | Frequency | Relative Adaptiveness |
|-------|-----------|----------------------|
| CUG | 0.47 | 1.00 (optimal) |
| UUA | 0.14 | 0.30 |
| UUG | 0.13 | 0.28 |
| CUU | 0.12 | 0.26 |
| CUC | 0.10 | 0.21 |
| CUA | 0.04 | 0.09 (rare) |

**E. coli K12 Arginine Codons:**
| Codon | Frequency | Relative Adaptiveness |
|-------|-----------|----------------------|
| CGU | 0.36 | 1.00 (tied optimal) |
| CGC | 0.36 | 1.00 (tied optimal) |
| CGG | 0.11 | 0.31 |
| CGA | 0.07 | 0.19 |
| AGA | 0.07 | 0.19 (rare) |
| AGG | 0.04 | 0.11 (rare) |

### Hand-Calculated Test Cases

**Test Case 1: Single Met (AUG)**
- w_AUG = 1.0 / 1.0 = 1.0
- CAI = 1.0^(1/1) = 1.0

**Test Case 2: CUG-CCG-ACC (E. coli)**
- CUG: w = 0.47/0.47 = 1.0 (Leu optimal)
- CCG: w = 0.49/0.49 = 1.0 (Pro optimal)
- ACC: w = 0.40/0.40 = 1.0 (Thr optimal)
- CAI = (1.0 × 1.0 × 1.0)^(1/3) = 1.0

**Test Case 3: CUA-CCA-ACA (E. coli rare)**
- CUA: w = 0.04/0.47 ≈ 0.085 (Leu rare)
- CCA: w = 0.20/0.49 ≈ 0.408 (Pro suboptimal)
- ACA: w = 0.17/0.40 ≈ 0.425 (Thr suboptimal)
- CAI = (0.085 × 0.408 × 0.425)^(1/3) ≈ 0.24

## Assumptions

1. **ASSUMPTION:** Empty sequence returns CAI = 0 (not undefined)
   - Rationale: Practical convention to avoid special handling

2. **ASSUMPTION:** Minimum w value clamped to small positive value (0.01) to avoid log(0)
   - Rationale: Implementation detail to handle codons not in reference table

## References

- Sharp, P.M. & Li, W.H. (1987). Nucleic Acids Res. 15(3):1281-1295
- Jansen, R. et al. (2003). Nucleic Acids Res. 31(8):2242-2251 (CAI refinements)
- Wikipedia: Codon Adaptation Index
- Kazusa Codon Usage Database: https://www.kazusa.or.jp/codon/
