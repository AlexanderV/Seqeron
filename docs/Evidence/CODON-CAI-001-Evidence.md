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

1. **Range:** 0 ≤ CAI ≤ 1
   - CAI = 1 when all codons are the most frequent for their amino acids
   - CAI = 0 when any codon has zero frequency in the reference set, or when input is empty
   - 0 < CAI < 1 for typical genes

2. **Geometric Mean:** CAI uses geometric mean, which is sensitive to low values
   - A single rare codon significantly lowers overall CAI

3. **Single-Codon Amino Acids:** Methionine (AUG) and Tryptophan (UGG) have w=1.0 always
   - Only one codon exists, so it's always the "most frequent"

4. **Stop Codons:** Excluded from CAI calculation
   - Source: Sharp & Li (1987) — stop codons do not encode amino acids

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

The implementation follows Sharp & Li (1987) with one deviation:
1. Converts T→U and handles case-insensitively
2. Splits sequence into codons
3. For each non-stop codon:
   - Finds the amino acid
   - Calculates relative adaptiveness: w = codon_freq / max_synonym_freq
   - If amino acid unknown or maxFreq = 0: returns NaN (skipped by caller)
   - If codon_freq = 0 but maxFreq > 0: clamps w to 1e-6 (incomplete table protection)
   - Accumulates ln(w) sum
4. Returns exp(sum / count)
5. **Deviation:** 1e-6 clamp for zero-frequency codons when amino acid has other codons in table (see CODON-CAI-001.md Deviations section for rationale)

## Test Datasets

### Reference Codon Tables (Kazusa MG1655, species=316407)

**E. coli K12 Leucine Codons:**
| Codon | Frequency | Relative Adaptiveness |
|-------|-----------|----------------------|
| CUG | 0.50 | 1.00 (optimal) |
| UUA | 0.13 | 0.26 |
| UUG | 0.13 | 0.26 |
| CUU | 0.10 | 0.20 |
| CUC | 0.10 | 0.20 |
| CUA | 0.04 | 0.08 (rare) |

**E. coli K12 Arginine Codons:**
| Codon | Frequency | Relative Adaptiveness |
|-------|-----------|----------------------|
| CGC | 0.40 | 1.00 (optimal) |
| CGU | 0.38 | 0.95 |
| CGG | 0.10 | 0.25 |
| CGA | 0.06 | 0.15 |
| AGA | 0.04 | 0.10 (rare) |
| AGG | 0.02 | 0.05 (rare) |

### Hand-Calculated Test Cases

**Test Case 1: Single Met (AUG)**
- w_AUG = 1.0 / 1.0 = 1.0
- CAI = 1.0^(1/1) = 1.0

**Test Case 2: CUG-CCG-ACC (E. coli)**
- CUG: w = 0.50/0.50 = 1.0 (Leu optimal)
- CCG: w = 0.53/0.53 = 1.0 (Pro optimal)
- ACC: w = 0.44/0.44 = 1.0 (Thr optimal)
- CAI = (1.0 × 1.0 × 1.0)^(1/3) = 1.0

**Test Case 3: CUA-CCA-ACA (E. coli rare)**
- CUA: w = 0.04/0.50 = 0.08 (Leu rare)
- CCA: w = 0.19/0.53 = 0.3585 (Pro suboptimal)
- ACA: w = 0.13/0.44 = 0.2955 (Thr suboptimal)
- CAI = (0.08 × 0.3585 × 0.2955)^(1/3) = 0.1980

## Assumptions

One documented deviation from strict Sharp & Li (1987):
- **1e-6 clamp:** When codon_freq = 0 but max_synonym_freq > 0, w is clamped to 1e-6 instead of 0. This protects against incomplete codon usage tables. See CODON-CAI-001.md Deviations section.
- Empty sequence returns 0 by convention (no codons to evaluate)
- All codon frequency tables verified against Kazusa database (March 2026)

## References

- Sharp, P.M. & Li, W.H. (1987). Nucleic Acids Res. 15(3):1281-1295
- Jansen, R. et al. (2003). Nucleic Acids Res. 31(8):2242-2251 (CAI refinements)
- Wikipedia: Codon Adaptation Index
- Kazusa Codon Usage Database: https://www.kazusa.or.jp/codon/
