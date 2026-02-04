# Codon Adaptation Index (CAI) Calculation

## Overview

The Codon Adaptation Index (CAI) is a quantitative measure of codon usage bias that predicts gene expression level based on codon composition. It was introduced by Sharp and Li in 1987 and remains the most widely used metric for analyzing synonymous codon usage bias.

## Algorithm

### Mathematical Definition

CAI is calculated as the geometric mean of relative adaptiveness values for all codons in a sequence.

**Step 1: Calculate Relative Adaptiveness (w_i)**

For each codon `i` that encodes amino acid `a`:

```
w_i = f_i / max(f_j)  for all j encoding amino acid a
```

Where:
- `f_i` = frequency of codon i in the reference set
- `max(f_j)` = maximum frequency among all synonymous codons for amino acid a

**Step 2: Calculate CAI**

```
CAI = (∏_{i=1}^{L} w_i)^{1/L}
```

Equivalently using logarithms:
```
CAI = exp((1/L) × Σ_{i=1}^{L} ln(w_i))
```

Where `L` = number of non-stop codons in the sequence.

### Complexity

- **Time:** O(n) where n = sequence length
- **Space:** O(1) for computation (reference tables are constant)

## Properties

### Range
- **Minimum:** 0 (theoretically; practically > 0 due to minimum frequency thresholds)
- **Maximum:** 1 (when all codons are optimal for their amino acids)

### Invariants
1. **Single-codon amino acids always contribute w=1.0**
   - Methionine (AUG) and Tryptophan (UGG) have no synonymous alternatives
   
2. **Geometric mean sensitivity**
   - A single rare codon significantly lowers overall CAI
   - Product of probabilities amplifies the effect of low values

3. **Organism specificity**
   - Same sequence has different CAI values for different organisms
   - Reflects organism-specific tRNA pools and codon preferences

## Implementation

### Current Implementation

```csharp
// From CodonOptimizer.CalculateCAI
public static double CalculateCAI(string codingSequence, CodonUsageTable table)
```

**Key behaviors:**
- Returns 0 for empty sequences (by convention)
- Converts T→U and normalizes to uppercase
- Excludes stop codons from calculation
- Uses log-sum-exp for numerical stability

### Edge Case Handling

| Case | Behavior |
|------|----------|
| Empty sequence | Returns 0 |
| Null sequence | Returns 0 (safe handling) |
| DNA input (with T) | Converts T→U automatically |
| Lowercase input | Converts to uppercase |
| Unknown codon | Uses minimum frequency (0.01) |
| Stop codon | Excluded from calculation |

## Reference Data

### Supported Organism Tables

1. **E. coli K12** (`CodonOptimizer.EColiK12`)
   - Highly biased toward specific codons
   - Example: CUG dominant for Leucine (0.47)

2. **S. cerevisiae** (`CodonOptimizer.Yeast`)
   - Different bias pattern
   - Example: UUA/UUG preferred for Leucine

3. **H. sapiens** (`CodonOptimizer.Human`)
   - Less extreme bias than bacteria
   - Example: CUG still preferred for Leucine (0.41)

### Example Calculations

**Optimal sequence for E. coli:**
```
Sequence: AUGCUGCCGACC (Met-Leu-Pro-Thr)
Codons: AUG(1.0), CUG(1.0), CCG(1.0), ACC(1.0)
CAI = (1.0 × 1.0 × 1.0 × 1.0)^(1/4) = 1.0
```

**Suboptimal sequence for E. coli:**
```
Sequence: AUGCUACCAACU (Met-Leu-Pro-Thr)
Codons: AUG(1.0), CUA(0.085), CCA(0.408), ACU(0.475)
CAI = (1.0 × 0.085 × 0.408 × 0.475)^(1/4) ≈ 0.37
```

## Applications

1. **Gene expression prediction:** Higher CAI correlates with higher expression
2. **Codon optimization:** Target CAI improvement during sequence design
3. **Horizontal gene transfer detection:** Anomalous CAI may indicate foreign genes
4. **Evolutionary analysis:** Compare codon adaptation across species

## Sources

- Sharp, P.M. & Li, W.H. (1987). "The codon adaptation index-a measure of directional synonymous codon usage bias, and its potential applications." Nucleic Acids Research, 15(3):1281-1295.
- Wikipedia: Codon Adaptation Index
- Plotkin, J.B. & Kudla, G. (2011). "Synonymous but not the same." Nature Reviews Genetics, 12(1):32-42.
