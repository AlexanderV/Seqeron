# Rare Codon Detection

## Overview

Rare codon detection identifies codons in a coding sequence that occur at unusually low frequencies in a target organism's genome. These codons correspond to low-abundance tRNAs and can cause translation slowdown or stalling during heterologous protein expression.

## Biological Context

### Why Rare Codons Matter

1. **Translation Efficiency**: Rare codons are decoded by low-abundance tRNAs, causing ribosome pausing
2. **Protein Folding**: Strategic rare codon placement can facilitate co-translational protein folding
3. **Expression Optimization**: Replacing rare codons improves heterologous protein yields
4. **Translational Regulation**: Some genes use rare codons for regulatory purposes

### Common Rare Codons in E. coli K12

| Codon | Amino Acid | Frequency | tRNA Gene |
|-------|------------|-----------|-----------|
| AGA   | Arginine   | 0.07      | argU      |
| AGG   | Arginine   | 0.04      | argW      |
| CGA   | Arginine   | 0.07      | argW      |
| CUA   | Leucine    | 0.04      | leuV      |

Source: Kazusa Codon Usage Database; Shu et al. (2006)

## Algorithm

### Input
- `codingSequence`: DNA or RNA coding sequence (string)
- `table`: Codon usage table for target organism
- `threshold`: Frequency cutoff (default 0.15)

### Output
Enumerable of tuples: `(Position, Codon, AminoAcid, Frequency)`

### Process

```
1. Convert DNA to RNA (T → U)
2. Split sequence into codons (triplets)
3. For each codon at position i:
   a. Look up frequency in codon table
   b. If frequency < threshold:
      - Translate codon to amino acid
      - Yield (i * 3, codon, amino_acid, frequency)
```

### Complexity
- Time: O(n) where n = sequence length / 3
- Space: O(1) for streaming output

## Implementation Details

### Seqeron Implementation

```csharp
public static IEnumerable<(int Position, string Codon, string AminoAcid, double Frequency)> 
    FindRareCodons(string codingSequence, CodonUsageTable table, double threshold = 0.15)
```

**Behavior Notes**:
- Converts T to U automatically
- Incomplete codons (trailing 1-2 nucleotides) are ignored
- Unknown codons return frequency 0
- Position is nucleotide position (codon_index * 3)

### Threshold Selection

| Threshold | Description | Use Case |
|-----------|-------------|----------|
| 0.10      | Very rare   | Critical optimization |
| 0.15      | Default     | Standard analysis |
| 0.20      | Moderately rare | Broad screening |

## Usage Examples

### Basic Detection

```csharp
string sequence = "AUGAGAAGGCGA"; // M-R-R-R
var rareList = CodonOptimizer.FindRareCodons(sequence, CodonOptimizer.EColiK12, 0.10);
// Returns: (3, "AGA", "R", 0.07), (6, "AGG", "R", 0.04), (9, "CGA", "R", 0.07)
```

### Optimization Workflow

```csharp
// 1. Find rare codons
var rare = CodonOptimizer.FindRareCodons(sequence, CodonOptimizer.EColiK12, 0.10).ToList();

// 2. If rare codons found, optimize
if (rare.Count > 0)
{
    var result = CodonOptimizer.OptimizeSequence(
        sequence, 
        CodonOptimizer.EColiK12,
        CodonOptimizer.OptimizationStrategy.AvoidRareCodeons);
}
```

## Related Test Unit

- **Test Unit ID**: CODON-RARE-001
- **TestSpec**: [CODON-RARE-001.md](../../tests/TestSpecs/CODON-RARE-001.md)

## References

1. Sharp PM, Li WH. The codon adaptation index-a measure of directional synonymous codon usage bias. Nucleic Acids Res. 1987;15(3):1281-1295.
2. Shu P, Dai H, Gao W, Goldman E. Inhibition of Translation by Consecutive Rare Leucine Codons in E. coli. Gene Expr. 2006;13(2):97-106.
3. Plotkin JB, Kudla G. Synonymous but not the same: causes and consequences of codon bias. Nat Rev Genet. 2011;12(1):32-42.
4. Kane JF. Effects of rare codon clusters on high-level expression of heterologous proteins in E. coli. Curr Opin Biotechnol. 1995;6(5):494-500.

## See Also

- [Sequence Optimization](Sequence_Optimization.md) - CODON-OPT-001
- [CAI Calculation](CAI_Calculation.md) - CODON-CAI-001
