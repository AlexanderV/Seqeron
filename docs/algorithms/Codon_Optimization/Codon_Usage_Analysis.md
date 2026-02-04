# Codon Usage Analysis

## Algorithm Overview

Codon usage analysis measures the frequency of synonymous codons in coding sequences to understand codon bias patterns and compare sequences across organisms.

## Theory

### Codon Redundancy

The genetic code contains 64 codons that encode 20 amino acids plus 3 stop signals. Most amino acids are encoded by 2-6 synonymous codons. Organisms exhibit preferences for certain codons, reflecting:

- tRNA abundance patterns
- GC content of the genome
- Selection pressures on translation efficiency
- Mutational biases

### Codon Usage Calculation

Given a coding sequence, codon usage is calculated by:

1. Splitting the sequence into non-overlapping triplets (codons)
2. Counting occurrences of each codon
3. Optionally normalizing to frequencies

**Formula** (raw counts):
$$\text{Count}(c) = |\{i : \text{seq}[3i:3i+3] = c, \forall i \in [0, \lfloor n/3 \rfloor)\}|$$

**Formula** (frequency):
$$f(c) = \frac{\text{Count}(c)}{\sum_{c'} \text{Count}(c')}$$

### Codon Usage Comparison

Comparing codon usage between two sequences quantifies how similar their codon preferences are. Common metrics include:

1. **Cosine similarity** - angle between frequency vectors
2. **Pearson correlation** - linear correlation of frequencies
3. **Manhattan distance** - sum of absolute differences

The implementation uses a normalized Manhattan distance-based similarity:

$$\text{Similarity} = 1 - \frac{\sum_{c} |f_1(c) - f_2(c)|}{2}$$

Where the division by 2 normalizes the range to [0, 1], since the maximum possible sum of absolute differences between two probability distributions is 2.

### Properties

| Property | Value |
|----------|-------|
| Range (similarity) | [0, 1] |
| Identity | Sim(s, s) = 1.0 for non-empty s |
| Symmetry | Sim(a, b) = Sim(b, a) |
| Empty sequences | Returns 0 (no data to compare) |

## Complexity

| Operation | Time | Space |
|-----------|------|-------|
| CalculateCodonUsage | O(n) | O(64) = O(1) |
| CompareCodonUsage | O(n + m) | O(64) = O(1) |

Where n and m are sequence lengths.

## Implementation Notes

### Current Implementation (Seqeron.Genomics.MolTools)

The `CodonOptimizer` class provides:

```csharp
public static Dictionary<string, int> CalculateCodonUsage(string codingSequence)
public static double CompareCodonUsage(string sequence1, string sequence2)
```

**Key behaviors**:
- Converts DNA (T) to RNA (U) internally
- Case-insensitive processing
- Incomplete trailing codons are ignored
- Empty sequences return empty dictionary / 0.0 similarity

### Edge Cases Handled

1. **Empty sequence**: Returns empty dictionary or 0.0 similarity
2. **Incomplete codons**: Trailing 1-2 nucleotides ignored
3. **Mixed T/U**: Converted to standard RNA (U)
4. **Case variation**: Normalized to uppercase

## Applications

1. **Codon optimization**: Selecting preferred codons for expression host
2. **Evolutionary analysis**: Comparing codon usage across species
3. **Gene expression prediction**: Highly expressed genes use preferred codons
4. **Horizontal gene transfer detection**: Atypical codon usage may indicate foreign genes

## References

1. Sharp PM, Li WH (1987). "The codon adaptation index-a measure of directional synonymous codon usage bias, and its potential applications." *Nucleic Acids Research* 15(3):1281-1295.
2. Plotkin JB, Kudla G (2011). "Synonymous but not the same: The causes and consequences of codon bias." *Nature Reviews Genetics* 12(1):32-42.
3. Wikipedia: Codon usage bias - https://en.wikipedia.org/wiki/Codon_usage_bias
4. Kazusa Codon Usage Database - https://www.kazusa.or.jp/codon/

## Related Algorithms

- [CAI_Calculation.md](CAI_Calculation.md) - Codon Adaptation Index uses codon frequencies
- [Rare_Codon_Detection.md](Rare_Codon_Detection.md) - Finding rare codons based on usage
- [Sequence_Optimization.md](Sequence_Optimization.md) - Optimizing codon usage for expression
