# Evidence: CODON-USAGE-001 - Codon Usage Analysis

## Test Unit ID
CODON-USAGE-001

## Methods Under Test
| Method | Class | Type |
|--------|-------|------|
| `CalculateCodonUsage(string)` | CodonOptimizer | Canonical |
| `CompareCodonUsage(string, string)` | CodonOptimizer | Canonical |

## Sources

### Primary Sources

1. **Wikipedia - Codon usage bias**
   - URL: https://en.wikipedia.org/wiki/Codon_usage_bias
   - Key information:
     - Codon usage bias refers to differences in the frequency of occurrence of synonymous codons in coding DNA
     - 64 codons encode 20 amino acids + 3 stop codons, leading to redundancy
     - Codon usage tables detail genomic codon usage bias for organisms
     - RSCU (Relative Synonymous Codon Usage) measures codon preference
     - Comparison metrics include cosine similarity, correlation coefficients
   - Relevant quote: "The overabundance in the number of codons allows many amino acids to be encoded by more than one codon. Because of such redundancy it is said that the genetic code is degenerate."

2. **Kazusa Codon Usage Database**
   - URL: https://www.kazusa.or.jp/codon/
   - Key information:
     - Standard format for codon usage tables
     - 35,799 organisms with codon frequency data
     - Codons are counted per 1000 bases (frequency per thousand)
     - Tables contain all 64 codons with observed counts

3. **Sharp & Li (1987) - CAI Paper**
   - Citation: Sharp PM, Li WH. "The codon adaptation index-a measure of directional synonymous codon usage bias, and its potential applications." Nucleic Acids Research. 1987;15(3):1281-1295
   - Key information:
     - Codon usage frequencies are typically normalized by amino acid
     - Relative adaptiveness (w) = frequency of codon / frequency of most common synonymous codon

### Supporting Sources

4. **Plotkin & Kudla (2011)** - Nature Reviews Genetics
   - "Synonymous but not the same: The causes and consequences of codon bias"
   - Discusses codon usage optimization and comparison methods

5. **Athey et al. (2017)** - BMC Bioinformatics
   - "A new and updated resource for codon usage tables"
   - Describes HIVE-CUTs database with comprehensive codon usage data

## Algorithm Theory

### Codon Usage Calculation

1. **Input**: Coding sequence (RNA or DNA)
2. **Process**:
   - Split sequence into codons (triplets)
   - Count occurrences of each codon
   - Return dictionary mapping codon → count
3. **Output**: Dictionary<string, int> of codon counts

**Mathematical definition**:
$$\text{Count}(c) = |\{i : \text{sequence}[i:i+3] = c\}|$$

### Codon Usage Comparison

The implementation uses a frequency-based similarity metric based on absolute differences:

1. **Input**: Two coding sequences
2. **Process**:
   - Calculate codon usage for both sequences
   - Normalize counts to frequencies (count / total codons)
   - Calculate Manhattan distance between frequency vectors
   - Convert to similarity: 1 - (distance / 2)
3. **Output**: Similarity value in [0, 1]

**Mathematical definition**:
$$\text{Similarity} = 1 - \frac{\sum_{c \in \text{AllCodons}} |f_1(c) - f_2(c)|}{2}$$

Where $f_i(c)$ is the frequency of codon $c$ in sequence $i$.

### Edge Cases (from sources)

| Case | Expected Behavior | Source |
|------|-------------------|--------|
| Empty sequence | Empty dictionary / 0 similarity | Standard practice |
| Incomplete final codon | Ignore trailing nucleotides | Kazusa, EMBOSS |
| Identical sequences | Similarity = 1.0 | Mathematical definition |
| No overlapping codons | Similarity = 0.0 | Mathematical definition |
| T/U conversion | Convert DNA T to RNA U internally | Implementation note |

## Test Data

### Reference Codon Tables
- **E. coli K12**: Well-characterized codon usage preferences
- **S. cerevisiae (Yeast)**: Different preferences than E. coli
- **H. sapiens (Human)**: GC-rich codon bias

### Manually Verified Test Cases

1. **Simple counting**:
   - Input: "AUGGCUGCU" (M-A-A)
   - Expected: {"AUG": 1, "GCU": 2}

2. **Identical sequence comparison**:
   - Input: seq1 = seq2 = "AUGGCUGCACUG"
   - Expected: Similarity = 1.0

3. **Completely different codons**:
   - Input: seq1 = "CUGCUGCUGCUG" (all CUG)
   - Input: seq2 = "CUACUACUACUA" (all CUA)
   - Expected: Similarity < 1.0 (no shared codons)

4. **Empty sequences**:
   - Input: "", ""
   - Expected: Similarity = 0 (no data)

## Corner Cases

| Corner Case | Expected Result | Rationale |
|-------------|-----------------|-----------|
| Empty string | {} / 0.0 | No codons to process |
| Null input | ArgumentNullException | Defensive programming |
| 1-2 nucleotides | {} | Not a complete codon |
| Non-standard characters | Ignore or handle gracefully | Robustness |
| Mixed case | Case-insensitive | Standard practice |
| DNA (T) vs RNA (U) | Treat equivalently | Biological equivalence |

## Implementation Notes

### Current Implementation Observations

1. `CalculateCodonUsage`:
   - Converts T→U internally (works with RNA representation)
   - Uses uppercase normalization
   - Returns counts, not frequencies

2. `CompareCodonUsage`:
   - Uses frequency-normalized comparison
   - Returns 0 for empty sequences (not NaN or exception)
   - Calculates Manhattan distance-based similarity

### Potential Testing Focus

1. **Invariant**: Sum of all codon counts = total codons in sequence
2. **Range**: Similarity always in [0, 1]
3. **Symmetry**: CompareCodonUsage(a, b) = CompareCodonUsage(b, a)
4. **Identity**: CompareCodonUsage(a, a) = 1.0 (for non-empty)

## Related Test Units

- **CODON-OPT-001**: Uses codon usage for optimization
- **CODON-CAI-001**: Uses codon frequencies for CAI calculation
- **CODON-RARE-001**: Identifies rare codons based on usage tables

## Last Updated
2026-02-04
