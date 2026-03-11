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

The implementation uses Total Variation Distance (TVD) similarity between codon frequency distributions:

1. **Input**: Two coding sequences
2. **Process**:
   - Calculate codon usage for both sequences
   - Normalize counts to frequencies (count / total codons)
   - Calculate TVD = (1/2) × L¹ distance between frequency vectors
   - Similarity = 1 - TVD
3. **Output**: Similarity value in [0, 1]

**Mathematical definition**:
$$\text{Similarity} = 1 - \frac{\sum_{c \in \text{AllCodons}} |f_1(c) - f_2(c)|}{2}$$

Where $f_i(c)$ is the frequency of codon $c$ in sequence $i$.

**Proven properties** (from TVD theory):
- **Identity**: sim(s,s) = 1.0 (zero distance for identical distributions)
- **Symmetry**: sim(a,b) = sim(b,a) (|x-y| = |y-x|)
- **Range**: [0,1] (TVD of probability distributions ∈ [0,1])
- **Disjoint → 0**: For non-overlapping codon sets, Σ|f₁-f₂| = 1+1 = 2, so sim = 0
- **Partial overlap**: Analytically derivable for any input; e.g. 2/3 shared codons → sim = 2/3

### Edge Cases (derived from TVD formula and standard practice)

| Case | Expected Behavior | Source |
|------|-------------------|--------|
| Empty sequence | Empty dictionary / 0 similarity | Convention: no data → 0 |
| Incomplete final codon | Ignore trailing nucleotides | Kazusa, EMBOSS |
| Identical sequences | Similarity = 1.0 | TVD = 0 for identical distributions |
| Disjoint codons | Similarity = 0.0 | TVD = 1 for orthogonal distributions |
| Partial overlap | Exact value derivable | TVD formula computation |
| T/U conversion | Convert DNA T to RNA U internally | Biological equivalence |

## Test Data

### Reference Codon Tables
- **E. coli K12**: Well-characterized codon usage preferences
- **S. cerevisiae (Yeast)**: Different preferences than E. coli
- **H. sapiens (Human)**: GC-rich codon bias

### Manually Verified Test Cases

1. **Simple counting**:
   - Input: "AUGGCUGCU" (M-A-A)
   - Expected: {"AUG": 1, "GCU": 2}

2. **All 64 codons**:
   - Input: concatenation of all 64 standard RNA codons
   - Expected: 64 distinct keys, each with count 1

3. **Identical sequence comparison**:
   - Input: seq1 = seq2 = "AUGGCUGCACUG"
   - Expected: Similarity = 1.0

4. **Partial overlap sim=0.5** (exact TVD derivation):
   - seq1 = "CUGCUGCUGCUA" → f(CUG)=3/4, f(CUA)=1/4
   - seq2 = "CUACUACUACUG" → f(CUA)=3/4, f(CUG)=1/4
   - Σ|f₁-f₂| = 1/2 + 1/2 = 1 → Similarity = 0.5

5. **Symmetry + exact sim=0.75** (TVD derivation):
   - seq1 = "AUGAUGCCCUUU" → f(AUG)=1/2, f(CCC)=1/4, f(UUU)=1/4
   - seq2 = "AUGUUUUUUCCC" → f(AUG)=1/4, f(UUU)=1/2, f(CCC)=1/4
   - Σ = 1/4+0+1/4 = 1/2 → sim(a,b) = sim(b,a) = 0.75

6. **High-difference sim=0.25** (TVD derivation):
   - seq1 = "AUGAUGAUGAUG" → f(AUG)=1
   - seq2 = "AUGCCCCCCCCC" → f(AUG)=1/4, f(CCC)=3/4
   - Σ = 3/4+3/4 = 3/2 → Similarity = 0.25

7. **Low-difference sim=0.75** (TVD derivation):
   - seq1 = "AUGAUGAUGCCC" → f(AUG)=3/4, f(CCC)=1/4
   - seq2 = "AUGAUGCCCCCC" → f(AUG)=1/2, f(CCC)=1/2
   - Σ = 1/4+1/4 = 1/2 → Similarity = 0.75

8. **Disjoint codons** (TVD derivation):
   - seq1 = "UUUUUUUUU" (all UUU)
   - seq2 = "GGGGGGGGG" (all GGG)
   - Σ|f₁-f₂| = 1 + 1 = 2 → Similarity = 0.0

9. **2/3 shared codons** (TVD derivation):
   - seq1 = "AUGGCUAUG" → f(AUG)=2/3, f(GCU)=1/3
   - seq2 = "AUGUUUAUG" → f(AUG)=2/3, f(UUU)=1/3
   - Σ|f₁-f₂| = 0 + 1/3 + 1/3 = 2/3 → Similarity = 2/3

10. **Empty sequences**:
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

## Kazusa Verification

All predefined codon usage tables verified against Kazusa Codon Usage Database (March 2026):

| Organism | Species ID | CDS Count | Codons | Status |
|----------|-----------|-----------|--------|--------|
| E. coli K12 (W3110) | 316407 | 4,332 | 1,372,057 | ✅ All 64 relative fractions match |
| S. cerevisiae | 4932 | 14,411 | 6,534,504 | ✅ All 64 relative fractions match |
| H. sapiens | 9606 | 93,487 | 40,662,582 | ✅ All 64 relative fractions match |

Verification method: per-thousand frequencies from Kazusa converted to relative fractions per amino acid and compared with implementation values (2 decimal places).

## Implementation Notes

### Current Implementation Observations

1. `CalculateCodonUsage`:
   - Converts T→U internally (works with RNA representation)
   - Uses uppercase normalization
   - Returns counts, not frequencies

2. `CompareCodonUsage`:
   - Uses TVD-based similarity: 1 - Σ|f₁-f₂|/2
   - Returns 0 for empty sequences (not NaN or exception)
   - All expected values for test cases analytically derivable

### Proven Test Properties

1. **Invariant**: Sum of all codon counts = total codons in sequence
2. **Range**: Similarity always in [0, 1] (TVD ∈ [0,1] for probability distributions)
3. **Symmetry**: CompareCodonUsage(a, b) = CompareCodonUsage(b, a) (|x-y| = |y-x|)
4. **Identity**: CompareCodonUsage(a, a) = 1.0 (zero distance for identical distributions)
5. **Disjoint → 0**: Proven via Σ|f₁-f₂| = 2 for orthogonal distributions

## Related Test Units

- **CODON-OPT-001**: Uses codon usage for optimization
- **CODON-CAI-001**: Uses codon frequencies for CAI calculation
- **CODON-RARE-001**: Identifies rare codons based on usage tables

## Last Updated
2026-03-11
