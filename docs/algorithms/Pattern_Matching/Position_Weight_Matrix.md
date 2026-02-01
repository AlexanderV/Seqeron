# Position Weight Matrix (PWM)

**Algorithm Group:** Pattern Matching
**Test Unit:** PAT-PWM-001
**Last Updated:** 2026-01-22

---

## 1. Overview

A Position Weight Matrix (PWM), also known as Position-Specific Scoring Matrix (PSSM), is a mathematical model representing sequence motifs in biological sequences. PWMs are widely used for identifying transcription factor binding sites and other conserved sequence patterns.

---

## 2. Theory

### 2.1 Definition (Wikipedia)

> A position weight matrix (PWM) is a commonly used representation of motifs
> (patterns) in biological sequences. PWMs are often derived from a set of aligned
> sequences that are thought to be functionally related.

### 2.2 Construction Process

**Step 1: Position Frequency Matrix (PFM)**

Count occurrences of each nucleotide at each position:

$$PFM_{k,j} = \sum_{i=1}^{N} \mathbf{1}(X_{i,j} = k)$$

Where:
- N = number of sequences
- j = position (1 to L, sequence length)
- k ∈ {A, C, G, T}

**Step 2: Position Probability Matrix (PPM)**

Normalize by number of sequences:

$$PPM_{k,j} = \frac{PFM_{k,j}}{N}$$

**Step 3: Position Weight Matrix (PWM) with Log-Odds**

Convert to log-odds using background model:

$$PWM_{k,j} = \log_2\left(\frac{PPM_{k,j} + p}{b_k}\right)$$

Where:
- p = pseudocount (typically 0.25)
- b_k = background frequency (typically 0.25 for uniform)

### 2.3 Pseudocounts (Nishida et al., 2008)

> Pseudocounts are often applied when calculating PPMs if based on a small dataset,
> in order to avoid matrix entries having a value of 0.

This prevents -∞ scores for unseen nucleotides. Common pseudocount values:
- 0.25: Equivalent to adding one observation distributed equally
- 0.5: More conservative smoothing
- 1.0: Strong smoothing (Laplace estimator)

### 2.4 Scoring (Wikipedia)

> The sequence score can be calculated by adding the relevant values at each
> position in the PWM. The score gives an indication of how different the sequence
> is from a random sequence.

$$Score(S) = \sum_{j=1}^{L} PWM_{S_j, j}$$

**Interpretation:**
- Score = 0: Equal probability as random
- Score > 0: More likely functional
- Score < 0: More likely random

### 2.5 Properties

| Property | Description |
|----------|-------------|
| Length | Number of positions (width of motif) |
| Consensus | String with highest-scoring base at each position |
| MaxScore | Sum of maximum values at each position |
| MinScore | Sum of minimum values at each position |

---

## 3. Complexity

| Operation | Time | Space |
|-----------|------|-------|
| CreatePwm | O(N × L) | O(L) |
| ScanWithPwm | O(S × L) | O(1) per match |

Where:
- N = number of training sequences
- L = motif length
- S = target sequence length

---

## 4. Implementation Notes

### 4.1 Current Implementation

The `MotifFinder.CreatePwm()` method implements the standard PWM construction:

```
1. Validate input sequences (non-empty, equal length)
2. Count nucleotide frequencies at each position
3. Apply pseudocounts (default 0.25)
4. Convert to log-odds with background 0.25
5. Generate consensus from maximum scores
```

### 4.2 Matrix Layout

```
Matrix[4, Length] where:
  Row 0 = A scores
  Row 1 = C scores
  Row 2 = G scores
  Row 3 = T scores
```

### 4.3 Scanning Behavior

The `MotifFinder.ScanWithPwm()` method:
- Slides PWM along sequence
- Computes sum of position scores
- Returns matches where score ≥ threshold
- Skips positions with non-ACGT characters (ASSUMPTION)

### 4.4 Edge Cases

| Case | Behavior |
|------|----------|
| Empty input | ArgumentException |
| Unequal lengths | ArgumentException |
| Single sequence | Valid PWM with high confidence |
| Non-ACGT in training | Counted as 0 at that position |
| Non-ACGT in scanning | Position skipped, match invalidated |

---

## 5. Sources

| Source | Reference |
|--------|-----------|
| Wikipedia | https://en.wikipedia.org/wiki/Position_weight_matrix |
| Kel et al. (2003) | MATCH: A tool for searching TF binding sites. Nucleic Acids Res. 31(13):3576-3579 |
| Nishida et al. (2008) | Pseudocounts for transcription factor binding sites. Nucleic Acids Res. 37(3):939-944 |
| Rosalind | https://rosalind.info/problems/cons/ (Consensus and Profile) |
| Stormo (2000) | DNA binding sites: representation and discovery. Bioinformatics Review |

---

## 6. Related Algorithms

- **Consensus Sequence:** Simplified motif representation (single character per position)
- **Hidden Markov Models:** Extension with insertion/deletion probabilities (Pfam)
- **IUPAC Degenerate Matching:** Pattern matching with ambiguity codes
- **Information Content:** Measure of PWM specificity

---

## 7. ASSUMPTIONS

| ID | Assumption | Rationale |
|----|------------|-----------|
| A1 | Background frequency 0.25 for all bases | Standard for DNA without GC bias |
| A2 | Default pseudocount 0.25 | Implementation default |
| A3 | Non-ACGT characters invalidate match | Implementation-specific |
