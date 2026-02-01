# Approximate Matching: Hamming Distance

## Overview

Hamming distance is a fundamental string metric used for approximate pattern matching when only substitutions (mismatches) are allowed. It is widely used in bioinformatics for identifying point mutations, SNP detection, and primer/probe binding analysis.

## Definition

The **Hamming distance** between two strings $s$ and $t$ of equal length $n$ is the number of positions at which the corresponding symbols differ:

$$d_H(s, t) = |\{i : s[i] \neq t[i], 0 \leq i < n\}|$$

### Constraint

Hamming distance is **only defined for strings of equal length**. For strings of unequal length, use Edit distance (Levenshtein distance) instead.

## Mathematical Properties

From Wikipedia and Robinson (2003):

| Property | Definition | Implication |
|----------|------------|-------------|
| Non-negativity | $d_H(s, t) \geq 0$ | Distance is always non-negative |
| Identity | $d_H(s, t) = 0 \iff s = t$ | Zero distance means identical strings |
| Symmetry | $d_H(s, t) = d_H(t, s)$ | Order of comparison doesn't matter |
| Triangle inequality | $d_H(s, u) \leq d_H(s, t) + d_H(t, u)$ | Forms a proper metric space |

## Algorithm

### Basic Hamming Distance (O(n))

```
function HammingDistance(s1, s2):
    if length(s1) ≠ length(s2):
        throw "Strings must have equal length"

    distance = 0
    for i = 0 to length(s1) - 1:
        if s1[i] ≠ s2[i]:
            distance++

    return distance
```

### Pattern Matching with k Mismatches (O(n × m))

The approximate pattern matching problem finds all positions in a text where a pattern occurs with at most $k$ mismatches:

```
function FindWithMismatches(sequence, pattern, maxMismatches):
    results = []

    for i = 0 to length(sequence) - length(pattern):
        window = sequence[i : i + length(pattern)]
        distance = HammingDistance(window, pattern)

        if distance ≤ maxMismatches:
            results.append((i, window, distance))

    return results
```

## Complexity

| Operation | Time | Space |
|-----------|------|-------|
| Hamming distance | O(n) | O(1) |
| Find with k mismatches (brute force) | O(n × m) | O(z) where z = matches |
| Find with k mismatches (optimized) | O(n√k log k) | O(n + m) |

**Reference:** Nicolae & Rajasekaran (2015) achieved O(n√k log k) for k-mismatch pattern matching.

## Test Cases from Literature

### Rosalind HAMM Problem

**Problem:** Given two DNA strings of equal length, return the Hamming distance.

**Sample Input:**
```
GAGCCTACTAACGGGAT
CATCGTAATGACGGCCT
```

**Sample Output:** 7

**Source:** https://rosalind.info/problems/hamm/

### Wikipedia Examples

| String 1 | String 2 | Hamming Distance |
|----------|----------|------------------|
| "karolin" | "kathrin" | 3 |
| "karolin" | "kerstin" | 3 |
| "kathrin" | "kerstin" | 4 |
| "0000" | "1111" | 4 |

## Applications in Bioinformatics

### 1. Point Mutation Detection

Hamming distance measures the number of point mutations (substitutions) between two aligned sequences.

### 2. SNP Analysis

Single Nucleotide Polymorphisms (SNPs) are detected by comparing sequence variants where Hamming distance = 1 indicates a single SNP.

### 3. Approximate Motif Finding

Finding regulatory motifs that may contain mutations from a consensus sequence.

### 4. Primer/Probe Binding

Evaluating primer specificity by finding potential binding sites with allowed mismatches.

## Implementation Notes

### Current Implementation

The `ApproximateMatcher` class provides:

1. **`HammingDistance(string s1, string s2)`**
   - Case-insensitive comparison (normalized to uppercase)
   - Throws `ArgumentNullException` for null inputs
   - Throws `ArgumentException` for unequal lengths

2. **`FindWithMismatches(string sequence, string pattern, int maxMismatches)`**
   - Returns `IEnumerable<ApproximateMatchResult>`
   - Each result includes position, matched sequence, distance, and mismatch positions
   - Throws `ArgumentOutOfRangeException` for negative maxMismatches

3. **`HammingDistance(ReadOnlySpan<char>)` extension** (SequenceExtensions)
   - High-performance span-based API
   - Case-insensitive comparison

### Edge Case Handling

| Case | Behavior |
|------|----------|
| Empty pattern | Returns empty (no matches) |
| Empty sequence | Returns empty (no matches) |
| Pattern longer than sequence | Returns empty |
| maxMismatches = 0 | Equivalent to exact matching |
| maxMismatches ≥ pattern.Length | All positions match |

## References

1. **Hamming, R.W. (1950).** "Error detecting and error correcting codes." *Bell System Technical Journal*, 29(2): 147-160.

2. **Wikipedia.** "Hamming distance." https://en.wikipedia.org/wiki/Hamming_distance

3. **Rosalind.** "Counting Point Mutations." https://rosalind.info/problems/hamm/

4. **Navarro, G. (2001).** "A guided tour to approximate string matching." *ACM Computing Surveys*, 33(1): 31-88.

5. **Gusfield, D. (1997).** *Algorithms on Strings, Trees and Sequences.* Cambridge University Press.

6. **Nicolae, M. & Rajasekaran, S. (2015).** "On string matching with mismatches." *Algorithms*, 8(2): 248-270.

---

*Document generated for Test Unit PAT-APPROX-001*
*Last updated: 2026-01-22*
