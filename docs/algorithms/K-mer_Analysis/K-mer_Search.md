# K-mer Search Algorithms

## Overview

K-mer search algorithms identify k-mers of interest within sequences, including the most frequent k-mers, unique (singleton) k-mers, and clumps (localized clusters of repeated k-mers).

## Theoretical Background

### Most Frequent K-mers

**Definition:** A pattern P is a most frequent k-mer in Text if it maximizes Count(Text, P) among all k-mers. Multiple k-mers may share the maximum count.

**Source:** Rosalind BA1B

**Example:** In "ACGTTGCATGTCGCATGATGCATGAGAGCT", the most frequent 4-mers are "CATG" and "GCAT" (each appearing 3 times).

### Unique K-mers

**Definition:** A k-mer that appears exactly once in a sequence. Unique k-mers are useful for:
- Identifying marker sequences
- Primer design (avoiding repetitive regions)
- Genomic fingerprinting

**Source:** Wikipedia (K-mer)

### Clump Finding

**Definition:** Given integers L and t, a pattern P forms an (L, t)-clump inside Genome if there exists an interval of Genome of length L in which P appears at least t times.

**Source:** Rosalind BA1E

**Biological significance:** Clumps suggest regulatory regions, origins of replication (DnaA boxes), or other functionally important sequences where multiple copies of a motif cluster together.

**Example:** "TGCA" forms a (25, 3)-clump in:
```
gatcagcataagggtcccTGCAATGCATGACAAGCCTGCAgttgttttac
```
Where "TGCA" appears 3 times within a 25-bp window.

## Algorithms

### FindMostFrequentKmers

**Input:** Sequence (string), k (k-mer length)  
**Output:** Set of k-mers with maximum count

**Algorithm:**
1. Count all k-mers in the sequence
2. Find the maximum count among all k-mers
3. Return all k-mers with count equal to maximum

**Complexity:** O(n) where n = sequence length

### FindUniqueKmers

**Input:** Sequence (string), k (k-mer length)  
**Output:** Set of k-mers appearing exactly once

**Algorithm:**
1. Count all k-mers in the sequence
2. Filter to retain only k-mers with count = 1

**Complexity:** O(n) where n = sequence length

### FindClumps

**Input:** Sequence (string), k (k-mer length), L (window size), t (minimum occurrences)  
**Output:** Set of k-mers forming (L, t)-clumps

**Algorithm (Sliding Window):**
1. Initialize counts for k-mers in first window of size L
2. Check if any k-mer meets threshold t; add to clumps
3. Slide window by 1 position:
   - Decrement count of k-mer leaving window
   - Increment count of k-mer entering window
   - Check for new clumps
4. Repeat until end of sequence

**Complexity:** O(n × (L - k + 1)) in worst case, typically O(n) with efficient data structures

## Implementation Notes

### Current Implementation

The `KmerAnalyzer` class provides:

- `FindMostFrequentKmers(string sequence, int k)`: Returns all k-mers with maximum count
- `FindUniqueKmers(string sequence, int k)`: Returns k-mers appearing exactly once
- `FindClumps(string sequence, int k, int windowSize, int minOccurrences)`: Returns k-mers forming clumps

**Implementation characteristics:**
- Case insensitive (sequences converted to uppercase)
- Returns empty collection for invalid parameters (empty sequence, k ≤ 0, k > sequence length, windowSize > sequence length)
- Clump results are returned as a set (no duplicates)

### Edge Cases

| Condition | Behavior |
|-----------|----------|
| Empty sequence | Returns empty collection |
| k ≤ 0 | Returns empty (FindClumps), throws for CountKmers |
| k > sequence length | Returns empty collection |
| windowSize > sequence length | Returns empty collection (FindClumps) |
| windowSize < k | Returns empty collection (FindClumps) |
| All k-mers equally frequent | All returned as most frequent |
| No clumps exist | Returns empty collection |

## References

1. **Rosalind BA1B** - Find the Most Frequent Words in a String
   - Problem definition and sample datasets for frequent k-mers
   - https://rosalind.info/problems/ba1b/

2. **Rosalind BA1E** - Find Patterns Forming Clumps in a String
   - (L, t)-clump definition and algorithm
   - Sample: "CGGACTCGACAGATGTGAAGAAATGTGAAGACTGAGTGAAGAGAAGAGGAAACACGACACGACATTGCGACATAATGTACGAATGTAATGTGCCTATGGC" with k=5, L=75, t=4
   - Expected clumps: CGACA, GAAGA, AATGT
   - https://rosalind.info/problems/ba1e/

3. **Wikipedia (K-mer)**
   - K-mer definition, k-mer spectrum, applications in bioinformatics
   - Pseudocode for k-mer extraction
   - https://en.wikipedia.org/wiki/K-mer

---

*Document generated: 2026-01-23*  
*Test Unit: KMER-FIND-001*
