# Edit Distance (Levenshtein Distance)

## Overview

Edit distance is a string metric for measuring the difference between two sequences. The Levenshtein distance, named after Soviet mathematician Vladimir Levenshtein (1965), counts the minimum number of single-character edits (insertions, deletions, substitutions) required to transform one string into the other.

**Test Unit:** PAT-APPROX-002  
**Source Files:** `ApproximateMatcher.cs`

---

## Mathematical Definition

The Levenshtein distance between two strings $a, b$ (of length $|a|$ and $|b|$ respectively) is given by:

$$
\text{lev}(a,b) = \begin{cases}
|a| & \text{if } |b| = 0 \\
|b| & \text{if } |a| = 0 \\
\text{lev}(\text{tail}(a), \text{tail}(b)) & \text{if } \text{head}(a) = \text{head}(b) \\
1 + \min \begin{cases}
\text{lev}(\text{tail}(a), b) & \text{(deletion)} \\
\text{lev}(a, \text{tail}(b)) & \text{(insertion)} \\
\text{lev}(\text{tail}(a), \text{tail}(b)) & \text{(substitution)}
\end{cases} & \text{otherwise}
\end{cases}
$$

**Source:** Wikipedia - Levenshtein Distance

---

## Properties and Invariants

### Metric Properties
Edit distance with non-negative cost satisfies the axioms of a metric space:

1. **Identity:** $d(a, b) = 0$ if and only if $a = b$
2. **Positivity:** $d(a, b) > 0$ when $a \neq b$  
3. **Symmetry:** $d(a, b) = d(b, a)$
4. **Triangle inequality:** $d(a, c) \leq d(a, b) + d(b, c)$

**Source:** Wikipedia - Edit Distance

### Upper and Lower Bounds

- Minimum: absolute value of the difference of string sizes
- Maximum: length of the longer string
- Zero if and only if strings are equal
- For equal-length strings: Hamming distance ≥ Levenshtein distance
- **Example:** "flaw" → "lawn" has Levenshtein distance 2 (delete 'f', insert 'n'), but Hamming distance 4

**Source:** Wikipedia - Levenshtein Distance

---

## Canonical Test Cases

From authoritative sources:

| String 1 | String 2 | Distance | Source |
|----------|----------|----------|--------|
| "kitten" | "sitting" | 3 | Wikipedia, Rosetta Code |
| "rosettacode" | "raisethysword" | 8 | Rosetta Code |
| "saturday" | "sunday" | 3 | Rosetta Code |
| "" | "abc" | 3 | Definition (insertion of all chars) |
| "abc" | "" | 3 | Definition (deletion of all chars) |
| "flaw" | "lawn" | 2 | Wikipedia |

### "kitten" → "sitting" (distance = 3)
1. kitten → sitten (substitution of 's' for 'k')
2. sitten → sittin (substitution of 'i' for 'e')
3. sittin → sitting (insertion of 'g' at end)

**Source:** Wikipedia - Levenshtein Distance

---

## Algorithm Complexity

- **Time:** O(m × n) where m and n are string lengths (Wagner-Fischer algorithm)
- **Space:** O(min(m, n)) with two-row optimization

The implementation uses the space-optimized two-row variant.

**Source:** Wikipedia - Wagner-Fischer algorithm

---

## Implementation Notes

### `ApproximateMatcher.EditDistance(string s1, string s2)`

The implementation in this library:
- Uses the two-row space-optimized Wagner-Fischer algorithm
- Normalizes input to uppercase (case-insensitive comparison)
- Returns `n` if `s1` is empty (length of `s2`)
- Returns `m` if `s2` is empty (length of `s1`)
- Throws `ArgumentNullException` if either input is null

### `ApproximateMatcher.FindWithEdits(string sequence, string pattern, int maxEdits)`

Finds all approximate matches using edit distance:
- Uses sliding window with variable length (pattern ± maxEdits)
- Returns matches with distance ≤ maxEdits
- Distinguishes substitution-only matches from edit matches (insertions/deletions)

---

## Bioinformatics Applications

In bioinformatics, Levenshtein distance measures the difference between biological sequences. The edits correspond to genetic mutations:
- **Insertion:** Addition of a nucleotide
- **Deletion:** Removal of a nucleotide  
- **Substitution:** Replacement of one nucleotide with another

A lower distance indicates closer evolutionary or functional relationship.

**Source:** Wikipedia - Levenshtein Distance (Bioinformatics section), Berger et al. (2021)

---

## Related Algorithms

| Algorithm | Operations | Note |
|-----------|------------|------|
| Hamming distance | Substitution only | Equal-length strings only |
| Damerau-Levenshtein | + Transposition | Adjacent character swap |
| LCS distance | Insertion, Deletion | No substitution |

**Source:** Wikipedia - Edit Distance

---

## References

1. Levenshtein, V.I. (1966). "Binary codes capable of correcting deletions, insertions, and reversals." Soviet Physics Doklady, 10(8): 707–710.
2. Wagner, R.A.; Fischer, M.J. (1974). "The String-to-String Correction Problem." Journal of the ACM, 21(1): 168–173.
3. Navarro, G. (2001). "A guided tour to approximate string matching." ACM Computing Surveys, 33(1): 31–88.
4. Berger, B.; Waterman, M.S.; Yu, Y.W. (2021). "Levenshtein Distance, Sequence Comparison and Biological Database Search." IEEE Transactions on Information Theory, 67(6): 3287–3294.
5. Rosetta Code - Levenshtein Distance: https://rosettacode.org/wiki/Levenshtein_distance
6. Wikipedia - Levenshtein Distance: https://en.wikipedia.org/wiki/Levenshtein_distance
7. Wikipedia - Edit Distance: https://en.wikipedia.org/wiki/Edit_distance
