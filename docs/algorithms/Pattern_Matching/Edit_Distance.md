# Edit Distance (Levenshtein Distance)

| Field | Value |
|-------|-------|
| Algorithm Group | Pattern Matching |
| Test Unit ID | PAT-APPROX-002 |
| Related Projects | N/A |
| Implementation Status | Simplified |
| Last Reviewed | 2026-04-30 |

## 1. Overview

Edit distance measures the minimum number of insertions, deletions, and substitutions required to transform one string into another. In this repository, the core Levenshtein distance implementation uses a two-row Wagner-Fischer dynamic program, and the approximate-search surface scans variable-length windows around the pattern length to find matches within a maximum edit threshold.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

Levenshtein distance is the standard edit-distance metric for strings when insertions, deletions, and substitutions all cost one unit. In bioinformatics it models sequence differences produced by insertion, deletion, and substitution events and provides a more appropriate notion of distance than Hamming distance when indels are possible. Sources: Levenshtein (1966), Wagner & Fischer (1974), Wikipedia (Levenshtein distance, Edit distance), Berger et al. (2021), Navarro (2001).

### 2.2 Core Model

The recursive definition preserved from the original document is:

$$
lev(a,b) = \begin{cases}
|a| & \text{if } |b| = 0 \\
|b| & \text{if } |a| = 0 \\
lev(tail(a), tail(b)) & \text{if } head(a) = head(b) \\
1 + \min\begin{cases}
lev(tail(a), b) \\
lev(a, tail(b)) \\
lev(tail(a), tail(b))
\end{cases} & \text{otherwise}
\end{cases}
$$

The implemented search routine `FindWithEdits(...)` compares the pattern against windows whose lengths range from `pattern.Length - maxEdits` to `pattern.Length + maxEdits`.

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | `d(a, b) = 0` iff the strings are identical under the exact comparison used by `EditDistance(...)` | Zero cost is incurred only when every aligned character matches and lengths agree |
| INV-02 | `d(a, b) >= |len(a) - len(b)|` | At least the length difference must be repaired by insertions or deletions |
| INV-03 | `d(a, b) <= max(len(a), len(b))` | One string can be transformed into the other by deletions plus substitutions |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `s1`, `s2` | `string` | required | Strings compared by `EditDistance(...)` | Null input throws `ArgumentNullException` |
| `[FindWithEdits(string)] sequence` | `string` | required | Sequence searched by `FindWithEdits(...)` | Null or empty input yields no matches |
| `[FindWithEdits(DnaSequence)] sequence` | `DnaSequence` | required | Sequence searched through the typed wrapper | Null input throws because the wrapper dereferences `sequence.Sequence` |
| `pattern` | `string` | required | Pattern compared against variable-length windows | Null or empty input yields no matches |
| `maxEdits` | `int` | required | Maximum allowed edit distance | Negative values throw `ArgumentOutOfRangeException` |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `distance` | `int` | Levenshtein distance between two strings |
| `Position` | `int` | Start of a matching window in `FindWithEdits(...)` |
| `MatchedSequence` | `string` | Window whose edit distance is within threshold |
| `Distance` | `int` | Observed edit distance |
| `MismatchType` | `MismatchType` | `Substitution` when the edit distance equals the Hamming distance on equal-length windows; otherwise `Edit` |

### 3.3 Preconditions and Validation

`EditDistance(...)` throws `ArgumentNullException` when either string is null. `FindWithEdits(string, ...)` returns no matches when the sequence or pattern is null or empty and throws `ArgumentOutOfRangeException` when `maxEdits < 0`. The `DnaSequence` wrapper overload dereferences `sequence.Sequence` directly and therefore throws when `sequence` is null. The approximate-search method uppercases both sequence and pattern before scanning, but the core `EditDistance(...)` routine compares characters as-is.

## 4. Algorithm

### 4.1 High-Level Steps

1. For direct distance, initialize a two-row dynamic-programming table.
2. Fill the current row from the previous row using insertion, deletion, and substitution costs.
3. Return the final value in the last DP row as the Levenshtein distance.
4. For approximate search, uppercase the sequence and pattern.
5. Compare the pattern against every window whose length is between `pattern.Length - maxEdits` and `pattern.Length + maxEdits`.
6. Yield windows whose edit distance is at most `maxEdits`.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `EditDistance` | `O(m × n)` | `O(n)` | Two-row Wagner-Fischer dynamic program for strings of lengths `m` and `n` |
| `FindWithEdits` | `O(s × (2e + 1) × p × (p + e))` | `O(p + e)` | `s` = sequence length, `p` = pattern length, `e` = `maxEdits`; each candidate window invokes `EditDistance(...)` on a window whose length ranges from `p - e` to `p + e` |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [ApproximateMatcher.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/ApproximateMatcher.cs)

- `ApproximateMatcher.EditDistance(string, string)`: Two-row Levenshtein distance.
- `ApproximateMatcher.FindWithEdits(string, string, int)`: Approximate search by variable-length windows.
- `ApproximateMatcher.FindWithEdits(DnaSequence, string, int)`: Typed wrapper over the string implementation.

### 5.2 Current Behavior

The core `EditDistance(...)` method is case-sensitive because it compares characters directly. `FindWithEdits(...)` uppercases both the sequence and pattern before scanning and distinguishes substitution-only matches from general edits by comparing the edit distance to a helper Hamming-like distance on equal-length windows. For edit matches that involve insertions or deletions, `MismatchPositions` is returned as an empty list. The `DnaSequence` overload is a thin wrapper over the string implementation and does not add its own null guard.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Levenshtein distance with insertion, deletion, and substitution costs of one.
- Space-optimized Wagner-Fischer dynamic programming.
- Approximate search by accepting windows with edit distance at most `maxEdits`.

**Intentionally simplified:**

- `FindWithEdits(...)` uses a brute-force variable-window scan; **consequence:** it is easy to reason about but not optimized for very large texts or very permissive edit thresholds.
- Edit-match results do not reconstruct insertion or deletion coordinates; **consequence:** callers receive the edit distance and matched window, but not a full alignment trace.

**Not implemented:**

- Damerau-style transpositions or other extended edit operations; **users should rely on:** specialized edit-distance variants if those operations are required.

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| `""` vs `"abc"` | Distance `3` | Three insertions or deletions are required |
| Empty string sequence or pattern in `FindWithEdits(...)` | Returns no matches | Explicit source guard |
| Null `DnaSequence` input in `FindWithEdits(...)` | Throws | The typed wrapper dereferences `sequence.Sequence` |
| `maxEdits < 0` | Throws `ArgumentOutOfRangeException` | Invalid threshold |
| Equal-length strings | Levenshtein distance is at most the Hamming distance | Substitutions alone can realize the Hamming path |

### 6.2 Limitations

The current search routine is a brute-force approximation layer over the core distance function and does not expose a full alignment traceback. The core distance method is also case-sensitive, so callers who need normalized comparisons must uppercase or otherwise normalize inputs before calling it directly.

## 7. Examples and Related Material

### 7.1 Worked Example

**Numerical / biological walk-through (optional):**

`kitten -> sitting` has edit distance `3`:

1. `kitten -> sitten`
2. `sitten -> sittin`
3. `sittin -> sitting`

## 8. References

1. Levenshtein, V.I. (1966). "Binary codes capable of correcting deletions, insertions, and reversals." Soviet Physics Doklady, 10(8): 707–710.
2. Wagner, R.A.; Fischer, M.J. (1974). "The String-to-String Correction Problem." Journal of the ACM, 21(1): 168–173.
3. Navarro, G. (2001). "A guided tour to approximate string matching." ACM Computing Surveys, 33(1): 31–88.
4. Berger, B.; Waterman, M.S.; Yu, Y.W. (2021). "Levenshtein Distance, Sequence Comparison and Biological Database Search." IEEE Transactions on Information Theory, 67(6): 3287–3294.
5. Rosetta Code - Levenshtein Distance: https://rosettacode.org/wiki/Levenshtein_distance
6. Wikipedia - Levenshtein Distance: https://en.wikipedia.org/wiki/Levenshtein_distance
7. Wikipedia - Edit Distance: https://en.wikipedia.org/wiki/Edit_distance
