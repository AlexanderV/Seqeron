# Approximate Matching: Hamming Distance

| Field | Value |
|-------|-------|
| Algorithm Group | Pattern Matching |
| Test Unit ID | PAT-APPROX-001 |
| Related Projects | N/A |
| Implementation Status | N/A |
| Last Reviewed | 2026-04-30 |

## 1. Overview

Hamming distance is a substitutions-only string metric for equal-length sequences. In this repository, it supports both direct distance computation and brute-force approximate pattern matching with a maximum mismatch threshold. The implementation is case-insensitive and records mismatch positions for each reported approximate match.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

Hamming distance is widely used for point-mutation analysis, SNP detection, approximate motif search, and primer or probe specificity checks when only substitutions are allowed. Because it forbids insertions and deletions, it applies to aligned strings of equal length. Sources: Hamming (1950), Wikipedia (Hamming distance), Rosalind HAMM, Navarro (2001), Gusfield (1997), Nicolae & Rajasekaran (2015).

### 2.2 Core Model

For equal-length strings `s` and `t`:

$$
d_H(s, t) = |\{ i : s[i] \neq t[i], 0 \le i < n \}|
$$

Approximate matching with at most `k` mismatches compares the pattern to every equal-length window in the target sequence and reports windows whose Hamming distance is at most `k`.

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | `d_H(s, t) >= 0` | Distance is a count of mismatching positions |
| INV-02 | `d_H(s, t) = 0` iff the compared strings are identical under case-insensitive comparison | The implementation increments only on unequal uppercased characters |
| INV-03 | `d_H(s, t) = d_H(t, s)` | Positionwise mismatch counting is symmetric |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `s1`, `s2` | `string` | required | Strings compared by `HammingDistance(...)` | Must be non-null and of equal length |
| `[FindWithMismatches(string)] sequence` | `string` | required | Sequence searched for approximate matches | Empty or null input yields no matches |
| `[FindWithMismatches(DnaSequence)] sequence` | `DnaSequence` | required | Sequence searched for approximate matches through the typed wrapper | Null input throws because the wrapper dereferences `sequence.Sequence` |
| `pattern` | `string` | required | Pattern compared against equal-length windows | Empty or null pattern yields no matches |
| `maxMismatches` | `int` | required | Maximum allowed substitutions | Negative values throw `ArgumentOutOfRangeException` |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `distance` | `int` | Hamming distance between two strings or spans |
| `Position` | `int` | Start position of an approximate match |
| `MatchedSequence` | `string` | Sequence window that matched within tolerance |
| `Distance` | `int` | Observed mismatch count |
| `MismatchPositions` | `IReadOnlyList<int>` | Zero-based mismatch positions |
| `MismatchType` | `MismatchType` | Always `Substitution` for Hamming-based matching |

### 3.3 Preconditions and Validation

`ApproximateMatcher.HammingDistance(string, string)` throws `ArgumentNullException` for null inputs and `ArgumentException` for unequal lengths. `FindWithMismatches(string, ...)` returns no matches when the sequence or pattern is null or empty or when the pattern is longer than the sequence. The `DnaSequence` wrapper overloads do not add a null guard and therefore throw when `sequence` is null because they dereference `sequence.Sequence`. Matching uppercases both sequence and pattern before comparison and throws `ArgumentOutOfRangeException` when `maxMismatches < 0`.

## 4. Algorithm

### 4.1 High-Level Steps

1. For direct distance, verify equal lengths and compare the strings character by character.
2. For approximate matching, uppercase the sequence and pattern.
3. Slide an equal-length window across the sequence.
4. Count mismatches between the window and pattern and collect mismatch positions.
5. Yield windows whose mismatch count is at most `maxMismatches`.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `HammingDistance` | `O(n)` | `O(1)` | Single pass over aligned strings |
| `FindWithMismatches` | `O(n × m)` | `O(z)` | `n` is sequence length, `m` is pattern length, and `z` is number of reported matches |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [ApproximateMatcher.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/ApproximateMatcher.cs), [SequenceExtensions.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Core/SequenceExtensions.cs)

- `ApproximateMatcher.HammingDistance(string, string)`: Case-insensitive string Hamming distance.
- `ApproximateMatcher.FindWithMismatches(...)`: Brute-force approximate pattern matching with substitutions only.
- `SequenceExtensions.HammingDistance(ReadOnlySpan<char>, ReadOnlySpan<char>)`: Case-insensitive span-based distance helper.

### 5.2 Current Behavior

The string and span Hamming-distance helpers both compare uppercased characters, so they are case-insensitive. `FindWithMismatches(...)` returns `ApproximateMatchResult` values whose `MismatchType` is always `Substitution`, and the `DnaSequence` overloads are thin wrappers over the string implementation that dereference `sequence.Sequence` without adding an extra null guard.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Hamming distance over equal-length strings.
- Approximate matching with a maximum mismatch threshold.
- Reporting of mismatch counts and mismatch positions for each match.

**Intentionally simplified:**

- Approximate matching uses a brute-force sliding-window scan; **consequence:** the implementation favors clarity over asymptotically faster mismatch-search algorithms.

**Not implemented:**

- Insertions and deletions in the Hamming-based path; **users should rely on:** the edit-distance workflow for those operations.

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Empty string sequence or pattern | Returns no matches | Explicit source guard |
| Null `DnaSequence` input | Throws | The typed wrapper dereferences `sequence.Sequence` |
| Pattern longer than sequence | Returns no matches | No equal-length windows exist |
| `maxMismatches = 0` | Equivalent to exact matching | Only windows with zero mismatches qualify |
| `maxMismatches >= pattern.Length` | All windows of equal length qualify | Every possible mismatch count is within threshold |

### 6.2 Limitations

The Hamming-based workflow models substitutions only and requires equal-length comparisons. It is not suitable when insertions or deletions are biologically relevant or when a more specialized approximate-matching algorithm is required for scale.

## 7. Examples and Related Material

### 7.1 Worked Example

**Numerical / biological walk-through (optional):**

Rosalind HAMM sample:

```text
GAGCCTACTAACGGGAT
CATCGTAATGACGGCCT
```

Expected Hamming distance: `7`.

## 8. References

1. Hamming, R.W. (1950). "Error detecting and error correcting codes." *Bell System Technical Journal*, 29(2): 147-160.
2. Wikipedia. "Hamming distance." https://en.wikipedia.org/wiki/Hamming_distance
3. Rosalind. "Counting Point Mutations." https://rosalind.info/problems/hamm/
4. Navarro, G. (2001). "A guided tour to approximate string matching." *ACM Computing Surveys*, 33(1): 31-88.
5. Gusfield, D. (1997). *Algorithms on Strings, Trees and Sequences.* Cambridge University Press.
6. Nicolae, M. & Rajasekaran, S. (2015). "On string matching with mismatches." *Algorithms*, 8(2): 248-270.
