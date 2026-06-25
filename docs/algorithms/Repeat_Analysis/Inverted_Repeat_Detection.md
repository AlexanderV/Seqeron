# Inverted Repeat Detection

| Field | Value |
|-------|-------|
| Algorithm Group | Repeat Analysis |
| Test Unit ID | REP-INV-001 |
| Related Projects | Seqeron.Genomics.Analysis |
| Implementation Status | Simplified |
| Last Reviewed | 2026-04-30 |

## 1. Overview

Inverted repeat detection identifies a sequence segment followed downstream by its reverse complement, optionally separated by a loop [1][2]. Such structures can form stem-loops or hairpins in single-stranded contexts and are closely related to palindromes, which are the special case with loop length zero [1]. The repository implements exact inverted-repeat detection in `RepeatFinder.FindInvertedRepeats`, returning explicit arm coordinates, sequences, loop sequence, and a `CanFormHairpin` flag. The implementation is deterministic and exact, but unlike score-based tools such as EMBOSS `einverted`, it does not tolerate mismatches or gaps.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

An inverted repeat has the form:

```text
5'---TTACG------nnnnnn------CGTAA---3'
  Left arm      Loop       Right arm
```

where the right arm is the reverse complement of the left arm [1][2]. Related terminology includes:

- Stem-loop or hairpin: a folded single-stranded structure with paired arms and an intervening loop [1].
- Palindrome: an inverted repeat with loop length `0` [1].
- Cruciform: a double-stranded extrusion formed by inverted-repeat regions [2].

The legacy reference set notes that hairpin loops are often optimal around 4-8 bases and that loops shorter than 3 bases are sterically unfavorable [1][2]. Inverted repeats occur at replication origins, transposon boundaries, rho-independent terminators, riboswitches, and tRNA structural elements, and long inverted repeats are associated with genomic instability, deletions, recombination, and mutation hotspots [2][3].

### 2.2 Core Model

For a left arm $L$, loop $X$, and right arm $R$, an inverted repeat satisfies:

$$
R = \operatorname{ReverseComplement}(L)
$$

with total structure length:

$$
\mathrm{TotalLength} = 2 \times \text{ArmLength} + \text{LoopLength}
$$

The implementation searches all possible left-arm starts and arm lengths, computes the reverse complement of each candidate left arm, and then checks downstream substrings within the configured loop-length window.

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | `ReverseComplement(LeftArm) = RightArm` for every result. | A result is emitted only when the candidate right arm equals the computed reverse complement. |
| INV-02 | `TotalLength = 2 × ArmLength + LoopLength`. | `InvertedRepeatResult.TotalLength` is defined directly from those fields. |
| INV-03 | `LoopLength = RightArmStart - (LeftArmStart + ArmLength)`. | Loop length is computed from the stored coordinates. |
| INV-04 | `CanFormHairpin` is true exactly when `LoopLength >= 3`. | The constructor sets `CanFormHairpin` from that Boolean test. |
| INV-05 | Each `(LeftArmStart, RightArmStart, ArmLength)` tuple is unique. | A hash set suppresses duplicate results. |

### 2.5 Comparison with Related Implementations

| Feature | Repository implementation | EMBOSS `einverted` |
|---------|---------------------------|--------------------|
| Matching model | Exact arm matching | Dynamic programming [4] |
| Mismatches | Not allowed | Allowed with penalty [4] |
| Gaps | Not allowed | Allowed with penalty [4] |
| Acceptance rule | Arm length and loop constraints | Score threshold [4] |
| Overlap handling | Duplicate coordinates suppressed | Tool-specific reporting |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `sequence` | `DnaSequence` or `string` | required | DNA sequence to search. | The `DnaSequence` overload throws on `null`; the raw-string overload yields no results for `null` or empty input. |
| `minArmLength` | `int` | `4` | Minimum length of each repeat arm. | Both overloads reject values below `2` with `ArgumentOutOfRangeException`. |
| `maxLoopLength` | `int` | `50` | Maximum loop length between arms. | Used as an upper bound when scanning downstream candidate starts. |
| `minLoopLength` | `int` | `3` | Minimum loop length between arms. | Both overloads reject negative values with `ArgumentOutOfRangeException`. |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `LeftArmStart` | `int` | 0-based start position of the left arm. |
| `RightArmStart` | `int` | 0-based start position of the right arm. |
| `ArmLength` | `int` | Length of each arm. |
| `LoopLength` | `int` | Number of intervening nucleotides between the two arms. |
| `LeftArm` | `string` | Left-arm sequence. |
| `RightArm` | `string` | Right-arm sequence, equal to the reverse complement of `LeftArm`. |
| `Loop` | `string` | Intervening sequence between the arms. |
| `CanFormHairpin` | `bool` | `true` when `LoopLength >= 3`. |

### 3.3 Preconditions and Validation

`FindInvertedRepeats(DnaSequence, ...)` throws `ArgumentNullException` when `sequence` is `null`, throws `ArgumentOutOfRangeException` when `minArmLength < 2`, and throws `ArgumentOutOfRangeException` when `minLoopLength < 0`. The raw-string overload uppercases non-empty input and yields no results for `null` or empty strings; it now enforces the SAME numeric validation as the `DnaSequence` overload (`minArmLength < 2` and `minLoopLength < 0` both throw `ArgumentOutOfRangeException`), so a degenerate `minArmLength = 0` can no longer emit nonsense zero-length-arm results on either surface. Coordinates are 0-based throughout the returned results.

## 4. Algorithm

### 4.1 High-Level Steps

1. Normalize the sequence to uppercase when the raw-string overload is used.
2. For each left-arm start position that leaves room for two minimum-length arms and the minimum loop, choose a candidate left-arm start.
3. For each arm length from `minArmLength` upward, extract the left arm and compute its reverse complement.
4. Search downstream start positions from `i + armLength + minLoopLength` through the configured loop-length bound.
5. When the downstream substring equals the reverse complement, compute loop length and loop sequence, suppress duplicate coordinate triples, and emit an `InvertedRepeatResult`.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| Inverted-repeat detection | `O(n × A^2 × L)` | `O(k + A)` | `A` is the maximum tested arm length, `L` is the effective loop-length search bound, and each candidate performs substring plus reverse-complement work proportional to arm length. |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [RepeatFinder.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/RepeatFinder.cs)

- `RepeatFinder.FindInvertedRepeats(DnaSequence, int, int, int)`: Validating overload for `DnaSequence` input.
- `RepeatFinder.FindInvertedRepeats(string, int, int, int)`: Uppercases raw string input and yields results for non-empty strings.

### 5.2 Current Behavior

The implementation uses `DnaSequence.GetReverseComplementString()` to derive the right-arm target from each candidate left arm, then scans downstream positions bounded by `minLoopLength` and `maxLoopLength`. Results are deduplicated with a hash set keyed by `(LeftArmStart, RightArmStart, ArmLength)`. `CanFormHairpin` is set from `loopLength >= 3`, and loops of length `0` are possible only if the caller explicitly lowers `minLoopLength` below its default. The raw-string overload normalizes case with `ToUpperInvariant()`.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Exact reverse-complement matching between left and right arms [1][2].
- Explicit reporting of arm positions, arm length, loop length, loop sequence, and total structure length.
- Hairpin-viability flag derived from the biologically motivated `loopLength >= 3` rule described in the legacy doc [1][2].

**Intentionally simplified:**

- Exact matching only, with no mismatch or gap penalties; **consequence:** approximate inverted repeats that would be scored by dynamic-programming tools are not reported.

**Not implemented:**

- Score-based inverted-repeat search with mismatch and gap tolerance; **users should rely on:** EMBOSS `einverted` when they need that richer search model [4].

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | The raw-string overload enforces the same numeric validation as the `DnaSequence` overload. | Resolved | Invalid `minArmLength` or `minLoopLength` values now throw `ArgumentOutOfRangeException` on both surfaces; no nonsense zero-length-arm output on `minArmLength = 0`. | resolved | Validation asymmetry removed during the REP-INV-001 fuzzing pass (a `minArmLength = 0` previously emitted spurious empty-arm results via the raw-string overload). |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Empty sequence | Returns empty enumerable. | No structure can be formed. |
| Sequence shorter than `2 × minArmLength + minLoopLength` | Returns empty enumerable. | There is not enough space for two arms and the required loop. |
| No complementary regions | Returns empty enumerable. | No downstream substring matches a reverse complement candidate. |
| Homopolymer input such as `AAAA` | Typically returns empty. | The reverse complement of `AAAA` is `TTTT`, which is absent from the same homopolymer. |
| Self-complementary sequence such as `GCGC` | Can match when loop-length settings permit it. | A palindrome is an inverted repeat with loop length `0`. |
| Loop length `0` | Not returned under default settings. | The default `minLoopLength` is `3`. |

### 6.2 Limitations

The algorithm is DNA-specific and uses DNA complement rules. It does not score approximate inverted repeats or RNA base-pairing variants. For RNA-specific inverted-repeat discovery, the repository provides a separate helper in [RnaSecondaryStructure.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/RnaSecondaryStructure.cs) with RNA complement rules and tuple-based output.

## 7. Examples and Related Material

### 7.3 Related Tests, Evidence, or Documents

- Tests: [RepeatFinder_InvertedRepeat_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/RepeatFinder_InvertedRepeat_Tests.cs)
- Test spec: [REP-INV-001.md](../../../tests/TestSpecs/REP-INV-001.md)
- Related RNA smoke tests: [RnaSecondaryStructureTests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/RnaSecondaryStructureTests.cs)

## 8. References

1. Wikipedia. 2026. Inverted repeat. Wikipedia. https://en.wikipedia.org/wiki/Inverted_repeat
2. Pearson CE, Zorbas H, Price GB, Zannis-Hadjopoulos M. 1996. Inverted repeats, stem-loops, and cruciforms: significance for initiation of DNA replication. Journal of Cellular Biochemistry. 63(1):1-22.
3. Bissler JJ. 1998. DNA inverted repeats and human disease. Frontiers in Bioscience. 3:d408-d418.
4. Rice P, Longden I, Bleasby A. 2000. EMBOSS: the European Molecular Biology Open Software Suite. Trends in Genetics. 16(6):276-277.
