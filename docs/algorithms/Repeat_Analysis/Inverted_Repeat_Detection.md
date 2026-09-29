# Inverted Repeat Detection

| Field | Value |
|-------|-------|
| Algorithm Group | Repeat Analysis |
| Test Unit ID | REP-INV-001 |
| Related Projects | Seqeron.Genomics.Analysis |
| Implementation Status | Complete (exact-stem finder; EMBOSS `palindrome` semantics with 0 mismatches) |
| Last Reviewed | 2026-09-29 |

## 1. Overview

Inverted repeat detection identifies a sequence segment followed downstream by its reverse complement, optionally separated by a loop [1][2]. Such structures can form stem-loops or hairpins in single-stranded contexts and are closely related to palindromes, which are the special case with loop length zero [1]. The repository implements exact inverted-repeat detection in `RepeatFinder.FindInvertedRepeats`, returning explicit arm coordinates, sequences, loop sequence, and a `CanFormHairpin` flag. Only **maximal** stems are reported, using the reporting rule of EMBOSS `palindrome` [5] run with 0 mismatches: a stem that lies inside another stem in both arms is dropped. It is an exact-stem finder by design; unlike the score-based EMBOSS `einverted` [4], it does not tolerate mismatches or gaps.

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

A stem is a triple $(i, j, A)$ with $s_{i+k}$ Watson–Crick-complementary to $s_{j+A-1-k}$ for $0 \le k < A$ (pairs among A, C, G, T only), $A \ge$ `minArmLength` and `minLoopLength` $\le j-(i+A) \le$ `maxLoopLength`. The reported set is the stems not contained in both arms of another stem (EMBOSS `palindrome_AInB`, `-overlap Y`). Equivalently: each reported stem cannot be extended outward, can be extended inward only by making the loop shorter than `minLoopLength`, and is not a "slipped" re-pairing lying inside a longer stem on a neighbouring diagonal.

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | `ReverseComplement(LeftArm) = RightArm` for every result. | A result is emitted only when the candidate right arm equals the computed reverse complement. |
| INV-02 | `TotalLength = 2 × ArmLength + LoopLength`. | `InvertedRepeatResult.TotalLength` is defined directly from those fields. |
| INV-03 | `LoopLength = RightArmStart - (LeftArmStart + ArmLength)`. | Loop length is computed from the stored coordinates. |
| INV-04 | `CanFormHairpin` is true exactly when `LoopLength >= 3`. | The constructor sets `CanFormHairpin` from that Boolean test. |
| INV-05 | No reported stem lies inside another reported stem in both arms; each tuple is unique. | Containment filter (EMBOSS `palindrome_AInB`). |
| INV-06 | Revcomp symmetry: the stems of revcomp(S) are the mirror images of the stems of S. | The definition is strand-symmetric. |

### 2.5 Comparison with Related Implementations

| Feature | Repository implementation | EMBOSS `einverted` |
|---------|---------------------------|--------------------|
| Matching model | Exact arm matching | Dynamic programming [4] |
| Mismatches | Not allowed | Allowed with penalty [4] |
| Gaps | Not allowed | Allowed with penalty [4] |
| Acceptance rule | Arm length and loop constraints | Score threshold [4] |
| Overlap handling | Stems contained in both arms of another stem dropped (= EMBOSS `palindrome -overlap Y`) | Best-scoring local alignments |
| Ambiguous bases | Only A/C/G/T pair | Only a/c/g/t score as match [4] |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `sequence` | `DnaSequence` or `string` | required | DNA sequence to search. | The `DnaSequence` overload throws on `null`; the raw-string overload yields no results for `null` or empty input. |
| `minArmLength` | `int` | `4` | Minimum length of each repeat arm. | Both overloads reject values below `2` with `ArgumentOutOfRangeException`. |
| `maxLoopLength` | `int` | `50` | Maximum loop length between arms. | Must be ≥ `minLoopLength` (`ArgumentOutOfRangeException`); values beyond the sequence length are capped. |
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

`FindInvertedRepeats(DnaSequence, ...)` throws `ArgumentNullException` when `sequence` is `null`. Both overloads throw `ArgumentOutOfRangeException` **eagerly (at the call, not on enumeration)** when `minArmLength < 2`, `minLoopLength < 0` or `maxLoopLength < minLoopLength`. The raw-string overload uppercases non-empty input and yields no results for `null` or empty strings; characters other than A/C/G/T (N, IUPAC codes, U) never pair. Coordinates are 0-based; results are ordered by `LeftArmStart`, then `RightArmStart`.

## 4. Algorithm

### 4.1 High-Level Steps

1. Normalize the sequence to uppercase when the raw-string overload is used.
2. For every innermost pair `(iIn, rIn)` with loop `L = rIn − iIn − 1` in `[minLoopLength, maxLoopLength]`: skip it unless the bases pair and the stem cannot move inward (the next inner pair fails, or it would leave a loop shorter than `minLoopLength`).
3. Extend the stem outward while the bases pair; discard it if the arm is shorter than `minArmLength`. These are EMBOSS `palindrome`'s candidate stems with 0 mismatches.
4. Drop a stem that lies inside another stem in both arms. Only two families of stems can contain it: for a shift `1 ≤ m ≤ L − minLoopLength`, the stem `(i, j, A + m)` (same arm starts) or `(i − m, j − m, A + m)` (same arm ends). Checking these is exact.
5. Emit arms, loop and `CanFormHairpin = L ≥ 3`, sorted by `(LeftArmStart, RightArmStart)`.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| Candidate scan | `O(n · W)` | `O(k)` | `W = maxLoopLength − minLoopLength + 1`, capped by `n`. |
| Stem extension + containment check | `O(A + W · A)` per candidate stem, early exit | — | Output-bound in periodic input; `int.MaxValue` loop bounds terminate (capped by `n`). |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [RepeatFinder.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/RepeatFinder.cs)

- `RepeatFinder.FindInvertedRepeats(DnaSequence, int, int, int)`: Validating overload for `DnaSequence` input.
- `RepeatFinder.FindInvertedRepeats(string, int, int, int)`: Uppercases raw string input and yields results for non-empty strings.

### 5.2 Current Behavior

Pairing uses the canonical `SequenceExtensions.GetComplementBase` restricted to A/C/G/T (no private complement table). The implementation scans innermost pairs within the loop window, extends each inward-maximal stem outward, and removes stems contained in both arms of another stem. With `minLoopLength = 0` its output equals the EMBOSS 6.6.0 `palindrome` binary (`-nummismatches 0 -overlap Y -gaplimit maxLoopLength -minpallen minArmLength`, unbounded `-maxpallen`) on 1000 random sequences, and it equals an independent Python brute force (enumerate all exact stems, keep the non-contained ones) on 3000 random cases with `minLoopLength` 0–6, N, U and lowercase ([REP-INV-001-Evidence.md](../../Evidence/REP-INV-001-Evidence.md)).

### 5.3 Conformance to Theory / Spec

**Implemented:**

- Exact reverse-complement matching between left and right arms [1][2].
- Maximal-stem reporting of EMBOSS `palindrome` with 0 mismatches [5], plus a minimum loop length.
- Hairpin-viability flag `loopLength >= 3` [1].

**Not implemented (by design):**

- Score-based inverted-repeat search with mismatch and gap tolerance; **users should rely on:** EMBOSS `einverted` [4] (or `palindrome -nummismatches k`) when they need that richer model.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | N / IUPAC codes never pair (EMBOSS `palindrome` lets `n` complement `n` and only rejects all-`n` stems). | Deviation (stricter) | A run such as `GGN…NCC` is not a 3-bp stem. | accepted | Follows EMBOSS `einverted`, which scores a match only for a/c/g/t. |
| 2 | `minLoopLength` (default 3) has no `palindrome` counterpart. | Extension | Inward extension stops at the minimum loop. | accepted | With `minLoopLength = 0` output equals `palindrome`. |
| 3 | No maximum arm length (`palindrome -maxpallen`). | Simplification of API | Long stems are reported whole. | accepted | Compared with unbounded `-maxpallen`. |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Empty sequence | Returns empty enumerable. | No structure can be formed. |
| Sequence shorter than `2 × minArmLength + minLoopLength` | Returns empty enumerable. | There is not enough space for two arms and the required loop. |
| No complementary regions | Returns empty enumerable. | No downstream substring matches a reverse complement candidate. |
| Homopolymer input such as `AAAA` | Typically returns empty. | The reverse complement of `AAAA` is `TTTT`, which is absent from the same homopolymer. |
| Self-complementary sequence such as `GCGC` | Can match when loop-length settings permit it. | A palindrome is an inverted repeat with loop length `0`. |
| Perfect stem longer than `minArmLength` | Reported once, as the whole stem. | Sub-stems lie inside it in both arms. |
| Run of `N` | Returns empty. | N never pairs. |
| Loop length `0` | Not returned under default settings. | The default `minLoopLength` is `3`. |

### 6.2 Limitations

The algorithm is DNA-specific and uses DNA complement rules. It does not score approximate inverted repeats or RNA base-pairing variants. For RNA-specific inverted-repeat discovery, the repository provides a separate helper in [RnaSecondaryStructure.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/RnaSecondaryStructure.cs) with RNA complement rules and tuple-based output.

## 7. Examples and Related Material

### 7.3 Related Tests, Evidence, or Documents

- Tests: [RepeatFinder_InvertedRepeat_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Analysis/RepeatFinder_InvertedRepeat_Tests.cs)
- Test spec: [REP-INV-001.md](../../../tests/TestSpecs/REP-INV-001.md)
- Related RNA smoke tests: [RnaSecondaryStructureTests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Analysis/RnaSecondaryStructureTests.cs)

## 8. References

1. Wikipedia. 2026. Inverted repeat. Wikipedia. https://en.wikipedia.org/wiki/Inverted_repeat
2. Pearson CE, Zorbas H, Price GB, Zannis-Hadjopoulos M. 1996. Inverted repeats, stem-loops, and cruciforms: significance for initiation of DNA replication. Journal of Cellular Biochemistry. 63(1):1-22.
3. Bissler JJ. 1998. DNA inverted repeats and human disease. Frontiers in Bioscience. 3:d408-d418.
4. Rice P, Longden I, Bleasby A. 2000. EMBOSS: the European Molecular Biology Open Software Suite. Trends in Genetics. 16(6):276-277. `einverted.c` (Durbin R, Thierry-Mieg J, 1993).
5. EMBOSS `palindrome` (Faller M), `emboss/palindrome.c` and application documentation https://emboss.sourceforge.net/apps/cvs/emboss/apps/palindrome.html.
