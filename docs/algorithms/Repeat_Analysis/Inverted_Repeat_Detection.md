# Inverted Repeat Detection

| Field | Value |
|-------|-------|
| Algorithm Group | Repeat Analysis |
| Test Unit ID | REP-INV-001 |
| Related Projects | Seqeron.Genomics.Analysis |
| Implementation Status | Complete — maximal stems = EMBOSS `palindrome` (exact by default; `-nummismatches` / `-maxpallen` options), G·U wobble option, scored gap-tolerant variant = EMBOSS `einverted` |
| Last Reviewed | 2026-09-30 |

## 1. Overview

Inverted repeat detection identifies a sequence segment followed downstream by its reverse complement, optionally separated by a loop [1][2]. Such structures can form stem-loops or hairpins in single-stranded contexts and are closely related to palindromes, which are the special case with loop length zero [1]. The repository implements exact inverted-repeat detection in `RepeatFinder.FindInvertedRepeats`, returning explicit arm coordinates, sequences, loop sequence, and a `CanFormHairpin` flag. Only **maximal** stems are reported, using the reporting rule of EMBOSS `palindrome` [5]: a stem that lies inside another stem in both arms is dropped. By default stems are exact; optional parameters reproduce `palindrome -nummismatches k` (mismatch-tolerant stems) and `-maxpallen` (maximum arm length), and `allowWobble` admits G·U pairs [6][7]. For scored inverted repeats with mismatches **and gaps** the repository provides `RepeatFinder.FindInvertedRepeatsScored`, a reimplementation of EMBOSS `einverted` [4] (Durbin & Thierry-Mieg dynamic programming) that returns the same repeats, scores and alignments as the EMBOSS 6.6.0 binary.

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

### 2.3 Options and the scored variant

- **Mismatches (`maxMismatches = k`, EMBOSS `palindrome -nummismatches k`).** Every stem starts at a *pairing* outer pair and is extended inward pair by pair until the (k+1)-th mismatch, or until the next pair would leave a loop shorter than `minLoopLength` (`palindrome`: until the arms meet); mismatches at the inner end are trimmed (`count -= mismatchAtEnd`), interior mismatches stay and are reported in `Mismatches`. Stems with arm ≥ `minArmLength` and loop ≤ `maxLoopLength` are candidates; a candidate inside another candidate in both arms is dropped (`palindrome_AInB`). Consequently both end pairs of a reported stem pair and it has at most k mismatches.
- **Maximum arm length (`maxArmLength`, EMBOSS `palindrome -maxpallen`).** `palindrome` starts stems only from outer pairs with `rev ≤ current + 2·maxpallen + gaplimit` and prints only stems of arm ≤ `maxpallen`; longer candidates stay in its list and still suppress the stems they contain. A longer stem is therefore **not split** into `maxpallen`-sized pieces — it disappears, and only stems not nested in a longer candidate surface (e.g. `G₁₀AAAC₁₀` gives nothing with `maxpallen` 4–6). This rule is reproduced exactly; `maxArmLength ≥ n/2` is unbounded.
- **G·U wobble (`allowWobble`).** G·U is the principal non-Watson–Crick pair in RNA helices, isosteric enough to substitute for Watson–Crick pairs in stems (Crick 1966 wobble hypothesis [6]; Varani & McClain 2000 [7]). With `allowWobble` pairing uses the canonical `RnaSecondaryStructure.CanPair` (A·U, G·C, G·U; T read as U — the ViennaRNA default pair set); the maximal-stem rule is unchanged. EMBOSS `palindrome`/`einverted` have no wobble option; the semantics are validated against a brute force (all exact stems under the RNA pair set, minus nested ones).
- **Scored variant (`FindInvertedRepeatsScored`, EMBOSS `einverted`).** Local alignment of the sequence with its reverse complement, grown outward from the loop: cell (i, k) pairs right base `i` with left base `i−1−k`, $H(i,k) = \maxig(s(i,i{-}1{-}k) + \max(0, H(i{-}1,k{-}2)),\ \max(H(i,k{-}1), H(i{-}1,k{-}1)) - gig)$ with `s` = match (3) for a Watson–Crick pair of a/c/g/t, mismatch (−4) otherwise, gap `g` = 12, threshold 50 and `maxrepeat` 2000 (einverted.acd, EMBOSS 6.6.0; Durbin's original compile-time value was 4000). The first best cell ≥ threshold of each row is remembered per left start; `maxrepeat` bases later the best end in the window is traced back (same-row gap, previous-row gap, else diagonal) and reported, and the window is cleared — so repeats do not overlap. Results carry 0-based inclusive arm coordinates, score, matches/mismatches/gaps and the three alignment rows exactly as einverted computes them.

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | Default (exact, Watson–Crick): `ReverseComplement(LeftArm) = RightArm`. With `maxMismatches = k`: both end pairs pair and `Mismatches ≤ k` equals the number of non-pairing positions. | Exact stems pair at every position; the mismatch walk stops before the (k+1)-th mismatch and trims inner mismatches. |
| INV-02 | `TotalLength = 2 × ArmLength + LoopLength`. | `InvertedRepeatResult.TotalLength` is defined directly from those fields. |
| INV-03 | `LoopLength = RightArmStart - (LeftArmStart + ArmLength)`. | Loop length is computed from the stored coordinates. |
| INV-04 | `CanFormHairpin` is true exactly when `LoopLength >= 3`. | The constructor sets `CanFormHairpin` from that Boolean test. |
| INV-05 | No reported stem lies inside another reported stem in both arms; each tuple is unique. | Containment filter (EMBOSS `palindrome_AInB`). |
| INV-06 | Revcomp symmetry (Watson–Crick modes, any k / `maxArmLength`): the stems of revcomp(S) are the mirror images `(n−R−A, n−L−A, A)` of the stems of S. With `allowWobble` the mirror symmetry holds for reverse(S) (G·U is not complement-symmetric). | The pairing rule and the candidate/nesting rules are strand-symmetric. |
| INV-07 | Scored: `LeftArmStart ≤ LeftArmEnd < RightArmStart ≤ RightArmEnd`, `Score ≥ threshold`, repeats do not overlap. | einverted reporting (`jstart > lastReported`). |

### 2.5 Comparison with Related Implementations

| Feature | `FindInvertedRepeats` | EMBOSS `palindrome` [5] | `FindInvertedRepeatsScored` = EMBOSS `einverted` [4] |
|---------|------------------------|--------------------------|------------------------------------------------------|
| Matching model | Maximal stems (exact by default) | Same | Dynamic programming, score threshold |
| Mismatches | `maxMismatches` (default 0) | `-nummismatches` | match/mismatch scores |
| Gaps (bulges) | No | No | Linear gap penalty |
| Arm limit | `maxArmLength` (= `-maxpallen` rule) | `-maxpallen` | `maxRepeatLength` (= `-maxrepeat`) |
| Minimum loop | `minLoopLength` (extension) | — | — |
| G·U wobble | `allowWobble` | — | — |
| Overlap handling | Nested (both arms) stems dropped | `-overlap Y` | Non-overlapping best local alignments |
| Ambiguous bases | Only A/C/G/T pair | `n` pairs `n` (all-`n` stems rejected), IUPAC complements | Only a/c/g/t score as match |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `sequence` | `DnaSequence` or `string` | required | DNA sequence to search. | The `DnaSequence` overload throws on `null`; the raw-string overload yields no results for `null` or empty input. |
| `minArmLength` | `int` | `4` | Minimum length of each repeat arm. | Both overloads reject values below `2` with `ArgumentOutOfRangeException`. |
| `maxLoopLength` | `int` | `50` | Maximum loop length between arms. | Must be ≥ `minLoopLength` (`ArgumentOutOfRangeException`); values beyond the sequence length are capped. |
| `minLoopLength` | `int` | `3` | Minimum loop length between arms. | Both overloads reject negative values with `ArgumentOutOfRangeException`. |
| `maxMismatches` | `int` | `0` | Maximum mismatched pairs inside a stem (`palindrome -nummismatches`). | ≥ 0. |
| `maxArmLength` | `int` | `int.MaxValue` | Maximum arm length (`palindrome -maxpallen` rule, §2.3). | ≥ `minArmLength`; ≥ n/2 means unbounded. |
| `allowWobble` | `bool` | `false` | Admit G·U (G·T) pairs via `RnaSecondaryStructure.CanPair` (U read as T). | — |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `LeftArmStart` | `int` | 0-based start position of the left arm. |
| `RightArmStart` | `int` | 0-based start position of the right arm. |
| `ArmLength` | `int` | Length of each arm. |
| `LoopLength` | `int` | Number of intervening nucleotides between the two arms. |
| `LeftArm` | `string` | Left-arm sequence. |
| `RightArm` | `string` | Right-arm sequence (the reverse complement of `LeftArm` for default exact stems). |
| `Loop` | `string` | Intervening sequence between the arms. |
| `CanFormHairpin` | `bool` | `true` when `LoopLength >= 3`. |
| `Mismatches` | `int` (init) | Mismatched positions inside the stem (0 for exact stems). |

**Scored variant** `FindInvertedRepeatsScored(sequence, gapPenalty = 12, threshold = 50, matchScore = 3, mismatchScore = −4, maxRepeatLength = 2000)` returns `ScoredInvertedRepeatResult(LeftArmStart, LeftArmEnd, RightArmStart, RightArmEnd, Score, Matches, Mismatches, Gaps, LeftArmAlignment, MatchLine, RightArmAlignment)` (0-based inclusive coordinates; `RightArmAlignment` reads the right arm from its outer end inward, column-aligned with `LeftArmAlignment`; '-' = gap, '|' = Watson–Crick pair) plus `LeftArmLength`, `RightArmLength`, `LoopLength`, `PercentMatches`. Validation (eager, einverted.acd ranges): `gapPenalty ≥ 0`, `threshold ≥ 0`, `matchScore ≥ 0`, `mismatchScore ≤ 0`, each below einverted's sentinel 1 000 000; `maxRepeatLength ≥ 2`. `null`/empty string → empty; `DnaSequence` overload throws on `null`. Results are in einverted report order.

### 3.3 Preconditions and Validation

`FindInvertedRepeats(DnaSequence, ...)` throws `ArgumentNullException` when `sequence` is `null`. Both overloads throw `ArgumentOutOfRangeException` **eagerly (at the call, not on enumeration)** when `minArmLength < 2`, `minLoopLength < 0`, `maxLoopLength < minLoopLength`, `maxMismatches < 0` or `maxArmLength < minArmLength`. The raw-string overload uppercases non-empty input and yields no results for `null` or empty strings; characters other than A/C/G/T (N, IUPAC codes, U) never pair. Coordinates are 0-based; results are ordered by `LeftArmStart`, then `RightArmStart`.

## 4. Algorithm

### 4.1 High-Level Steps

1. Normalize the sequence to uppercase when the raw-string overload is used.
2. For every innermost pair `(iIn, rIn)` with loop `L = rIn − iIn − 1` in `[minLoopLength, maxLoopLength]`: skip it unless the bases pair and the stem cannot move inward (the next inner pair fails, or it would leave a loop shorter than `minLoopLength`).
3. Extend the stem outward while the bases pair; discard it if the arm is shorter than `minArmLength`. These are EMBOSS `palindrome`'s candidate stems with 0 mismatches.
4. Drop a stem that lies inside another stem in both arms. Only two families of stems can contain it: for a shift `1 ≤ m ≤ L − minLoopLength`, the stem `(i, j, A + m)` (same arm starts) or `(i − m, j − m, A + m)` (same arm ends). Checking these is exact.
5. Emit arms, loop and `CanFormHairpin = L ≥ 3`, sorted by `(LeftArmStart, RightArmStart)`.

**Mismatch / bounded-arm path** (`maxMismatches > 0` or `maxArmLength < n/2`): all pairs of a stem share the anti-diagonal `D = left + right`. Per diagonal the pairs are indexed outward from the innermost admissible pair (loop ≥ `minLoopLength`); walking outward, the mismatch positions are recorded, and every pairing index `u` yields the `palindrome` candidate whose inner end is the first pairing index above the (k+1)-th mismatch below `u` (O(1) per index). The walk stops once that stop lies beyond the largest inner index allowed by `maxLoopLength`, or beyond the `2·maxArmLength + maxLoopLength + 1` start span. Nesting is tested with per-diagonal prefix maxima: a container on diagonal `D ± m` must be ≥ m bases longer and `m ≤ loop − minLoopLength`, so each query is a binary search on at most `min(loop − minLoop, longestArm − arm)` diagonals each side. Candidates longer than `maxArmLength` are then not emitted.

### 4.2 Scored variant (einverted)

A ring of `W = min(maxRepeatLength, n + 2)` DP rows (row `i` in slot `i mod W`), the per-row best score (`localMax`) and the left-start → right-end table (`back`), exactly as `einverted.c`; for `W ≥ n + 2` no slot is ever reused and every row is scored in full, so the cap does not change the output. Each row costs O(W); reporting and trace-back cost O(W) per reported repeat. Degenerate parameters (threshold ≤ match score, gap 0) can make einverted 6.6.0 itself abort with SIGFPE (its percentage divides by matches + mismatches): a trace-back that records no column is consumed without a result, a repeat made only of gap columns is returned with `PercentMatches = 0`, and the scan otherwise continues as einverted's code (verified against an einverted build with only that division guarded).

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| Candidate scan | `O(n · W)` | `O(k)` | `W = maxLoopLength − minLoopLength + 1`, capped by `n`. |
| Stem extension + containment check | `O(A + W · A)` per candidate stem, early exit | — | Output-bound in periodic input; `int.MaxValue` loop bounds terminate (capped by `n`). |
| Mismatch / bounded-arm path | `O(n · (W/2 + k + A))` diagonal walk + `O(C · min(W, A_max) · log)` nesting | `O(C)` | `C` = candidates; 200 kb random, k = 1: 0.46 s; k = 3: 0.86 s (Release). |
| Scored (einverted) | `O(n · W)` | `O(W²)` | `W = min(maxRepeatLength, n + 2)`; 200 kb, W = 2000: 3.2 s; 1 Mb: 16 s (Release). |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [RepeatFinder.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/RepeatFinder.cs)

- `RepeatFinder.FindInvertedRepeats(DnaSequence, int, int, int, int maxMismatches = 0, int maxArmLength = int.MaxValue, bool allowWobble = false)`: Validating overload for `DnaSequence` input.
- `RepeatFinder.FindInvertedRepeats(string, …same options…)`: Uppercases raw string input and yields results for non-empty strings.
- `RepeatFinder.FindInvertedRepeatsScored(string | DnaSequence, int gapPenalty = 12, int threshold = 50, int matchScore = 3, int mismatchScore = −4, int maxRepeatLength = 2000)`: EMBOSS `einverted`.

### 5.2 Current Behavior

Pairing uses the canonical `SequenceExtensions.GetComplementBase` restricted to A/C/G/T (no private complement table). The implementation scans innermost pairs within the loop window, extends each inward-maximal stem outward, and removes stems contained in both arms of another stem. With `minLoopLength = 0` its output equals the EMBOSS 6.6.0 `palindrome` binary (`-nummismatches 0 -overlap Y -gaplimit maxLoopLength -minpallen minArmLength`, unbounded `-maxpallen`) on 1000 random sequences, and it equals an independent Python brute force (enumerate all exact stems, keep the non-contained ones) on 3000 random cases with `minLoopLength` 0–6, N, U and lowercase. The options were added on 2026-09-30: `maxMismatches`/`maxArmLength` equal the `palindrome` binary on 5100 random cases (144 727 stems, 113 676 with mismatches, k 0–6, about half with `-maxpallen` < len/2, 0 mismatches) and a literal transcription of `palindrome.c` generalised to `minLoopLength`/wobble on 6000 cases (0 mismatches); wobble stems equal the all-stems brute force; `FindInvertedRepeatsScored` equals the `einverted` binary on 4700 random cases with planted imperfect/gapped repeats (11 100 repeats, 0 mismatches in coordinates, score, counts and alignment rows) and an einverted build with only the SIGFPE division guarded on 3400 further cases incl. the degenerate ones (8 560 repeats, 0 mismatches) ([REP-INV-001-Evidence.md](../../Evidence/REP-INV-001-Evidence.md)).

### 5.3 Conformance to Theory / Spec

**Implemented:**

- Exact reverse-complement matching between left and right arms [1][2].
- Maximal-stem reporting of EMBOSS `palindrome` [5] incl. `-nummismatches` and `-maxpallen`, plus a minimum loop length.
- G·U wobble pairing option [6][7].
- Score-based inverted-repeat search with mismatches and gaps: EMBOSS `einverted` [4] (`FindInvertedRepeatsScored`).
- Hairpin-viability flag `loopLength >= 3` [1].

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | N / IUPAC codes never pair (EMBOSS `palindrome` lets `n` complement `n` and only rejects all-`n` stems). | Deviation (stricter) | A run such as `GGN…NCC` is not a 3-bp stem. | accepted | Follows EMBOSS `einverted`, which scores a match only for a/c/g/t. |
| 2 | `minLoopLength` (default 3) has no `palindrome` counterpart. | Extension | Inward extension stops at the minimum loop. | accepted | With `minLoopLength = 0` output equals `palindrome`. |
| 3 | Scored variant prints upper-case input characters (einverted: lower case, non-ACGT shown as '-'); coordinates 0-based inclusive. | Output format | None on repeats/scores. | accepted | Normalised in the cross-check. |
| 4 | einverted 6.6.0 aborts (SIGFPE) on some degenerate parameter sets; the port continues (§4.2). | Robustness | Only threshold ≤ match or gap 0. | accepted | Verified against a division-guarded build. |

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

Default pairing is DNA Watson–Crick; RNA stems with G·U pairs need `allowWobble`. `FindInvertedRepeats` allows mismatches but no bulges (as `palindrome`); gapped repeats come from `FindInvertedRepeatsScored` (as `einverted`, which does not score G·U). The scored variant needs O(min(maxRepeatLength, n)²) memory (einverted's DP ring). For RNA-specific inverted-repeat discovery, the repository also provides a separate helper in [RnaSecondaryStructure.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/RnaSecondaryStructure.cs) with tuple-based output.

## 7. Examples and Related Material

### 7.3 Related Tests, Evidence, or Documents

- Tests: [RepeatFinder_InvertedRepeat_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Analysis/RepeatFinder_InvertedRepeat_Tests.cs), [RepeatFinder_InvertedRepeatOptions_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Analysis/RepeatFinder_InvertedRepeatOptions_Tests.cs)
- Test spec: [REP-INV-001.md](../../../tests/TestSpecs/REP-INV-001.md)
- Related RNA smoke tests: [RnaSecondaryStructureTests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Analysis/RnaSecondaryStructureTests.cs)

## 8. References

1. Wikipedia. 2026. Inverted repeat. Wikipedia. https://en.wikipedia.org/wiki/Inverted_repeat
2. Pearson CE, Zorbas H, Price GB, Zannis-Hadjopoulos M. 1996. Inverted repeats, stem-loops, and cruciforms: significance for initiation of DNA replication. Journal of Cellular Biochemistry. 63(1):1-22.
3. Bissler JJ. 1998. DNA inverted repeats and human disease. Frontiers in Bioscience. 3:d408-d418.
4. Rice P, Longden I, Bleasby A. 2000. EMBOSS: the European Molecular Biology Open Software Suite. Trends in Genetics. 16(6):276-277. `einverted.c` (Durbin R, Thierry-Mieg J, 1993).
5. EMBOSS `palindrome` (Faller M), `emboss/palindrome.c` and application documentation https://emboss.sourceforge.net/apps/cvs/emboss/apps/palindrome.html.
6. Crick FHC. 1966. Codon–anticodon pairing: the wobble hypothesis. Journal of Molecular Biology. 19(2):548-555.
7. Varani G, McClain WH. 2000. The G·U wobble base pair: a fundamental building block of RNA structure crucial to RNA function in diverse biological systems. EMBO Reports. 1(1):18-23.
