# RNA Inverted Repeats (Potential Stem Regions)

| Field | Value |
|-------|-------|
| Algorithm Group | RNA Secondary Structure |
| Test Unit ID | RNA-INVERT-001 |
| Related Projects | Seqeron.Genomics.Analysis |
| Implementation Status | Complete (delegates to `RepeatFinder`; gapped einverted scoring = `FindInvertedRepeatsScored`) |
| Last Reviewed | 2026-09-25 (B12 review F22) |

## 1. Overview

Finds the **maximal stems** of a nucleotide sequence: a left arm `W`, a loop `G` and a right arm that pairs antiparallel with `W` (`W G W̄ᴿ`) [1][2][3]. They are the sequence signature of potential RNA hairpin stems. This method is the RNA entry point of the one stem engine of the repository (`RepeatFinder.FindInvertedRepeats`, REP-INV-001, EMBOSS `palindrome` model [4]); it only fixes the **RNA pairing rule** (ViennaRNA `CanPair`: A·U, G·C, G·U wobble; `T` read as `U`; ambiguity codes never pair [5]).

## 2. Scientific / Formal Basis

### 2.1 Domain Context

An inverted repeat is "a single stranded sequence of nucleotides followed downstream by its reverse complement"; the loop may have any length, including zero [2]. In RNA the two arms fold back to a double-stranded stem and the spacer becomes the loop.

### 2.2 Core Model

`W G W̄ᴿ` with `|G| ≥ 0`; *within k mismatches* when `δH(W, W̄ᴿ) ≤ k` [1]. For an arm of length `m` at `s1` / `s2`: pair `j` is `(seq[s1+j], seq[s2+m-1-j])`. A **stem** is a valid `(s1, s2, m)` (all pairs pair, `m ≥ minLength`, loop `s2 − (s1+m) ∈ [minSpacing, maxSpacing]`). The **reported set** is the set of stems that are not contained, in both arms, in another valid stem (EMBOSS `palindrome` nesting filter [4]); consequently every reported stem is outward- and inward-maximal, and different stems may share bases (alternative pairings are all kept; there is no unsourced "non-overlapping" selection).

### 2.3 Modeling Assumptions

| ID | Assumption | Consequence if Violated |
|----|------------|--------------------------|
| ASM-01 | Pair set = ViennaRNA `CanPair` (A·U/T, G·C, G·U); `allowWobble:false` = A·U/T, G·C. N/R/Y/… never pair. | G·U-only stems appear/disappear with the option. |
| ASM-02 | Loop length bounded by `[minSpacing, maxSpacing]` (defaults 3 / 100, the repository 3-nt hairpin convention). | Stems outside the window are not reported. |
| ASM-03 | Ungapped. Mismatches (option) are interior only (EMBOSS `palindrome -nummismatches`). Gapped/scored (EMBOSS `einverted`: match +3, mismatch −4, gap −12) is `RepeatFinder.FindInvertedRepeatsScored`. | Bulged stems not reported here. |

### 2.4 Properties and Invariants

| ID | Invariant |
|----|-----------|
| INV-01 | `maxMismatches = 0`: `CanPair(seq[Start1+k], seq[End2-k])` for all `k`; with k > 0 at most k interior mismatches, outer and inner pairs pair. |
| INV-02 | `End1-Start1+1 == Length == End2-Start2+1`. |
| INV-03 | Loop `Start2 - End1 - 1 ∈ [minSpacing, maxSpacing]`. |
| INV-04 | `Length ≥ minLength`. |
| INV-05 | `Start1 ≤ End1 < Start2 ≤ End2`. |
| INV-06 | Completeness / no extras / no duplicates: the result equals the brute-force set of §2.2 (verified on 12 000 random cases, Evidence §F22). |
| INV-07 | Ordered by `Start1`, then `Start2`. |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| sequence | string | required | DNA or RNA, case-insensitive | null/empty → empty result |
| minLength | int | 4 | Minimum arm length (pairs) | `< 1` → empty result |
| minSpacing | int | 3 | Minimum loop length | `< 0` → empty result |
| maxSpacing | int | 100 | Maximum loop length | `< minSpacing` → empty result |
| maxMismatches | int | 0 | Interior mismatches per stem | `< 0` → `ArgumentOutOfRangeException` |
| allowWobble | bool | true | G·U counts as a pair | — |

### 3.2 Output

`IEnumerable<(int Start1, int End1, int Start2, int End2, int Length)>`; 0-based inclusive coordinates, ordered by `Start1`, then `Start2` (materialised list).

### 3.3 Preconditions

Input is upper-cased; `U` is folded to `T` for the shared engine (positions unchanged). The legacy no-throw contract for out-of-range size parameters is kept; only the new `maxMismatches` is validated.

## 4. Algorithm

### 4.1 Steps

1. Normalise (upper case, U→T) and delegate to `RepeatFinder.FindMaximalStemsUpper` (shared with the DNA API, which additionally validates `minArmLength ≥ 2`).
2. The engine enumerates inward-maximal stems by their innermost pair and loop, extends outward, and removes stems contained in both arms of a stem on a shifted diagonal (exact test, `IsContainedInShiftedStem`); with k mismatches it uses the EMBOSS `palindrome` candidate walk.
3. Map `(LeftArmStart, RightArmStart, ArmLength)` to the tuple.

### 4.2 Complexity

O(n · maxSpacing · m) time for k = 0 (m = arm length), O(#stems) space; see REP-INV-001.

## 5. Implementation Notes

### 5.1 Entry points

[RnaSecondaryStructure.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/RnaSecondaryStructure.cs) `FindInvertedRepeats(string, int, int, int, int, bool)` → [RepeatFinder.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/RepeatFinder.cs) `FindMaximalStemsUpper` (internal). MCP: `find_rna_inverted_repeats`.

### 5.2 Conformance, duplication, changes (review 2026-09, F22)

- **Duplication removed:** the former private search (per-`p` longest arm + greedy longest-first selection) was a second, weaker implementation of the stem model in `RepeatFinder`. Differences proven by example: (a) `GGGGAAACCCCAAAGGGG` — two stems `(0,3,7,10,4)`, `(7,10,14,17,4)` share `CCCC`; the greedy rule kept one. (b) per loop-start it kept a single `q` (`UUUUAAAAAAAA` has two maximal stems `(0,3,7,10,4)`, `(0,3,8,11,4)`). (c) the strict IUPAC complement made `N·N`, `S·S`, `W·W` pairs and any non-IUPAC character self-pair, and had no G·U.
- **Canonical:** `RepeatFinder` (B04 finished; REP-INV-001 has its own EMBOSS `palindrome` cross-check, 0 mismatches). Semantics are identical to RepeatFinder's with `allowWobble` (pairs via the same `CanPair`) once `U` is folded to `T`; the RNA wrapper adds only normalisation and the RNA default `allowWobble = true`. Test `FindInvertedRepeats_EqualsRepeatFinderCanonicalEngine` locks equality.
- **Simplifications removed:** "perfect stems only" → `maxMismatches`; "no G·U" → `allowWobble`. The gapped Gotoh-style einverted scoring is not re-implemented (no duplication): use `RepeatFinder.FindInvertedRepeatsScored`.

### 5.3 Search reuse (suffix tree)

Not used: the task is antiparallel reverse-complement extension over a bounded loop window, not exact-substring enumeration; the established tools (`palindrome`, `einverted`) use diagonal scans / local DP.

## 6. Edge Cases and Limitations

| Case | Behavior |
|------|----------|
| null / empty / shorter than `2·minLength + minSpacing` | empty |
| parallel direct repeat (`AAGG…AAGG`) | empty (right arm must pair antiparallel) |
| loop outside window, arm < minLength | not reported |
| N, IUPAC codes, other characters | never pair |
| lower case, `T` | read as upper case / `U` |
| G·U-only stem | reported by default; empty with `allowWobble:false` |

Limitations: loops longer than `maxSpacing` not searched; no bulges/gaps (see `FindInvertedRepeatsScored`); the reported set includes overlapping alternative stems, so it is not a non-overlapping tiling or a fold.

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
// UUACG ... loop AAAAAA ... CGUAA  (CGUAA = reverse complement of UUACG)
var irs = RnaSecondaryStructure.FindInvertedRepeats("UUACGAAAAAACGUAA").ToList();
// irs[0] == (Start1: 0, End1: 4, Start2: 11, End2: 15, Length: 5)
// G·U wobble: FindInvertedRepeats("GGGGAAAAUUUU") == (0,3,8,11,4), (1,4,8,11,4); allowWobble:false → empty
```

**Numerical walk-through:** innermost pair `p = 4` (G), `q = 11` (C); loop length `11-4-1 = 6 ∈ [3,100]`. Extend outward: (G,C)✓ (C,G)✓ (A,U)✓ (U,A)✓ (U,A)✓ → arm length 5; left arm `[0,4]` = UUACG, right arm `[11,15]` = CGUAA = reverse complement of UUACG.

### 7.3 Related Tests, Evidence, or Documents

- Tests: [RnaSecondaryStructure_FindInvertedRepeats_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Analysis/RnaSecondaryStructure_FindInvertedRepeats_Tests.cs) — covers INV-01..INV-05
- Evidence: [RNA-INVERT-001-Evidence.md](../../../docs/Evidence/RNA-INVERT-001-Evidence.md)

## 8. References

1. Alamro H, Alzamel M, Iliopoulos CS, Pissis SP, Watts S. 2021. IUPACpal: efficient identification of inverted repeats in IUPAC-encoded DNA sequences. BMC Bioinformatics 22:51. https://doi.org/10.1186/s12859-021-03983-2
2. Ussery DW, Wassenaar TM, Borini S. 2008. Computing for Comparative Microbial Genomics: Bioinformatics for Microbiologists. Springer. https://doi.org/10.1007/978-1-84800-255-5 (definition via Wikipedia "Inverted repeat": https://en.wikipedia.org/wiki/Inverted_repeat)
3. Rice P, Longden I, Bleasby A. 2000. EMBOSS: The European Molecular Biology Open Software Suite (einverted). Trends in Genetics 16(6):276–277. https://doi.org/10.1016/S0168-9525(00)02024-2 (manual: https://emboss.bioinformatics.nl/cgi-bin/emboss/help/einverted)
4. EMBOSS `palindrome` (Rice et al. 2000), maximal stems / `-nummismatches` / nesting filter; see docs/algorithms/Repeats REP-INV-001.
5. ViennaRNA `pair[][]` / `vrna_ptype` (Lorenz et al. 2011, doi:10.1186/1748-7188-6-26); `RnaSecondaryStructure.CanPair` (RNA-PAIR-001).
