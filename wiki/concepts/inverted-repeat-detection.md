---
type: concept
title: "Inverted repeat detection (reverse-complement arms with a loop)"
tags: [annotation, algorithm]
sources:
  - docs/algorithms/Repeat_Analysis/Inverted_Repeat_Detection.md
source_commit: 8767a4acf4f5df25167ee59d10e9e94725d5c0cf
created: 2026-07-16
updated: 2026-07-16
---

# Inverted repeat detection (reverse-complement arms with a loop)

Finding a sequence segment (the **left arm**) followed downstream by its **reverse
complement** (the **right arm**), optionally separated by an intervening **loop** — the
*reverse-complement* member of the Seqeron `RepeatFinder` repeat family (test unit
**REP-INV-001**, `RepeatFinder.FindInvertedRepeats` in
`src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/RepeatFinder.cs`). Such structures can
extrude as stem-loops / hairpins in single-stranded contexts or cruciforms in duplex DNA,
and occur at replication origins, transposon boundaries, rho-independent terminators,
riboswitches, and tRNA elements; long inverted repeats associate with genomic instability.

```text
5'---TTACG------nnnnnn------CGTAA---3'
   Left arm       Loop       Right arm   (Right = ReverseComplement(Left))
```

This is a *distinct* operation from its neighbouring repeat concepts and should be linked,
not re-derived:

- [[direct-repeat-detection]] — the `RepeatFinder` sibling (`FindDirectRepeats`,
  REP-DIRECT-001). There the downstream copy preserves the original 5'→3' sequence; here it
  is the **reverse complement**. Both are exact, brute-force, hash-set-deduped, and share the
  same validation-symmetry fix pattern.
- [[repetitive-element-detection]] — the repeats/tandem family anchor. Its own inverted-repeat
  section covers the **annotation `RepeatAnalyzer`** IUPACpal `W·G·W̄ᴿ` model (ANNOT-REPEAT-001,
  which permits imperfect arms up to `k` mismatches) and the **RNA** `W·G·W̄ᴿ` stem-region
  model (RNA-INVERT-001, RNA base-pairing). `RepeatFinder.FindInvertedRepeats` is a *third,
  distinct* implementation: DNA-only, **exact** reverse-complement arms, no mismatch/gap
  tolerance, with an explicit `CanFormHairpin` viability flag.

See [[test-unit-registry]] for how the unit is tracked and [[algorithm-validation-evidence]]
for the artifact pattern.

## Core model

For a left arm `L`, loop `X`, and right arm `R`, an inverted repeat satisfies:

```text
R = ReverseComplement(L)
TotalLength = 2 × ArmLength + LoopLength
LoopLength  = RightArmStart − (LeftArmStart + ArmLength)
```

A **palindrome** is the special case with `LoopLength = 0` (even-length reverse-complement
palindrome, e.g. self-complementary `GCGC`); it is not returned under default settings because
the default `minLoopLength` is `3`. The right arm is derived by
`DnaSequence.GetReverseComplementString()` on each candidate left arm.

## API contract

`FindInvertedRepeats` has two overloads over the same exact-match model:

| Parameter | Type | Default | Notes |
|---|---|---|---|
| `sequence` | `DnaSequence` or `string` | required | `DnaSequence` overload throws `ArgumentNullException` on `null`; `string` overload uppercases (`ToUpperInvariant`) and yields nothing for `null`/empty. |
| `minArmLength` | `int` | `4` | Rejected below `2` (`ArgumentOutOfRangeException`). |
| `maxLoopLength` | `int` | `50` | Upper bound when scanning downstream candidate starts. |
| `minLoopLength` | `int` | `3` | Rejected when negative (`ArgumentOutOfRangeException`). |

Each result (`InvertedRepeatResult`) carries `LeftArmStart`, `RightArmStart` (both
**0-based**), `ArmLength`, `LoopLength`, `LeftArm`, `RightArm` (= reverse complement of
`LeftArm`), `Loop`, and `CanFormHairpin`.

**Validation is now symmetric across both overloads.** The raw-string overload was a
documented deviation (REP-INV-001 fuzzing): it did **not** mirror the `DnaSequence` overload's
range checks, so a degenerate `minArmLength = 0` emitted spurious zero-length-arm results. The
fix hoists the `minArmLength < 2` / `minLoopLength < 0` guards into an eager wrapper on the
raw-string overload; both surfaces now throw identically, and a
[[research-grade-limitations|degenerate]] `minArmLength = 0` can no longer produce nonsense
output.

## Algorithm and complexity

1. Normalize case (raw-string overload only, `ToUpperInvariant`).
2. For each left-arm start with room for two minimum-length arms plus the minimum loop,
3. for each arm length from `minArmLength` upward, extract the left arm and compute its reverse
   complement,
4. search downstream starts from `i + armLength + minLoopLength` through the `maxLoopLength`
   bound, and
5. when the downstream substring equals the reverse complement, compute the loop, suppress
   duplicate `(LeftArmStart, RightArmStart, ArmLength)` triples with a hash set, and emit an
   `InvertedRepeatResult`.

Detection is `O(n × A² × L)` time / `O(k + A)` space, where `A` is the maximum tested arm
length and `L` the effective loop-length search bound (each candidate does substring +
reverse-complement work proportional to arm length). No suffix tree is built — unlike its
`FindDirectRepeats` sibling this is a direct brute-force scan.

## Invariants (test oracles)

- `ReverseComplement(LeftArm) = RightArm` for every result (emitted only on match).
- `TotalLength = 2 × ArmLength + LoopLength`; `LoopLength = RightArmStart − (LeftArmStart + ArmLength)`.
- `CanFormHairpin` is `true` exactly when `LoopLength ≥ 3` — the biologically motivated
  sterically-viable-loop threshold (loops < 3 bases are unfavorable; hairpin loops are often
  optimal around 4–8 bases).
- Each `(LeftArmStart, RightArmStart, ArmLength)` tuple is unique (hash-set dedup).

## Edge cases and scope

Empty sequence, a sequence shorter than `2 × minArmLength + minLoopLength`, or no complementary
downstream region → empty enumerable. A homopolymer such as `AAAA` typically returns nothing
(its reverse complement `TTTT` is absent). A self-complementary `GCGC` can match when the
loop-length window permits it. Loop length `0` (pure palindrome) is returned only if the caller
lowers `minLoopLength` below its default of `3`.

**Exact matching only** — no mismatch or gap penalties, so approximate inverted repeats that a
dynamic-programming tool would score are not reported; for that richer model users should rely
on **EMBOSS `einverted`** (DP alignment, mismatch/gap penalties, score threshold). A documented
Framework/Simplified [[research-grade-limitations|limitation]], not an invented constraint. The
algorithm is DNA-specific (DNA complement rules); RNA-specific inverted-repeat / stem discovery
lives in `RnaSecondaryStructure.cs` (RNA-INVERT-001, see [[repetitive-element-detection]]).
