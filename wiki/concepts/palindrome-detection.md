---
type: concept
title: "DNA palindrome detection (zero-loop reverse-complement, FindPalindromes)"
tags: [annotation, algorithm]
sources:
  - docs/algorithms/Repeat_Analysis/Palindrome_Detection.md
source_commit: 96bfcfe823892dd86e97504488ef753c9b93dc01
created: 2026-07-16
updated: 2026-07-16
---

# DNA palindrome detection (zero-loop reverse-complement, FindPalindromes)

Finding DNA segments that equal their own **reverse complement** — the *biological*
palindrome, distinct from textual single-strand symmetry — the **zero-loop special case** of
an [[inverted-repeat-detection|inverted repeat]] given its own dedicated `RepeatFinder` entry
point (test unit **REP-PALIN-001**, `RepeatFinder.FindPalindromes` in
`src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/RepeatFinder.cs`). Such sites are central to
restriction-enzyme recognition, can extrude as cruciforms/hairpins, and mark chromosome
fragility where long self-complementary tracts occur.

```text
5'-GAATTC-3'
3'-CTTAAG-5'
     ↓ read the reverse strand 5'→3'
5'-GAATTC-3'      GAATTC = ReverseComplement(GAATTC)   (EcoRI site)
```

This is a *dedicated implementation*, not a re-derivation of its neighbours — link, don't copy:

- [[inverted-repeat-detection]] — the general reverse-complement-arm case (`FindInvertedRepeats`,
  REP-INV-001). A palindrome is the **`LoopLength = 0`** member of that family. `FindInvertedRepeats`
  will not report it under its default `minLoopLength = 3`; `FindPalindromes` is the specialised
  even-length, zero-loop scanner that does.
- [[repetitive-element-detection]] — the repeats/tandem family anchor; its inverted-repeat section
  covers the imperfect IUPACpal `W·G·W̄ᴿ` annotation model and the RNA stem model. Palindrome
  detection here is **exact, DNA-only** and even-length-restricted.
- [[restriction-site-detection]] — palindrome discovery is a structural prerequisite for many
  Type II restriction-enzyme recognition sites; this unit finds the self-complementary window,
  not the enzyme catalog or cleavage semantics (use `RestrictionAnalyzer` for that).

See [[test-unit-registry]] for how the unit is tracked and [[algorithm-validation-evidence]] for
the artifact pattern.

## Core model

For a candidate window `S`:

```text
Palindrome(S) ⇔ S = ReverseComplement(S)
```

Basewise for an even-length `S = s₁…s₂ₙ`: `sᵢ = complement(s₂ₙ₊₁₋ᵢ)` for all `i`, with
`A↔T`, `G↔C`. Biological palindromes therefore **require even length** so every base pairs with a
complementary partner across the symmetry axis — the scan steps candidate lengths in increments
of `2`.

## Two entry points, one algorithmic core

The repository exposes **two** `FindPalindromes` implementations over the same exact
reverse-complement / even-length model, differing only in return type and validation surface:

| Entry point | Returns | Validation |
|---|---|---|
| `RepeatFinder.FindPalindromes(DnaSequence, int, int)` | `PalindromeResult` | Validating canonical API. Throws `ArgumentNullException` on null; `ArgumentOutOfRangeException` when `minLength < 4`, `minLength` odd, or `maxLength < minLength`. |
| `RepeatFinder.FindPalindromes(string, int, int)` | `PalindromeResult` | Same length validation; yields nothing for null/empty; uppercases non-empty input before scanning. |
| `GenomicAnalyzer.FindPalindromes(DnaSequence, int, int)` | `PalindromeInfo` | Same scan, but **no** explicit null/range checks — dereferences `sequence.Sequence` directly. |

`minLength` defaults to `4` (must be even, `≥ 4` — the validating API excludes trivial 2-bp
matches), `maxLength` defaults to `12` (must be `≥ minLength`). Each result carries a 0-based
`Position`, the palindromic `Sequence`, and its `Length`.

## Algorithm and complexity

1. Validate `minLength`/`maxLength` (RepeatFinder only); normalize raw-string input to uppercase.
2. For each **even** candidate length from `minLength` through `maxLength`, enumerate every start
   where a full window fits (`start ≤ seq.Length − len`).
3. Extract the window, compute its reverse complement, and test exact equality.
4. Emit a result for every matching window — **overlapping** palindromes of different even lengths
   at the same or neighbouring positions are all reported (each even window is checked
   independently).

Cost `O(n · r · m)` time / `O(m)` space, where `r` is the number of even candidate lengths and `m`
the reverse-complement comparison cost per window. This is a direct brute-force scan — no suffix
tree.

## Invariants (test oracles)

- **INV-01** Every reported palindrome satisfies `Sequence = ReverseComplement(Sequence)` (emitted
  only after exact comparison).
- **INV-02** All reported lengths are **even** (scan steps by `2`).
- **INV-03** Reported positions are within bounds (`start ≤ seq.Length − len`).
- **INV-04** `Length` equals the reported sequence's actual length.

Reference restriction-palindrome examples the detector recovers when present: `GAATTC` (EcoRI,6),
`GGATCC` (BamHI,6), `AAGCTT` (HindIII,6), `TCGA` (TaqI,4), `AGCT` (AluI,4), `GCGGCCGC` (NotI,8),
`CCCGGG` (SmaI,6), `GATATC` (EcoRV,6).

## Edge cases and scope

Empty sequence, a sequence shorter than `minLength`, or no self-complementary window → empty
enumerable. The full sequence is reported if it is palindromic and length constraints permit. Odd
`minLength`, `minLength < 4`, or `maxLength < minLength` throw from `RepeatFinder` (but
`GenomicAnalyzer` may fail differently — the one accepted deviation: equivalent core, divergent
validation surface).

**Exact matching only** — no IUPAC degeneracy or ambiguity codes, so degenerate recognition motifs
are not reported unless the sequence resolves to an exact palindrome. Restriction-enzyme catalog
lookup and cleavage semantics are out of scope; users needing those rely on `RestrictionAnalyzer`
(see [[restriction-site-detection]]). Long palindromic tracts have secondary-structure
consequences, but only local exact windows are reported — no cruciform energetics or digestion
modelling. A documented Framework/Simplified [[research-grade-limitations|limitation]], not an
invented constraint.
