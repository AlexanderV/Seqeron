---
type: source
title: "Evidence: COMPGEN-REVERSAL-001 (Reversal distance, unsigned breakpoint lower bound)"
tags: [validation, comparative-genomics]
doc_path: docs/Evidence/COMPGEN-REVERSAL-001-Evidence.md
sources:
  - docs/Evidence/COMPGEN-REVERSAL-001-Evidence.md
source_commit: c6c3b0169735a83d79bdf659368b539f39fc6995
ingested: 2026-07-09
created: 2026-07-09
updated: 2026-07-09
---

# Evidence: COMPGEN-REVERSAL-001

The validation-evidence artifact for test unit **COMPGEN-REVERSAL-001** — **reversal (inversion)
distance** returned as the **breakpoint-based integer lower bound** `⌈b/2⌉`. This is a
**Comparative-genomics** family Evidence file and one instance of the templated per-algorithm
[[algorithm-validation-evidence|evidence artifact]] pattern.

It is **not a distinct algorithm** from the signed-permutation breakpoint theory already documented:
`CalculateReversalDistance` computes the *same* reversal-distance lower bound `d ≥ b(π)/2` as
[[genome-rearrangement-breakpoint-distance]] (test unit COMPGEN-REARR-001), but on **unsigned**
gene-order indices and returning the tightest integer `⌈b/2⌉` rather than a raw breakpoint count. The
theory, invariants, worked oracles, and the unsigned-vs-signed distinction are summarized in that
shared concept. Its sibling COMPGEN units are [[average-nucleotide-identity]],
[[ortholog-detection-reciprocal-best-hits]], [[conserved-gene-clusters-common-intervals]],
[[dot-plot-word-match]], and [[genome-comparison-core-dispensable]]. See [[test-unit-registry]] for
how units are tracked.

## What this file records

- **Online sources:**
  - **Bafna & Pevzner (1998)** ("Sorting by Transpositions", *SIAM J. Discrete Math.* 11(2):224–240;
    authority rank 1) — the **breakpoint definition** "`(π_i, π_{i+1})` is a breakpoint if
    `π_{i+1} ≠ π_i + 1`" on the extended permutation, the fact that the **identity is the only
    permutation with 0 breakpoints** (so sorting = driving `b(π)` to 0), and the
    operation/breakpoint lower-bound construction (a transposition removes ≤ 3 breakpoints ⇒
    `d ≥ #breakpoints/3`; the reversal analogue removes ≤ 2 ⇒ `d ≥ #breakpoints/2`).
  - **Hunter College CompBio, Lecture 16** (authority rank 2, restating Bafna–Pevzner /
    Hannenhalli–Pevzner) — the **signed** breakpoint criterion, the verbatim "a reversal `ρ=[i,j]`
    can reduce the number of breakpoints by at most two", the **lower-bound derivation**
    `b(α) ≤ 2t ⇒ d(α) ≥ b(α)/2` ("this lower bound is not very tight"), and **symmetry**
    `d_β(α) = d_α(β)`. Its signed worked example `α=(−2,−3,+1,+6,−5,−4)` → `b=6`.
  - **Hübotter (2020)** survey (authority rank 4, corroboration only) — the **unsigned** breakpoint
    "`|π_{i+1} − π_i| ≠ 1`" and unsigned bound `b(π)/2 ≤ d_r(π)`, the specialization of
    Bafna–Pevzner §2 this unit implements.
  - **Bergeron, Mixtacki & Stoye (2009)** ("The Inversion Distance Problem", textbook chapter,
    authority rank 1) — **adjacency vs breakpoint** ("P1 has two adjacencies … all other points are
    breakpoints") and the attribution that the signed breakpoint graph is Bafna–Pevzner's extension
    of the common unsigned one.

- **Documented corner cases / failure modes:** identity → 0 breakpoints → distance 0 (only
  breakpoint-free permutation); the bound is a **lower bound, not exact** ("a permutation with few
  breakpoints may be more distant … than one with more") — a reversal removing two breakpoints need
  not make progress, so the exact distance needs the Hannenhalli–Pevzner cycle/hurdle refinement,
  **not implemented here**.

- **Datasets (documented oracles), unsigned model:**
  - *Hunter example, unsigned specialization* — perm1 `[2,3,1,6,5,4]` vs identity `[1,2,3,4,5,6]`,
    extended `[0,2,3,1,6,5,4,7]` → breakpoints at `0→2`, `3→1`, `1→6`, `4→7` ⇒ **`b = 4`**,
    `⌈b/2⌉ = ` **2**.
  - *Fully reversed* — perm1 `[4,3,2,1]` vs `[1,2,3,4]`, extended `[0,4,3,2,1,5]` → breakpoints
    `0→4`, `1→5` ⇒ **`b = 2`**, `⌈b/2⌉ = ` **1** (one full reversal sorts it).
  - *Identity* — `[1,2,3,4,5] = [1,2,3,4,5]` → `b = 0`, distance **0**.

- **Test-coverage recommendations:** MUST — identity → 0; Hunter unsigned example → 2; fully reversed
  `[4,3,2,1]` → 1; lower-bound property `0 ≤ d ≤ true distance`. SHOULD — single-element / empty → 0
  (no internal adjacency); unequal lengths throw `ArgumentException`. COULD — symmetry
  `d(α,β) = d(β,α)`.

## Deviations and assumptions

**Two ASSUMPTIONS**, both source-backed, neither an open correctness gap:

- **Integer rounding of the real-valued bound.** The theorem gives `d ≥ b/2` (real). The
  implementation returns the smallest integer satisfying it, `⌈b/2⌉ = (b + 1) / 2` in integer
  arithmetic — the canonical integer breakpoint lower bound, the tightest integer the theorem
  guarantees, not an invented value.
- **Unequal-length inputs throw.** Reversal distance is defined only between two permutations of the
  same marker set (same `n`); the implementation throws `ArgumentException` when lengths differ — the
  only well-defined behaviour, not separately specified by the sources.

No source contradictions — Bafna–Pevzner §2 (breakpoint definition + lower-bound construction),
Hunter (reversal removes ≤ 2, `d ≥ b/2`, symmetry), Hübotter (unsigned specialization), and
Bergeron–Mixtacki–Stoye (adjacency/breakpoint) are mutually consistent; the unsigned statement is the
`|Δ|`-magnitude specialization of the signed criterion.
