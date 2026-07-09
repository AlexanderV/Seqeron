---
type: concept
title: "Genome rearrangement detection by breakpoints (signed permutations)"
tags: [comparative-genomics, algorithm]
sources:
  - docs/Evidence/COMPGEN-REARR-001-Evidence.md
  - docs/Evidence/COMPGEN-REVERSAL-001-Evidence.md
source_commit: c6c3b0169735a83d79bdf659368b539f39fc6995
created: 2026-07-09
updated: 2026-07-09
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: compgen-rearr-001-evidence
      evidence: "Test Unit ID: COMPGEN-REARR-001 ... Genome Rearrangement Detection by Breakpoints (signed gene-order comparison)"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:test-unit-registry
      source: compgen-reversal-001-evidence
      evidence: "Test Unit ID: COMPGEN-REVERSAL-001 ... Reversal Distance (breakpoint-based lower bound); CalculateReversalDistance returns ceil(b/2) on the unsigned specialization of the same breakpoint criterion"
      confidence: high
      status: current
    - predicate: alternative_to
      object: concept:synteny-and-rearrangement-detection
      source: compgen-rearr-001-evidence
      evidence: "Both detect the rearrangements separating two genomes: this unit counts breakpoints on a signed gene-order permutation (b(α), d≥b/2, Inversion/Transposition per Hunter/Bafna-Pevzner); CHROM-SYNT-001's DetectRearrangements classifies from adjacent synteny-block coordinate signals (Translocation/Inversion/Deletion/Duplication)."
      confidence: high
      status: current
---

# Genome rearrangement detection by breakpoints (signed permutations)

Quantifying and classifying the **rearrangements** that separate two genomes by comparing their
**signed gene-order permutations** — the classical Hannenhalli–Pevzner / Bafna–Pevzner breakpoint
theory. This is a **Comparative-genomics** family unit and the **signed-permutation / breakpoint**
formulation of rearrangement detection; it is deliberately kept distinct from — and modelled as an
`alternative_to` — the block-signal [[synteny-and-rearrangement-detection]] concept, which solves the
same "what rearrangements separate two genomes?" problem from adjacent **synteny-block coordinates**
rather than a permutation. Its sibling COMPGEN units are [[average-nucleotide-identity]],
[[ortholog-detection-reciprocal-best-hits]], [[conserved-gene-clusters-common-intervals]],
[[dot-plot-word-match]], and [[genome-comparison-core-dispensable]]. Validated under test unit
**COMPGEN-REARR-001** (signed breakpoint count + rearrangement classification) and
**COMPGEN-REVERSAL-001** (the unsigned reversal-distance lower bound `⌈b/2⌉`, see below); the
validation records are [[compgen-rearr-001-evidence]] and [[compgen-reversal-001-evidence]],
[[test-unit-registry]] tracks the units, and [[algorithm-validation-evidence]] describes the artifact
pattern.

## Signed permutation model

A genome's gene order over markers `1..n` is a **signed permutation** `α`, where each element is
`±a` — the sign is the gene's strand (`'+'` → positive, `'-'` → negative). The comparison target `β`
is usually the identity `(+1,+2,…,+n)`. The permutation is **extended** by prepending `0` and
appending `n+1`:

```
extended α = (0, α(1), …, α(n), n+1)
```

The endpoints `0` and `n+1` model the **telomeres** (ends of a linear chromosome). A **reversal**
`α[i,j]` reverses a contiguous block *and negates the sign* of every element in it — this is the
**Inversion** operation. A **transposition** `ρ(i,j,k)` instead *moves* a block to a new location
**preserving its internal orientation** (no sign change) — the discriminator between the two events.

## Breakpoints

A consecutive pair `(x, y)` in extended α is a **breakpoint** of α with respect to β when neither
`(x, y)` nor `(−y, −x)` appears as a consecutive pair in extended β. Intuitively an *adjacency*
(non-breakpoint) is a pair that survives in the target under either orientation; a breakpoint is a
disrupted adjacency that some reversal must separate. `b(α)` = the number of breakpoints.

Two equivalent quantities:

```
b(α)          = number of breakpoint pairs in extended α vs β       (breakpoint count)
d_BP(π₁,π₂)   = n − sim(π₁,π₂)                                       (breakpoint distance, Tannier)
```

where `sim` = the number of **common adjacencies** (Tannier et al.; plus half the common telomeric
adjacencies in models with telomeres). The identity target has `b(β) = 0` and `sim = n` common
adjacencies; the fully-disrupted case has `sim = 0`, `d_BP = n`.

**Reversal-distance lower bound.** A single reversal changes the breakpoint count by at most 2
(`b(α) − b(αρ) ≤ 2`), so over `t` reversals `b(α) ≤ 2t`; since the reversal distance `d(α) ≥ t`, this
gives the classic bound:

```
d(α) ≥ b(α) / 2
```

`b(α)` is therefore a computable lower bound on how many inversions must separate the two genomes.

### Unsigned reversal distance `⌈b/2⌉` (`CalculateReversalDistance`, COMPGEN-REVERSAL-001)

The **same** lower bound applied to **unsigned** gene-order indices (strand ignored), returned as an
integer distance estimate rather than a raw breakpoint count. On unsigned indices the breakpoint
criterion is the magnitude specialization of the signed one (Bafna–Pevzner §2 / Hübotter 2020):

```
(π_i, π_{i+1}) is a breakpoint  ⇔  |π_{i+1} − π_i| ≠ 1     (extended with 0 … n+1)
```

`CalculateReversalDistance` counts these unsigned breakpoints `b` and returns the tightest integer
satisfying `d ≥ b/2`, i.e. `⌈b/2⌉ = (b + 1) / 2` (integer arithmetic). It is a **lower bound, not the
exact reversal distance** — the exact value needs the Hannenhalli–Pevzner cycle/hurdle refinement,
**not implemented** here. Documented oracles (validation record [[compgen-reversal-001-evidence]]):

| perm1 (vs identity) | extended | unsigned breakpoints `b` | returned `⌈b/2⌉` |
|---------------------|----------|--------------------------|-------------------|
| `[2,3,1,6,5,4]` (Hunter unsigned) | `[0,2,3,1,6,5,4,7]` | 4 (`0→2`,`3→1`,`1→6`,`4→7`) | **2** |
| `[4,3,2,1]` (fully reversed)      | `[0,4,3,2,1,5]`     | 2 (`0→4`,`1→5`)             | **1** |
| `[1,2,3,4,5]` (identity)          | —                   | 0                           | **0** |

Contracts: single-element / empty input → 0 (no internal adjacency); unequal-length inputs throw
`ArgumentException` (distance is defined only within one marker set); distance is symmetric
(`d(α,β) = d(β,α)`, Hunter).

## Classification (`DetectRearrangements` / `ClassifyRearrangement`)

From a single in-order signed permutation, only the two operations definable that way are classified:

| Signal in the permutation | Event |
|---------------------------|-------|
| Sign-flip + reversed local order of a block | **Inversion** (reversal) |
| Block relocated, internal orientation preserved (same strand) | **Transposition** |

Translocation / Deletion / Insertion / Duplication are **not implemented** by these two methods:
they need chromosome ids (translocation) or gene-set differences (indel/duplication) that a single
in-order permutation cannot express, and no authoritative single-permutation rule assigns them. (The
block-signal [[synteny-and-rearrangement-detection]] `DetectRearrangements` *does* classify those from
adjacent-block coordinates — the two units are complementary formulations.)

## Documented oracles

- **Hunter worked example** — `α=(−2,−3,+1,+6,−5,−4)` vs identity, extended
  `(0,−2,−3,+1,+6,−5,−4,+7)` → breakpoints `(0,−2)(−2,−3)(−3,+1)(+1,+6)(+6,−5)(−4,+7)`, **`b(α)=6`**,
  reversal lower bound `d≥3`. `(−5,−4)` is **not** a breakpoint because `(4,5)∈β` (via `(−y,−x)`).
- **Identity / collinear** — `α=β=(+1,+2,+3,+4,+5)` → **`b(α)=0`** (no events), for any size.
- **Single inversion** — reverse block [2,4] of identity → `α=(+1,−4,−3,−2,+5)`, extended
  `(0,+1,−4,−3,−2,+5,+6)` → **`b(α)=2`** (`(+1,−4)` and `(−2,+5)`; the internal `(−4,−3)`,`(−3,−2)`
  map back to identity adjacencies), `d≥1` (one reversal suffices).

## Edge cases and invariants

- `b(α) ∈ [0, n+1]` (there are `n+1` internal pairs in the extended permutation).
- Identical inputs → `b=0` regardless of size (idempotence of the identity).
- Fewer than 2 mappable orthologs → no events (a permutation of `<2` markers has no internal
  adjacency).
- An ortholog whose target is absent in genome2 is skipped, not a crash.
- Null inputs → `ArgumentNullException` on every public method.

## Assumptions

Three source-backed scoping decisions (detailed in [[compgen-rearr-001-evidence]]): the permutation is
derived from two gene lists plus an `orthologMap` (anchor *generation* delegated to the
[[ortholog-detection-reciprocal-best-hits|ortholog]]/synteny units); strand `'+'/'-'` encodes the
sign; and only Inversion/Transposition are classified (the rest a documented "Not implemented"
limitation). No deviations from the sources are recorded — Hunter (signed-permutation / breakpoint /
reversal bound), Tannier et al. (`d=n−sim` adjacency vocabulary), and Bafna–Pevzner (transposition vs
inversion) are mutually consistent.
