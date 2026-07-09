---
type: concept
title: "IUPAC-degenerate consensus generation (threshold → ambiguity code)"
tags: [motif, algorithm]
sources:
  - docs/Evidence/MOTIF-GENERATE-001-Evidence.md
  - docs/algorithms/Pattern_Matching/IUPAC_Degenerate_Consensus.md
source_commit: d36905351108ae77101357e168b9823952ca6dec
created: 2026-07-09
updated: 2026-07-09
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: motif-generate-001-evidence
      evidence: "Test Unit ID: MOTIF-GENERATE-001 ... Algorithm: IUPAC-Degenerate Consensus Generation (MotifFinder.GenerateConsensus)"
      confidence: high
      status: current
    - predicate: alternative_to
      object: concept:consensus-from-alignment
      source: motif-generate-001-evidence
      evidence: "Both collapse equal-length aligned columns to one consensus in the same MotifFinder area, but differ in decision rule: GenerateConsensus keeps every base above a 25% frequency threshold and emits the NC-IUB IUPAC symbol for that base set; CreateConsensusFromAlignment emits the single most-frequent base with an alphabetical tie-break. The evidence file names GenerateConsensus as the separate degenerate method alongside the plain consensus."
      confidence: high
      status: current
---

# IUPAC-degenerate consensus generation

Collapsing **equal-length aligned DNA sequences** into a single consensus in which each
column that carries appreciable variability is encoded with an **IUPAC ambiguity code**
(R, Y, B, N, …) rather than one plurality base. Where a plain consensus discards the minority
bases, the degenerate consensus **names the set of bases present**, the standard way to report
position variability, mixed probes, and consensus motifs. Seqeron exposes it as
`MotifFinder.GenerateConsensus`. This is the **degenerate** end of the MotifFinder consensus
family — the sibling of the plurality [[consensus-from-alignment]] and the exact
[[known-motif-search]]. Validated under test unit **MOTIF-GENERATE-001**; the validation record
is [[motif-generate-001-evidence]], [[test-unit-registry]] tracks the unit, and
[[algorithm-validation-evidence]] describes the artifact pattern.

## Decision rule (per column)

1. **Count the standard bases** A/C/G/T in column *j* over *n* sequences (case-insensitive;
   non-ACGT characters ignored).
2. **Threshold-filter.** Keep the base set `B_j = { b : count(b) > θ·n }` with **θ = 0.25**
   (named constant `threshold = total * 0.25`), strict `>`. Low-frequency bases are dropped
   *before* encoding. This is the DECIPHER threshold-consensus family (drop characters below a
   frequency threshold, then encode the survivors).
3. **Encode with the IUPAC symbol** for `B_j` — the single NC-IUB 1984 code the (bijective)
   nomenclature assigns to that base set.
4. **No-pass fallback.** If **no** base exceeds θ (e.g. four bases each at exactly 25 %), emit
   the single **most-frequent** base, ties broken **alphabetically** (A<C<G<T).

## The NC-IUB 1984 set→symbol table (authoritative)

The mapping is bijective over the 15 non-empty subsets of {A,C,G,T}; identical across
Cornish-Bowden/NC-IUB 1984, UCSC, and Wikipedia:

| Base set | Symbol | Base set | Symbol |
|---|---|---|---|
| {A} A · {C} C · {G} G · {T} T | (self) | {A,C,G,T} | **N** |
| {A,G} | **R** | {C,T} | **Y** |
| {C,G} | **S** | {A,T} | **W** |
| {G,T} | **K** | {A,C} | **M** |
| {C,G,T} | **B** (not-A) | {A,G,T} | **D** (not-C) |
| {A,C,T} | **H** (not-G) | {A,C,G} | **V** (not-T) |

Any single missing base yields a three-base **not-X** symbol (B/D/H/V), **not** N; N is only
the full four-base set.

## Worked oracles

- **Two-base columns** → R/Y/S/W/K/M; **three-base columns** → B/D/H/V (each base 1 > θ).
- **Threshold filtering precedes encoding**: `A,A,G,G,C` (n=5, θ=1.25) → C (count 1) dropped →
  {A,G} → **R** (not a three-base code).
- **Strict-`>` boundary**: `A,A,A,G` (n=4, θ=1.0) → G at count 1 is **not** > 1.0 → dropped →
  {A} → **A**.
- **No-pass fallback**: `A,C,G,T` (n=4, θ=1.0) → each count 1, none > 1.0 → fallback to
  most-frequent → **A** (alphabetical tie).

## The contrast with the plurality consensus (why this is a separate unit)

Same MotifFinder area, same aligned-column input as [[consensus-from-alignment]] (test unit
MOTIF-CONS-001), but a different decision rule — they are `alternative_to` each other:

| | This unit (MOTIF-GENERATE-001) | [[consensus-from-alignment]] (MOTIF-CONS-001) |
|---|---|---|
| Variable column | keep **all** bases > 25 % → **IUPAC code** (R/Y/B/N…) | single **most-frequent** base only |
| Encoding | NC-IUB ambiguity symbol for the base set | one standard base per column |
| Tie / no-pass | fallback to most-frequent, alphabetical tie | alphabetical tie-break |
| Threshold | **25 %** inclusion cut, strict `>` | **none** (pure most-frequent) |

The [[consensus-sequence]] assembly consensus (ASSEMBLY-CONSENSUS-001, Biopython
`dumb_consensus`) is a third relative: it uses a *plurality* cut-off but emits a single
ambiguous `N`/`X` on ties rather than the specific IUPAC set-symbol this unit produces.

## Scope and siblings

This is the **degenerate consensus** — it summarises an already-aligned instance set into an
ambiguity-coded string. It sits between the plain plurality [[consensus-from-alignment]] and the
richer position-weight-matrix model (`CreatePwm`, the input to motif *scanning*). It is
**generation/encoding**, distinct from the matching direction (scanning a subject for a
degenerate/IUPAC pattern) and from de novo [[overrepresented-kmer-discovery]] (finding
**unknown** motifs in unaligned sequence). See [[known-motif-search]] for the exact-match end of
the motif family, and [[regulatory-element-detection]] for the matching direction of an IUPAC
pattern in practice (its E-box `CANNTG` uses the `N` = any-base ambiguity code from the table above). **No source contradictions**; the only assumptions are the 25 % threshold
value, the strict-`>` boundary, and the no-pass fallback — all documented in
[[motif-generate-001-evidence]].
