---
type: concept
title: "IUPAC-degenerate motif matching (scanning a sequence for an ambiguity-coded pattern)"
tags: [motif, algorithm]
sources:
  - docs/algorithms/Pattern_Matching/IUPAC_Degenerate_Matching.md
source_commit: 19070d6ba2f6b3d30d50a67db9183a714db89787
created: 2026-07-15
updated: 2026-07-15
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: iupac-degenerate-matching
      evidence: "Test Unit ID: PAT-IUPAC-001 ... Algorithm Group: Pattern Matching; IUPAC Degenerate Motif Matching (MotifFinder.FindDegenerateMotif / IupacHelper.MatchesIupac)."
      confidence: high
      status: current
    - predicate: alternative_to
      object: concept:known-motif-search
      source: iupac-degenerate-matching
      evidence: "The spec Overview frames this as 'extends exact pattern search by allowing ambiguity codes at each motif position' — same find-all-occurrences-of-a-pattern problem as the exact multi-motif scan, but the per-position equality test is relaxed to allowed-base-set membership (INV-01) instead of character equality."
      confidence: high
      status: current
---

# IUPAC-degenerate motif matching

**IUPAC-degenerate motif matching** is the *matching* (scanning) direction of the IUPAC
ambiguity vocabulary: given a DNA subject and a **motif written with IUPAC ambiguity codes**
(N, R, Y, S, W, K, M, B, D, H, V plus the four bases), it slides the motif across the
sequence and reports every window whose bases all fall inside the base set each motif
position allows. It is the deliberate counterpart to the *generation* direction
[[iupac-degenerate-consensus]] (which **builds** an ambiguity code from an alignment): this
unit **consumes** an ambiguity-coded pattern and tests a subject against it. Seqeron exposes
it as `MotifFinder.FindDegenerateMotif(...)` over the shared `IupacHelper.MatchesIupac(...)`
decision table. Validated under test unit **PAT-IUPAC-001**; [[test-unit-registry]] tracks the
unit and [[algorithm-validation-evidence]] describes the per-unit artifact pattern.

## The matching rule (per window)

A window `S[i..i+m-1]` matches motif `P` iff, for **every** position `j`, the subject base
`S[i+j]` belongs to the base set encoded by the motif code `P[j]`
(`IupacHelper.MatchesIupac(P[j], S[i+j])`). The scan is a brute-force slide over all
`n - m + 1` start positions; a match emits a `MotifMatch` with a fixed `Score = 1.0` — this
mode distinguishes **match / non-match only**, it does not rank partial or probabilistic
agreement. The inner loop exits on the **first** disallowed base (INV-01: every reported hit
satisfies the allowed-base constraint at every position).

## The 15 standard IUPAC DNA codes (the decision table)

`MotifFinder` maps each code to its string of allowed bases; `IupacHelper.MatchesIupac`
exposes the same table as a switch expression. Membership in this dictionary is what pattern
validation checks (INV-02) — the matcher accepts **only** these 15 codes.

| Code | Bases | Mnemonic | | Code | Bases | Mnemonic |
|---|---|---|---|---|---|---|
| A/C/G/T | self | — | | K | G,T | Keto |
| R | A,G | puRine | | M | A,C | aMino |
| Y | C,T | pYrimidine | | B | C,G,T | not-A |
| S | G,C | Strong | | D | A,G,T | not-C |
| W | A,T | Weak | | H | A,C,T | not-G |
| N | A,C,G,T | aNy | | V | A,C,G | not-T |

This is the same NC-IUB 1984 bijective set↔symbol vocabulary tabulated in
[[iupac-degenerate-consensus]]; here it is read in the matching direction (code → accepted
base set) rather than the generation direction (base set → code).

## Contract and invariants

| Aspect | Behaviour |
|---|---|
| `Position` | **0-based** start index of the matching window |
| `MatchedSequence` | the subject substring that satisfied the motif (original window preserved, INV-03) |
| `Pattern` | the **uppercased** motif |
| `Score` | always **`1.0`** (acceptance-set match, no weighting) |
| Motif normalization | `ToUpperInvariant()` applied **before** validation and scanning |
| Cancellation | the cancellation-aware core checks the token every **1000** start positions |
| Complexity | single-window check `O(m)` / full scan `O(n·m)`, `O(1)` auxiliary — brute-force, not sublinear |

## Edge cases and guards

| Case | Behaviour | Why |
|---|---|---|
| Null `DnaSequence` | `ArgumentNullException` | explicit guard in `FindDegenerateMotif` |
| Empty motif | no matches | explicit guard (both standard and cancellation-aware overloads) |
| Pattern longer than sequence | no matches | the outer window loop never runs |
| Lowercase motif | normalized to uppercase | `ToUpperInvariant()` before validation |
| Invalid IUPAC code | `ArgumentException` (`MotifFinder`) / `ArgumentOutOfRangeException` (`IupacHelper`) | invalid symbols are **rejected**, never treated as literals or silently matched |

The one implementation-specific assumption confirmed in source: an invalid motif character is
**rejected up front** rather than falling back to literal matching. No deviation from the
standard IUPAC DNA code table was found.

## Worked motifs (from the spec)

| Motif | Pattern | Role |
|---|---|---|
| E-box | `CANNTG` | bHLH transcription-factor binding site (any central dinucleotide) |
| TATA box | `TATAAA` | core-promoter element (all-exact codes still scan through this engine) |
| Kozak | `GCCGCCRCCATG` | translation-initiation context with the −3 purine (`R`) degeneracy |

Note the E-box `CANNTG` is exactly the single degenerate entry the fixed catalog
[[regulatory-element-detection]] scans through this matcher — that unit is the
**canonical-catalog** application of this engine (a cited constants table), whereas this unit
is the **general** caller-supplied degenerate-pattern scan. (The Kozak example here uses the
degenerate `...RCC...` form; the regulatory catalog instead stores the exact preferred-base
string `GCCGCCACCATGG` — two representations of the same signal, not a contradiction.)

## Scope and siblings

This is degenerate **matching**, one branch of the motif family:

- **Exact** end — [[known-motif-search]] (`GenomicAnalyzer.FindMotif`) and the underlying
  [[exact-pattern-search]] suffix-tree engine: character equality, no ambiguity. This unit is
  the `alternative_to` that relaxes the per-position equality test to allowed-base-set
  membership.
- **Mismatch-tolerant** end — [[approximate-pattern-matching-mismatches]] (Hamming ≤ d):
  tolerates a *bounded number* of substitutions anywhere, whereas degenerate matching tolerates
  *specified* alternatives at *specified* positions (an exact acceptance-set test, zero
  mismatches allowed).
- **Generation** direction — [[iupac-degenerate-consensus]] (`MotifFinder.GenerateConsensus`)
  produces the ambiguity code that this matcher consumes.

Declared out of scope (deliberate simplifications): partial/probabilistic scoring,
PWM-style weighting combined with ambiguity ([[position-weight-matrix]] scanning is a separate
unit), and any index for faster repeated searches. **No source contradictions.**
</content>
</invoke>
