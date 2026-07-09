---
type: concept
title: "Consensus from a multiple alignment (most-frequent residue, alphabetical tie-break)"
tags: [motif, algorithm]
sources:
  - docs/Evidence/MOTIF-CONS-001-Evidence.md
  - docs/algorithms/Pattern_Matching/Consensus_From_Alignment.md
source_commit: de59ece45cd0b9e5969d6589c1c935e8522d4e4c
created: 2026-07-09
updated: 2026-07-09
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: motif-cons-001-evidence
      evidence: "Test Unit ID: MOTIF-CONS-001 ... Algorithm: Consensus Sequence from a Multiple Alignment (plurality / most-frequent residue)"
      confidence: high
      status: current
    - predicate: alternative_to
      object: concept:consensus-sequence
      source: motif-cons-001-evidence
      evidence: "Both collapse aligned sequences to one consensus but differ in decision rule: MOTIF-CONS-001 is pure most-frequent with a deterministic alphabetical tie-break and no threshold; ASSEMBLY-CONSENSUS-001 uses a plurality threshold (default 0.5) and emits an ambiguous symbol on ties / sub-threshold columns (Biopython dumb_consensus)."
      confidence: high
      status: current
---

# Consensus from a multiple alignment (most-frequent residue)

Collapsing a set of **equal-length aligned sequences** into a single **consensus** — the
column-wise most-frequent (plurality) residue, *"the calculated sequence of most frequent
residues found at each position in a sequence alignment"* (Wikipedia). Seqeron exposes it as
`MotifFinder.CreateConsensusFromAlignment(alignedSequences)`. This is the **motif-family**
consensus: pure most-frequent, deterministic, no threshold — the sibling of the exact
[[known-motif-search]] in the same MotifFinder area. Validated under test unit
**MOTIF-CONS-001**; the validation record is [[motif-cons-001-evidence]], [[test-unit-registry]]
tracks the unit, and [[algorithm-validation-evidence]] describes the artifact pattern.

## Decision rule (per column)

1. **Build the profile.** For each column `j`, count occurrences of each symbol (the 4×n
   **profile matrix** of Rosalind CONS — `P[base, j]` = number of times `base` appears in
   column `j`).
2. **Pick the maximum-count symbol** in the column as the consensus character.
3. **Tie-break — alphabetical (A<C<G<T).** When two symbols share the maximum count, the one
   earlier in the alphabet is chosen (Geneious/LANL rule), making the method **deterministic**.
4. **No threshold.** There is always a most-frequent symbol, so a residue is **always**
   emitted — there is no "no-consensus" position and no ambiguous symbol.

## The contrast with the assembly consensus (why this is a separate unit)

This is **not** the same method as [[consensus-sequence]] (test unit ASSEMBLY-CONSENSUS-001,
Biopython `dumb_consensus`). They `alternative_to` each other — same goal, different rule:

| | This unit (MOTIF-CONS-001) | [[consensus-sequence]] (ASSEMBLY-CONSENSUS-001) |
|---|---|---|
| Threshold | **none** — pure most-frequent | plurality cut-off (default **0.5**), strict `>=` |
| Tie | **alphabetical** (A<C<G<T), deterministic pick | → **ambiguous** symbol (`N`/`X`), never a pick |
| Sub-threshold column | n/a (always emits) | → **ambiguous** symbol |
| Gaps / ragged reads | equal-length precondition | gaps skipped; consensus length = full alignment |
| Area | `MotifFinder` (motif family) | assembly / OLC **C** step + MSA consensus |

EMBOSS `cons` is the *parameterised plurality* reference that ASSEMBLY-CONSENSUS-001 follows
and MOTIF-CONS-001 deliberately does not (its `n`/`x` no-consensus output is out of scope here).

## Worked oracles

- **Rosalind CONS** (rank-5 dataset): the 7×8 sample → profile A=`5 1 0 0 5 5 0 0`,
  C=`0 0 1 4 2 0 6 1`, G=`1 1 6 3 0 1 0 0`, T=`1 5 0 0 0 1 1 6` → consensus **`ATGCAACT`**
  (no ties affect this result).
- **Alphabetical tie-break**: `AT`, `GT` → column 1 A/G tied → **A**; column 2 unanimous T →
  **`AT`**.
- **Identical sequences** → that exact sequence (every column unanimous). **Single sequence**
  → returned unchanged. **Lowercase input** → normalised via `ToUpperInvariant`.

## Contract and preconditions

| Aspect | Behaviour |
|--------|-----------|
| Input | a collection of **equal-length** strings (Rosalind precondition); unequal length → `ArgumentException` |
| Alphabet | ACGT; non-ACGT → `ArgumentException` (consistent with `CreatePwm`) |
| Case | normalised to upper-case before counting |
| Null / empty | null → `ArgumentNullException`; empty collection → empty string |

## Scope and siblings

This is the **plain most-frequent** consensus. The richer motif representations live in the
same MotifFinder area but are **separate** methods (not this unit): an **IUPAC-degenerate**
consensus via `GenerateConsensus` (collapses column variability to an ambiguity code) and a
**position-weight matrix** via `CreatePwm` (probabilistic per-column model — the input to
motif *scanning/discovery*). It shares the aligned-input model with the
[[multiple-sequence-alignment|MSA]] consensus step but applies the strict no-threshold /
alphabetical-tie-break rule rather than MSA's own voting. See [[known-motif-search]] for the
exact-match end of the motif family. **No source contradictions**; the only assumptions are
the alphabetical tie-break and the no-threshold scope, both documented in [[motif-cons-001-evidence]].
