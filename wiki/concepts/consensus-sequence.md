---
type: concept
title: "Consensus sequence computation (column-wise majority/threshold)"
tags: [assembly, algorithm]
sources:
  - docs/Evidence/ASSEMBLY-CONSENSUS-001-Evidence.md
  - docs/algorithms/Assembly/Consensus_Computation.md
  - docs/Validation/reports/ASSEMBLY-CONSENSUS-001.md
source_commit: ac0a26ca9923867f41c95e2d4a7046a9712ccca6
created: 2026-07-09
updated: 2026-07-10
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: assembly-consensus-001-evidence
      evidence: "Test Unit ID: ASSEMBLY-CONSENSUS-001 ... Consensus Computation (column-wise majority/threshold consensus from aligned reads)"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:multiple-sequence-alignment
      source: assembly-consensus-001-evidence
      evidence: "Both compute a per-column majority-vote consensus; MSA's consensus step and the assembly 'C' of Overlap-Layout-Consensus share the same column-wise voting model"
      confidence: high
      status: current
---

# Consensus sequence computation

Collapsing a set of aligned reads/sequences into a single **consensus** — "the calculated
sequence of most frequent residues found at each position in a sequence alignment"
(Wikipedia). This is the **C** in Overlap-Layout-**Consensus** assembly and the same operation
the [[multiple-sequence-alignment|MSA]] consensus step performs. Validated under test unit
**ASSEMBLY-CONSENSUS-001**; the validation record is [[assembly-consensus-001-evidence]], and
[[test-unit-registry]] tracks the unit. See [[algorithm-validation-evidence]] for the artifact
pattern.

The motif-family sibling [[consensus-from-alignment]] (test unit MOTIF-CONS-001,
`MotifFinder.CreateConsensusFromAlignment`) is the `alternative_to` this rule: same goal, but
**pure most-frequent with no threshold** and a **deterministic alphabetical tie-break**
(A<C<G<T) instead of this unit's plurality cut-off + tie→ambiguous behaviour.

## Decision rule (per column)

Traced verbatim to Biopython `Bio.Align.AlignInfo.SummaryInfo.dumb_consensus`
(`threshold=0.7, ambiguous="X"`, v1.79):

1. **Tally non-gap residues only** — `-` and `.` are skipped; `num_atoms` = count of non-gap
   residues contributing to the column.
2. **Max-count set** — collect every residue sharing the maximum count (`max_atoms`).
3. **Emit a residue only when** exactly one residue holds the max (`len(max_atoms) == 1`)
   **and** its frequency among non-gap residues meets the cut-off
   (`max_size / num_atoms >= threshold`, strict `>=`). Otherwise emit the **ambiguous** symbol.

Consequences:

- **Tie → ambiguous**, never an arbitrary pick of one tied residue.
- **Sub-threshold single majority → ambiguous** (e.g. A,A,T = 2/3 ≈ 0.667 < 0.7 ⇒ ambiguous).
- **All-gap / empty column** (`num_atoms == 0`) → ambiguous; Python short-circuits the `and`
  so `max_size/num_atoms` is never evaluated — **no division by zero**.
- **Consensus length = the full alignment length** (longest read), not the first record's
  length; a read shorter than a column contributes nothing there (**ragged reads** handled).

EMBOSS `cons` corroborates the semantics: its *plurality* cut-off "sets the number of positive
matches below which there is no consensus" — insufficient support ⇒ no committed residue.

## Ambiguity symbol and threshold (repository defaults)

Two presentation defaults deviate from Biopython but are fully parameterized, so every
source-defined value is reachable:

- **Ambiguous symbol `N`, not Biopython's `X`.** For DNA/RNA assembly the IUPAC "any base"
  symbol is `N` (Wikipedia); `X` (protein) or any symbol can be requested via the parameter.
  Presentation-only — it does not change the decision rule.
- **Default threshold `0.5`, not Biopython's `0.7`.** A simple-majority ("plurality") cut-off
  is 0.5; exact Biopython behaviour is reproduced with `threshold: 0.7`. The rule itself
  (strict `>=`, tie→ambiguous, gaps skipped) is source-backed, not assumed.

Both are the only deviations; both are documented assumptions, not invented constants. The
independent two-stage re-validation ([[assembly-consensus-001-report]], Stage A
PASS-WITH-NOTES / Stage B PASS / ✅ CLEAN) re-ran the **Biopython 1.85** `dumb_consensus`
reference on the datasets — 10/10 match — and confirmed both divergences are parameter-reachable,
so no defect and no code change.
