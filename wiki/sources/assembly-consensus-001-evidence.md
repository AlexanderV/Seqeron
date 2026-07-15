---
type: source
title: "Evidence: ASSEMBLY-CONSENSUS-001 (Consensus Computation)"
tags: [validation, assembly]
doc_path: docs/Evidence/ASSEMBLY-CONSENSUS-001-Evidence.md
sources:
  - docs/Evidence/ASSEMBLY-CONSENSUS-001-Evidence.md
source_commit: 6a70d38db6e348e18fa4be3f6bcd6b2c8ba0c9f8
ingested: 2026-07-09
created: 2026-07-09
updated: 2026-07-09
---

# Evidence: ASSEMBLY-CONSENSUS-001

The validation-evidence artifact for test unit **ASSEMBLY-CONSENSUS-001** (Consensus
Computation — column-wise majority/threshold consensus from aligned reads). One instance of the
templated per-algorithm [[algorithm-validation-evidence|evidence artifact]] pattern; the
algorithm's decision rule and defaults are summarized in [[consensus-sequence]], the anchor for
the assembly CONSENSUS family (and shared with the [[multiple-sequence-alignment|MSA]] consensus
step). See [[test-unit-registry]] for how units are tracked.

## What this file records

- **Online sources** (accessed 2026-06-13):
  - **Biopython `AlignInfo.SummaryInfo.dumb_consensus`** (rank 3, v1.79 tag `biopython-179`,
    quoted verbatim) — the reference decision rule: signature `dumb_consensus(threshold=0.7,
    ambiguous="X", require_multiple=False)`; consensus over the full alignment length; per-column
    tally skipping `-`/`.`; `max_atoms` tie set; emit residue iff `len(max_atoms) == 1 and
    max_size/num_atoms >= threshold`, else ambiguous; tie ⇒ ambiguous; all-gap column ⇒ ambiguous
    with no division by zero (short-circuit `and`).
  - **EMBOSS `cons`** (rank 3) — the "plurality" cut-off concept: support below the cut-off ⇒ no
    consensus residue, corroborating Biopython's threshold semantics.
  - **Wikipedia "Consensus sequence"** (rank 4) — definition ("most frequent residues ... at each
    position") and IUPAC degenerate notation (`N` any base, `Y`, `R`) for uncommitted positions.
- **Datasets** — sub-threshold `A,A,T` (2/3 ≈ 0.667 < 0.7 ⇒ ambiguous); threshold-met `A,A,A,T`
  (3/4 = 0.75 ≥ 0.7 ⇒ `A`); tie `A,G` (`max_atoms`={A,G} ⇒ ambiguous); gaps+ragged
  `A-GT`/`ACGT`/`ACG` ⇒ `ACGT` (gap skipped, absent-read column decided on the 2 present reads).
- **Corner cases / failure modes** — gap chars skipped; residue tie ⇒ ambiguous not arbitrary
  winner; sub-threshold single majority ⇒ ambiguous; all-gap/empty column ⇒ ambiguous (no
  div-by-zero); ragged reads ⇒ consensus spans the longest read.
- **Recommended coverage** — nine MUST tests (unanimous, majority-above-threshold,
  sub-threshold, tie, gap-skipping, ragged/longest-length, all-gap, empty-read-list ⇒ empty
  string, configurable-threshold reproduces Biopython 0.7), plus SHOULD (null ⇒ throw, per repo
  contract) and COULD (custom ambiguous symbol `X`).

## Assumptions (from the artifact)

Two presentation-only defaults deviate from Biopython while remaining parameterized:

1. **Ambiguous symbol default `N`, not Biopython's `X`.** DNA/RNA assembly uses the IUPAC
   "any base" `N` (Wikipedia; existing `ComputeConsensus` convention); `X`/other symbols are
   selectable via parameter. Does not change the decision rule — only the emitted character.
2. **Default threshold `0.5`, not Biopython's documented `0.7`.** A true simple-majority /
   plurality cut-off is 0.5; passing `threshold: 0.7` reproduces exact Biopython behaviour. The
   decision rule itself (strict `>=`, tie→ambiguous, gaps skipped) is source-backed, not assumed.

No contradictions among the three sources — EMBOSS plurality and Wikipedia's definition both
corroborate Biopython's threshold semantics. The only deviations are the two parameterized
default choices above, neither of which alters the source-backed decision rule.
