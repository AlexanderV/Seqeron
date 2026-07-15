---
type: source
title: "Evidence: ALIGN-MULTI-001 (Multiple sequence alignment)"
tags: [validation, alignment]
doc_path: docs/Evidence/ALIGN-MULTI-001-Evidence.md
sources:
  - docs/Evidence/ALIGN-MULTI-001-Evidence.md
source_commit: 4fc6a948f66d23331a4fe87fd3f5176d56789e13
ingested: 2026-07-09
created: 2026-07-09
updated: 2026-07-09
---

# Evidence: ALIGN-MULTI-001

The validation-evidence artifact for test unit **ALIGN-MULTI-001** (Multiple Sequence
Alignment, `SequenceAligner.MultipleAlign`). One instance of the templated per-algorithm
[[algorithm-validation-evidence|evidence artifact]] pattern; the algorithm itself and its two
sibling aligners are summarized in [[multiple-sequence-alignment]]. See [[test-unit-registry]]
for how units are tracked.

## What this file records

Unusually rich for an Evidence file: a main document plus **two dated addenda** collected
2026-06-23, one per additional aligner.

- **Online sources** — Wikipedia "Multiple sequence alignment", "Clustal", "Consensus
  sequence" (main doc); Edgar 2004 (MUSCLE), Barton & Sternberg 1987, Wallace et al. 2005,
  and Notredame et al. 2000 (T-Coffee) for the addenda.
- **Test datasets** — MSA edge cases (empty / single / identical / variable-length); a
  hand-derived gap-relocation case for iterative refinement (`CGA, GAGAT, CGC, GAC`,
  SP −8 → −6); and the T-Coffee GARFIELD worked example (primary weight 88 → extended 165).
- **Deviations and assumptions** — several, all flagged as design choices rather than spec
  departures (see below).

## Implementation notes (from the "Implementation Details" / "Assumptions" sections)

- `MultipleAlign` is an **anchor-based star alignment**, not classical progressive alignment:
  center chosen by 4-mer **cosine** similarity (ClustalV uses k-tuple counts), suffix-tree
  exact-match anchors, Needleman-Wunsch between anchors, gap reconciliation — steps 2–4 are
  implementation-specific with no external source.
- **SP score** convention: gap-gap pairs = 0 (standard, not stated in Wikipedia).
- **Consensus** design choices: gaps vote; gap-vs-nucleotide ties go to the nucleotide.
- `MultipleAlignIterative` uses MUSCLE **edge-partition** refinement rather than
  Barton-Sternberg remove-one; the accept rule (non-decreasing SP) is identical and is what
  the tests verify.
- `MultipleAlignConsistency` reuses the repo's existing **UPGMA** guide tree (T-Coffee's paper
  uses NJ — changes merge order only) and sets progressive-DP gap penalties to 0 per the paper.

No contradictions; all deviations are documented as API/structure choices that preserve the
verified properties.
