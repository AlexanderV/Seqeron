---
type: concept
title: "Multiple sequence alignment (MSA)"
tags: [alignment, algorithm]
sources:
  - docs/Evidence/ALIGN-MULTI-001-Evidence.md
source_commit: 4fc6a948f66d23331a4fe87fd3f5176d56789e13
created: 2026-07-09
updated: 2026-07-09
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: align-multi-001-evidence
      evidence: "Test Unit ID: ALIGN-MULTI-001 ... Multiple Sequence Alignment ... SequenceAligner.MultipleAlign"
      confidence: high
      status: current
    - predicate: alternative_to
      object: concept:global-alignment-needleman-wunsch
      source: align-multi-001-evidence
      evidence: "MSA aligns three or more sequences (pairwise alignment handles 2); the star aligner uses Needleman-Wunsch only to fill gaps between exact-match anchors against a chosen center sequence"
      confidence: high
      status: current
---

# Multiple sequence alignment (MSA)

Aligning **three or more** sequences at once so that every column expresses a shared
position across all of them. Finding the *optimal* MSA is NP-complete, so all practical
methods are heuristics. Seqeron exposes three, all validated under test unit
**ALIGN-MULTI-001**; the validation record is [[align-multi-001-evidence]], and pairwise
[[global-alignment-needleman-wunsch|global alignment]] is the two-sequence sibling.

## Shared model

- **Objective — sum-of-pairs (SP):** score each column by summing over all C(k,2)
  sequence pairs (match / mismatch / residue-gap); gap-gap pairs score 0 (standard
  convention, an implementation choice not stated verbatim in the sources).
- **Consensus:** majority vote of the most-frequent residue per column. Seqeron's design
  choices: gaps participate in the vote, and a gap-vs-nucleotide tie resolves to the
  nucleotide.
- **Invariants:** all rows equal length (L ≥ max input length); removing gaps recovers each
  input exactly (reversibility); no column is all-gaps; consensus drawn from {A,C,G,T,-}.

## Three Seqeron implementations

1. **`MultipleAlign` — anchor-based star alignment.** Picks a **center** sequence by highest
   4-mer cosine similarity to all others, builds a suffix tree on it, aligns every other
   sequence to the center using exact-match anchors with Needleman-Wunsch filling the gaps
   between anchors, then reconciles gap columns into one coordinate space. O(k²·m) time.
   The cosine-similarity center pick and suffix-tree anchoring are implementation-specific
   (ClustalV uses k-tuple match counts; there is no phylogenetic guide tree here).
2. **`MultipleAlignIterative` — tree-dependent restricted partitioning** (MUSCLE Stage 3,
   Edgar 2004). Splits the alignment along each guide-tree edge (edges visited by decreasing
   distance from the root), re-aligns the two sub-profiles, and **keeps the result only if the
   SP score does not decrease** — repeating to convergence or a cap. Removes the progressive
   "once a gap, always a gap" limitation. Barton-Sternberg "remove-one-sequence" is the
   equivalent classical scheme; the edge-partition variant is an API/structure choice with an
   identical accept rule.
3. **`MultipleAlignConsistency` — T-Coffee consistency** (Notredame et al. 2000). Optimises a
   **different objective class**: a primary library of weighted residue pairs (weight = percent
   identity of the pairwise alignment) is **extended** by triplets — each intermediate sequence
   contributes min(W₁,W₂), summed across all triplets — then progressively aligned with the
   extended-library weights as substitution scores and **gap penalties set to 0**. Uninformative
   triplets and never-aligned pairs contribute 0, so extension never lowers a weight (worked
   GARFIELD example: primary 88 → extended 165 = 88 + 77).

## Known limitations

Progressive methods propagate early alignment errors; the star result depends on which
sequence is chosen as center; accuracy degrades with evolutionary distance; and being
heuristics none guarantees the optimum (MSA is NP-complete). See
[[algorithm-validation-evidence]] for the evidence-artifact pattern behind this unit.
