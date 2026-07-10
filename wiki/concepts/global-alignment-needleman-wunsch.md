---
type: concept
title: "Global alignment (Needleman–Wunsch)"
tags: [alignment, algorithm]
mcp_tools:
  - global_align
sources:
  - docs/Evidence/ALIGN-GLOBAL-001-Evidence.md
  - docs/Validation/reports/ALIGN-GLOBAL-001.md
  - docs/algorithms/Alignment/Global_Alignment_Needleman_Wunsch.md
source_commit: 4bcdb5e14c3f69dc14f370a3bdd05044a7951032
created: 2026-07-09
updated: 2026-07-10
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: align-global-001-evidence
      evidence: "Test Unit ID: ALIGN-GLOBAL-001 ... Algorithm: Global Alignment (Needleman–Wunsch)"
      confidence: high
      status: current
---

# Global alignment (Needleman–Wunsch)

The canonical dynamic-programming method for the **optimal global alignment** of two
sequences end-to-end — one of the first applications of DP to biological sequences.
Seqeron implements it as `GlobalAlign`; its validation record is
[[align-global-001-evidence|ALIGN-GLOBAL-001]].

## Model

- **Scoring:** a similarity function S(a,b) (match positive, mismatch negative) plus a
  single **linear gap penalty** *d*. The alignment score is the sum of all pairing scores.
- **Border:** F(0,j) = d·j, F(i,0) = d·i.
- **Recurrence:** F(i,j) = max( F(i−1,j−1) + S(Aᵢ,Bⱼ), F(i,j−1) + d, F(i−1,j) + d ).
- **Result:** F(n,m) is the maximum score over all alignments; the alignment is recovered
  by **traceback** — diagonal = match/mismatch, horizontal = gap in one sequence, vertical
  = gap in the other.
- **Complexity:** O(nm) time and O(nm) space for lengths n, m.

Global alignment is most useful when the sequences are similar and of roughly equal length,
and may start or end in gaps. Multiple optimal tracebacks can exist for a single optimal
score (common with repeats / low-complexity regions); Seqeron returns one deterministically.
The [[semi-global-alignment-fitting]] mode reuses this same recurrence but frees the
reference end gaps to fit a short query inside a long reference.

## Seqeron implementation note

`ScoringMatrix.GapExtend` is the linear penalty *d*; `ScoringMatrix.GapOpen` is unused by
`GlobalAlign`. **Affine** gap penalties (gap-open + gap-extend) are a documented *extension*
to the basic model, not part of the standard NW recurrence used here.

Once an alignment is produced, [[alignment-statistics]] summarizes it as percent
identity / similarity / gaps and renders the EMBOSS/BLAST three-line display.

This algorithm is one [[test-unit-registry|test unit]] (ALIGN-GLOBAL-001), validated against
the Wikipedia worked example; see [[algorithm-validation-evidence]] for the evidence-artifact
pattern behind that validation. The 2026-06-24 validation write-up
([[align-global-001-report|Stage A/B report]]) closed it **CLEAN** — Stage A and Stage B both
PASS, 13/13 tests green, one documented (non-defect) empty-`DnaSequence`-input nuance.
