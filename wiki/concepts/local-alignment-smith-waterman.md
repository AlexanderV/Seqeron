---
type: concept
title: "Local alignment (Smith–Waterman)"
tags: [alignment, algorithm]
sources:
  - docs/algorithms/Alignment/Local_Alignment_Smith_Waterman.md
  - docs/Validation/reports/ALIGN-LOCAL-001.md
source_commit: b7e2c1eeb773db02af22541c7e80e6eb7019780c
created: 2026-07-10
updated: 2026-07-10
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: align-local-001-report
      evidence: "Test Unit ID: ALIGN-LOCAL-001 ... Algorithm: Local Alignment (Smith–Waterman)"
      confidence: high
      status: current
    - predicate: alternative_to
      object: concept:global-alignment-needleman-wunsch
      source: align-local-001-report
      evidence: "Same DP recurrence family; SW adds a zero floor and takes the optimum from the max cell (not F(m,n)), so it aligns best local subsequences rather than end-to-end."
      confidence: high
      status: current
---

# Local alignment (Smith–Waterman)

The canonical dynamic-programming method for the **optimal local alignment** of two
sequences: it returns the highest-scoring pair of aligned *subsequences* rather than
forcing both inputs end-to-end. Seqeron implements it as `SequenceAligner.LocalAlign`;
its validation record is [[align-local-001-report|ALIGN-LOCAL-001]].

## Model

- **Scoring:** a substitution score s(a,b) (match positive, mismatch negative) plus a single
  **linear gap penalty** W₁.
- **Border:** first row and first column initialized to **0** (no end-gap penalty).
- **Recurrence (linear gap):** H(i,j) = max( **0**, H(i−1,j−1)+s(aᵢ,bⱼ), H(i−1,j)−W₁, H(i,j−1)−W₁ ).
- **Zero floor** — the distinguishing feature vs Needleman–Wunsch: negative running scores
  reset to 0, so the alignment can *restart* after a low-similarity region.
- **Result:** the optimal score is the **maximum cell anywhere in the matrix** (not the
  bottom-right corner). Traceback starts at that max cell and **stops at the first cell
  with score 0**, yielding 0-based inclusive start/end coordinates for each subsequence.
- **Complexity:** O(mn) time and O(mn) space for lengths m, n.

## Where it sits in the alignment family

Local, [[global-alignment-needleman-wunsch|global (NW)]], and
[[semi-global-alignment-fitting|semi-global (fitting)]] share one DP recurrence and differ
only in **matrix border initialization** and **traceback start**:

| Mode | First row / column | Traceback starts at | Zero floor? |
|------|--------------------|---------------------|-------------|
| Global (NW) | d·j / d·i | bottom-right F(m,n) | no |
| **Local (SW)** — this page | **0 / 0** | **global max cell** | **yes** |
| Fitting | 0 / d·i | max of last row | no |

## Seqeron implementation note

`ScoringMatrix.GapExtend` (stored negative) is the linear penalty W₁; with the test's
`GapExtend = −2` this equals Wikipedia's W₁ = 2. `GapOpen` is **intentionally unused** —
**affine** (gap-open + gap-extend) and arbitrary Wₖ gap-cost variants are *not implemented*;
`LocalAlign` is the linear-gap form only. Traceback resolves ties deterministically
diagonal → up → left, returning one optimal alignment (not the full set of co-optimal paths).

The `string` overload uppercases its inputs and returns `AlignmentResult.Empty` on null/empty;
the `DnaSequence` overload null-guards (`ArgumentNullException`) but does **not** short-circuit
empty input — an empty `DnaSequence` runs the core and yields a `Local` result with empty
aligned strings and `−1` coordinates (a documented API-contract nuance, not a defect).

Once an alignment is produced, [[alignment-statistics]] summarizes it as percent
identity / similarity / gaps.

## Validation

This algorithm is one [[test-unit-registry|test unit]] (ALIGN-LOCAL-001). Its 2026-06-24
validation write-up ([[align-local-001-report|Stage A/B report]]) closed it **✅ CLEAN** —
both stages PASS, 7/7 tests green, no code changes. The validator hand-recomputed the full
Wikipedia DP matrix (`TGTTACGG` / `GGTTGACTA`, +3/−3/−2 → score **13**, alignment
`GTT-AC` / `GTTGAC`) and confirmed the code reproduces it. See
[[algorithm-validation-evidence]] for the evidence-artifact pattern behind the unit.
