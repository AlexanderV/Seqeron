---
type: source
title: "Validation report: ALIGN-GLOBAL-001 (Needleman–Wunsch global alignment)"
tags: [validation, alignment, governance]
doc_path: docs/Validation/reports/ALIGN-GLOBAL-001.md
sources:
  - docs/Validation/reports/ALIGN-GLOBAL-001.md
source_commit: 4bcdb5e14c3f69dc14f370a3bdd05044a7951032
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Validation report: ALIGN-GLOBAL-001

The two-stage **validation write-up** for test unit **ALIGN-GLOBAL-001** (Global
Alignment, Needleman–Wunsch), validated 2026-06-24. This is the *report* artifact
that feeds one row of the [[validation-ledger]]; it is distinct from the pre-implementation
[[align-global-001-evidence|evidence artifact]] (which traces the spec from external
sources) — this report is the validator's **verdict** on both the algorithm description
and the shipped code. The algorithm itself is summarized in
[[global-alignment-needleman-wunsch]]; the methodology behind the two stages is the
[[validation-protocol]].

## Verdict

**Stage A: PASS · Stage B: PASS · State: CLEAN.** No code changes, no defects logged.
Full `GlobalAlign` test class 13/13 green plus corroborating property and cancellation tests.

## Stage A — description (algorithm faithfulness)

- Canonical method: `SequenceAligner.GlobalAlign(DnaSequence, DnaSequence, ScoringMatrix?)`
  → `GlobalAlignCore` + `Traceback`, with string / cancellation / progress delegates.
- Checked against live Wikipedia "Needleman–Wunsch algorithm": border `F(0,j)=d·j`,
  `F(i,0)=d·i`; recurrence `F(i,j)=max(diag+S, left+d, up+d)`; optimum at `F(n,m)`.
- Worked example `GCATGCG` vs `GATTACA` (+1/−1/−1) → **optimal score 0**, hand-recomputed
  DP matrix reproduced the canonical Wikipedia table.
- Single linear gap `d` maps to `ScoringMatrix.GapExtend`; `GapOpen` is unused by NW
  (affine gaps are a documented extension). No divergences from the primary source.

## Stage B — implementation (code review + cross-check)

- Code path: `SequenceAligner.cs` `GlobalAlignCore` (L220–273) + `Traceback` (L468–538),
  pooled flat DP buffer; the cancellation overload (L96–202) carries its own equivalent
  inline fill+traceback and returns the identical `matrix[m,n]`.
- Cross-verification table recomputed vs code, **13/13 pass** — includes M1 (score 0),
  unequal lengths M3/M4 (−2), all-mismatch M5 (−4), identical M6 (8, no gaps), single
  deletion M7 (2), match-weight M8 (20), stats M10, and symmetry M11.
- Integer DP, no division; scores bounded by `len·max(|Match|,|GapExtend|)`, no overflow
  on stated DNA ranges. Property tests (FsCheck) cover symmetry, equal aligned length,
  aligned length ≥ max(input), identity→max, determinism.

## Findings

- **Minor (documented, not a defect):** the `DnaSequence` empty-input overload runs the
  core and returns a score-0 empty-coordinate alignment, whereas the **string** overload
  returns `AlignmentResult.Empty`. Both are defensible (NW of empty vs empty = 0); framed
  as an API-contract nuance, not an algorithm-spec issue.
- **Optional future tidy (non-blocking):** make the `DnaSequence` empty overload mirror the
  string overload's `AlignmentResult.Empty` guard for symmetry.
