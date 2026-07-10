---
type: source
title: "Validation report: ALIGN-SEMI-001 (semi-global / fitting alignment)"
tags: [validation, alignment, governance]
doc_path: docs/Validation/reports/ALIGN-SEMI-001.md
sources:
  - docs/Validation/reports/ALIGN-SEMI-001.md
source_commit: a8a874422c1e77ba29bb0c230252ed0e96a82660
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Validation report: ALIGN-SEMI-001

The two-stage **validation write-up** for test unit **ALIGN-SEMI-001** (Semi-Global
Alignment, fitting / query-in-reference variant), validated 2026-06-24. This is the
*report* artifact that feeds one row of the [[validation-ledger]]; it records the
validator's **verdict** on both the algorithm description and the shipped code. The
algorithm itself is summarized in [[semi-global-alignment-fitting]]; the two-stage
methodology is the [[validation-protocol]]. Distinct from the pre-implementation
[[align-semi-001-evidence]] artifact.

## Verdict

**Stage A: PASS · Stage B: PASS · State: ✅ CLEAN.** No defects found, no code changes.
All 17 canonical + property tests green (`SequenceAligner_SemiGlobalAlign_Tests` 17/0;
`AlignmentProperties` 22/0).

## Stage A — description (algorithm faithfulness)

- Canonical method: `SequenceAligner.SemiGlobalAlign(DnaSequence, DnaSequence, ScoringMatrix?)`
  → `SemiGlobalAlignCore` + shared `Traceback`.
- Checked against live **Rosalind SIMS** (fitting alignment), **Wikipedia** "Sequence
  alignment" (semi-global / glocal use case) and "Needleman–Wunsch". Confirms the fitting
  init: first row `F(0,j)=0` (free leading reference gaps), first column `F(i,0)=i·d`
  (query fully aligned), NW recurrence with **no zero floor**, optimum = `max_j F(m,j)`
  (free trailing reference gaps).
- **Role-labelling note** (not a defect): Rosalind labels the long string *s*; Seqeron
  labels seq1 = query (short, fully aligned), seq2 = reference (long, free end gaps). The
  transpose is verified consistent — semantics identical.
- Validator **hand-recomputed** four cases: M1 `ATGC`/`AAAATGCAAA` → 4; GAP `ACGT`/`AGT`
  → 2; MIX `AGT`/`AAACTAAA` → 1; MAX `ATG`/`ATGCCC` → 3 (confirms max-of-last-row, not
  bottom-right which is 0). No divergences.

## Stage B — implementation (code review + cross-check)

- Code path in `SequenceAligner.cs`: `SemiGlobalAlignCore` (L423–462) + `Traceback` (L468–538).
- Formula realised: first column `score[i,0]=i*GapExtend`; first row left at default 0;
  recurrence `Max(diag, up, left)` with no zero floor; `maxJ = argmax_j score[m,j]` with
  returned score `score[m, maxJ]` = INV-5; trailing reference suffix appended as gaps in
  seq1 (reversed for post-reverse output), leading reference handled by the `i==0, j>0`
  left-move branch.
- Cross-verification table **11/11 recomputed vs code** (M1 4, S1 8, S2 3, S3 3, S4 20 with
  5/−3/−2, NEG −4, MAX 3, OFS 3, MIX 1, GAP 2, INV 7).
- MCP wrapper `AlignmentTools.SemiGlobalAlign` delegates directly to the canonical method
  (no separate logic). `GapOpen` intentionally unused (linear gap only).

## Findings

- **None.** No defects, no code changes required. State ✅ CLEAN.
- Optional (non-blocking) doc nicety flagged: the unit is precisely the **fitting** member
  of the semi-global family (both reference ends free, query fully aligned), not the
  four-ends-free SMGB variant — already correctly stated in the spec/evidence. No action.
