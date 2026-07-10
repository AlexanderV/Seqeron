---
type: source
title: "Validation report: ALIGN-MULTI-001 (multiple sequence alignment — star + progressive + consistency)"
tags: [validation, alignment, governance]
doc_path: docs/Validation/reports/ALIGN-MULTI-001.md
sources:
  - docs/Validation/reports/ALIGN-MULTI-001.md
source_commit: 5035604e5bf429aebca354e08cf9e246b1c0079e
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Validation report: ALIGN-MULTI-001

The two-stage **validation write-up** for test unit **ALIGN-MULTI-001** (Multiple Sequence
Alignment), a full re-validation done 2026-06-24 in fresh context. This is the *report*
artifact that feeds one row of the [[validation-ledger]]; it records the validator's
**verdict** on both the algorithm description and the shipped code. The algorithm family is
summarized in [[multiple-sequence-alignment]]; the two-stage methodology is the
[[validation-protocol]]. Distinct from [[align-multi-001-evidence]] (the pre-impl evidence
artifact) — this is the independent re-validation verdict.

## Verdict

**Stage A: PASS · Stage B: PASS · State: ✅ CLEAN.** No defects, no code changed (only a
temporary probe test, added then removed). Re-confirms **all three** named variants — `MultipleAlign`
(center-star, SP), `MultipleAlignProgressive` (UPGMA guide-tree progressive), and the newest
`MultipleAlignConsistency` (T-Coffee, commit `5b7a9f37`, not covered by the prior 2026-06-17 report);
`MultipleAlignIterative` exercised for additivity via its own addendum. All 96 MultipleAlign-family
tests pass; full suite **18208 passed, 0 failed**, build clean. All variants **NOT LIMITED**.

## Stage A — description (algorithm faithfulness)

Theory checked against primary sources opened this session (not the repo's own assertions):
the **T-Coffee PDF** (Notredame, Higgins & Heringa 2000, J Mol Biol 302:205–217, pp.208–210
read page-by-page), Wikipedia "UPGMA" and "Multiple sequence alignment"/"Clustal", and
Feng & Doolittle (1987).

- **Star:** center = max total similarity → pairwise NW to center → gap reconciliation into one
  MSA coordinate space → majority consensus → SP score (match/mismatch from matrix, residue-gap =
  GapExtend, **gap-gap = 0**). ✔
- **Progressive:** distance d = 1 − fractional identity; **UPGMA** with smallest-distance-first,
  lowest-index tie-break, **proportional size-weighted** averaging `d((A∪B),X)=(|A|·d(A,X)+|B|·d(B,X))/(|A|+|B|)`
  matching Wikipedia symbol-for-symbol; **once-a-gap, always-a-gap** (existing columns copied
  verbatim, merges insert only whole all-gap columns) per Feng-Doolittle. ✔
- **Consistency:** primary library weight = pairwise **percent identity** (×100 integer), global (NW)
  + local (SW) combined by **signal addition**; triplet **extension** = direct primary + Σ over
  intermediates of **min(W₁,W₂)** (uninformative triplets contribute 0, so extension never lowers a
  weight — GARFIELD 88 → 165); progressive DP is Gotoh with **zero gap penalty**, once-a-gap enforced.
  Matches T-Coffee pp.207–210. ✔

**Documented minor divergence (carried, not a defect):** the paper averages group-vs-group column
library scores; the code **sums** them (`ColumnLibraryScore`) — same objective up to a per-merge constant,
so the DP argmax is unchanged; reported `TotalScore` is deliberately the SP score for cross-aligner
comparability, not the consistency objective. Already recorded as an assumption in the evidence doc.

## Stage B — implementation (code review + probe cross-check)

Code path in `src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.cs`:
`MultipleAlign` L702 / `SelectCenterSequence` L791 / `BuildConsensus` L1046 / `ComputeSumOfPairsScore`
L1078; `MultipleAlignProgressive` L1170 / `PairwiseIdentityDistance` L1232 / `BuildProgressiveGuideTree`
L1275 (UPGMA proportional averaging L1331) / `AlignProfiles` L1385; `MultipleAlignConsistency` L1918 /
`BuildExtendedLibrary` L1990 (`extended[key] += min(w1,w2)` L2036–2046) / `PercentIdentity` L2116 /
`ColumnLibraryScore` L2299 / `AlignConsistencyProfiles` L2203.

Independent probe derivations verified against the **live code** (probe since removed):
- Star-vs-progressive discriminator `["AAGAA","AACAA","GGTGG","GGTGG"]`: STAR len 6, SP **−13**;
  PROG len 5, gap-free, SP **−12** — the methods genuinely differ, progressive is better here. ✔
- Progressive `["ACGT","ACGT","AGT"]` → `ACGT/ACGT/A-GT`, SP **8** (hand-derivation matches). ✔
- Consistency library `["ACGT","ACGT","ACGA"]` via internal `GetLibraryWeights`: pair (S0.0,S1.0)
  primary **200** → extended **375** = 200 + min-triplet 175; supported pair (275→375) > less-supported —
  independently reproduces `extended = primary + Σ min-triplet` over the DNA alphabet. ✔

Correctness traps all hold: equal-length rows, degap recovers each input, no all-gap column, byte-identical
on repeat (deterministic), once-a-gap at leaf and profile level, k=2 consistency = pairwise global. All four
aligners reuse the shared cores; the consistency aligner perturbs nothing (TS01 + full green suite).

**Test-quality audit:** star/progressive tests pin exact rows + hand-computed SP (8, 24, the discriminator);
the 12 consistency tests recompute the objective independently (TM04 GARFIELD relation, TM08 via
`GetLibraryWeights`, not echoed from the DP) — strict, deterministic, not tautological.

## Findings

- **None.** No defects, no code changes. State ✅ CLEAN. One documented divergence (summed vs averaged
  group columns; SP as `TotalScore`) carried as a design assumption, not a defect.
