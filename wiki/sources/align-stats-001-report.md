---
type: source
title: "Validation report: ALIGN-STATS-001 (alignment statistics — identity/similarity/gaps + formatting)"
tags: [validation, alignment, governance]
doc_path: docs/Validation/reports/ALIGN-STATS-001.md
sources:
  - docs/Validation/reports/ALIGN-STATS-001.md
source_commit: 9d11fab0e13fcd2856375d2db84d917b0806f9ec
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Validation report: ALIGN-STATS-001

The two-stage **validation write-up** for test unit **ALIGN-STATS-001** (Pairwise
Alignment Statistics — Identity / Similarity / Gaps — and alignment formatting),
validated 2026-06-15. This is the *report* artifact that feeds one row of the
[[validation-ledger]]; it records the validator's **verdict** on both the algorithm
description and the shipped code. The metrics themselves are summarized in
[[alignment-statistics]]; the two-stage methodology is the [[validation-protocol]].
Distinct from the pre-implementation [[align-stats-001-evidence]] artifact.

## Verdict

**Stage A: PASS-WITH-NOTES · Stage B: PASS-WITH-NOTES · State: ✅ CLEAN.** No algorithm
defect, no code or test change required. Full unfiltered suite **6536 passed / 0 failed**,
build 0 errors. Test-quality gate PASS. The sole "note" on both stages is a documented,
rendering-only display simplification (see NOTE-1) that affects no counted statistic.

## Stage A — description (algorithm faithfulness)

- Canonical methods: `SequenceAligner.CalculateStatistics(AlignmentResult, ScoringMatrix?)`
  and `SequenceAligner.FormatAlignment(AlignmentResult, int, ScoringMatrix?)`.
- Checked against live **EMBOSS `needle`** (rel 6.6) worked example, **EMBOSS AlignFormats**
  (srspair markup legend), and the **pseqsid** reference implementation. Confirms the formula:
  Identity% = M/L×100, Similarity% = (M+Sim⁺)/L×100 (a non-identical column is similar iff its
  substitution score is **positive** — the EMBOSS/BLAST "positives" rule), Gaps% = G/L×100, with
  denominator **L = alignment length including gap columns** and partition **M + X + G = L**.
- Validator **independently re-derived** (Python) the published EMBOSS numbers: Identity 65/149
  = 43.6%, Similarity 90/149 = 60.4%, Gaps 9/149 = 6.0% ✓; and re-classified the hand alignment
  `ACGT-ACGT` / `ACCTAAC-T` → 6 match, 1 mismatch, 2 gap, L=9 (Identity 6/9, Gap 2/9).
- Edge-case semantics defined and consistent with EMBOSS: empty → `AlignmentStatistics.Empty`/`""`;
  null → `ArgumentNullException`; `lineWidth ≤ 0` → `ArgumentOutOfRangeException`; all-gap → Gaps==L,
  Identity 0%; perfect identity → 100/100/0.

## Stage B — implementation (code review + cross-check)

- Code path: `SequenceAligner.cs:566-618` (`CalculateStatistics`), `:630-673` (`FormatAlignment`);
  `AlignmentTypes.cs:8-53` (`ScoringMatrix`, `AlignmentResult.Empty`, `AlignmentStatistics.Empty`).
- Formula realised: single O(L) pass classifies each column (gap if either char is `-`; else
  identical if equal; else mismatch), incrementing the similar counter iff `score.Mismatch > 0`.
  Denominator is `alignmentLength` (includes gaps). For the scalar DNA model every non-identical
  pair scores exactly `Mismatch`, so `score.Mismatch > 0` is **mathematically exact** (not an
  approximation) for "is this column's substitution score positive". SimpleDna / BlastDna /
  HighIdentityDna all have `Mismatch < 0` ⇒ Similarity == Identity, as documented.
- Cross-verification table recomputed vs code via the suite (7/7 rows PASS): EMBOSS 149-col
  (43.6% / 60.4%, L=149, Matches=65); SimpleDna ACGT/ACCT (75%, Sim==Id); PositiveMismatch
  ACGT/ACCT (75% → 100%); hand 9-col (M6 X1 G2 L9, 66.67% / 66.67% / 22.22%); all-gap
  (Gaps==L); FormatAlignment SimpleDna markup `|| |`; FormatAlignment PositiveMismatch `||:|`.
- Test-quality audit (HARD gate) PASS: assertions are sourced (M1 reproduces published 43.6/60.4,
  FormatAlignment markup strings are literal srspair glyphs), exact `.Within(1e-10)`/`1e-3`
  tolerances, no green-washing, no skips; `Similarity Is.EqualTo(Identity)` asserts anchored by
  independent exact-value Identity asserts (not tautologies). All 7 MUST, 6 SHOULD, 2 COULD present.

## Findings

- **No algorithm defect. No code or test change.** State ✅ CLEAN. A BY-DESIGN note is recorded
  in FINDINGS_REGISTER §D; no logged defects.
- **NOTE-1 (rendering-only, documented on both stages):** EMBOSS srspair marks `:` for a
  substitution score **> 1.0** and `.` for a small positive score (0 < score ≤ 1.0). The scalar
  integer Match/Mismatch model has no graded positive tier, so `FormatAlignment` renders any
  positive non-identical column as `:` and never emits `.`. This is a display simplification
  explicitly disclosed in algorithm doc §5.3 / TestSpec Assumption #1; the **counted** statistics
  are unaffected (the "positive substitution score ⇒ similar" count rule is fully correct).
