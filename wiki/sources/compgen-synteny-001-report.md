---
type: source
title: "Validation report: COMPGEN-SYNTENY-001 (synteny / collinearity block detection — MCScanX DP scoring, ComparativeGenomics.FindSyntenicBlocks)"
tags: [validation, comparative-genomics, governance]
doc_path: docs/Validation/reports/COMPGEN-SYNTENY-001.md
sources:
  - docs/Validation/reports/COMPGEN-SYNTENY-001.md
source_commit: 3d86b2b7c044235f2082bf78748c355fefbb6176
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Validation report: COMPGEN-SYNTENY-001

The two-stage **validation write-up** for test unit **COMPGEN-SYNTENY-001** — whole-genome
syntenic-block detection using the **MCScanX collinearity DP scoring** model (Wang et al. 2012),
validated 2026-06-15. This is the *report* artifact that feeds one row of the [[validation-ledger]];
it records the validator's independent **verdict** on both the algorithm description (Stage A) and
the shipped code (Stage B), and the wider campaign is [[validation-and-testing]]. The DP scoring
model, its constants, oracles and edge cases are synthesized in the shared synteny anchor
[[synteny-and-rearrangement-detection]] (see its *MCScanX collinearity DP model* section);
[[test-unit-registry]] defines the unit.

Distinct from two sibling pages: the pre-implementation evidence artifact
[[compgen-synteny-001-evidence]] (sourced from `docs/Evidence/`, records the sourcing decisions and
assumptions) — this page is the independent two-stage re-validation **verdict**; and the
chromosome-scale [[chrom-synt-001-report]] (a different code path,
`ChromosomeAnalyzer.FindSyntenyBlocks` / `DetectRearrangements`, using the simpler
consecutive-consistent-ordering rule). This unit is the **comparative-genomics-scale** synteny unit:
the same collinearity concept realised with MCScanX's explicit DP scoring constants over an
ortholog/anchor map.

## Verdict

**Stage A: PASS-WITH-NOTES · Stage B: PASS · End state: ✅ CLEAN.** No code defect; no code change
required. Tests were **strengthened (+3)** this session to close real contract/branch gaps. Full
unfiltered suite = **6504 passed, 0 failed** (1 pre-existing benchmark skip); the changed test file
builds warning-free.

## Canonical methods & source under test

- **Canonical:** `ComparativeGenomics.FindSyntenicBlocks(genome1Genes, genome2Genes, orthologMap,
  minAnchors=5, maxGap=25)`; delegate `ComparativeGenomics.VisualizeSynteny(blocks)`.
- **Code path reviewed:** `src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/ComparativeGenomics.cs:84–299`
  — constants (`:88–109`), `FindSyntenicBlocks` (`:128–163`), `ChainCollinearAnchors` greedy chaining
  with direction/gap guards (`:171–224`), `IsReportable` score check (`:230–242`), `BuildBlocks`
  coordinate / inverted assembly (`:244–274`), `VisualizeSynteny` (`:280–299`).

## Stage A — description (algorithm faithfulness)

Confirmed against **MCScanX** (Wang et al. 2012, *Nucleic Acids Res.* 40(7):e49) — retrieved this
session via WebFetch from the Oxford Academic HTML — plus the **wyp1125/MCScanX** reference README
and **Wikipedia "Synteny"**. Verbatim from the paper:

- DP recurrence `Score(v) = max[Score(u) + MatchScore(v) + GapPenalty × NumberofGaps(u,v)]`, u
  preceding v; `MatchScore = 50` per gene pair; `GapPenalty = −1`; `NumberofGaps(u,v)` = max
  intervening genes between u and v, "should be fewer than 25".
- Report rule: "Non-overlapping chains with scores **over 250** (i.e. involving **at least 5
  collinear gene pairs**) are reported"; matches sorted "in both transcriptional directions" →
  forward and inverted blocks.

**Formula check:** the implementation's `score = n × MatchScore + GapPenalty × Σ NumberofGaps` (with
`NumberofGaps = |Δpos2| − 1`) is the **closed form** of the cited DP recurrence for a single monotone
chain; constants match the source exactly. The operative rule `score ≥ 250 AND count ≥ minAnchors(5)`
operationalises the paper's "over 250 (i.e. at least 5 collinear gene pairs)". All edge-case
semantics are sourced (sub-threshold <5 → not reported; gap ≥ cutoff → chain breaks; non-overlapping
reporting; reversed order → inverted block; single consistent direction per block; null args →
`ArgumentNullException`).

**Independent hand cross-check** (from sourced constants): 5 adjacent forward anchors Σgaps=0 →
50×5 − 0 = **250** → reported (1 forward block of 5); 5 reversed → 1 block, `IsInverted=true`;
4 adjacent → 50×4 = **200 < 250** → not reported; 6 anchors with one 1-gene gap (positions
0,1,2,4,5,6) → 50×6 − 1 = **299 ≥ 250** → 1 block of 6.

**Stage A notes (documentation-only → PASS-WITH-NOTES):**

1. **MAX_GAPS = 25 vs 20.** The cited 2012 *paper* says "fewer than 25"; the current MCScanX *tool*
   README defaults to 20. The repo follows its cited primary source (paper, `<25`), so the
   description is faithful to what it cites. Recorded as a documented note (FINDINGS_REGISTER
   F-SYNTENY-001), not a defect.
2. **"over 250" vs exactly 250.** The paper's two phrasings conflict at the boundary (5 pairs = 250,
   not strictly "over"). The paper itself equates "over 250" with "at least 5 collinear gene pairs",
   and the tool's `MATCH_SIZE` default is a count of 5 — so a 5-pair zero-gap block is intended to be
   reported. The repo's `≥ 250` is the correct reading; documented as Assumption #1.

## Stage B — implementation

**Formula realised correctly?** Yes. `IsReportable` computes `chain.Count*50 + (−1)*totalGaps` with
`totalGaps = Σ(|Δpos2|−1)` and requires `≥ 250` **and** `count ≥ minAnchors` — the cited recurrence's
closed form plus the sourced report rule. Gap guard `numberOfGaps < maxGap` (default 25) matches
"fewer than 25"; direction guard `direction == currentDir` enforces single-orientation chains;
`IsInverted = chain[^1].pos2 < chain[0].pos2`; coordinates use min/max parent `Start/End` with
Min/Max clamping. Inputs null-checked; empty genome / `< minAnchors` → empty (no throw).

**Cross-verification (recomputed vs code, tests executed):** 5 forward (score 250) → 1 fwd block,
count 5; 5 reversed → 1 inverted block; 4 adjacent (score 200) → 0 blocks; gap 47 ≥ maxGap → 0 blocks
(two 3-runs); 6 anchors 1 gap (score 299) → 1 block of 6; two 5-runs → 2 non-overlap blocks (Σ=10);
coords g0..g4 → Start1=0/End1=450/Start2=0/End2=450 exact; direction switch mid-run → 1 fwd block of
5 (new test). All ✓. `VisualizeSynteny` renders one line per block (ids, span, orientation, count) —
smoke-confirmed. `CompareGenomes` reuses `FindSyntenicBlocks` (out of scope here).

**Test-quality audit (HARD gate → PASS):** sourced not code-echoes (M1 locks the score-250 boundary —
would fail on strict `>250`; M3 locks the 200<250 negative; S1 locks the score-299 gap case; values
trace to the MCScanX constants retrieved this session). No green-washing — no weakened MUST
assertions; `GreaterThanOrEqualTo` appears only in the C2 invariant *property* test (appropriate); no
skipped/ignored tests. **Added this session** to close real gaps: `S2b` null genome2, `S2c` null
orthologMap (contract validates all three args), and `S4` direction-switch mid-run (exercises the
`direction == currentDir` flush branch / INV-3, previously untested). Honest green: full suite **6504
passed / 0 failed**, 1 pre-existing benchmark skip; the 4 NUnit2007 build warnings are pre-existing in
an unrelated `ApproximateMatcher_EditDistance_Tests.cs`.

## Findings

- **No code defect; no code change (State CLEAN).** Every worked example and cross-check value
  reproduced exactly. Tests strengthened (+3) to cover all contract validation and the
  direction-consistency branch.
- **Documented simplification (not a defect):** greedy single-pass chaining is a simplification of the
  full predecessor-DP; for direction-consistent, gap-bounded, non-interleaved inputs (all test inputs)
  it is provably identical to the DP segmentation. Adversarial interleaved inputs where greedy may
  differ are intentionally avoided and disclosed in the algorithm doc §5.2/§5.3.
- **Follow-up (non-blocking):** if alignment with the *current tool default* is later desired,
  expose/adjust `maxGap` default to 20; today the implementation correctly matches its cited 2012
  paper (`<25`).
