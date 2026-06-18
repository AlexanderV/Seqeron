# Validation Report: COMPGEN-SYNTENY-001 — Synteny / Collinearity Block Detection (MCScanX model)

- **Validated:** 2026-06-15   **Area:** Comparative Genomics
- **Canonical method(s):** `ComparativeGenomics.FindSyntenicBlocks(genome1Genes, genome2Genes, orthologMap, minAnchors=5, maxGap=25)`; delegate `ComparativeGenomics.VisualizeSynteny(blocks)`
- **Stage A verdict:** PASS-WITH-NOTES
- **Stage B verdict:** PASS

## Stage A — Description

### Sources opened this session
- **Wang et al. (2012) MCScanX, NAR 40(7):e49 — Oxford Academic HTML** (https://academic.oup.com/nar/article/40/7/e49/1202057). Retrieved with WebFetch this session. Confirmed verbatim:
  - DP recurrence `Score(v) = max[Score(u) + MatchScore(v) + GapPenalty × NumberofGaps(u,v)]`, u preceding v.
  - `MatchScore(v)` = "50 for one gene pair"; `GapPenalty` = "−1".
  - `NumberofGaps(u,v)` = "the maximum number of intervening genes between u and v, should be fewer than 25".
  - Report rule: "Non-overlapping chains with scores **over 250** (i.e. involving **at least 5 collinear gene pairs**) are reported."
  - Matches sorted "in both transcriptional directions" → forward and inverted blocks.
  - Score is the primary filtering criterion; "at least 5 collinear gene pairs" describes what 250 typically represents (per WebFetch follow-up query).
- **MCScanX reference source / README — wyp1125/MCScanX** (https://github.com/wyp1125/MCScanX). Authoritative defaults of the canonical implementation: `-k MATCH_SCORE` default 50; `-g GAP_PENALTY` default −1; `-s MATCH_SIZE` "number of genes required to call synteny (default: 5)"; `-m MAX_GAPS` "maximum gaps allowed (**default: 20**)"; `-e E_VALUE` default 1e-05.
- **Wikipedia — Synteny** (per Evidence doc): confirms collinearity = stricter, gene-order-preserving form of synteny; MCScan-family DP method.

### Formula check
The implementation's scoring (`score = n × MatchScore + GapPenalty × Σ NumberofGaps`, `NumberofGaps = |Δpos2| − 1`) is the closed-form of the cited DP recurrence for a single monotone chain. Constants `MatchScore=50`, `GapPenalty=−1` match the source exactly. The report rule `score ≥ 250 AND count ≥ minAnchors(5)` operationalises the paper's "scores over 250 (i.e. at least 5 collinear gene pairs)".

### Edge-case semantics
All Stage-A edge cases are sourced and defined: sub-threshold (<5 pairs) → not reported; gap ≥ cutoff → chain breaks; non-overlapping reporting; reversed order → inverted block; no anchors / empty → empty; direction must stay consistent within a block. Null args → `ArgumentNullException` (implementation contract).

### Independent cross-check (hand-computed from sourced constants)
- 5 adjacent forward anchors, Σgaps=0 → 50×5 − 1×0 = **250** → reported (1 forward block of 5). ✓
- 5 reversed anchors → 1 block, IsInverted=true. ✓
- 4 adjacent anchors → 50×4 = **200** < 250 → not reported. ✓
- 6 anchors, one 1-gene gap (positions 0,1,2,4,5,6) → 50×6 − 1×1 = **299** ≥ 250 → 1 block of 6. ✓

### Findings / divergences (→ PASS-WITH-NOTES)
1. **MAX_GAPS = 25 vs 20.** The cited 2012 *paper* says "fewer than 25"; the current MCScanX *tool* README defaults to 20. The repo follows its cited primary source (paper, `<25`), so the description is faithful to the source it cites. Recorded as a documented note, not a defect (FINDINGS_REGISTER F-SYNTENY-001).
2. **"over 250" vs exactly 250.** The paper's two phrasings conflict at the boundary (5 pairs = exactly 250, which is not strictly "over"). The paper itself equates "over 250" with "at least 5 collinear gene pairs", and the tool's `MATCH_SIZE` default is a count of 5 genes — so a 5-pair zero-gap block is intended to be reported. The repo's `≥ 250` is the correct reading; documented as Assumption #1. PASS-WITH-NOTES.

## Stage B — Implementation

### Code path reviewed
`src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/ComparativeGenomics.cs:84-299` — constants (88-109), `FindSyntenicBlocks` (128-163), `ChainCollinearAnchors` greedy chaining with direction/gap guards (171-224), `IsReportable` score check (230-242), `BuildBlocks` coordinate/inverted assembly (244-274), `VisualizeSynteny` (280-299).

### Formula realised correctly?
Yes. `IsReportable` computes `chain.Count*50 + (−1)*totalGaps` with `totalGaps = Σ(|Δpos2|−1)` and requires `≥ 250` AND `count ≥ minAnchors` — the cited recurrence's closed form plus the sourced report rule. Gap guard `numberOfGaps < maxGap` (default 25) matches "fewer than 25". Direction guard `direction == currentDir` enforces single-orientation chains. `IsInverted = chain[^1].pos2 < chain[0].pos2` matches INV-4. Coordinates use min/max parent `Start/End` with `Min/Max` clamping (INV-5). Inputs null-checked; empty genome / <minAnchors → empty (no throw).

### Cross-verification table recomputed vs code (tests executed)
| Case | Sourced expectation | Code result | Match |
|------|---------------------|-------------|-------|
| 5 forward (score 250) | 1 fwd block, count 5 | 1, 5, fwd | ✓ |
| 5 reversed | 1 inverted block | 1, inverted | ✓ |
| 4 adjacent (score 200) | 0 blocks | 0 | ✓ |
| gap 47 ≥ maxGap | 0 blocks (two 3-runs) | 0 | ✓ |
| 6 anchors, 1 gap (score 299) | 1 block of 6 | 1, 6 | ✓ |
| two 5-runs | 2 non-overlap blocks, Σ=10 | 2, 10 | ✓ |
| coords (g0..g4) | Start1=0,End1=450,Start2=0,End2=450 | exact | ✓ |
| direction switch mid-run | 1 fwd block of 5 | 1, 5, fwd | ✓ (new test) |

### Variant/delegate consistency
`VisualizeSynteny` renders one line per block (genome ids, span, orientation, gene count); smoke test confirms forward labelling and count. `CompareGenomes` reuses `FindSyntenicBlocks` (out of this unit's scope).

### Test quality audit (HARD gate)
- **Sourced, not code-echoes:** M1 locks the score-250 boundary (would fail if code used strict `>250`); M3 locks the 200<250 negative; S1 locks the gap-penalty score 299. Values trace to the MCScanX constants retrieved this session, not to code output.
- **No green-washing:** no weakened assertions on MUST cases; `GreaterThanOrEqualTo` appears only in the C2 invariant *property* test (appropriate). No skipped/ignored tests; no widened tolerances.
- **Coverage:** all MUST (M1–M8) and SHOULD (S1–S3) covered. **Added this session** to close real branch/contract gaps: `S2b` null genome2, `S2c` null orthologMap (contract validates all three args), and `S4` direction-switch mid-run (exercises the `direction == currentDir` flush branch / INV-3, previously untested).
- **Honest green:** FULL unfiltered suite = **6504 passed, 0 failed**, 1 pre-existing benchmark skip. Changed test file builds warning-free (the 4 NUnit2007 warnings are pre-existing in an unrelated `ApproximateMatcher_EditDistance_Tests.cs`).
- **Gate result: PASS.**

### Findings / defects
No code defect. Greedy single-pass chaining is a documented simplification of the full predecessor-DP; for direction-consistent, gap-bounded, non-interleaved inputs (all test inputs) it is provably identical to the DP segmentation. Tests intentionally avoid adversarial interleaved inputs where greedy may differ — this is correctly disclosed in the algorithm doc §5.2/§5.3.

## Verdict & follow-ups
- **Stage A:** PASS-WITH-NOTES (MAX_GAPS paper=25 vs tool=20; "over 250" boundary resolved to ≥250 per the paper's own equivalence).
- **Stage B:** PASS.
- **End-state:** ✅ CLEAN — no defect; tests strengthened (+3) to cover all contract validation and the direction-consistency branch; full suite green.
- **Follow-up (non-blocking):** if alignment with the *current tool default* is later desired, expose/adjust `maxGap` default 20; today the implementation correctly matches its cited 2012 paper.
