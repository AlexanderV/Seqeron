---
type: source
title: "Validation report: CRISPR-OFF-001 (CRISPR off-target scoring — MIT/Hsu-Zhang 2013 hit + aggregate specificity and CFD Doench-2016 off-target score, CrisprDesigner.CalculateMitHitScore / CalculateMitSpecificityScore / CalculateCfdScore)"
tags: [validation, primer, governance]
doc_path: docs/Validation/reports/CRISPR-OFF-001.md
sources:
  - docs/Validation/reports/CRISPR-OFF-001.md
source_commit: c763b50cd147a56f659298a8626d2ede865297f1
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Validation report: CRISPR-OFF-001

The two-stage **validation write-up** for test unit **CRISPR-OFF-001** (CRISPR **off-target scoring** —
scanning a genome for near-matches to a gRNA and scoring off-target risk), area **MolTools**. This is
the *report* artifact that feeds one row of the [[validation-ledger]]; it records the validator's
independent **verdict** on the algorithm description (Stage A) and the shipped code (Stage B), inside
the wider [[validation-and-testing]] campaign. It is the **off-target sibling** of the on-target report
[[crispr-guide-001-report]]; both scoring surfaces live on the same `CrisprDesigner` class and are
synthesized together in the concept [[crispr-guide-rna-design]] (the MolTools/CrisprDesigner anchor).
[[test-unit-registry]] defines the unit. Source under test:
`src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/CrisprDesigner.cs`.

The report is in **two parts**, each an independent re-validation from a fresh validator context: (1) the
**MIT/Hsu 2013** single-hit + aggregate specificity model (added commit `129c2ca`; re-confirmed
2026-06-24), and (2) the **CFD (Doench 2016)** off-target score (implemented this session — the last
off-target residual of finding **C7**).

## Verdict

**Both parts: Stage A PASS · Stage B PASS · End state ✅ CLEAN.** No defect found; no production-code
change required. One in-tree **orientation-guard** test (`CalculateMitHitScore_WeightOrientation_
PamProximalIsHighPenalty`) was verified source-correct and committed. Off-target + CFD suite ran
**71/71 passed**; the full unfiltered suite **6812 passed / 0 failed**, build **0 warnings / 0 errors**.

## Canonical methods & source under test

In `CrisprDesigner.cs`:

- `CalculateMitHitScore(string guide20, string offTarget20)` (`:672–717`) — the **single-hit** MIT/Hsu
  score in [0,100] for one 20-nt off-target against a 20-nt guide.
- `CalculateMitSpecificityScore(IEnumerable<double>)` (`:688–693`) — the **aggregate** guide specificity
  from a set of per-hit scores; plus the genome-scanning overload `CalculateMitSpecificityScore(string,
  DnaSequence, int, CrisprSystemType)` (`:705–721`) that reuses `FindOffTargets` (exact on-target
  excluded) and combines per-hit scores.
- `CalculateCfdScore(string sgRna20, string offTarget20, string offTargetPam)` — the **CFD** off-target
  score in [0,1]; matrices `CfdMismatchScores` (12 keys × 20 positions = 240 entries), `CfdPamScores`
  (16 keys), helper `CfdComplement`.
- The pre-existing honest-heuristic `FindOffTargets` / `CalculateSpecificityScore` /
  `CalculateOffTargetScore` are **unchanged**; this validation is scoped to the MIT/Hsu + CFD additions.
- Tests: 11 MIT `[Test]` methods + the orientation guard; new fixture `CrisprDesigner_Cfd_Tests` (32
  methods) → 71 off-target+CFD tests total.

## Stage A — description (algorithm faithfulness)

### MIT / Hsu-Zhang 2013 model
Grounded against **Hsu et al. 2013** (*Nat Biotechnol* 31:827, PMID 23873081, doi:10.1038/nbt.2647 — the
MIT "Scores of single hits" scheme) and the canonical reference **CRISPOR `crispor.py`** (`calcHitScore`
/ `calcMitGuideScore`, `hitScoreM` weight list), re-downloaded raw this session and read directly.

- **W mismatch-penalty vector — 20/20 byte-identical** to the reference `hitScoreM`. Code
  `MitHitScoreWeights` = `0,0,0.014,0,0,0.395,0.317,0,0.389,0.079,0.445,0.508,0.613,0.851,0.732,0.828,
  0.615,0.804,0.685,0.583` (index 0..19). No transcription error.
- **Orientation pinned from the source** (not the code): the loop indexes `hitScoreM[pos]` over the
  protospacer 5'→3', so **index 0 = PAM-distal 5' (low/zero weight), index 19 = PAM-proximal / seed**;
  the high-weight positions (max `W[13]=0.851`) form the PAM-proximal seed — the biologically expected
  least-tolerated end (Hsu 2013 / Wikipedia off-target editing).
- **Formula (verbatim from `calcHitScore`):** `score1 = Π over mismatched positions i of (1 − W[i])`;
  `score2 = 1` if `nmm < 2` else `1/(((19 − meanInterMismatchDist)/19)·4 + 1)` (`maxDist = 19`);
  `score3 = 1` if `nmm == 0` else `1/nmm²`; `hitScore = score1·score2·score3·100`. Aggregate
  (`calcMitGuideScore`): `100/(100 + Σ hitScores)·100` — CRISPOR rounds to int; the C# method returns
  the unrounded value (documented, benign difference).
- **Independent cross-check** (Python re-derived from the source formula, **not** the C# array) for guide
  `GACGCATAAAGATGAGACGC`: perfect→**100.0**; mm@0 (W=0)→**100.0**; mm@5→**60.5**; mm@13 (W=0.851)→**14.9**;
  mm@19→**41.7**; two-mm {5,15}→**0.8987**; aggregate {60.5}→**62.30529595**; aggregate {60.5,41.7}→
  **49.45598417** — all match the test assertions exactly.
- **Orientation counterfactual** (reversed W): mm@13→68.3, mm@0→41.7 — both differ from the correct
  14.9 / 100.0, so the committed orientation-guard test is a genuine reversal-detector, not a code-echo.

### CFD (Doench 2016) off-target score
CFD's mismatch + PAM matrices originate from **Doench et al. 2016** (*Nat Biotechnol* 34:184, PMID
26780180) and ship as binary pickles. The faithfulness boundary (verbatim numbers cross-checked across
**two independent sources**) was cleared this session:

- **Source 1** CRISPOR `CFD_Scoring/mismatch_score.pkl` + `pam_scores.pkl` + `cfd-score-calculator.py`
  (the canonical John Doench calculator); **Source 2** `bm2-lab/iGWOS` `CFD/*.pkl` + `otscore.py`
  (independent repo, documented `calcCfdScore` doctest oracles).
- **Cross-check:** both pickles decoded to text and diffed element-by-element — **mismatch matrix
  240/240 identical, PAM table 16/16 identical, zero diffs**; reproduced into C# at full `double`
  precision.
- **Algorithm (verbatim from `calc_cfd`):** `score = 1`; guide/off `T→U`; for each position `i` (0..19,
  5'→3'): matched → ×1, else key `r{guide[i]}:d{complement(off[i])},{i+1}` → `mm_scores[key]`; finally
  ×`pam_scores[last-2-PAM-nt]`. **Orientation** index 0 = PAM-distal, position 20 = PAM-proximal seed
  (getting this backwards is the classic CFD bug). **Key convention:** `rX` = guide (RNA) base T→U,
  `dY` = complement of the off-target base.
- **PAM table** (16 NGG-region dinucleotides): `GG=1.0`, `AG=0.259259`, `CG=0.107143`, `GA=0.069444`,
  `TG=0.038961`, `GC=0.022222`, `GT=0.016129`, all others `=0.0`.
- **Independent cross-checks** (Python re-derived from the decoded pickle): perfect+GG→**1.0**; iGWOS
  doctests **0.4635989007074176** and **0.5140384614450001**; perfect+GA/AG/TG→0.069444/0.259259/0.038961;
  perfect+AA→0.0; single mm rG:dT,1→0.9, rC:dT,5→0.571428571, rU:dG,7→0.6875, rG:dA,16→0.0, rC:dT,20→0.5;
  product rG:dT,1·rC:dT,20→0.45 — all match test assertions. Orientation counterfactual (guide C×20:
  rC:dT,1=1.0 but rC:dT,20=0.5) makes the guard a real reversal-detector.

## Stage B — implementation

- **MIT/Hsu code path** (`:607–721`) realises `score1/score2/score3` exactly as `calcHitScore`
  (`maxDist=19`, `mmCount<2`→score2=1, `mmCount==0`→score3=1, ×100); the aggregate implements
  `100/(100+Σ)·100`; the genome overload combines per-hit off-target scores (exact on-target excluded).
  Edge cases: null/empty guide or off-target → throw; non-20-nt → `ArgumentException`; empty hit set → 100;
  genome with no off-targets → 100. 39/39 off-target tests + the 11 MIT tests pass.
- **CFD code path** is additive (no existing method/signature/test changed); loops i=0..19, skips matches,
  builds the key, multiplies `CfdMismatchScores[key][i]` then × `CfdPamScores[pam]`, returns [0,1]. Edge
  cases: null guide/off/PAM → `ArgumentNullException`; wrong length (≠20) → `ArgumentException`; PAM length
  ∉{2,3} → `ArgumentException`; non-ACGT → `ArgumentException`. All tested.
- **Test-quality audit (HARD gate): PASS, no green-washing.** Every expected value traces to an
  independent Python re-derivation or a published doctest — none read off the C# arrays; the computed value
  is the NUnit `actual`. The MIT orientation guard asserts pos13(W=0.851)→14.9 **and** pos0(W=0)→100 (fails
  on a reversed W → 68.3/41.7); the CFD guard asserts pos1→1.0 **and** pos20→0.5 (fails on reversal). The
  CFD guard caught a real off-base-complementation distinction during authoring (off C → dG, not dC),
  confirming the suite is non-tautological.
- **Defects:** none in either part.

## Findings & follow-ups

- **No code defect, no production-code change (State CLEAN)** in either part. W vector 20/20 vs CRISPOR
  `hitScoreM`; CFD matrices 240/240 + 16/16 cross-source identical (CRISPOR + iGWOS); orientation and key
  conventions pinned from the primary sources; all worked numbers reproduced independently.
- The in-tree **MIT orientation-guard** test was verified source-correct and **committed** with this work.
- **Finding C7 fully resolved:** CFD was the last off-target residual (implemented here); Doench Rule Set 2 /
  Azimuth — the earlier on-target residual — was cleared separately (2026-06-18, see
  [[crispr-guide-001-report]]).
- Not-yet-ingested MolTools algorithm doc `docs/algorithms/MolTools/Off_Target_Analysis.md` (backlog slug
  `off-target-analysis`) describes the underlying `FindOffTargets` scan surface; left pending — this report
  validates the MIT/CFD scoring layer, not that doc. Research-grade, not for clinical use.
