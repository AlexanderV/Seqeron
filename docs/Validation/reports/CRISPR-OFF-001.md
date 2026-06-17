# Validation Report: CRISPR-OFF-001 — CRISPR Off-Target Analysis (MIT/Hsu 2013 model)

- **Validated:** 2026-06-17   **Area:** MolTools
- **Canonical method(s):** `CrisprDesigner.CalculateMitHitScore(string guide20, string offTarget20)` (single-hit),
  `CrisprDesigner.CalculateMitSpecificityScore(IEnumerable<double>)` and the genome-scanning overload
  `CalculateMitSpecificityScore(string, DnaSequence, int, CrisprSystemType)` (aggregate).
  (Pre-existing honest-heuristic `FindOffTargets` / `CalculateSpecificityScore` / `CalculateOffTargetScore`
  unchanged; this Phase-3 re-validation is scoped to the MIT/Hsu enhancement added in commit 129c2ca.)
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS
- **End-state:** ✅ CLEAN (no defect; uncommitted in-tree orientation test verified source-correct and committed)

## Scope

Independent re-validation (fresh validator context, different from the implementer session) of the
MIT/Hsu 2013 off-target scoring model added to `CrisprDesigner`, plus the uncommitted test
`CalculateMitHitScore_WeightOrientation_PamProximalIsHighPenalty` left in the working tree. Governing
rule: tests follow the primary source; the code obeys the tests; every expected value re-derived from
the source, never read off the code's own array.

## Stage A — Description

### Primary source retrieved this session
- **Hsu, Scott, Weinstein, Ran, Konermann, et al. "DNA targeting specificity of RNA-guided Cas9
  nucleases." Nat Biotechnol 31, 827–832 (2013).** PMID 23873081, doi:10.1038/nbt.2647. The MIT
  "Scores of single hits" scheme (originally on crispr.mit.edu/about, now retired).
- **Canonical reference implementation:** CRISPOR `crispor.py` — `calcHitScore` / `calcMitGuideScore`
  and the `hitScoreM` weight list.
  `https://github.com/maximilianh/crisporWebsite/blob/master/crispor.py`
  (raw read this session: `https://raw.githubusercontent.com/maximilianh/crisporWebsite/master/crispor.py`).
  The `calcHitScore` docstring reads: `" see 'Scores of single hits' on http://crispr.mit.edu/about "`
  and the comment `# The Patrick Hsu weighting scheme`; the list is preceded by `# aka Matrix "M"`.
- **Orientation corroboration:** Wikipedia "Off-target genome editing" — PAM-proximal seed (10–12 nt
  adjacent to PAM) is critical for specificity; "mismatches in the 5' end … are more tolerated".

### W mismatch-penalty vector — element-by-element diff vs PRIMARY source
Reference `hitScoreM` (CRISPOR, transcribing crispr.mit.edu / Hsu 2013):

```
[0, 0, 0.014, 0, 0, 0.395, 0.317, 0, 0.389, 0.079, 0.445, 0.508, 0.613, 0.851, 0.732, 0.828, 0.615, 0.804, 0.685, 0.583]
```

Code `MitHitScoreWeights` (`CrisprDesigner.cs:611-615`): **identical, all 20 elements** —
index 0..19 = `0,0,0.014,0,0,0.395,0.317,0,0.389,0.079,0.445,0.508,0.613,0.851,0.732,0.828,0.615,0.804,0.685,0.583`.
No transcription error.

### Orientation — pinned from the source (not from the code)
In `calcHitScore`, the loop indexes `hitScoreM[pos]` for `pos = 0 .. 19` over the protospacer given
5'→3'. The protospacer 5' end is **PAM-distal**, the 3' end (index 19) is **PAM-proximal / seed**.
Therefore **index 0 = PAM-distal (low/zero weight), index 19 = PAM-proximal**, and the high-weight
positions (max `W[13]=0.851`, plus `0.732/0.828/0.804/0.685/0.583` toward the 3' end) form the
PAM-proximal seed — exactly the biologically expected high-penalty end (Wikipedia/Hsu: seed mismatches
least tolerated). The code's doc comment and array orientation match the source.

### Formula (verbatim from `calcHitScore`)
- `score1 = Π over mismatched positions i of (1 − W[i])`
- `score2 = 1` if `nmm < 2`, else `1 / ( ((19 − meanInterMismatchDist)/19) · 4 + 1 )`, where
  `meanInterMismatchDist` = mean of consecutive inter-mismatch gaps (`dists`), `maxDist = 19`.
- `score3 = 1` if `nmm == 0`, else `1 / nmm²`.
- `hitScore = score1 · score2 · score3 · 100`.
- Aggregate (`calcMitGuideScore`): `100 / (100 + Σ hitScores) · 100` (CRISPOR rounds to int; the C#
  method returns the unrounded value — a documented, benign difference).

### Independent cross-check (Python re-derivation from the source, NOT the C# array)
Guide `GACGCATAAAGATGAGACGC`:

| Case | Re-derived (Python from source formula) | Test expects |
|------|------------------------------------------|--------------|
| perfect match | 100.0 | 100.0 ✓ |
| 1 mm @ pos 0 (W=0) | 100.0 | 100.0 ✓ |
| 1 mm @ pos 5 (W=0.395) | 60.5 | 60.5 ✓ |
| 1 mm @ pos 19 (W=0.583) | 41.7 | 41.7 ✓ |
| 2 mm @ {5,15} | 0.10406·0.3454545·0.25·100 = 0.8987 | 0.8987 ✓ |
| 1 mm @ pos 13 (W=0.851) | 14.9 | 14.9 ✓ |
| aggregate {60.5} | 100/(100+60.5)·100 = 62.30529595 | 62.30529595 ✓ |
| aggregate {60.5,41.7} | 100/(100+102.2)·100 = 49.45598417 | 49.45598417 ✓ |

**Orientation counterfactual** (reversed W vector): pos 13 → `1−W_rev[13]=1−0.317` → 68.3;
pos 0 → `1−W_rev[0]=1−0.583` → 41.7. Both differ from the correct 14.9 / 100.0, confirming the
orientation guard test is a genuine reversal-detector.

**Stage A verdict: PASS** — formula, W values (20/20), orientation, ×100, and aggregate all match the
primary source and reference implementation; all worked numbers reproduced independently.

## Stage B — Implementation

- **Code path:** `CrisprDesigner.cs:607-721` (`MitHitScoreWeights`, `CalculateMitHitScore`,
  `CalculateMitSpecificityScore` ×2).
- **Formula realised correctly:** `CalculateMitHitScore` (lines 629-674) implements score1/score2/score3
  exactly as `calcHitScore`, with `maxDist=19`, the `mmCount<2`→score2=1 and `mmCount==0`→score3=1
  special-cases, and `×100`. `CalculateMitSpecificityScore(IEnumerable<double>)` (688-693) implements
  `100/(100+Σ)·100`. The genome overload (705-721) reuses `FindOffTargets` (off-targets only; exact
  on-target excluded) and combines per-hit scores. All recomputed values match the C# output
  (39/39 off-target tests pass).
- **Edge/error cases:** null/empty guide or off-target → throw (`ArgumentNullException`); non-20-nt →
  `ArgumentException` (MIT-011); empty hit set → 100 (MIT-006); genome with no off-targets → 100
  (MIT-010). All present and tested.
- **Test quality audit:** the 11 MIT tests assert exact source-derived values (perfect→100 boundary,
  zero-weight position, two known single-mismatch penalties, two-mismatch all-three-terms,
  two aggregate cases, two genome-scan cases, length/null guards). Not tautological, deterministic.
  The **uncommitted** test `CalculateMitHitScore_WeightOrientation_PamProximalIsHighPenalty` is
  **source-correct**: it asserts pos13(W=0.851)→14.9 AND pos0(W=0)→100, both of which I re-derived
  from the source, and it fails if the W vector is reversed (would give 68.3 / 41.7). It is a true
  orientation guard, not a code-echo. **Kept as-is and committed** with this work.
- **Defects:** none. No transcription or orientation defect found; no code or test change required
  beyond committing the in-tree orientation test.

**Stage B verdict: PASS.**

## Verdict & follow-ups

- **Stage A: PASS · Stage B: PASS · End-state: ✅ CLEAN.**
- W vector confirmed 20/20 against CRISPOR `hitScoreM` / crispr.mit.edu (Hsu 2013); orientation
  confirmed index 0 = PAM-distal, 19 = PAM-proximal/seed.
- Uncommitted orientation guard test verified source-correct and committed.
- Residual (documented, not a defect of this unit; see FINDINGS_REGISTER C7): Doench Rule Set 2 /
  Azimuth (GBT, no coefficient table) and CFD (Doench 2016, binary pickle) intentionally not
  implemented. These are on-target / alternative off-target models, out of scope for the MIT/Hsu unit.
- Full unfiltered suite: **6773 passed, 0 failed**; build 0 errors.
