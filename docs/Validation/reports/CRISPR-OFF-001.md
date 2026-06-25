# Validation Report: CRISPR-OFF-001 — CRISPR Off-Target Analysis (MIT/Hsu 2013 model)

> **Independent re-confirmation 2026-06-24 (fresh validator context).** Re-downloaded `crispor.py`
> from `maximilianh/crisporWebsite` (master, 8956 lines) and read `hitScoreM`/`calcHitScore`/
> `calcMitGuideScore` directly. W vector confirmed 20/20 identical to the C# `MitHitScoreWeights`
> (`CrisprDesigner.cs:654-658`). C# `CalculateMitHitScore` (lines 672-717) realises `calcHitScore`
> verbatim (maxDist=19; mmCount<2→score2=1; mmCount==0→score3=1; ×100). Re-derived all cross-check
> values in Python from the freshly-downloaded reference (not the C# array): perfect→100, mm@0(W=0)→100,
> mm@5→60.5, mm@13(W=0.851)→14.9, mm@19→41.7, mm@{5,15}→0.8987, agg[60.5]→62.30529595,
> agg[60.5,41.7]→49.45598417 — all match the test assertions exactly. Orientation confirmed
> index 0 = PAM-distal 5', index 19 = PAM-proximal seed. Off-target + CFD suite: 71/71 passed,
> build 0 warnings/0 errors. **No code change — re-confirmed ✅ CLEAN.**

- **Validated:** 2026-06-17 (re-confirmed 2026-06-24)   **Area:** MolTools
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
  Azimuth (GBT, no coefficient table) and — at the time of the MIT/Hsu validation — CFD (Doench 2016,
  binary pickle) intentionally not implemented. **CFD is now implemented (see the CFD section below); Rule Set 2 / Azimuth is now also implemented (2026-06-18) — see `reports/CRISPR-GUIDE-001.md`. C7 fully resolved.**
- Full unfiltered suite: **6773 passed, 0 failed**; build 0 errors.

---

# Validation Report: CRISPR-OFF-001 (continued) — CFD off-target score (Doench 2016)

- **Validated/Implemented:** 2026-06-17   **Area:** MolTools
- **Canonical method added:** `CrisprDesigner.CalculateCfdScore(string sgRna20, string offTarget20, string offTargetPam)` → CFD score in [0,1].
- **Stage A verdict:** PASS · **Stage B verdict:** PASS · **End-state:** ✅ CLEAN (CFD residual cleared).

## Scope

Implementation of the deferred CFD (Cutting Frequency Determination) off-target score — the last
off-target residual of CRISPR-OFF-001 / C7. Governing rule: tests follow the primary source/spec; the
code obeys the tests; every matrix value and expected number re-derived independently from the external
source, never read off the C# code.

## Stage A — Description

### Faithfulness boundary — CLEARED this session
CFD's mismatch + PAM matrices originate from Doench et al. 2016 (Nat Biotechnol 34:184, PMID 26780180)
and are shipped by the reference tools as binary pickles. The boundary requires obtaining them as
verbatim numbers cross-checked across **two independent sources**. Both conditions were met:

- **Authoritative source 1:** `maximilianh/crisporWebsite`, `CFD_Scoring/mismatch_score.pkl` +
  `pam_scores.pkl` + `cfd-score-calculator.py` (the canonical John Doench reference calculator).
- **Authoritative source 2:** `bm2-lab/iGWOS`, `CFD/mismatch_score.pkl` + `pam_scores.pkl` +
  `otscore.py` (independent repository shipping the same matrices, with documented `calcCfdScore`
  doctest oracles).
- **Cross-check (this session):** both pickles decoded to text and diffed element-by-element —
  **mismatch matrix 240/240 entries identical, PAM table 16/16 entries identical, ZERO diffs.**
  The decoded values were reproduced into C# at full `double` precision (exact decoded bit patterns).

Decoding the authoritative pickle to text IS obtaining the verbatim numbers (the pickle is the
canonical distribution); nothing was fabricated or approximated.

### Algorithm (verbatim from `calc_cfd`, cfd-score-calculator.py / otscore.py — identical in both)
```
score = 1
sg = offTarget.replace('T','U'); wt = guide.replace('T','U')
for i, off_base in enumerate(sg):            # i = 0..19, 5'->3'
    if wt[i] == off_base: score *= 1
    else: key = 'r'+wt[i]+':d'+complement(off_base)+','+str(i+1); score *= mm_scores[key]
score *= pam_scores[ pam ]                    # pam = off[-2:]  (last two PAM nt)
```

### Orientation — pinned from the SOURCE (not the code)
In the reference, `sg = off[:-3]` and `pam = off[-2:]`, so the 20-nt protospacer precedes the PAM and
the loop enumerates it 5'→3' with key position `i+1`. Therefore **position 1 (index 0) = 5' / PAM-DISTAL
end; position 20 = 3' / PAM-PROXIMAL (seed) end.** Getting this backwards is the classic CFD bug; the
orientation guard test detects reversal.

### Key convention — pinned from the source
`rX` = the **guide (RNA)** base, T written as U. `dY` = the **complement of the off-target base** (the
base on the off-target's non-target DNA strand that pairs the guide). A position contributes a penalty
only when guide[i] ≠ offTarget[i]; matched positions contribute 1.0. Perfect match + GG PAM → 1.0.

### PAM table (16 NGG-region dinucleotides; the N of NGG contributes 1)
`GG=1.0` (canonical), `AG=0.259259`, `CG=0.107143`, `GA=0.069444`, `TG=0.038961`, `GC=0.022222`,
`GT=0.016129`, all others (AA/AC/AT/CA/CC/CT/TA/TC/TT) `=0.0`. (Verbatim from both pickles.)

### Contract implemented
20-nt guide vs 20-nt off-target protospacer + the off-target PAM (2-nt or 3-nt; only the last two nt
scored), A/C/G/T only (case-insensitive; guide T treated as U). Insertions/deletions and non-ACGT bases
are unsupported (CFD undefined) and throw.

### Independent cross-checks (Python re-derivation from the decoded pickle, NOT from the C# arrays)
| Case | Re-derived (Python from decoded matrices) | Test expects |
|------|-------------------------------------------|--------------|
| perfect match + GG | 1.0 | 1.0 ✓ |
| published iGWOS doctest: G×20 vs G…AAA + GG | 0.4635989007074176 (= rG:dT,18·rG:dT,19·rG:dT,20) | 0.4635989007074176 ✓ |
| published iGWOS doctest: G×20 vs aaaaGaGaG… +gg | 0.5140384614450001 | 0.5140384614450001 ✓ |
| perfect + GA / AG / TG PAM | 0.069444 / 0.259259 / 0.038961 | same ✓ |
| perfect + AA PAM | 0.0 | 0.0 ✓ |
| single mm rG:dT,1 (pos 1) | 0.9 | 0.9 ✓ |
| single mm rC:dT,5 (pos 5) | 0.571428571 | 0.571428571 ✓ |
| single mm rU:dG,7 (pos 7, guide T→U, off C) | 0.6875 | 0.6875 ✓ |
| single mm rG:dA,16 (pos 16) | 0.0 | 0.0 ✓ |
| single mm rC:dT,20 (pos 20) | 0.5 | 0.5 ✓ |
| two mm rG:dT,1 · rC:dT,20 (product) | 0.45 | 0.45 ✓ |

**Orientation counterfactual** (guide C×20): rC:dT,1 = **1.0** but rC:dT,20 = **0.5**; if the position
axis were reversed the two would swap, so the guard test (asserting 1.0 at pos 1 AND 0.5 at pos 20)
fails on reversal — a genuine reversal-detector.

**Stage A verdict: PASS** — matrices 240/240 + 16/16 cross-source identical, algorithm verbatim,
orientation + key convention pinned from the source, all worked numbers reproduced independently.

## Stage B — Implementation

- **Code path:** `CrisprDesigner.cs` — `CfdMismatchScores` (12 keys × 20 positions), `CfdPamScores`
  (16 keys), `CfdComplement`, `CalculateCfdScore`. Additive; no existing method/signature/test changed.
- **Realised correctly:** loops i=0..19, skips matches (×1), builds key `r{guide,T→U}:d{complement(off)}`
  and multiplies `CfdMismatchScores[key][i]`, then × `CfdPamScores[pam last-2-nt]`. Returns [0,1].
- **Edge/error cases:** null guide/off/PAM → `ArgumentNullException`; empty guide → `ArgumentNullException`;
  wrong-length guide/off (≠20) → `ArgumentException`; PAM length ∉{2,3} → `ArgumentException`; non-ACGT in
  guide/off/PAM → `ArgumentException`. All tested.
- **Tests:** new fixture `CrisprDesigner_Cfd_Tests` (32 `[Test]` methods): perfect→1.0 (×2: 2-nt and 3-nt
  PAM), two published doctest oracles, six single-mismatch exact-matrix-entry, four PAM-application
  (GA/AG/TG/AA), product-of-penalties, mismatch×PAM combined, **orientation guard**
  `Cfd_OrientationGuard_Position1VsPosition20_NotReversed`, unit-interval / lowercase / determinism
  invariants, and the full edge/error set. Every expected value traces to the independent Python
  re-derivation or a published doctest — none read off the C# arrays; computed value is the NUnit
  `actual` (no NUnit2007). The orientation guard caught a real off-base-complementation distinction
  during authoring (off C → dG, not dC), confirming the suite is non-tautological.
- **Defects:** none.

**Stage B verdict: PASS.**

## Verdict & follow-ups

- **Stage A: PASS · Stage B: PASS · End-state: ✅ CLEAN.**
- CFD matrices confirmed 240/240 + 16/16 across CRISPOR and iGWOS; orientation index 0 = PAM-distal,
  19 = PAM-proximal; key `rX` = guide(T→U), `dY` = complement(off-target base).
- **Remaining C7 residual (as of this 2026-06-17 session): ONLY Doench Rule Set 2 / Azimuth** (gradient-boosted-tree, no coefficient
  table — not reproducible from published numbers without the trained model). CFD is no longer a residual.
  **UPDATE 2026-06-18: Rule Set 2 / Azimuth now implemented (sklearn-free reconstruction of the trained model); C7 fully resolved — see `reports/CRISPR-GUIDE-001.md`.**
- Full unfiltered suite: **6812 passed, 0 failed**; build 0 errors AND 0 warnings.
