# Validation Report: CRISPR-OFF-001 — Off-Target Analysis

- **Validated:** 2026-06-12   **Area:** MolTools (CRISPR)
- **Canonical method(s):**
  - `CrisprDesigner.FindOffTargets(string, DnaSequence, int, CrisprSystemType)`
  - `CrisprDesigner.CalculateSpecificityScore(string, DnaSequence, CrisprSystemType)`
  - (private) `CalculateOffTargetScore`, `CountMismatches`, `GetMismatchPositions`
- **Source file:** `src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/CrisprDesigner.cs:331-440`
- **Test file:** `tests/Seqeron/Seqeron.Genomics.Tests/CrisprDesigner_OffTarget_Tests.cs`
- **Stage A verdict:** PASS-WITH-NOTES
- **Stage B verdict:** PASS
- **End state:** CLEAN (honestly declared heuristic; one cosmetic comment fix)

---

## Stage A — Description

### Sources opened & what they confirm
- **Hsu et al. (2013), *Nat. Biotechnol.* (DNA targeting specificity of RNA-guided Cas9 nucleases)** — SpCas9 tolerates guide:DNA mismatches in a position- and number-dependent manner; the PAM-proximal **seed** region (8–12 nt at the 3′ end of the 20-nt guide) is most sensitive, so seed mismatches abolish/reduce cleavage more than PAM-distal mismatches. Off-target activity observed broadly with up to ~3–5 mismatches.
- **Fu et al. (2013), *Nat. Biotechnol.*** — substantial off-target mutagenesis can occur; truncated (17–18 nt) gRNAs improve specificity.
- **Wikipedia "Off-target genome editing"** — off-targets are near-matches with mismatches at PAM-supported loci; PAM-distal mismatches more tolerated than PAM-proximal (seed).

### The published MIT/Hsu (2013) off-target score (canonical reference formula)
For one candidate off-target with mismatch set `M` (per Hsu 2013 / the MIT CRISPR Design tool, as reimplemented in CRISPOR / crisprScore):

```
hit_score = ( ∏_{i ∈ M} (1 − W[i]) )  ×  1 / ( ((19 − d̄)/19) × 4 + 1 )  ×  1 / (nmm²)
```

- `W[i]` = published 20-element position-specific mismatch penalty weight vector
  (position 1 = PAM-distal 5′ end … position 20 = PAM-proximal 3′ end):
  `[0.000, 0.000, 0.014, 0.000, 0.000, 0.395, 0.317, 0.000, 0.389, 0.079, 0.445, 0.508, 0.613, 0.851, 0.732, 0.828, 0.615, 0.804, 0.685, 0.583]`
  (verified weight values, e.g. W[3]=0.014, W[6]=0.395, W[14]=0.851 — larger toward the PAM-proximal end).
- `d̄` = mean pairwise distance (bp) between mismatches (1 if a single mismatch).
- `nmm` = number of mismatches.
The genome-level **MIT specificity score** = `100 / (100 + Σ hit_scores)`, scaled ×100.

### Edge-case semantics (sourced)
- 0 mismatches ⇒ on-target (an *exact* match), **not** an off-target — highest cleavage risk but excluded from the off-target list by definition.
- Mismatch in seed (PAM-proximal) ⇒ greater specificity penalty than a PAM-distal mismatch.
- More mismatches / wider spread ⇒ lower off-target hit score (less likely to cut).
- PAM (NGG, also NAG at ~5× lower efficiency) is required for cleavage; no PAM ⇒ no off-target.

### Independent hand-check of the *intended biology direction*
A single seed-region mismatch must score worse (more penalty / lower specificity) than a single PAM-distal mismatch. Confirmed against Hsu 2013 and the W-vector (PAM-proximal W[i] ≫ PAM-distal W[i]).

### Findings / divergences (Stage A → PASS-WITH-NOTES)
The TestSpec and algorithm doc **do not claim to implement the published MIT/Hsu numeric score or the 20 W[i] weights.** The doc (`docs/algorithms/MolTools/Off_Target_Analysis.md`) is explicit: *Implementation Status = "Simplified"*, scores are *"a simple position-dependent mismatch penalty"* and *"intended for heuristic ranking rather than experimental prediction."* It tabulates its own weights (seed = 5, non-seed = 2). Because the description is presented as a **declared heuristic** consistent with the Hsu seed-region biology — not as the MIT/CFD model — the abstract description is correct for what it claims. Note (the reason for PASS-WITH-NOTES): the heuristic is qualitatively, not quantitatively, faithful to Hsu 2013 (no W-vector, no mean-pairwise-distance term, no `1/nmm²` term, additive instead of multiplicative). This is acceptable per protocol because it is **declared**, not advertised as the MIT score.

---

## Stage B — Implementation

### Code path reviewed
- `FindOffTargets` (`CrisprDesigner.cs:331-373`): validates inputs (empty guide → `ArgumentNullException`; null genome → `ArgumentNullException`; `maxMismatches` ∉ [0,5] → `ArgumentOutOfRangeException`; guide length ≠ system guide length → `ArgumentException`), scans PAM sites on both strands, yields a site only when `0 < mismatches ≤ maxMismatches`.
- `CalculateOffTargetScore` (`:421-440`): seed = PAM-proximal 12 bp (`PamAfterTarget` ⇒ last 12; else first 12); per mismatch adds `inSeed ? 5 : 2`.
- `CalculateSpecificityScore` (`:378-391`): finds off-targets at a fixed cap of 4 mismatches, returns `100` if none, else `max(0, 100 − Σ OffTargetScore)`.

### Formula realised correctly (vs the *declared heuristic*)? Yes.
The code computes exactly the seed/distal additive penalty the doc specifies. Seed orientation is correct: for SpCas9 (`PamAfterTarget = true`, guide length 20) seed = positions 8..19 (3′ PAM-proximal); for Cas12a (`PamAfterTarget = false`) seed = positions 0..11 (5′ PAM-proximal). This matches Hsu's PAM-proximal seed definition.

### Cross-verification table (recomputed by hand, confirmed by tests, all green)
| Scenario | Expected (hand) | Test asserts | Result |
|---|---|---|---|
| 1 distal mismatch @ pos 0 | OffTargetScore = 2 | 2.0 | ✓ |
| 2 distal mismatches @ [0,1] | 2+2 = 4 | 4.0 | ✓ |
| 3 distal mismatches @ [0,1,2] | 6 | 6.0 | ✓ |
| 1 seed mismatch @ pos 19 (SpCas9) | 5 | 5.0 | ✓ |
| seed vs distal | 5 > 2 (seed worse) | distal 2.0, seed 5.0 | ✓ |
| Cas12a 1 mismatch @ pos 0 (seed, 5′) | 5 | 5.0 | ✓ |
| Specificity, 1 distal off-target | 100 − 2 = 98 | 98.0 | ✓ |
| Specificity, seed off-target | 100 − 5 = 95 | 95.0 | ✓ |
| No off-targets | 100 | 100 | ✓ |
| Exact match | excluded (mismatches = 0) | Is.Empty | ✓ |
| No PAM | excluded | Is.Empty | ✓ |
| Reverse strand | found, IsForwardStrand=false | ✓ | ✓ |

### Edge cases
0 mismatches excluded (`mismatches > 0` guard); `mismatches ≤ maxMismatches` enforced; `MismatchPositions.Count == Mismatches` (both walk identical loop); PAM required (only PAM sites scanned); score clamped to [0,100]; deterministic. All covered by tests M-001…M-012, S-001…S-004, plus edge/invariant tests.

### Defect fixed
- **Misleading comment** at `CrisprDesigner.cs:432` previously read *"Mismatches in seed region are more tolerated"* — this is **backwards** (Hsu 2013: seed mismatches are *less* tolerated). The numeric behaviour was already correct (seed = 5 > distal = 2). Replaced with a correct comment citing Hsu 2013. No behavioural change; all tests still pass.

### Test quality audit
Tests assert exact sourced/heuristic values (not "no-throw" tautologies), are deterministic, and cover every Stage-A edge case. Strong.

---

## Verdict & follow-ups

- **Stage A: PASS-WITH-NOTES** — description is biologically correct and, importantly, **honestly declared as a simplified heuristic** (doc status "Simplified", "heuristic ranking rather than experimental prediction"). It does **not** falsely advertise the MIT/Hsu numeric score, so this is not a LIMITED scientific-correctness gap.
- **Stage B: PASS** — code faithfully realises the declared heuristic; all worked examples recomputed and match; one misleading comment corrected.
- **End state: CLEAN.**
- **Files changed:** `src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/CrisprDesigner.cs` (comment-only fix at line 432).
- **Tests:** filter `~OffTarget` = 33 passed / 0 failed; full suite = **4461 passed, 0 failed**.
- **Optional future enhancement (not a defect):** an `MitOffTargetScore`/CFD variant implementing the published W-vector + mean-pairwise-distance + `1/nmm²` terms could be added as an explicitly-named experimental-grade scorer; the current heuristic remains valid for ranking.
