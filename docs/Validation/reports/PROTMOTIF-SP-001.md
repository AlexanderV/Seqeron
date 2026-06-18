# Validation Report: PROTMOTIF-SP-001 — Signal Peptide Cleavage-Site Prediction (von Heijne 1986 weight matrix)

- **Validated:** 2026-06-16   **Area:** ProteinMotif
- **Canonical method(s):** `ProteinMotifFinder.PredictSignalPeptide(string, bool, double)` (+ private `BuildWeightMatrix(int[][], double[])`)
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS-WITH-NOTES (one test-quality defect found and fixed in-session; no code defect)
- **End-state:** ✅ CLEAN

## Stage A — Description

### Sources opened this session (retrieved, not trusted from citation label)

| Source | URL | What it confirmed |
|--------|-----|-------------------|
| EMBOSS 6.6.0 `data/Esig.euk` | https://raw.githubusercontent.com/lauringlab/CodonShuffle/master/lib/EMBOSS-6.6.0/emboss/data/Esig.euk | The full 20×16 eukaryotic count matrix (161 sequences, von Heijne 1986) + Expect column. Header: "The cleavage site is between +1 and -1". |
| EMBOSS 6.6.0 `data/Esig.pro` | https://raw.githubusercontent.com/lauringlab/CodonShuffle/master/lib/EMBOSS-6.6.0/emboss/data/Esig.pro | The full 20×16 prokaryotic count matrix (36 sequences) + Expect column. |
| EMBOSS 6.6.0 `sigcleave.c` | https://raw.githubusercontent.com/lauringlab/CodonShuffle/master/lib/EMBOSS-6.6.0/emboss/sigcleave.c | Scoring loop (`pval=-13`, `nval=2`, `ic=13+pval=0`, window `j=i-13..i+1`); `maxweight`/`maxsite` argmax; matrix transform `mat=log(mat/expected)`, expected=last column; zero-count pseudocount `1.0e-10` iff `j==10 || j==12` else `1.0`; signal `[maxsite+pval, maxsite-1]`, mature start `maxsite`. |
| EMBOSS `sigcleave` application doc | https://emboss.sourceforge.net/apps/release/6.6/emboss/apps/sigcleave.html | ACH2_DROME worked example (9 hits); `-minweight ≥ 3.5` (≈95% sens/spec, 75-80% cleavage accuracy); eukaryotic default, `-prokaryote` selects the other matrix. |

### Formula check
- Score `S(i) = Σ_{col=0..14} ln( C(residue, col) / E(residue) )` over the 15-column window (positions −13..+2), natural log. **Matches `sigcleave.c` verbatim** (`log()` is C natural log; expected read from the last column `mat[i][d2-1]`).
- Pseudocount: zero count → `1.0e-10` at columns 10 (−3) and 12 (−1), else `1.0`, applied **before** the log. Matches `if(j==10||j==12) … 1.0e-10; else … 1.0`.
- Argmax single-site selection; cleavage between −1 and +1; mature protein starts at the +1 residue (`maxsite`). Matches source.

### Edge-case semantics
- No intrinsic cutoff: a best site is always returned for an in-window sequence; `-minweight` only flags likelihood. Consistent with EMBOSS.
- Min-length: EMBOSS scores any length by skipping off-window columns; the implementation returns `null` below 15 residues. This is a *documented, accepted* simplification (Evidence ASSUMPTION-1, doc §5.3): a partial-window score is not meaningful for cleavage prediction, and all in-scope signal peptides are ≥ 15 aa. Does not change any in-window candidate score. Recorded as a minor divergence, not a defect.
- Non-standard residues (X/B/Z/*) are not in the matrix → contribute 0, no throw. Consistent with EMBOSS handling of unmapped residues.

### Independent cross-check (numbers retrieved/derived this session)
Independent Python re-implementation (separate from the repo) using the **verbatim Esig.euk matrix** and natural-log transform, run on UniProt P17644 (ACH2_DROME):

| Quantity | External re-derivation | EMBOSS doc | Match |
|----------|-----------------------|------------|-------|
| Best score | 13.7390400704164 | "Score 13.739" | ✅ |
| Mature start (1-based) | 42 | 29→41 signal ⇒ mature 42 | ✅ |
| Signal residues −13..−1 | `LLVLLLLCETVQA` | `LLVLLLLCETVQA` | ✅ |
| Runner-up score / site | 12.135281025149874 @ mature 39 | "(2) Score 12.135 … 26->38" | ✅ |
| Hits 3/4/5 | 10.465 @41, 7.360 @541, 6.981 @343 | 10.465, 7.360, 6.981 | ✅ |

The 161-column-sum / Expect values and the full eukaryotic and prokaryotic matrices in the code were compared cell-by-cell against the EMBOSS data files: **identical**.

### Findings / divergences
- Description is biologically and mathematically correct and faithful to von Heijne (1986) via the EMBOSS reference. One documented simplification (min length = 15). **Stage A: PASS.**

## Stage B — Implementation

### Code path reviewed
`src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/ProteinMotifFinder.cs`
- `PredictSignalPeptide` (lines ~554-608): window loop `j = i + (-13) + col`, col 0..14; `MatrixResidues.IndexOf` skips unmapped; argmax via `weight > bestScore`; `CleavagePosition = bestSite+1`; window `[bestSite-13, bestSite+1]`.
- `BuildWeightMatrix` (lines ~616-637): `ln(count/expect)` with `1.0e-10` at columns `MinusThreeColumn=10`/`MinusOneColumn=12`, else `1.0`.
- Matrices `EukaryoticCounts`/`EukaryoticExpect`/`ProkaryoticCounts`/`ProkaryoticExpect` (lines ~451-520).

### Formula realised correctly?
Yes — verbatim match to `sigcleave.c`: column-0 = position −13, columns 10/12 = −3/−1, window spans −13..+2 (15 columns), natural log, argmax. Confirmed by reproducing the worked example to full double precision (13.7390400704164) and the runner-up (12.135281025149874).

### Cross-verification table recomputed vs code
| Case | External value (this session) | Code output | Match |
|------|------|------|------|
| ACH2_DROME best score | 13.7390400704164 | 13.7390400704164 (M2, 1e-3) | ✅ |
| ACH2_DROME mature start | 42 | 42 (M1) | ✅ |
| ACH2_DROME window −13..+2 | `LLVLLLLCETVQANP` | `LLVLLLLCETVQANP` (M3) | ✅ |
| `AAAAAAAAAAAAGAN` best score | 3.2425554865825688 | 3.2425554865825688 (M10) | ✅ |
| Snapshot `TestProtein` | score 11.152834574562595, pos 17 | 11.152834574562595, 17 | ✅ |

### Variant/delegate consistency
Single public method; eukaryotic vs prokaryotic matrices both built by the same `BuildWeightMatrix`; verified distinct (S3). Case-insensitivity (S1) and `minWeight` flag independence (S2) confirmed.

### Test quality audit (HARD gate)
- **Sourced expectations:** M1/M2/M3 assert the externally-sourced EMBOSS values (42, 13.739, `LLVLLLLCETVQA`/window). M2 tolerance 1e-3 is anchored to the EMBOSS-printed precision and the code matches to full precision anyway.
- **Defect found & fixed (test-only):** **M10** previously asserted only `Score == HandScoreEukaryotic(seq, …)` — a helper that re-implements the formula from its *own* copy of the matrix. While it cross-checks an independent code path, it is partly a code-echo (a shared blind spot in a copied constant would pass). **Fix:** added an assertion against the externally re-derived literal `3.2425554865825688` (independent Python, verbatim Esig.euk) at `Within(1e-12)`, keeping the helper cross-check as a secondary check. No production code change; full suite re-run green.
- **M4 argmax:** uses `Is.GreaterThan(runnerUp=12.135281025149874)`. The runner-up exact value is sourced, but the public API returns only the best site, so a strict-inequality against the sourced runner-up is the strongest assertion the API permits. Acceptable (not green-washing).
- **Edge/error cases covered:** null (M6), empty (M7), 14-aa short (M8), exactly-15 boundary (M9), case-insensitivity (S1), threshold flip (S2), prokaryotic matrix (S3), non-standard residue X (C1). All Stage-A branches exercised.
- **Snapshot** (`ProteinMotifSnapshotTests`) value verified correct against external re-derivation (not a fabricated lock).
- **Honest green:** full unfiltered suite **6579 passed, 0 failed, 0 skipped** (the one "Skipped" is an `[Explicit]` MFE benchmark, unrelated); `dotnet build` 0 errors. The 4 NUnit2007 warnings are pre-existing in `ApproximateMatcher_EditDistance_Tests.cs` (unrelated unit), unchanged by this session.

### Findings / defects
- **No code defect.** Implementation is a faithful, full-precision realisation of EMBOSS `sigcleave` / von Heijne (1986).
- **One Stage-B test-quality defect (M10 code-echo) — fixed in-session** by locking the externally-sourced score literal.

## Verdict & follow-ups
- **Stage A: PASS.** **Stage B: PASS-WITH-NOTES** (M10 strengthened). **End-state: ✅ CLEAN.**
- Test-quality gate: **PASS** after the M10 fix; full suite honest-green (6579/0).
- Follow-up: none. The min-length=15 simplification is documented and accepted.
