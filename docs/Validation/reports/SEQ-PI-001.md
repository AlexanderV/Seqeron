# Validation Report: SEQ-PI-001 — Isoelectric Point (pI) Calculation

- **Validated:** 2026-06-15   **Area:** Statistics
- **Canonical method(s):** `SequenceStatistics.CalculateIsoelectricPoint(string)` (private helper `NetCharge`)
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS

## Stage A — Description

### Sources opened this session (independent, not the repo's own artifacts)

1. **EMBOSS `iep` documentation** — fetched https://emboss.sourceforge.net/apps/release/6.6/emboss/apps/iep.html
   - Confirms pI = "the pH at which the numbers of positive and negative charges on the protein
     are equal" (net charge = 0), computed from amino-acid composition assuming no electrostatic
     coupling.
   - **Epk.dat pKa values, extracted verbatim:** Amino/N-terminus **8.6**, Carboxyl/C-terminus **3.6**,
     C **8.5**, D **3.9**, E **4.1**, H **6.5**, K **10.8**, R **12.5**, Y **10.1**.
2. **Peptides R package `charge_pI.cpp`** — fetched
   https://raw.githubusercontent.com/cran/Peptides/master/src/charge_pI.cpp
   - **Net-charge formula extracted:** basic groups (N-term, R, H, K) contribute
     `+1/(1+10^(pH−pKa))`; acidic groups (C-term, D, E, C, Y) contribute
     `−1/(1+10^(−(pH−pKa)))` = `−1/(1+10^(pKa−pH))`. Termini added once per chain.
   - EMBOSS-scale pKa values in this file match the EMBOSS doc exactly.

### Formula check
The Henderson–Hasselbalch net-charge model in the description matches the Peptides reference
source character-for-character (acidic `−1/(1+10^(pKa−pH))` is algebraically identical to the
source's `−1/(1+10^(−(pH−pKa)))`). The EMBOSS pKa scale, the [0,14] bisection window, and the
"net charge = 0" pI definition all match the EMBOSS `iep` documentation.

### Edge-case semantics
- **Empty/null → 7.0:** correctly flagged as a non-sourced input-guard *convention* (ASM-03);
  no literature defines pI for a zero-length protein. Acceptable as a documented sentinel.
- **Non-ionizable chars (gaps/whitespace/non-standard residues) ignored:** documented in §3.3;
  independently confirmed correct (such a string reduces to the termini-only pI).
- **Bounds 0 ≤ pI ≤ 14:** guaranteed by the bisection window (INV-01).
- **Composition-only (order-independent):** correct given the no-coupling assumption (INV-02).

### Independent cross-check (numbers retrieved/derived this session)
A standalone Python reference was written from the externally-confirmed EMBOSS pKa scale +
charge formula (not from the repo code). It **reproduces the Peptides published worked example
exactly**: net charge of `FLPVLAGLTPSIVPKLVCLLTKKC` = **3.037398 / 2.914112 / 0.7184524** at
pH 5/7/9 (6 dp). Using that validated charge function, the independent pI values are:

| Sequence | Independent pI |
|----------|----------------|
| FLPVLAGLTPSIVPKLVCLLTKKC | 9.67 |
| A / AG (termini only) | 6.10 |
| D | 3.75 |
| K | 9.70 |
| DDDD | 3.23 |
| KKKK | 11.27 |
| ACDEFGHIKLMNPQRSTVWY | 7.36 |
| DDDDDDDD | 2.96 |
| RRRRRRRR | 13.35 |
| DKDK / KDDK (permutation) | 6.45 (equal) |

INV-04 midpoint (8.6+3.6)/2 = 6.10 confirmed for termini-only sequences.

### Findings / divergences
None. The description's biology, formula, pKa scale, and edge semantics all trace to the
retrieved authoritative sources. **Stage A: PASS.**

## Stage B — Implementation

### Code path reviewed
`src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/SequenceStatistics.cs`
- `CalculateIsoelectricPoint` lines 295–326 (guard, residue counting, bisection over [0,14] to
  precision 0.01, round to 2 dp).
- `NetCharge` lines 334–350 (N-term `+1/(1+10^(pH−8.6))`, C-term `−1/(1+10^(3.6−pH))`, side
  chains by sign).
- pKa constants lines 256–269 — match EMBOSS Epk.dat exactly.

### Formula realised correctly?
Yes. The code computes the exact Henderson–Hasselbalch terms confirmed in Stage A (not an
approximation), with the EMBOSS pKa values verbatim. Bisection direction (charge>0 ⇒ search
higher pH) is correct for a monotonically non-increasing charge curve.

### Cross-verification table recomputed vs code
Every expected value asserted in the tests matches the independent reference table above
(A 6.10, D 3.75, K 9.70, DDDD 3.23, KKKK 11.27, all-20 7.36, ref-peptide 9.67, DDDDDDDD 2.96,
RRRRRRRR 13.35). Full suite re-run after edits: **6517 passed, 0 failed.**

### Variant/delegate consistency
Single public entry point; `NetCharge` is the only delegate and is exercised transitively by
every pI assertion. No `*Fast`/instance duplicate exists.

### Test quality audit (HARD gate)
- **Sourced expectations:** all exact pI values now trace to the independent reference built
  from external sources, not to code output.
- **Defect found & fixed:** the original M2 test (`...StaysWithinPhBounds`) asserted only
  `[0,14]` bounds on `DDDDDDDD`/`RRRRRRRR` — it would pass against any implementation returning
  anything in range (a code-echo weakness per the gate). **Fixed:** rewritten to assert the
  exact sourced values **2.96** and **13.35** (bounds kept as additional invariant checks).
- **Coverage gap fixed:** the documented "non-ionizable characters ignored / no throw" path
  (§3.3) had no test. **Added** `..._NonIonizableCharacters_IgnoredEqualsTerminiOnly` asserting
  `"A B!G"` and `"XZ"` → 6.10.
- **No green-washing:** no assertions weakened, no tolerances widened, no tests skipped. The
  one remaining range assertion (INV-01 bounds) is a genuine invariant, retained *alongside*
  the new exact-value checks.
- **Coverage:** the single public method and all Stage-A branches/edge cases (basic, acidic,
  termini-only, all-20, empty, null, lowercase, order-independence, bounds, non-ionizable
  noise) are exercised with exact sourced values.
- **Honest green:** full unfiltered suite `Failed: 0, Passed: 6517`; build `0 Error(s)`,
  no new warnings.

### Findings / defects
- Stage-B defect (test quality): weak bounds-only M2 assertion — **FIXED this session**.
- Coverage gap: untested non-ionizable-character path — **FIXED this session** (new test).
- No implementation defect found; the code is correct against the validated description.

## Verdict & follow-ups
- **Stage A: PASS.** **Stage B: PASS.**
- **End-state: CLEAN** — algorithm fully functional; the two test-quality issues found were
  completely fixed (exact sourced values locked in, missing edge-case covered), build + full
  suite green.
- Test-quality gate: **PASS** (after fixes).
