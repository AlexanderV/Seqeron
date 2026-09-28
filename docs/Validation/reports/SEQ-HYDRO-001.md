# Validation Report: SEQ-HYDRO-001 — Hydrophobicity Analysis (Kyte-Doolittle GRAVY + sliding-window profile)

- **Validated:** 2026-06-15   **Area:** Statistics
- **Canonical method(s):** `SequenceStatistics.CalculateHydrophobicity(string)`, `SequenceStatistics.CalculateHydrophobicityProfile(string, int windowSize=9)`
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS-WITH-NOTES

## Stage A — Description

### Sources opened this session (independent of the repo)

| # | Source | URL | What it confirmed |
|---|--------|-----|-------------------|
| 1 | Biopython `ProtParamData.py` (kd dict) | https://raw.githubusercontent.com/biopython/biopython/master/Bio/SeqUtils/ProtParamData.py | Full kd dict verbatim; attributed to "J. Mol. Biol. 157:105-132 (1982)" |
| 2 | Expasy ProtParam doc (GRAVY) | https://web.expasy.org/protparam/protparam-doc.html | GRAVY = "sum of hydropathy values of all the amino acids, divided by the number of residues"; scale = Kyte & Doolittle 1982 |
| 3 | Expasy ProtScale `Hphob.Doolittle.html` | https://web.expasy.org/protscale/pscale/Hphob.Doolittle.html | All 20 KD values, independent rendering of the 1982 table |
| 4 | Biopython `ProtParam.py` (gravy, protein_scale) | https://raw.githubusercontent.com/biopython/biopython/master/Bio/SeqUtils/ProtParam.py | `gravy = total_gravy/length`; `protein_scale` loops `range(N-W+1)`; `edge=1.0` ⇒ equal per-position weight |

### Formula check

- **GRAVY = Σ kd(s_i) / n** — confirmed verbatim by Expasy ProtParam doc (src 2) and Biopython `gravy()` `total_gravy / self.length` (src 4).
- **Profile** returns **N − W + 1** values; with default `edge=1.0` Biopython weights every window position equally, i.e. an unweighted window mean (src 4). For odd windows `sum_of_weights = window`, so the divisor is exactly W — matching the description.

### Kyte-Doolittle scale — two-source cross-check (every value retrieved this session)

| AA | Impl | Biopython (src 1) | Expasy ProtScale (src 3) | AA | Impl | src 1 | src 3 |
|----|------|------|------|----|------|------|------|
| A | 1.8 | 1.8 | 1.800 | M | 1.9 | 1.9 | 1.900 |
| R | -4.5 | -4.5 | -4.500 | F | 2.8 | 2.8 | 2.800 |
| N | -3.5 | -3.5 | -3.500 | P | -1.6 | -1.6 | -1.600 |
| D | -3.5 | -3.5 | -3.500 | S | -0.8 | -0.8 | -0.800 |
| C | 2.5 | 2.5 | 2.500 | T | -0.7 | -0.7 | -0.700 |
| Q | -3.5 | -3.5 | -3.500 | W | -0.9 | -0.9 | -0.900 |
| E | -3.5 | -3.5 | -3.500 | Y | -1.3 | -1.3 | -1.300 |
| G | -0.4 | -0.4 | -0.400 | V | 4.2 | 4.2 | 4.200 |
| H | -3.2 | -3.2 | -3.200 | I | 4.5 | 4.5 | 4.500 |
| L | 3.8 | 3.8 | 3.800 | K | -3.9 | -3.9 | -3.900 |

All 20 values match both independent sources exactly.

### Edge-case semantics

- **W > N → empty profile** — sourced from Biopython `range(N-W+1)` = 0 iterations (src 4). Correct.
- **Empty/null → GRAVY 0, empty profile** — a defined library contract (sources are silent; not a scoring constant). Reasonable.
- **Non-standard residues (B/Z/X/gaps)** — undefined by KD/Expasy (only 20 canonical). Biopython raises `KeyError`; this library instead **skips** them in GRAVY (divides by recognized count) and treats them as **0** in profile windows. Documented deviation; does not alter any canonical value. Acceptable per protocol (undefined input, sourced behavior absent).

### Independent cross-check (hand-computed numbers)

- FLIV GRAVY = (2.8+3.8+4.5+4.2)/4 = 15.3/4 = **3.825** ✓
- RKDE GRAVY = (−4.5−3.9−3.5−3.5)/4 = −15.4/4 = **−3.85** ✓
- FLIV profile W=3 = [(2.8+3.8+4.5)/3, (3.8+4.5+4.2)/3] = [**3.7**, **4.16666…**] ✓
- I×19 W=19 = 4.5 (> 1.6 transmembrane threshold) ✓

### Stage A findings

Description is biologically and mathematically correct and fully sourced. **PASS.**

## Stage B — Implementation

- **Code path:** `src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/SequenceStatistics.cs:354-425`.
- **GRAVY realised correctly:** uppercases input, sums kd over recognized residues, divides by recognized count, returns 0 when count==0 or input null/empty. Matches Σkd/n for canonical input.
- **Profile realised correctly:** guards `windowSize > length` and null/empty → `yield break`; otherwise loops `i = 0 … N-W` (N−W+1 iterations) and yields `sum/windowSize` per window. Lazy via `yield return`.

### Cross-verification table recomputed vs code (full suite run)

| Case | Sourced expected | Code result |
|------|------------------|-------------|
| GRAVY "A" | 1.8 | 1.8 ✓ |
| GRAVY FLIV | 3.825 | 3.825 ✓ |
| GRAVY RKDE | −3.85 | −3.85 ✓ |
| GRAVY "fliv" (case) | 3.825 | 3.825 ✓ |
| GRAVY ""/null | 0 | 0 ✓ |
| Profile FLIV W=3 | [3.7, 4.16667] | match ✓ |
| Profile AG W=3 (W>N) | empty | empty ✓ |
| Profile I×19 W=19 | 4.5 | 4.5 ✓ |

### Variant/delegate consistency

`CalculateHydrophobicity` is also called by the aggregate stats record (`SequenceStatistics.cs:118`) which simply forwards the GRAVY value — consistent.

### Test quality audit (gate)

- **Sourced, not code-echoed:** every expected value (kd constants, GRAVY sums, window means) traces to Biopython/Expasy retrieved this session, not to code output. A deliberately-wrong kd value or a divide-by-length-vs-count bug would fail M2/M3/S1.
- **No green-washing:** exact `Is.EqualTo(...).Within(1e-10)` assertions; no ranges/Greater/Contains where an exact value exists; no widened tolerances; nothing skipped.
- **Coverage gaps found and fixed this session (test-only, no code change):**
  - `CalculateHydrophobicity_AllUnknownResidues_ReturnsZero` ("XXZB" → 0) — exercises the `count==0` branch for a non-empty all-unknown string (distinct from the null/empty short-circuit); guards against NaN (0/0).
  - `CalculateHydrophobicityProfile_UnknownResidueInWindow_ContributesZero` ("FXIV" W=3 → [7.3/3, 8.7/3]) — exercises the profile's `TryGetValue` skip-to-0 branch and confirms the divisor is W (not recognized count), locking the documented deviation 5.4.
- **Result:** full unfiltered suite **6519 passed, 0 failed, 0 skipped-as-defect**; `dotnet build` 0 errors; changed test file builds warning-free (the 4 build warnings are pre-existing NUnit2007 in unrelated files).

### Stage B findings / notes

1. **Even-window divisor divergence (minor, out of contract).** Biopython `protein_scale` divides by `sum_of_weights = 2·(W//2)+1`, which for **even** W equals W+1 (slight down-weight), whereas this library always divides by W. For all standard/odd windows (3, 9, 19 — the only ones used in KD analysis and in every test) the two agree exactly. The library docs only ever recommend odd windows, so this is outside the documented contract. Recorded as a note, not a defect.
2. **Non-standard residue handling** diverges from Biopython `KeyError` (skips / treats as 0) — already documented deviation 5.4; now locked by the two added tests.

## Verdict & follow-ups

- **Stage A: PASS** — formula and all 20 scale constants confirmed against two independent authoritative sources.
- **Stage B: PASS-WITH-NOTES** — code faithfully realises GRAVY and the odd-window profile; two minor documented divergences (even-window divisor; non-standard residue handling) do not affect any canonical value. Test coverage was strengthened (+2 tests) to lock the previously-untested `count==0` and profile-unknown branches.
- **End-state: CLEAN** — no correctness defect; coverage gaps fully fixed in-session; full suite green.
- No FINDINGS_REGISTER defect logged (no behavioural defect; only test-coverage strengthening).

---

## Re-review 2026-09 (campaign review-2026-09, batch B03)

- **Stage A:** PASS-WITH-NOTES — GRAVY and the KD scale re-confirmed against Biopython 1.88 (installed `Bio.SeqUtils.ProtParam` source inspected: `gravy`, `protein_scale`, `_weight_list`). The ExPASy ProtScale page itself is egress-blocked; its window/edge conventions (odd window centred on residue i, centre weight 100 %, edges set 0–100 %, linear or exponential variation) were taken from the Biopython docstring ("Similar to expasy's ProtScale") and WebSearch snippets of the ProtScale documentation.
- **Stage B:** FAIL → FIXED.
  1. `CalculateHydrophobicityProfile(seq, 0)` returned `NaN` for every position (0/0); negative windows returned garbage (N+|W|+1 values of −0). Now `ArgumentOutOfRangeException` for W < 1, validated eagerly (the method was a lazy iterator). Biopython raises `IndexError` for W = 0.
  2. Simplification removed: the profile had no edge weighting although ProtScale/Biopython `protein_scale(window, edge)` define it. Added `edgeWeight` (default 1.0 → bit-identical to the previous unweighted mean) implementing the linear weight model; confirmed against Biopython 1.88 on ubiquitin (W = 5/9/19, edge 0.4/0.4/0.1), a signal peptide, FLIVAG (edge 0 and 0.4) — all ≤ 1e-10. Edge < 1 with an even window is rejected (no centre residue).
  3. Notes (by design, documented): even window with edge 1 = plain mean (Biopython double-counts residue W/2 and divides by W+1: "AV" W=2 → 3.4 vs 3.0); off-centre unknown residue contributes 0 without dropping its symmetric partner (Biopython drops both).
- **Duplication:** `SequenceStatistics` holds the canonical KD table/GRAVY (DisorderPredictor already delegates). Remaining copies: `ProteinSequence.Gravy()` (Core, rounds to 3 dp, cannot reference Analysis) and `ProteinMotifFinder.HydropathyScale`/`CalculateHydropathyProfile` (divides by recognised-residue count — not Biopython behaviour) — cross-batch requests (B02, B14).
- **Tests:** +13 cases in `SequenceStatistics_CalculateHydrophobicity_Tests.cs` (Biopython-sourced).
- **End-state:** CLEAN.
