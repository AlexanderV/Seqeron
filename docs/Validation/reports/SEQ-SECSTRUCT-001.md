# Validation Report: SEQ-SECSTRUCT-001 — Protein Secondary Structure Prediction (Chou-Fasman propensity profile)

- **Validated:** 2026-06-15   **Area:** Statistics
- **Canonical method(s):** `SequenceStatistics.PredictSecondaryStructure(string proteinSequence, int windowSize = 7)`
- **Stage A verdict:** PASS-WITH-NOTES
- **Stage B verdict:** PASS

## Stage A — Description

### Sources opened this session (not the repo's own claims)

| # | Source | Retrieved | What it confirmed |
|---|--------|-----------|-------------------|
| 1 | ravihansa3000/ChouFasman `ChouFasman.py` (raw) | curl, this session | Full 20-residue Pa/Pb/Pt integer table (×100). Verbatim values below. |
| 2 | Wikipedia — Chou–Fasman method | WebFetch | Method definition; helix formers (Ala, Glu, Leu, Met); nucleation rules 4-of-6 (helix, cutoff 1.03) and 3-of-5 (sheet, cutoff 1.00). No numeric table on the page. |
| 3 | CSB\|SJU (Jakubowski) propensity table | WebFetch | Pa/Pb for all 20 residues; independently confirms 19 of 20 helix/sheet values. Lists **Lys Pα = 1.16**. |
| 4 | SlideShare "Chou-Fasman algorithm" (partial table) | WebFetch | Independently confirms E 1.51/0.37/0.74, M 1.45/1.05/0.60, A 1.42/0.83/0.66, V 1.06/1.70/0.50, Y 0.69/1.47/1.14, P 0.57/0.55/1.52, G 0.57/0.75/1.56. |

### Formula check

The method computes, per sliding window of size *w* stepping by 1 (N→C), the per-component
arithmetic mean of the Chou-Fasman propensities Pα/Pβ/Pt over the known residues in the
window. This is the windowed mean-propensity used during Chou-Fasman nucleation-region
evaluation (Kelley lecture STEP 3: "compute mean P(a) and P(b) over the region"). The
implementation is explicitly a *continuous profile*, not the full nucleation/extension/turn
state machine — correctly scoped and documented as such (algorithm doc §5.3).

### Propensity table — independently verified

Raw values keyed by **one-letter code** from source #1 (the `.py` file mislabels the long
names "Asparagine"/"Aspartic Acid" but the one-letter codes are correct):

```
A 142/83/66   R 98/93/95    N 67/89/156   D 101/54/146  C 70/119/119
E 151/37/74   Q 111/110/98  G 57/75/156   H 100/87/95   I 108/160/47
L 121/130/59  K 114/74/101  M 145/105/60  F 113/138/60  P 57/55/152
S 77/75/143   T 83/119/96   W 108/137/96  Y 69/147/114  V 106/170/50
```

Dividing by 100, this matches the implementation's `SecondaryStructurePropensity` dictionary
**exactly for all 20 residues**. Source #3 independently confirms 19/20 helix+sheet values;
source #4 independently confirms 7 residues across all three columns.

### Edge-case semantics

- Non-standard residues (X/B/Z/gaps): no defined propensity (table covers exactly the 20
  standard residues) → excluded from window mean (sources #1, #3). This is an ASSUMPTION for
  the *averaging* behaviour (no source mandates skip-vs-error), documented and deterministic.
- Window > sequence length → no scan positions (window-vs-length, Kelley lecture). Sourced.
- Null/empty, w < 1 → empty result. Precondition contract; reasonable.

### Independent cross-check (numbers)

- A → (1.42, 0.83, 0.66); E → (1.51, 0.37, 0.74); V → (1.06, 1.70, 0.50): all from sources #1, #4.
- K → (1.14, 0.74, 1.01): source #1 lists `['K', 114, 74, 101]` verbatim.
- "AE" window 2: helix (1.42+1.51)/2 = 1.465; sheet (0.83+0.37)/2 = 0.60; turn (0.66+0.74)/2 = 0.70. Hand-computed.

### Findings / divergences (PASS-WITH-NOTES)

1. **Lysine Pα = 1.14 vs 1.16 (documented literature conflict).** Source #1 (reference impl)
   gives 1.14 and the Przytycka NCBI lecture gives 1.14; source #3 (CSB\|SJU) gives 1.16. The
   implementation and tests adopt **1.14** on a two-source majority. This is a genuine
   divergence in published reproductions of the 1978 Annu Rev Biochem table, not an error;
   it shifts only window means containing K. Documented in Evidence Assumption 1 and doc §5.4.
   Verdict reflects this as a note, not a fail.
2. **Default windowSize = 7** is an API convenience, not a Chou-Fasman constant (nucleation
   windows are 6/5). Documented; tests pass the window explicitly. Acceptable.
3. Turn (Pt) values vary slightly across reproductions (e.g. Ile Pt = 0.47 in source #1 vs
   0.50 in source #4); the implementation uses the source-#1 set consistently. Minor.

## Stage B — Implementation

### Code path reviewed

`src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/SequenceStatistics.cs:806-870`.

### Formula realised correctly?

Yes. `PredictSecondaryStructure` (lines 840-870):
- Guards null/empty, `windowSize < 1`, `windowSize > length` → `yield break` (line 844). INV-06.
- Uppercases input (line 847). INV-04.
- For each start i in `[0, n-w]`, sums Pα/Pβ/Pt over the *known* residues (TryGetValue skips
  unknowns, line 856), counts known ones, and emits `(sum/count)` only when `count > 0`
  (lines 865-868). Realises INV-01/02/03/05 exactly. No precision loss (doubles; counts small).

### Cross-verification table recomputed vs code (full suite run)

| Input (w) | Expected (sourced) | Source |
|-----------|--------------------|--------|
| "A" (1) | (1.42, 0.83, 0.66) | #1, #4 |
| "E" (1) | (1.51, 0.37, 0.74) | #1, #4 |
| "V" (1) | (1.06, 1.70, 0.50) | #1, #4 |
| "K" (1) | (1.14, 0.74, 1.01) — NOT 1.16 | #1 (114/74/101) |
| "AE" (2) | (1.465, 0.60, 0.70) | hand-computed from #1 |
| "AXE" (3) | (1.465, 0.60, 0.70) — X excluded | #1 + skip rule |
| "XBZ" (3) | empty | skip rule |

All pass against the actual code (full suite: Failed 0).

### Variant/delegate consistency

Single public method; no `*Fast`/instance variants. N/A.

### Test quality audit (HARD gate)

- **Sourced expectations, not code echoes:** M1-M4 assert exact literal tuples taken from the
  external propensity table (A 1.42/0.83/0.66, E, V, K 1.14/0.74/1.01). A wrong impl returning
  Lys Pα = 1.16 would fail M4; a wrong table entry fails M1-M4. Not tautological.
- **No green-washing:** all means asserted as exact equalities `Within(1e-10)`. C1/C2 also
  carry a `GreaterThan` sanity check but each is *anchored by an exact-mean equality*, so the
  inequality is additive, not a weakening. No skipped/ignored/widened tests.
- **Cover all logic:** every Stage-A branch exercised — single-residue (INV-01), multi-residue
  mean (INV-02), slide+count (INV-03), case-insensitivity (INV-04), unknown-exclude + all-unknown
  (INV-05), and all four INV-06 error/edge paths (null, empty, window>len, window<1).
- **Honest green:** FULL unfiltered suite = **Failed: 0, Passed: 6523**; `dotnet build` 0 errors.
  The SECSTRUCT source/test files are warning-free (the 4 build warnings are pre-existing in the
  unrelated `ApproximateMatcher_EditDistance_Tests.cs`, untouched this session).

**Gate result: PASS.** No changes to code or tests were required; the existing tests already
lock the externally-sourced values (including the conflict-resolved Lys 1.14). Minor coverage
note: residues C/D/F/G/H/N/P/Q/R/S/T/W are not individually unit-tested, but the table was
verified verbatim against the external reference for all 20, and the formula paths are fully
covered — not a defect.

## Verdict & follow-ups

- **Stage A:** PASS-WITH-NOTES (Lys Pα 1.14-vs-1.16 literature conflict, resolved & documented;
  default window 7 is API-only; minor Pt-column variation across reproductions).
- **Stage B:** PASS (code realises the windowed-mean formula exactly; table matches the external
  reference verbatim for all 20 residues; all edge cases handled and tested).
- **End-state: CLEAN** — no defect found; algorithm fully functional. No code/test fixes needed.
- **Test-quality gate: PASS** — exact sourced expectations, all invariants/edges covered, full
  unfiltered suite green (6523/0).
