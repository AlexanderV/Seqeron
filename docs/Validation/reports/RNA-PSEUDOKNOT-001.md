# Validation Report: RNA-PSEUDOKNOT-001 — Pseudoknot Detection (crossing base pairs)

- **Validated:** 2026-06-16   **Area:** RnaStructure
- **Canonical method(s):** `RnaSecondaryStructure.DetectPseudoknots(IReadOnlyList<BasePair>)`
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS

## Stage A — Description

### Sources opened this session
1. **Wikipedia "Pseudoknot"** (https://en.wikipedia.org/wiki/Pseudoknot) — retrieved via WebFetch.
   Verbatim: *"A pseudoknot is a nucleic acid secondary structure containing at least two stem-loop
   structures in which half of one stem is intercalated between the two halves of another stem."* and
   *"The base pairing in pseudoknots is not well nested; that is, base pairs occur that 'overlap' one
   another in sequence position."* H-type confirmed as the simplest class.
2. **Staple & Butcher, "RNA pseudoknots: folding and finding"** (PMC2873773) — retrieved via WebFetch.
   Verbatim formal condition: *"for any base pairs i-j and k-l, i<j, k<l, and i<k, there are cases in
   which i<k<j<l"*; this crossing condition is *avoided* in regular nested secondary structure.
   Confirms the H-type / classical three-loop (Stem 1, Stem 2; Loop 1/2/3) configuration.
3. **WebSearch cross-check** (DP-algorithm and counting papers, arXiv/PMC) independently restated the
   nesting convention: a non-pseudoknotted structure obeys *"either i < k < l < j or k < i < j < l"*,
   i.e. crossing (`i<k<j<l`) is exactly the complement of nested/disjoint.

### Formula check
The crossing condition in the TestSpec/Evidence/doc — two base pairs (i,j) and (k,l) with open<close
cross iff **i < k < j < l** — matches **verbatim** the formal definitions from sources 2 and 3.
- Nested: `i < k < l < j` (fails `j < l`) — confirmed by source 3's nesting convention.
- Disjoint: `j < k` (fails `k < j`) — derived from overlap requirement, consistent with all sources.

### Independent cross-check (numbers, all hand-computed against sourced condition)
| Case | Pairs | i,k,j,l after open<close + open-first ordering | i<k<j<l? | Expected |
|------|-------|-----------------------------------------------|----------|----------|
| H-type `([)]` | (0,2),(1,3) | 0<1<2<3 | yes | 1 pseudoknot |
| Nested | (0,5),(1,4) | 0<1, j=5, l=4 → 5<4 false | no | 0 |
| Disjoint | (0,2),(3,5) | k=3, j=2 → 3<2 false | no | 0 |
| Order-2 `([{)]}` | (0,3),(1,4),(2,5) | each pair: 0<1<3<4, 0<2<3<5, 1<2<4<5 | all yes | 3 (pairwise) |

The `([)]` H-type and the nested/disjoint controls match the Antczak (2018) DBL example and the
Wikipedia "not well nested" wording. The condition is symmetric (crossing is a relation between two
pairs), so input/storage order independence is a genuine mathematical invariant.

### Edge-case semantics
- Empty / single / `null` → no pseudoknot: a crossing needs two pairs. This is a *derived* consequence
  of the definition (not an invented free parameter); sound.
- Degenerate pair (i=j) cannot satisfy strict `<` chain → never part of a pseudoknot. Sound.

### Findings / divergences
None. The description (TestSpec, Evidence, algorithm doc) is biologically and mathematically correct
and every non-trivial claim traces to a source retrieved this session. **Stage A = PASS.**

## Stage B — Implementation

### Code path reviewed
`src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/RnaSecondaryStructure.cs:1562–1600`
(`DetectPseudoknots`), plus the `BasePair` (line 21) and `Pseudoknot` (line 81) record structs.

### Formula realised correctly?
Yes. The code:
- returns empty for `null` or `< 2` pairs (line 1565);
- O(n²) all-pairs scan (lines 1568–1570);
- normalizes each pair to (open<close) via `Math.Min/Max` (lines 1576–1579) — handles reversed storage;
- reorders the two pairs so the first-opening pair plays (i,j) (lines 1583–1586) — makes the test
  direction/order-independent;
- emits when `i < k && k < j && j < l` (line 1589), i.e. exactly `i<k<j<l`;
- reports `Start1<End1`, `Start2<End2`, `Start1<Start2` with the two original `BasePair`s.

This is the validated formula verbatim — no approximation, no thermodynamics. Lazy `yield`, O(1) extra.

### Cross-verification table recomputed vs code (full suite run)
All four hand-computed cases above were exercised by tests and pass (M1, M2, M3, and the new S5
order-2 case). Reversed-endpoint (M4) and reported-coordinate (M5) cases confirm normalization and
output shape. Mixed (S3) and order-independence (S4) confirm INV-01..05.

### Variant/delegate consistency
Single public canonical method; no `*Fast`/delegate variants. `AnalyzeStructure` (line 1503) calls the
same `DetectPseudoknots`, so the structure-level path is consistent.

### Test quality audit (HARD gate)
- **Sourced, not code-echoed:** expected counts/coordinates come from the hand-computed crossing
  condition, not from running the code. The `([)]`, nested, disjoint, and `([{)]}` numbers all trace
  to sources 1–3.
- **No green-washing:** assertions use exact `EqualTo` counts and exact coordinate tuples; no
  `Greater`/`AtLeast`/range/tolerance weakening; no skips/ignores.
- **Coverage:** every documented branch is hit — crossing (M1), nested (M2), disjoint (M3),
  normalization (M4), output coords (M5), empty (S1), single (S2), null (S2b), mixed (S3), order
  independence (S4), and the O(n²) invariant property (C1).
- **Gap found & fixed:** the documented contract "each crossing pair-of-pairs reported separately /
  no order grouping" (algorithm doc Deviation 1, TestSpec Assumption 2) was **not locked by any test**
  — C1 only checked the per-result invariant, never that N mutually-crossing pairs yield C(N,2)
  results. Added **S5** (`DetectPseudoknots_ThreeMutuallyCrossingPairs_ReportsEachPairwiseCrossing`)
  asserting `([{)]}` → exactly 3 reported crossings with the exact tuple set
  `{(0,3,1,4),(0,3,2,5),(1,4,2,5)}`. This expectation is sourced (Antczak crossing condition applied
  pairwise), not a code echo.
- **Honest green:** full unfiltered suite = **6592 passed, 0 failed, 0 skipped of interest**
  (MFE benchmark intentionally skipped, unrelated). Changed test file builds warning-free.

### Findings / defects
No implementation defect. One test-coverage gap (documented separate-reporting contract) — **fixed**
in this session by adding test S5.

## Verdict & follow-ups
- **Stage A: PASS** — description matches authoritative sources verbatim (crossing `i<k<j<l`).
- **Stage B: PASS** — code realises the validated condition exactly; tests strengthened to lock the
  documented multi-crossing reporting contract.
- **End-state: CLEAN.** No code change needed; one test added; build 0 errors, full suite Failed: 0.
- Out of scope (documented Not Implemented, not a defect): pseudoknot *order* / DBL layer grouping
  (Antczak 2018) — this unit reports pairwise crossings only.
