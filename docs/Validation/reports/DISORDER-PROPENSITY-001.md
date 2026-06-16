# Validation Report: DISORDER-PROPENSITY-001 — Disorder Propensity (TOP-IDP lookup + Dunker classification)

- **Validated:** 2026-06-16   **Area:** ProteinPred
- **Canonical method(s):** `DisorderPredictor.GetDisorderPropensity(char)`, `IsDisorderPromoting(char)`,
  `DisorderPromotingAminoAcids`, `OrderPromotingAminoAcids`, `AmbiguousAminoAcids`
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS
- **End-state:** ✅ CLEAN

## Scope note

The session prompt described this unit as a per-residue disorder propensity with **sliding-window
smoothing** producing a profile / mean disorder. That windowing/profile logic is **not** part of this
unit — it lives in `DisorderPredictor.PredictDisorder(...)` under **DISORDER-PRED-001**. Per the
Registry (`ALGORITHMS_CHECKLIST_V2.md` §DISORDER-PROPENSITY-001, lines 3658-3676) and the Method Index
(lines 5070-5073), DISORDER-PROPENSITY-001 is scoped to the **O(1)** primitives: the single-character
TOP-IDP scale lookup (Campen 2008) and the Dunker (2001) order/disorder/ambiguous amino-acid
classification. This report validates the unit **as registered** (the windowing belongs to its own
unit and is out of scope here).

## Stage A — Description

### Sources opened this session

1. **Campen et al. (2008), TOP-IDP-Scale, PMC2676888 full text** — `WebFetch`
   (https://pmc.ncbi.nlm.nih.gov/articles/PMC2676888/). Extracted Table 2 verbatim (all 20 values),
   the rank string, the prediction cut-off, and the extreme anchors.
2. **Wikipedia — Intrinsically disordered proteins** (cites Dunker et al. 2001 primary) — `WebFetch`
   (https://en.wikipedia.org/wiki/Intrinsically_disordered_proteins). Extracted the three Dunker
   classification sets.
3. **localCIDER (Pappu lab)** — `git clone` of https://github.com/Pappulab/localCIDER; confirmed it
   cites Campen 2008 for the TOP-IDP scale (`sequenceParameters.py:221`) but does not embed the raw
   numeric table (it computes derived sequence parameters), so it corroborates the citation, not the
   per-residue numbers.
4. **PubMed 18991772 / USF DigitalCommons / Bentham** — confirmed the rank string
   `W,F,Y,I,M,L,V,N,C,T,A,G,R,D,H,Q,K,S,E,P` independently.

### Formula / data check (vs Campen 2008 Table 2, fetched this session)

| AA | Source value | Code value | AA | Source value | Code value |
|----|------|------|----|------|------|
| W | -0.884 | -0.884 | A | 0.06  | 0.060 |
| F | -0.697 | -0.697 | G | 0.166 | 0.166 |
| Y | -0.510 | -0.510 | R | 0.180 | 0.180 |
| I | -0.486 | -0.486 | D | 0.192 | 0.192 |
| M | -0.397 | -0.397 | H | 0.303 | 0.303 |
| L | -0.326 | -0.326 | Q | 0.318 | 0.318 |
| V | -0.121 | -0.121 | S | 0.341 | 0.341 |
| N | 0.007  | 0.007  | K | 0.586 | 0.586 |
| C | 0.02   | 0.020  | E | 0.736 | 0.736 |
| T | 0.059  | 0.059  | P | 0.987 | 0.987 |

All 20 match exactly. Min = **W = -0.884**, Max = **P = 0.987** — confirmed by the fetched source.

### Dunker (2001) classification check (vs Wikipedia, fetched this session)

| Class | Source | Code |
|-------|--------|------|
| Disorder-promoting (8) | A,R,G,Q,S,P,E,K | {A,R,G,Q,S,P,E,K} ✓ |
| Order-promoting (8) | W,C,F,I,Y,V,L,N | {W,C,F,I,Y,V,L,N} ✓ |
| Ambiguous (4) | H,M,T,D | {H,M,T,D} ✓ |

Union = 20 distinct standard residues; the three sets are pairwise disjoint (verified by hand from the
fetched residues) — invariant INV-4 holds.

### Documented pitfall — confirmed and correctly handled

The Campen/Wikipedia **rendered rank string** places "…Q, **K, S**, E, P", but the Table 2 **numeric
values** give S = 0.341 < K = 0.586 (so by value: …Q, **S, K**, E, P). This is a presentation-order
artifact in the paper's rank string; the per-residue Table 2 values are authoritative. The
implementation stores the numeric values (S=0.341, K=0.586) and the tests lock both — correct.

### Edge-case semantics

- Scale defined only for the 20 standard residues (Table 2); no value for B/J/O/U/X/Z or gaps.
- Unknown residue → 0.0 is an **implementation contract** (`GetValueOrDefault(...,0)`), not a
  source-defined value — honestly registered as an assumption (TestSpec §6 #1), tested as a contract.
- Case-insensitivity (upper-case before lookup) is an implementation contract.

**Stage A verdict: PASS.** Every non-trivial expected value traces to a source fetched this session
(Campen 2008 PMC2676888 Table 2; Wikipedia/Dunker 2001). No biological or mathematical error.

## Stage B — Implementation

### Code path reviewed

`src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/DisorderPredictor.cs`:
- `DisorderPropensity` dictionary (lines 85-107) — matches Table 2 verbatim (table above).
- `DisorderPromotingSet`/`OrderPromotingSet`/`AmbiguousSet` (lines 111-121) — match Dunker sets.
- `GetDisorderPropensity` (680-683) — `GetValueOrDefault(ToUpperInvariant(aa), 0)`; realises lookup
  + case-insensitivity + unknown→0.0 contract.
- `IsDisorderPromoting` (691-694) — `DisorderPromotingSet.Contains(ToUpperInvariant(aa))`.
- `DisorderPromotingAminoAcids`/`OrderPromotingAminoAcids`/`AmbiguousAminoAcids` (700-712) — return
  pre-sorted cached lists (`OrderBy(c => c)`), so they are ascending and stable.

### Cross-verification vs code

Ran the full unfiltered suite (below): the 14 tests in the canonical fixture recompute the 20 Table-2
values, the W/P anchors, the three Dunker sets (exact equivalence + counts), disjoint/cover-20, the
membership⇔predicate invariant, unknown→0.0, and case-insensitivity — all green against the actual
code.

### Variant / delegate consistency

`IsDisorderPromoting` ⇔ membership in `DisorderPromotingAminoAcids` (M10) — consistent. The property
getters expose the same underlying sets as the predicate/lookup.

### Test-quality audit (HARD gate)

- **Sourced, not code-echoes:** the fixture hardcodes the expected scale (`TopIdp`) and the three sets
  (`DisorderPromoting`/`OrderPromoting`/`Ambiguous`) as **independent literals copied from the
  published sources** (verified this session to match Campen Table 2 / Dunker). M1/M2 assert exact
  values; M6/M7/M8 assert exact set equivalence + exact counts. A wrong scale or a swapped set would
  fail. M4 (order→false) and M5 (ambiguous→false) plus M3 (disorder→true) exercise **both** branches
  of the predicate, so neither an always-true nor an always-false impl survives.
- **No green-washing:** exact `Is.EqualTo(...).Within(1e-10)` on values; `Is.EquivalentTo(set)` +
  `Count` on sets; no `Greater/AtLeast/Contains` where an exact value is known; no widened tolerance;
  no skipped/ignored tests.
- **Coverage:** all 5 public members exercised; both predicate branches; the documented edge/contract
  cases (unknown residue X/Z/B/`*`→0.0; lowercase value & predicate; sorted property order). The S<K
  rank-vs-value pitfall is locked because M1 asserts S=0.341 and K=0.586 exactly.
- **Honest green:** full **unfiltered** suite run, not a `--filter` subset.

### Build & full suite (this session)

- `dotnet build …Tests.csproj -c Debug` → **0 errors**, 4 warnings (pre-existing NUnit2007 in the
  unrelated `ApproximateMatcher_EditDistance_Tests.cs`; no file changed this session).
- `dotnet test … --no-build` → **Failed: 0, Passed: 6609, Skipped: 0**.

**Stage B verdict: PASS.** Code faithfully realises the validated description; tests assert exact
sourced values and cover all branches/edge cases.

## Verdict & follow-ups

- **Stage A: PASS · Stage B: PASS · End-state: ✅ CLEAN.** No defect found; no code/test/spec change
  required. The unit is fully functional.
- **Test-quality gate: PASS** — sourced expectations (not echoes), no green-washing, full coverage of
  all 5 members + branches + documented edges, honest green on the full unfiltered suite (6609 passed,
  0 failed).
- Note recorded: the prompt's "sliding-window smoothing / profile / mean disorder" description applies
  to **DISORDER-PRED-001** (`PredictDisorder`), not this unit; validated as registered (O(1)
  lookup + classification).
