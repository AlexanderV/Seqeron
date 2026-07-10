---
type: source
title: "Validation report: DISORDER-PROPENSITY-001 (per-residue disorder propensity — raw TOP-IDP Table-2 lookup + Dunker classification primitives)"
tags: [validation, analysis]
doc_path: docs/Validation/reports/DISORDER-PROPENSITY-001.md
sources:
  - docs/Validation/reports/DISORDER-PROPENSITY-001.md
source_commit: 540cb0dff42f0ca3f526ff0743d2c31e144ff4d6
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Validation report: DISORDER-PROPENSITY-001

The two-stage **validation write-up** for test unit **DISORDER-PROPENSITY-001** — the **O(1)
per-residue propensity primitives** beneath the sliding-window predictor: the single-character
**raw TOP-IDP scale lookup** (Campen et al. 2008) and the **Dunker (2001) order/disorder/ambiguous
classification** — validated 2026-06-16. This is the *report* artifact that feeds one row of the
[[validation-ledger]]; it records the validator's independent **verdict** on both the algorithm
description (Stage A) and the shipped code (Stage B), and the wider campaign is
[[validation-and-testing]]. The scale, the Dunker sets, the raw-vs-normalized value-space distinction
and the S/K ranking pitfall are synthesized on the concept
[[intrinsic-disorder-prediction-top-idp]] (the shared `PredictDisorder` anchor these primitives
underlie); [[test-unit-registry]] defines the unit. Distinct from
[[disorder-propensity-001-evidence]] — the pre-implementation evidence artifact sourced from
`docs/Evidence/` that records the source trace and recommended coverage — this page is the independent
two-stage re-validation verdict. Sibling reports [[disorder-pred-001-report]] (the windowed
predictor), [[disorder-lc-001-report]] (SEG low-complexity) and [[disorder-morf-001-report]] (MoRF)
cover different units of the same protein-disorder family.

## Verdict

**Stage A: PASS · Stage B: PASS · End state: CLEAN.** No algorithm defect, no code defect, **no
code/test/spec change this session**. Full unfiltered suite **6609 passed / 0 failed / 0 skipped**;
the 14-test canonical fixture is green against the actual code. No divergences at either stage.

## Scope note (as registered)

The session prompt framed this unit as a per-residue propensity with **sliding-window smoothing /
mean-disorder profile** — but that windowing logic is **not** part of this unit. It belongs to
`DisorderPredictor.PredictDisorder` under **DISORDER-PRED-001** ([[disorder-pred-001-report]]). Per
the Registry (`ALGORITHMS_CHECKLIST_V2.md` §DISORDER-PROPENSITY-001, lines 3658-3676) and the Method
Index (lines 5070-5073), DISORDER-PROPENSITY-001 is scoped to the **O(1) primitives** only: the raw
Table-2 lookup and the Dunker classification. The report validates the unit **as registered**.

## Canonical methods & source under test

Five public members on `DisorderPredictor.cs`
(`src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/`):

- `GetDisorderPropensity(char)` (`:680-683`) — `GetValueOrDefault(ToUpperInvariant(aa), 0)` over the
  `DisorderPropensity` dictionary (`:85-107`); realises lookup + case-folding + the **unknown→0.0**
  contract. Returns the **raw, un-normalized** Table-2 value, **not** the `[0,1]` normalized `p(c)`
  the windowed `Sᵢ` uses.
- `IsDisorderPromoting(char)` (`:691-694`) — `DisorderPromotingSet.Contains(ToUpperInvariant(aa))`.
- `DisorderPromotingAminoAcids` / `OrderPromotingAminoAcids` / `AmbiguousAminoAcids` (`:700-712`) —
  pre-sorted cached lists (`OrderBy(c => c)`), ascending and stable, over the three sets (`:111-121`).

## Stage A — description (algorithm faithfulness)

Sources re-opened live this session: **Campen et al. 2008 (TOP-IDP, PMC2676888)** via WebFetch of the
PMC full text — Table 2 extracted verbatim (all 20 values), the rank string, the prediction cut-off
and the extreme anchors; **Wikipedia "Intrinsically disordered proteins"** (citing Dunker 2001
primary) for the three classification sets; **localCIDER (Pappu lab)** `git clone` corroborating the
Campen citation (`sequenceParameters.py:221`) but not the raw numbers; and **PubMed 18991772 / USF
DigitalCommons / Bentham** independently confirming the rank string.

- **Scale check — all 20 values match exactly** vs the fetched Table 2 (W −0.884, F −0.697, Y −0.510,
  I −0.486, M −0.397, L −0.326, V −0.121, N 0.007, C 0.020, T 0.059, A 0.060, G 0.166, R 0.180,
  D 0.192, H 0.303, Q 0.318, S 0.341, K 0.586, E 0.736, P 0.987). Min = **W −0.884**, Max = **P
  0.987**, both source-confirmed.
- **Dunker classification check** — disorder-promoting {A,R,G,Q,S,P,E,K}, order-promoting
  {W,C,F,I,Y,V,L,N}, ambiguous {H,M,T,D}, all 8/8/4 matching the code; union = 20 distinct residues,
  three sets pairwise disjoint (verified by hand) — invariant INV-4 holds.
- **Documented pitfall confirmed & correctly handled** — the rendered rank string places
  "…Q, **K, S**, E, P", but the Table 2 **values** give S = 0.341 < K = 0.586 (so by value
  "…Q, **S, K**, E, P"). A presentation-order artifact; the per-residue values are authoritative and
  the implementation/tests lock both S=0.341 and K=0.586 — correct.
- **Edge-case semantics** (all defined) — scale defined only for the 20 standard residues (no value
  for B/J/O/U/X/Z or gaps); **unknown→0.0** is an implementation contract (`GetValueOrDefault(...,0)`),
  honestly registered as an assumption (TestSpec §6 #1), not a source value; case-insensitivity
  (upper-case before lookup) is an implementation contract.
- **Findings: None.** Every non-trivial expected value traces to a source fetched this session.

## Stage B — implementation

Code faithfully realises the validated description — the `DisorderPropensity` dictionary matches
Table 2 verbatim, the three sets match Dunker, and the getters return the pre-sorted cached lists.
Cross-verification was **run, not merely traced**: the full **unfiltered** suite (6609 passed) with
the 14-test canonical fixture recomputing the 20 Table-2 values, the W/P anchors, the three Dunker
sets (exact equivalence + counts), disjoint/cover-20, the membership⇔predicate invariant, unknown→0.0
and case-insensitivity — all green against the actual code.

- **Variant/delegate consistency:** `IsDisorderPromoting` ⇔ membership in `DisorderPromotingAminoAcids`
  (M10); the property getters expose the same underlying sets as the predicate/lookup.
- **Test-quality audit (HARD gate) — PASS.** The fixture hardcodes the expected scale (`TopIdp`) and
  the three sets as **independent literals copied from the published sources** (not code echoes); M1/M2
  assert exact values, M6/M7/M8 exact set equivalence + counts. M3 (disorder→true), M4 (order→false)
  and M5 (ambiguous→false) exercise **both** predicate branches, so neither an always-true nor an
  always-false impl survives. No green-washing (`Is.EqualTo(...).Within(1e-10)`, `Is.EquivalentTo` +
  `Count`; no `Greater/AtLeast/Contains` where an exact value is known; no widened tolerance; no
  skipped tests). Coverage: all 5 public members, both branches, the documented edge/contract cases
  (unknown X/Z/B/`*`→0.0; lowercase value & predicate; sorted property order); the S<K rank-vs-value
  pitfall is locked (M1 asserts S=0.341 and K=0.586). Honest green on the full unfiltered suite.
- **Build & suite:** `dotnet build …Tests.csproj -c Debug` → 0 errors, 4 pre-existing NUnit2007
  warnings (unrelated `ApproximateMatcher_EditDistance_Tests.cs`, no file changed); `dotnet test
  --no-build` → Failed 0, Passed 6609, Skipped 0.
- **Findings: None.**

## Findings & follow-ups

- **No algorithm defect and no code defect (State CLEAN).** No code/test/spec changed this session; the
  unit is fully functional.
- All 20 TOP-IDP Table-2 values, the anchors (W −0.884 / P 0.987), and the Dunker 8/8/4 classification
  were independently re-confirmed against the fetched Campen 2008 PMC text and Wikipedia/Dunker 2001,
  with localCIDER + PubMed/Bentham corroborating the citation and rank string. Test-quality gate PASS.
- Scope note recorded: the prompt's sliding-window / mean-disorder description applies to
  **DISORDER-PRED-001** (`PredictDisorder`, [[disorder-pred-001-report]]), not this O(1) unit. **No
  follow-ups.**
