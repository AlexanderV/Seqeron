---
type: source
title: "Validation report: COMPGEN-REVERSAL-001 (reversal / inversion distance — unsigned breakpoint lower bound ⌈b/2⌉, ComparativeGenomics.CalculateReversalDistance)"
tags: [validation, comparative-genomics, governance]
doc_path: docs/Validation/reports/COMPGEN-REVERSAL-001.md
sources:
  - docs/Validation/reports/COMPGEN-REVERSAL-001.md
source_commit: e4a1444b69f5b25d8a9f776d0c7f7c36746d8425
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Validation report: COMPGEN-REVERSAL-001

The two-stage **validation write-up** for test unit **COMPGEN-REVERSAL-001** — **reversal (inversion)
distance** returned as the **unsigned breakpoint lower bound** `⌈b/2⌉`
(`ComparativeGenomics.CalculateReversalDistance(permutation1, permutation2) → int`), validated
2026-06-16. This is the *report* artifact that feeds one row of the [[validation-ledger]]; it records
the validator's independent **verdict** on both the algorithm description (Stage A) and the shipped code
(Stage B), within the wider [[validation-and-testing]] campaign. The algorithm itself — the unsigned
breakpoint criterion `|Δ|≠1`, the `d ≥ b/2` bound, the extended-permutation sentinels, the oracle table
and the "lower bound, not exact" scoping — is synthesized in the shared concept
[[genome-rearrangement-breakpoint-distance]] (§ "Unsigned reversal distance `⌈b/2⌉`"), where this unit
is the **unsigned specialization** of the same breakpoint theory that backs COMPGEN-REARR-001.
[[test-unit-registry]] defines the unit. Distinct from [[compgen-reversal-001-evidence]] — the
pre-implementation Evidence artifact sourced from `docs/Evidence/` — this page is the independent
two-stage re-validation verdict. It is the sibling of [[compgen-rearr-001-report]] (the signed
`DetectRearrangements`/`ClassifyRearrangement` counterpart).

## Verdict

**Stage A: PASS · Stage B: PASS · End state: ✅ CLEAN.** No defect found; **no code and no test change
required**. The method is a correct, honestly-labelled *lower bound* — it returns the unsigned
breakpoint bound `⌈b/2⌉`, explicitly **not** the exact signed Hannenhalli–Pevzner reversal distance
(no cycle/hurdle refinement), and the docs say so. That limitation is by design, not a defect. Full
unfiltered suite **6605 passed, 0 failed**; `dotnet build` 0 errors; the validated file builds
warning-free (the 4 build warnings are pre-existing NUnit2007 warnings in the unrelated
`ApproximateMatcher_EditDistance_Tests.cs`).

## Canonical method & source under test

`ComparativeGenomics.CalculateReversalDistance(IReadOnlyList<int> permutation1, IReadOnlyList<int> permutation2)`
→ `int`, in `src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/ComparativeGenomics.cs:840–880`:

- Equal-length guard → `ArgumentException` (`:844`).
- `n ≤ 1` → `0` (`:849`, no internal adjacency).
- Relative permutation: target relabelled to identity `0..n−1` via a position map (`:853–857`).
- Left sentinel: `relative[0] != 0` ⇒ breakpoint (`:864`) — equivalent to sentinel `−1`
  (`|relative[0]−(−1)|=1` iff `relative[0]=0`).
- Interior: `|Δ| != 1` ⇒ breakpoint (`:870`) — **matches Hübotter Def 2.1 exactly.**
- Right sentinel: `relative[n−1] != n−1` ⇒ breakpoint (`:875`) — equivalent to sentinel `n`.
- Return `(b+1)/2` integer = `⌈b/2⌉` (`:879`); the constant `MaxBreakpointsRemovedPerReversal = 2` is
  the divisor.

Single public static method; no overloads / `*Fast` / instance variants (variant-consistency check
N/A).

## Stage A — description (algorithm faithfulness)

Confirmed against sources **independently retrieved this session** (see
[[compgen-reversal-001-evidence]] for the wider provenance): **Hübotter J (2020), "On Sorting by
Reversals"** (PDF fetched + `pdftotext`) and **Hunter College CompBio Lecture 16** (PDF fetched +
`pdftotext`), with WebSearch corroboration (Hübotter; Miranda et al. arXiv:1602.00778) attributing the
`b(π)/2 ≤ d_rev(π)` bound to Kececioglu–Sankoff / Bafna–Pevzner.

Verbatim confirmations from Hübotter:

- **Def 2.1** — adjacency vs breakpoint via `|i − j| = 1`; i.e. the *unsigned* test
  `|π_{i+1} − π_i| ≠ 1` the code implements.
- **Extended permutation** — expand with a leading `0` and trailing `n+1`, so `b(π)=0 iff π=id` (this
  disambiguates the reversed identity `(n … 2 1)`, otherwise breakpoint-free).
- **Corollary 2.1.1** — Kececioglu et al.'s first fundamental lower bound `d(π) ≥ b(π)/2`, because "a
  reversal can at most add or eliminate 2 breakpoints". Hunter restates `b(α) − b(αρ) ≤ 2` ⇒
  `d(α) ≥ b(α)/2`, and explicitly "this lower bound is not very tight" and "a reversal that decreases
  breakpoints by two is not necessarily one that makes progress" — confirming the value is a *lower
  bound*, not the exact distance.

**Formula check.** Unsigned breakpoint `|Δ|≠1` = Hübotter Def 2.1; sentinel extension makes identity the
unique 0-breakpoint permutation; `⌈b/2⌉ = (b+1)/2` is the tightest integer satisfying the real-valued
`d ≥ b/2` — not an invented value.

**Edge-case semantics** (all defined and sourced): identity ⇒ `b=0` ⇒ `0`; `n ≤ 1` ⇒ `0` (only the
sentinel adjacency); unequal lengths ⇒ `ArgumentException` (reversal distance is defined only within one
marker set — a documented assumption, the only well-defined behaviour, acceptable).

**Independent cross-check** (computed this session from the *definition*, not the C# code — a standalone
Python reimplementation of `|Δ|≠1` + sentinels + `⌈b/2⌉` reproduced every expected value):

| Case | Input → target | `b` | `⌈b/2⌉` |
|------|----------------|:---:|:---:|
| M1 identity | `[1,2,3,4,5]` → same | 0 | 0 |
| M2 Hunter unsigned | `[2,3,1,6,5,4]` → `[1..6]` | 4 | 2 |
| M3 fully reversed | `[4,3,2,1]` → `[1..4]` | 2 | 1 |
| M4 adjacent swap | `[1,2,4,3]` → `[1..4]` | 2 | 1 |
| C1 arbitrary labels | `[30,10,20]` → `[10,20,30]` | 3 | 2 |
| M5 (3 reversals on `[1..8]`) | `[4,5,1,3,8,7,6,2]` → `[1..8]` | 6 | 3 (= applied 3) |

Hübotter Figure 1 example `(0 7 5 6 3 2 4 1 8)` hand-checked: 2 adjacencies (5 6, 3 2) of 8 pairs ⇒
`b=6`, consistent. **No divergences.** The "Simplified" status is honest. → **Stage A: PASS.**

## Stage B — implementation

Formula realised correctly. The sentinel encoding (`relative[0] != 0` left, `|Δ| != 1` interior,
`relative[n−1] != n−1` right) is exactly equivalent to the extended-permutation `(−1, rel…, n)` form,
verified by the independent reimplementation reproducing all six values. The divisor
`MaxBreakpointsRemovedPerReversal = 2` is justified by "a reversal removes ≤ 2 breakpoints".

**Numerical robustness.** O(n) single pass, one `Dictionary<int,int>`; no overflow (counts ≤ n+1), no
div-by-zero (constant divisor 2). Distinct-label assumption: a marker in `permutation1` absent from
`permutation2` makes `positionMap[x]` throw `KeyNotFoundException` — acceptable, the input contract
requires the same marker set (documented in §3.3 of the algorithm doc).

**Test-quality audit (HARD gate) — PASS, no change.**
`tests/Seqeron/Seqeron.Genomics.Tests/ComparativeGenomics_CalculateReversalDistance_Tests.cs`, 10 tests:

- **Sourced expectations, not code echoes:** M1–M4, C1 assert *exact* integers (0/2/1/1/2)
  independently recomputed from Hübotter Def 2.1 + the `d ≥ b/2` bound, not from the implementation; a
  deliberately-wrong implementation would fail them.
- **No green-washing:** exact `Is.EqualTo` wherever an exact value is known. The only range assertion
  (M5, `0 ≤ d ≤ applied_reversals`, INV-05) is the genuine lower-bound *invariant* — a range is correct
  there because the invariant itself is an inequality, not a weakened known value. No skips / ignores /
  widened tolerances.
- **Branch coverage:** zero-breakpoint (M1); interior + both sentinels mixed (M2, M4, C1);
  end-boundary-only breakpoints (M3); lower-bound property (M5); empty (S1), single (S2),
  unequal-length throw (S3), symmetry INV-03 (S4), arbitrary non-`1..n` labels / relabeling (C1). Every
  branch (left sentinel, interior, right sentinel, `⌈b/2⌉` rounding) is exercised.
- **Honest green:** full unfiltered suite Failed 0 / Passed 6605; `dotnet build` 0 errors.

→ **Stage B: PASS.** No test was weakened or added — the existing tests already lock the
externally-sourced values and cover all branches / edge cases.

## Findings

- **No defect.** Description and code agree with Hübotter Def 2.1 / Corollary 2.1.1 and Hunter Lecture
  16; the sentinel encoding is provably equivalent to the extended permutation; the full suite is green.
- **End state: ✅ CLEAN — no code or test change.** The method is a correct, honestly-labelled unsigned
  reversal-distance **lower bound** `⌈b/2⌉`, not the exact signed Hannenhalli–Pevzner distance
  (cycle/hurdle refinement **not implemented**) — a by-design, documented limitation, not a defect.
- **No open follow-ups.**
