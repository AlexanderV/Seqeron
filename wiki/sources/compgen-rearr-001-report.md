---
type: source
title: "Validation report: COMPGEN-REARR-001 (genome rearrangement detection by breakpoints — signed gene-order permutation, ComparativeGenomics.DetectRearrangements / ClassifyRearrangement)"
tags: [validation, comparative-genomics, governance]
doc_path: docs/Validation/reports/COMPGEN-REARR-001.md
sources:
  - docs/Validation/reports/COMPGEN-REARR-001.md
source_commit: 4c3caf900067a440f88ab2a5d4addc3dac8cb20f
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Validation report: COMPGEN-REARR-001

The two-stage **validation write-up** for test unit **COMPGEN-REARR-001** — genome rearrangement
detection by **breakpoints** on a *signed gene-order permutation* (`ComparativeGenomics.DetectRearrangements`
+ `ComparativeGenomics.ClassifyRearrangement`), validated 2026-06-15. This is the *report* artifact
that feeds one row of the [[validation-ledger]]; it records the validator's independent **verdict** on
both the algorithm description (Stage A) and the shipped code (Stage B), and the wider campaign is
[[validation-and-testing]]. The algorithm, its signed-permutation model, breakpoint criterion,
invariants, documented oracles and edge cases are synthesized in the concept
[[genome-rearrangement-breakpoint-distance]]; [[test-unit-registry]] defines the unit. Distinct from
[[compgen-rearr-001-evidence]] — the pre-implementation evidence artifact sourced from
`docs/Evidence/` — this page is the independent two-stage re-validation verdict. It is the
signed-permutation counterpart to the block-signal [[synteny-and-rearrangement-detection]] (the two are
complementary `alternative_to` formulations of "what rearrangements separate two genomes?").

## Verdict

**Stage A: PASS-WITH-NOTES · Stage B: PASS-WITH-NOTES · End state: CLEAN.** No correctness defect in
`DetectRearrangements`; the two notes are the honestly-documented **simplified per-boundary classifier**
and a telomere-convention nuance in the `d_BP` equality — both scoped and labelled, not code defects.
The only defects found were **two test-coverage gaps**, both **fixed in-session**. Full unfiltered
suite **6508 passed, 0 failed, 0 skipped**(*); `dotnet build` 0 errors; the REARR test file + implementation
build warning-free. (*) the reported `Skipped` is an unrelated `[Explicit]` MFE benchmark.

## Canonical methods & source under test

In `src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/ComparativeGenomics.cs`:

- `DetectRearrangements(genome1Genes, genome2Genes, orthologMap)` — eager null checks `:581–591`
  (fire before enumeration, separated from the `yield`-iterator) → `DetectRearrangementsIterator`
  `:593–669`. Builds `signedRank = sign × (1-based genome-2 rank)`, relabels markers to genome-2 order
  positions `1..n` (sort by `|rank|`, preserve sign), then walks `[0, relabelled…, n+1]` emitting a
  breakpoint whenever `curr != prev + 1` (`:651`) — the verified exact criterion. `<2` markers →
  `yield break` (`:619–620`, the sourced "no internal adjacency" case); dangling ortholog guarded by
  `TryGetValue` (`:610–611`, skip unmapped/absent anchors).
- `ClassifyBoundary` `:680–694`; `ClassifyRearrangement(RearrangementEvent)` `:708–723` — re-parses
  `TargetPosition "x->y"` and delegates to the same `ClassifyBoundary` used at emission time, so the
  stored `Type` and a re-classification agree. `CompareGenomes` (COMPGEN-COMPARE-001) consumes
  `DetectRearrangements` (`:785`) unchanged; no other overloads.
- Tests: 16 REARR tests after the in-session additions.

## Stage A — description (algorithm faithfulness)

Confirmed against three authority-rank-1 sources (see [[compgen-rearr-001-evidence]] for the full
provenance): **Hunter College CompBio Lecture 16** (Hannenhalli–Pevzner / Bafna–Pevzner; PDF retrieved
and decoded with `pdftotext`), **Tannier/Sankoff breakpoint distance** (PMC3887456), and **Bafna &
Pevzner 1998**. Verbatim confirmations from Hunter: signed-permutation definition `α(i)=±a`; reversal
`ρ=[i,j]` reverses a contiguous block; extended permutation `(0, …, n+1)`; the breakpoint criterion
"if `(x,y)` appears in extended α but neither `(x,y)` nor `(−y,−x)` appears in extended β"; the worked
example `α=(0,−2,−3,+1,+6,−5,−4,+7)` → **b(α)=6** (with `(−5,−4)` *not* a breakpoint since `(4,5)∈β`);
`b_β(β)=0`; lower bound `d(α) ≥ b(α)/2`.

**Formula reduction proven exact (not an approximation).** The doc reduces the breakpoint test to
"breakpoint iff `y ≠ x+1`". Against β=identity (adjacencies are exactly `(i,i+1)`), both clauses of the
full Hunter criterion collapse to the same condition: clause A `(x,y)=(i,i+1) ⟺ y=x+1`; clause B
`(−y,−x)=(i,i+1) ⟺ −x=−y+1 ⟺ y=x+1`. So "non-breakpoint ⇔ `y=x+1`" is exact — the documented
non-breakpoint `(−5,−4)` confirms it (`−4=−5+1`), validating INV-02's "`y=x+1` subsumes both clauses".

**Edge-case semantics** (all defined and sourced): identity → 0; `<2` markers → 0 (no internal
adjacency); extended endpoints `(0,π₁)`/`(πₙ,n+1)` may themselves be breakpoints (sourced from the
worked example `(0,−2)`,`(−4,+7)`); null args → `ArgumentNullException`.

**Independent hand cross-check** (this session, computed *before* running the code) over the extended
permutations using `y ≠ x+1`:

| Case | Relative perm | b (independent) | Match |
|------|---------------|-----------------|-------|
| M1 Hunter | (−2,−3,+1,+6,−5,−4) | 6 — `(0,−2)(−2,−3)(−3,1)(1,6)(6,−5)(−4,7)` | 6 ✓ |
| M2 identity | (1,2,3,4,5) | 0 | 0 ✓ |
| M3/S2 reversal | (+1,−4,−3,−2,+5) | 2 — `(1,−4)(−2,5)` | 2 ✓ |
| M4 descending pair | (−2,−1,+3,+4,+5) | 2 — `(0,−2)(−1,3)`; `(−2,−1)` excluded | 2 ✓ |
| M6 transposition | (+1,+2,+4,+5,+3) | 3 (all positive) | all Transposition ✓ |
| C1 full-reverse | (4,3,2,1) | 5 | ∈[0,n+1] ✓ |

M1 reproduces the Hunter published count of 6 exactly, including the `(−5,−4)`-style exclusion.

**Stage A notes (documented divergences, not errors):**

1. **Classification is a documented simplification, not a formal result.** Assigning Inversion vs
   Transposition from a *single boundary's local sign signature* (INV-05) has no formal theorem behind
   it (a minimal scenario needs global cycle/graph analysis). Evidence Assumption 3 and doc §5.3 mark
   it "intentionally simplified" and scope it to Inversion/Transposition. It is a reasonable heuristic
   consistent with the sourced operation definitions (reversal negates signs; transposition preserves
   orientation), not a proven invariant.
2. **`d_BP` bookkeeping is convention-dependent.** The doc equates the extended-permutation breakpoint
   *count* `b(α)` with Tannier's `d = n − sim`; these differ by telomere-accounting conventions. The
   unit's reported quantity is `b(α)` (the Hunter quantity, which every test checks) — not a
   correctness defect, but the equality as written depends on the telomere convention.

→ **Stage A: PASS-WITH-NOTES.**

## Stage B — implementation

Formula realised correctly. The iterator's `curr != prev + 1` walk over the relabelled extended
permutation is the verified exact criterion; the relabelling is a no-op when ranks are contiguous (all
tests), so the relative permutation equals the intended signed perm — confirmed by tracing M1/M3/M4.
Null checks eager and separated from the `yield`-iterator; `<2` markers `yield break`; dangling
orthologs `TryGetValue`-guarded.

**Cross-verification recomputed vs code** — the suite was run and all counts match the independent hand
computation above (M1=6, M2=0, M3=2, M4=2, S2=2, C1∈[0,n+1], C2=0). Classification: M5 boundaries
`(1,−4)`,`(−2,5)` → **Inversion** (sign flip); M6 all-positive boundaries → **Transposition**. Matches
code. `ClassifyRearrangement` re-parses `TargetPosition` and delegates to the same `ClassifyBoundary`
used at emission → stored `Type` and re-classification agree.

**Test-quality audit (HARD gate) — PASS after fix.** M1=6 is the Hunter published count; M3/M4/S2=2 and
the C1 bounds were independently hand-computed *before* running the code; assertions are exact
`Has.Count.EqualTo(...)` (identity `Is.Empty`/`Is.Zero`, nulls `Assert.Throws<ArgumentNullException>`) —
none would survive a deliberately-wrong count. No green-washing (the legacy permissive
`Any(...)`/`GreaterThanOrEqualTo(0)` tests were already removed per the spec).

**Two coverage gaps found & fixed this session** (the only Stage-B defects):

- **M9b** — **null genome2** was unexercised even though the code validates it (added: throws
  `ArgumentNullException`).
- **M10** — the **`ClassifyRearrangement` fallback branch** (null / unparsable `TargetPosition` →
  return the stored `Type`) was unexercised (added: fallback returns stored `Type` for null and
  malformed `TargetPosition`).

Both are sourced to the contract §3.3 and the documented fallback; the suite is now 16 REARR tests.
→ **Stage B: PASS-WITH-NOTES** (the only notes are the documented simplified classifier).

## Findings

- **No correctness defect in `DetectRearrangements`** — the `y≠x+1` reduction is provably exact and
  reproduces the Hunter published example (b=6) plus every independent hand computation.
- **Two test-coverage gaps fixed in-session** (M9b null genome2; M10 `ClassifyRearrangement` fallback) —
  the only defects, both closed. **End state: CLEAN.**
- **Documented (not defects):** the per-boundary Inversion/Transposition classifier is an intentionally
  simplified heuristic with no formal single-permutation basis (scoped by doc §5.3 / Evidence
  Assumption 3); the `d_BP = n − sim` equality is telomere-convention-dependent while the unit reports
  and tests the extended-permutation count `b(α)`.
</content>
</invoke>
