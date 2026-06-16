# Validation Report: COMPGEN-REVERSAL-001 — Reversal Distance (unsigned breakpoint lower bound)

- **Validated:** 2026-06-16   **Area:** Comparative
- **Canonical method(s):** `ComparativeGenomics.CalculateReversalDistance(IReadOnlyList<int> permutation1, IReadOnlyList<int> permutation2)` → `int` (unsigned breakpoint lower bound ⌈b/2⌉)
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS

## Stage A — Description

### Sources opened this session (independently retrieved)

1. **Hübotter J (2020). On Sorting by Reversals.** PDF fetched and text-extracted via `pdftotext` this session.
   - **Definition 2.1 (verbatim):** "Let i ∼ j if |i − j| = 1. A pair of consecutive elements π_i and π_j forms an **adjacency** if π_i ∼ π_j and a **breakpoint** if π_i ̸∼ π_j." → This is exactly the *unsigned* breakpoint test `|π_{i+1} − π_i| ≠ 1` that the code implements.
   - **Extended permutation (verbatim):** "we extend our definition … by expanding permutations with one initial 0 and one trailing n + 1 … π = (0 π₁ … πₙ n + 1) … b(π) = 0 if and only if π = id." (This expansion is what disambiguates the reversed identity (n …2 1), which is otherwise breakpoint-free.)
   - **Corollary 2.1.1 (verbatim):** "Kececioglu et al. found a first fundamental lower bound … d(π) ≥ b(π)/2." Plus: "a reversal … k ∈ [−2, 2] as a reversal can at most add or eliminate 2 breakpoints."
2. **Hunter College CompBio Lecture 16.** PDF fetched and text-extracted via `pdftotext` this session.
   - Extended version definition; "a reversal ρ = [i, j] can reduce the number of breakpoints by at most two … b(α) − b(αρ) ≤ 2"; "Therefore b(α) ≤ 2t. But d(α) ≥ t; therefore, d(α) ≥ b(α)/2. This lower bound is not very tight."
   - "a reversal ρ that decreases the number of breakpoints by two is not necessarily one that makes progress" (confirms the value is a *lower bound*, not exact).
3. **WebSearch corroboration** (Hübotter, Miranda et al., arXiv:1602.00778): breakpoint lower bound b(π)/2 ≤ d_rev(π) attributed to Kececioglu–Sankoff / Bafna–Pevzner; "a reversal can at most add or eliminate 2 breakpoints"; sorting unsigned permutations by reversals is NP-hard, so an O(n) breakpoint bound is necessarily only a *lower bound*, never the exact distance.

### Formula check
- Unsigned breakpoint = consecutive pair with `|Δ| ≠ 1` — **matches Hübotter Def 2.1 exactly.**
- Extended permutation prepends/appends sentinels so identity is the unique 0-breakpoint permutation — **matches Hübotter/Hunter.**
- Lower bound `d ≥ b/2`; smallest satisfying integer `⌈b/2⌉ = (b+1)/2` — **matches Corollary 2.1.1 / Hunter.** The integer ceiling is the tightest integer guaranteed by the real-valued theorem; it is not an invented value.

### Edge-case semantics check
- Identity ⇒ b=0 ⇒ 0 (Hübotter: "b(π)=0 iff π=id"). Defined and sourced.
- n ≤ 1 ⇒ 0 (no internal adjacency). Consistent with the extended definition (only sentinel pair, which is an adjacency for the identity-mapped single element).
- Unequal lengths ⇒ `ArgumentException`. Not separately specified by sources (reversal distance is defined only between permutations of the same marker set); the throw is the only well-defined behaviour — documented assumption, acceptable.

### Independent cross-check (numbers I computed this session, from the definition — not from the C# code)
A standalone Python reimplementation of the sourced definition (`|Δ|≠1`, extended with sentinels, `⌈b/2⌉`) reproduced every expected test value:

| Case | Input → target | b (recomputed) | ⌈b/2⌉ |
|------|----------------|:---:|:---:|
| M1 identity | [1,2,3,4,5]→same | 0 | 0 |
| M2 Hunter unsigned | [2,3,1,6,5,4]→[1..6] | 4 | 2 |
| M3 fully reversed | [4,3,2,1]→[1..4] | 2 | 1 |
| M4 adjacent swap | [1,2,4,3]→[1..4] | 2 | 1 |
| C1 arbitrary labels | [30,10,20]→[10,20,30] | 3 | 2 |
| M5 (3 reversals on [1..8]) | [4,5,1,3,8,7,6,2]→[1..8] | 6 | 3 (= applied 3) |

Hübotter Figure 1 example `(0 7 5 6 3 2 4 1 8)` independently hand-checked: 2 adjacencies (5 6, 3 2) of 8 pairs ⇒ b=6, consistent with the definition used.

### Findings / divergences
None. The description is mathematically correct. The "Simplified" status is honest: it returns the breakpoint *lower bound* (unsigned), explicitly not the exact signed Hannenhalli–Pevzner distance, and the docs say so. Stage A **PASS**.

## Stage B — Implementation

### Code path reviewed
`src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/ComparativeGenomics.cs:840–880`.

- Equal-length guard → `ArgumentException` (line 844).
- `n ≤ 1` → 0 (line 849).
- Relative permutation: target → identity 0..n−1 via position map (lines 853–857).
- Left sentinel: `relative[0] != 0` ⇒ BP (line 864) — equivalent to sentinel −1 since `|relative[0]−(−1)|=1` iff `relative[0]=0`.
- Interior: `|Δ| != 1` ⇒ BP (line 870) — **matches Hübotter Def 2.1.**
- Right sentinel: `relative[n−1] != n−1` ⇒ BP (line 875) — equivalent to sentinel n.
- Return `(b+1)/2` integer = `⌈b/2⌉` (line 879).

### Formula realised correctly?
Yes. The sentinel encoding is exactly equivalent to the extended-permutation `(−1, rel…, n)` form (verified by independent reimplementation reproducing all values). `MaxBreakpointsRemovedPerReversal = 2` is the divisor justified by "a reversal removes ≤ 2 breakpoints" (Hübotter/Hunter).

### Cross-verification table recomputed vs code
All six cases above were produced by an independent (Python) implementation derived from the sourced definitions, then confirmed against the C# via the test suite. Every value matches.

### Variant/delegate consistency
Single public static method; no overloads/`*Fast`/instance variants. N/A.

### Numerical robustness
O(n) single pass; one `Dictionary<int,int>`. No overflow (counts ≤ n+1), no div-by-zero (constant divisor 2). Distinct-label assumption: if `permutation1` contains a marker absent from `permutation2`, `positionMap[x]` throws `KeyNotFoundException` — acceptable (input contract requires the same marker set; documented in §3.3 of the algorithm doc).

### Test quality audit (HARD gate)
File: `tests/Seqeron/Seqeron.Genomics.Tests/ComparativeGenomics_CalculateReversalDistance_Tests.cs` (10 tests).

- **Sourced expectations, not code echoes:** M1–M4, C1 assert *exact* integers (0/2/1/1/2) that I independently recomputed from Hübotter Def 2.1 + the d≥b/2 bound, **not** from the implementation. A deliberately-wrong implementation would fail them.
- **No green-washing:** exact `Is.EqualTo` used wherever an exact value is known. The only range assertion (M5) is the genuine *lower-bound invariant* test (`0 ≤ d ≤ applied_reversals`, INV-05) — a range is correct there because the invariant itself is an inequality, not because a known value was weakened. No skips/ignores/widened tolerances.
- **Coverage of all branches:** zero-breakpoint (M1), interior + both sentinels mixed (M2, M4, C1), end-boundary-only breakpoints (M3), lower-bound property (M5), empty (S1), single (S2), unequal-length throw (S3), symmetry INV-03 (S4), arbitrary non-1..n labels / relabeling (C1). Every code branch (left sentinel, interior, right sentinel, ⌈b/2⌉ rounding) is exercised.
- **Honest green:** FULL unfiltered suite = **Failed: 0, Passed: 6605**. `dotnet build` = 0 errors. The 4 build warnings are pre-existing NUnit2007 warnings in the unrelated `ApproximateMatcher_EditDistance_Tests.cs`; the validated file builds warning-free.

**Gate result: PASS.** No test was weakened or added — the existing tests already lock the externally-sourced values and cover all branches/edge cases.

### Findings / defects
None.

## Verdict & follow-ups
- **Stage A: PASS.** Description (unsigned breakpoint model, extended permutation, d ≥ b/2 ⇒ ⌈b/2⌉) is mathematically correct and matches Hübotter Def 2.1 / Corollary 2.1.1 and Hunter Lecture 16, independently retrieved and hand-verified this session.
- **Stage B: PASS.** Code faithfully realises the validated formula; tests assert exact sourced values, cover all branches and documented edge cases, and the full suite is green.
- **End-state: ✅ CLEAN.** No defect found; no code or test change required. The method is a correct, honestly-labelled *lower bound* (not exact reversal distance) — this limitation is documented in the algorithm doc/Evidence and is by design, not a defect.
