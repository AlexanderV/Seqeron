# Validation Report: COMPGEN-REARR-001 — Genome Rearrangement Detection by Breakpoints

- **Validated:** 2026-06-15   **Area:** Comparative
- **Canonical method(s):** `ComparativeGenomics.DetectRearrangements(genome1Genes, genome2Genes, orthologMap)`, `ComparativeGenomics.ClassifyRearrangement(RearrangementEvent)`
- **Stage A verdict:** PASS-WITH-NOTES
- **Stage B verdict:** PASS-WITH-NOTES

## Stage A — Description

### Sources opened & what they confirm

1. **Hunter College CompBio Lecture 16** (PDF retrieved this session and decoded with `pdftotext` —
   the binary saved by WebFetch, lines 252–312 of the extracted text). Confirmed **verbatim**:
   - Signed permutation definition (line 253): `α(i) = +a or −a, a ∈ L`.
   - Reversal `ρ=[i,j]` reverses a contiguous block (line 255–257) — captures reversal of contiguous blocks.
   - Extended permutation with `0` prepended and `n+1` appended (lines 282, 286–287).
   - **Breakpoint definition (line 284–285):** "If `(x, y)` appear in (extended) α but neither `(x, y)`
     nor `(−y, −x)` appear in (extended) β, then `(x, y)` is a breakpoint of α with respect to β."
   - **Worked example (lines 286–289):** `α=(0,−2,−3,+1,+6,−5,−4,+7)`, breakpoints listed as
     `(0,−2),(−2,−3),(−3,+1),(+1,+6),(+6,−5),(−4,+7)` → **b(α)=6**; "`(−5,−4)` is not a breakpoint
     since `(4,5)` appears in β."
   - `b_β(β)=0` (line 293); lower bound `d(α) ≥ b(α)/2` (lines 311–312).
2. **Tannier/Sankoff breakpoint distance** (PMC3887456, per evidence doc): `d(π₁,π₂) = n − sim`,
   sim = common adjacencies; telomeric adjacencies = the extended endpoints.
3. **Bafna & Pevzner (1998)** (evidence doc): operation classes (inversion, transposition, deletion,
   insertion, duplication); transposition moves a block **preserving orientation** (no sign change);
   extension with `π₀=0`, `π_{n+1}=n+1`.

### Formula check

The doc's reduction is "breakpoint iff `y ≠ x+1`". I verified this is **mathematically equivalent**
to the full Hunter criterion, not an approximation:

- β = identity ⇒ its adjacencies are exactly the pairs `(i, i+1)`, i = 0..n.
- Clause A `(x,y) = (i,i+1)` ⟺ `y = x+1`.
- Clause B `(−y,−x) = (i,i+1)` ⟺ `−y = i ∧ −x = i+1` ⟺ `−x = −y+1` ⟺ `x = y−1` ⟺ **`y = x+1`**.

Both clauses collapse to the same condition `y = x+1`, so "non-breakpoint ⇔ `y = x+1`" is exact.
The documented non-breakpoint `(−5,−4)` confirms it: `−4 = −5+1`. INV-02's claim that `y=x+1`
"subsumes the `(x,y)` and `(−y,−x)` clauses" is therefore correct.

### Edge-case semantics check

Identity → 0 (sourced `b(β)=0`); <2 markers → no internal adjacency → 0 (sourced); extended
endpoints `(0,π₁)`/`(πₙ,n+1)` may themselves be breakpoints (sourced from worked example
`(0,−2)`, `(−4,+7)`); null args → `ArgumentNullException` (contract). All defined and sourced.

### Independent cross-check (hand computation, this session)

Ran the verified criterion `y ≠ x+1` over the extended permutations:

| Case | Relative perm | b (independent) | Spec/test expects |
|------|---------------|-----------------|-------------------|
| M1 Hunter | (−2,−3,+1,+6,−5,−4) | **6** — `(0,−2),(−2,−3),(−3,1),(1,6),(6,−5),(−4,7)` | 6 ✓ |
| M2 identity | (1,2,3,4,5) | 0 | 0 ✓ |
| M3/S2 reversal | (+1,−4,−3,−2,+5) | **2** — `(1,−4),(−2,5)` | 2 ✓ |
| M4 descending pair | (−2,−1,+3,+4,+5) | **2** — `(0,−2),(−1,3)`; `(−2,−1)` excluded | 2 ✓ |
| M6 transposition | (+1,+2,+4,+5,+3) | 3 — all positive | non-empty, all Transposition ✓ |
| C1 full-reverse | (4,3,2,1) | 5 | ∈[0,n+1] ✓ |

M1 reproduces the Hunter published count of 6 exactly, including the `(−5,−4)`-style exclusion.

### Findings / divergences (Stage A)

- **NOTE 1 (classification is a documented simplification, not a formal result).** Assigning
  Inversion vs Transposition from a *single boundary's local sign signature* has no formal
  theorem behind it (a minimal scenario requires global cycle/graph analysis). The spec is honest:
  Evidence Assumption 3 and doc §5.3 mark it "intentionally simplified" and limit scope to
  Inversion/Transposition. INV-05 is a reasonable heuristic consistent with the sourced operation
  definitions (reversal negates signs; transposition preserves orientation), not a proven invariant.
- **NOTE 2 (`d_BP` formula bookkeeping).** The doc equates the breakpoint *count* (over n+1 internal
  pairs of the extended permutation) with Tannier's `d = n − sim`. These differ by telomere
  accounting conventions; the unit's reported quantity is the extended-permutation breakpoint count
  `b(α)`, which is the Hunter quantity and is what every test checks. Not a correctness defect, but
  the equality as written is convention-dependent.

These are documented divergences, not errors → **Stage A: PASS-WITH-NOTES.**

## Stage B — Implementation

### Code path reviewed

`src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/ComparativeGenomics.cs`
- `DetectRearrangements` L581–591 (eager null checks) → `DetectRearrangementsIterator` L593–669.
- `ClassifyBoundary` L680–694; `ClassifyRearrangement` L708–723.

### Formula realised correctly?

Yes. The iterator builds `signedRank = sign × (1-based genome-2 rank)`, relabels markers to their
genome-2 order positions `1..n` (sort by |rank|, preserve sign), then walks `[0, relabelled…, n+1]`
emitting a breakpoint whenever `curr != prev + 1` (L651). This is the verified exact criterion.
- Null checks eager (L586–588), separated from the `yield`-iterator so they fire before enumeration.
- `<2` markers → `yield break` (L619–620): the sourced "no internal adjacency" case.
- Dangling ortholog: `TryGetValue` guards (L610–611) skip unmapped/absent anchors.
- The relabelling is a no-op when ranks are contiguous (all tests), so the relative permutation
  equals the intended signed perm — confirmed by tracing M1/M3/M4.

### Cross-verification table recomputed vs code

Ran the actual code via the test suite; all counts match the independent hand computation above
(M1=6, M2=0, M3=2, M4=2, S2=2, C1∈[0,n+1], C2=0). Classification: M5 boundaries `(1,−4)`,`(−2,5)`
→ Inversion (sign flip / negative involved); M6 boundaries all-positive → Transposition. Matches code.

### Variant/delegate consistency

`ClassifyRearrangement` re-parses `TargetPosition "x->y"` and delegates to the same `ClassifyBoundary`
used at emission time → the event's stored `Type` and a re-classification agree. `CompareGenomes`
consumes `DetectRearrangements` (L785) unchanged. No other overloads.

### Test quality audit (TEST-QUALITY GATE)

- **Sourced, not code-echoes:** M1=6 is the Hunter published count; M3/M4/S2=2 and the C1 bounds
  were independently hand-computed this session from the verified criterion, *before* running the
  code. Exact `Has.Count.EqualTo(...)` assertions — none would survive a deliberately-wrong count.
- **No green-washing:** counts asserted with exact equality; identity asserted `Is.Empty`/`Is.Zero`;
  nulls with `Assert.Throws<ArgumentNullException>`. No weakened/ranged/skipped assertions found.
  (Legacy permissive `Any(...)`/`GreaterThanOrEqualTo(0)` tests were already removed per the spec.)
- **Coverage gaps found & fixed this session:** the original 14 tests omitted (a) **null genome2**
  even though the code validates it, and (b) the **`ClassifyRearrangement` fallback branch** (null /
  unparsable `TargetPosition` → return stored `Type`). Added **M9b** (null genome2 throws) and
  **M10** (fallback returns stored Type for null and malformed `TargetPosition`). Both sourced to the
  contract §3.3 and the documented fallback. Suite now 16 REARR tests.
- **Honest green:** full unfiltered suite **6508 passed, 0 failed, 0 skipped(*)**; `dotnet build`
  0 errors; the REARR test file and implementation build warning-free (the 4 NUnit2007 warnings are
  pre-existing in `ApproximateMatcher_EditDistance_Tests.cs`, unrelated to this unit).
  (*) `MFE_Benchmark_AllScenarios` is an `[Explicit]` benchmark, reported Skipped — unrelated.

### Findings / defects (Stage B)

- No correctness defect in `DetectRearrangements`. The `y≠x+1` reduction is provably exact.
- The two missing-branch test gaps were the only defects; both **fixed** this session.
- Classification limitation is inherent to the documented simplified model (Stage A NOTE 1); it is
  scoped and labelled, not a code defect.

→ **Stage B: PASS-WITH-NOTES** (passes; the only notes are the documented simplified classifier).

## Verdict & follow-ups

- **Stage A: PASS-WITH-NOTES** — biology/maths correct and verified against Hunter/Tannier/Bafna–Pevzner;
  `y≠x+1` reduction proven equivalent to the full breakpoint criterion. Notes: simplified
  per-boundary classification (no formal basis, scoped) and telomere-convention nuance in the
  `d_BP` equality.
- **Stage B: PASS-WITH-NOTES** — code realises the verified criterion exactly; all counts reproduce
  the Hunter published example and independent hand computations. Two test-coverage gaps (null
  genome2; `ClassifyRearrangement` fallback) found and fixed (M9b, M10).
- **Test-quality gate: PASS** (after fix) — sourced exact assertions, all public methods and Stage-A
  branches/edge cases now exercised, full unfiltered suite green (6508/0).
- **End-state: CLEAN** — no correctness defect; the test gaps found were completely fixed in-session;
  build + full suite green.
