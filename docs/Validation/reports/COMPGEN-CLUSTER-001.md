# Validation Report: COMPGEN-CLUSTER-001 — Conserved Gene Clusters (common intervals of permutations)

- **Validated:** 2026-06-16   **Area:** Comparative
- **Canonical method(s):** `ComparativeGenomics.FindConservedClusters(genomes, orthologGroups, minClusterSize=3, maxGap=2)`
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS-WITH-NOTES (tests strengthened: 3 weak `Contains`/`Does.Not.Contain` assertions rewritten to exact sourced sets)
- **End-state:** ✅ CLEAN

## Scope clarification

The orchestrator prompt described this unit generically as "ortholog/gene clustering … COG/OrthoMCL …
connected-components / MCL". The unit's Registry row, TestSpec, and Evidence all unambiguously define
COMPGEN-CLUSTER-001 as the **common-interval gene-cluster model** (a set of ortholog-group labels that
forms a contiguous *interval* in every genome — Uno & Yagiura 2000; Heber & Stoye 2001; Bui-Xuan, Habib
& Paul 2013; sequence generalisation Didier et al. 2013). The single public method under test is
`FindConservedClusters`. COG/OrthoMCL ortholog-grouping is a *different* concept and is covered by other
units (COMPGEN-ORTHO-001, COMPGEN-RBH-001). This validation therefore validates the common-interval
model, which is what the code actually computes.

## Stage A — Description

### Sources opened & what they confirm

All sources fetched fresh this session as PDFs and extracted with `pdftotext`:

1. **Bui-Xuan, Habib, Paul (2013), arXiv:1304.5140** — fetched `https://arxiv.org/pdf/1304.5140`.
   - "The interval [i, j] of Pk, **defined only for 1 ≤ i < j ≤ n**, is the set of elements located
     between position i (included) and position j (included) in Pk." → intervals have size ≥ 2;
     singletons excluded. (Generalities §2.)
   - **Definition 1 [23]:** "A common interval of P is a set of integers that is an interval of each Pk,
     k ∈ [K]." → contiguous in *every* permutation.
   - **Example 1 (verbatim):** "Let P1 = Id7 and P2 = (7 2 1 3 6 4 5). Then the common intervals of
     P = {P1, P2} are (1..2), (1..3), (1..6), (1..7), (3..6), (4..5) and (4..6)." → sets
     {1,2},{1,2,3},{1..6},{1..7},{3,4,5,6},{4,5},{4,5,6}.

2. **Didier, Schmidt, Stoye, Tsur (2013), arXiv:1310.4290** — fetched `https://arxiv.org/pdf/1310.4290`.
   - Interval of a sequence: "any set I of integers from Σ such that there exist i, j with 1 ≤ i ≤ j ≤ n
     and **I = Set(T[i..j])**." (§ Preliminaries, line 132.)
   - **Definition 1 [8,17]:** a common interval of two sequences T and S is a set I that is an interval of
     both. (line 143.)
   - **Example 1 (verbatim):** T = 1 2 5 2 1 4 3 1 2 6 5, S = 5 6 4 2 3 4 1 5. "{1, 2} is an interval of
     T … but is **not** a common interval of T and S. An example of common interval is **{1, 2, 3, 4}**."

3. **Uno & Yagiura (2000)** and **Heber & Stoye (2001)** — used for provenance of the K=2 originating
   model and the k-permutation generalisation; cited as Definition 1 source [23] and the k-permutation
   algorithm respectively. (Full text behind Springer auth; bibliographic + definitional facts confirmed
   via the two arXiv papers that cite them.)

### Formula / definition check

- A conserved cluster = a **common interval** = a set of ortholog-group labels occupying a contiguous
  window (interval) in *every* genome, where an interval of a (possibly repeated-label) sequence is the
  *set of all elements* of some window `Set(T[i..j])`. Matches Source 1 Def. 1 and Source 2 §Prelim/Def.1
  exactly.
- Size ≥ 2 (interval defined only for `i < j`); whole set always trivially common. Matches Source 1 §2.

### Edge-case semantics

- **< 2 genomes:** common interval is a family notion; vacuous → empty. (Source 1 §2 family definition.)
- **Repeated labels (paralogs):** any matching window in each genome counts (Source 2 Example 1).
- **Foreign group inside window:** breaks the interval — an interval is the set of *all* window elements
  (Source 2 §Prelim).

### Independent cross-check (numbers brute-forced this session)

Brute force over all subsets (`/tmp/ci.py`, `/tmp/seq.py`), contiguity = positions span exactly |S|−1:

| Dataset | Result (independent) | Matches |
|---|---|---|
| P1=Id7, P2=(7 2 1 3 6 4 5), min 2 | {1,2},{1,2,3},{3,4,5,6},{4,5},{4,5,6},{1..6},{1..7} (7 sets) | **= paper Example 1** |
| same, min 4 | {3,4,5,6},{1..6},{1..7} | M4 |
| {2,3} interval of P2? | No (positions 1,3) | M2 |
| Didier T,S: {1,2} / {1,2,3,4} | No / Yes | matches Source 2 Example 1 |
| M3 A=(a,b,c,x),B=(x,c,b,a,x), min 3 | {a,b,c},{a,b,c,x},{b,c,x} | M3 (strengthened) |
| M5 A=(a,b,c),B=(a,z,b,c), min 3 | ∅ | M5 (strengthened) |
| S2 identical 1..5 ×3, min 3 | {1,2,3},{2,3,4},{3,4,5},{1,2,3,4},{2,3,4,5},{1,2,3,4,5} | S2 |
| S3 A=B=(1,2,3,4),C=(1,9,2,3,4), min 2 | {2,3},{3,4},{2,3,4} | S3 (strengthened) |
| C2 A=(1,2,3,4),B=(2,4,1,3), min 3 | {1,2,3,4} only | C2 |

Every expected value traces to either the paper's verbatim Example 1 or an independent brute force run
this session — none is a code echo.

### Findings / divergences (Stage A)

None. Description, formulas, definitions, edge cases, and worked example all match the cited primary
literature verbatim. The `maxGap` parameter is documented (Assumption 1) as API-shape only; the validated
contract is the strict gap-free common-interval model — consistent and correct.

**Stage A verdict: PASS.**

## Stage B — Implementation

### Code path reviewed

`src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/ComparativeGenomics.cs:914–1021`
(`FindConservedClusters` + helper `IsIntervalOf`).

### Formula realised correctly?

- Maps each genome to its ordered ortholog-group label sequence (genes with no group → window-breaking
  sentinel). ✓
- **Candidate generation (`:948–965`):** the group-set of every foreign-free contiguous window of genome 0.
  This is *complete*: every common interval is in particular an interval of genome 0, hence equals
  `Set` of some genome-0 window, and that window cannot contain a foreign (non-member) group. ✓
- **Filter (`:967–978`):** keep only candidates that are `IsIntervalOf` *every* genome. `IsIntervalOf`
  (`:1000–1021`) finds a window whose label set equals the target exactly (breaks on any foreign element),
  i.e. the `Set(T[i..j])` sequence-interval test — matches Didier et al. exactly, including paralogs. ✓
- **Size threshold:** `effectiveMin = max(minClusterSize, 2)` (`:929`, `MinCommonIntervalSize=2`),
  enforcing the `i < j` ≥2 rule. ✓
- **< 2 genomes → empty** (`:926`). ✓
- **Determinism:** sorted by size then lexicographic joined labels (`:980–986`). ✓
- **Null guards:** `ArgumentNullException.ThrowIfNull` on both args (`:920–921`). ✓

### Cross-verification vs code (via the test suite)

The strengthened exact-set tests now assert the independently brute-forced sets above; the full suite
passes, so the code produces **exactly** those sets for every dataset (not a superset/subset). This is the
strongest available confirmation that candidate completeness + the interval filter realise the definition.

### Variant/delegate consistency

Single public method; no `*Fast` / instance variant. MCP wrapper (`AnalysisTools.cs`) delegates to the
same method — not a separate algorithm, out of scope.

### Test quality audit (HARD gate)

Canonical file `tests/Seqeron/Seqeron.Genomics.Tests/ComparativeGenomics_FindConservedClusters_Tests.cs`
(12 tests). Old weak/duplicate `FindConservedClusters` tests were already removed from
`ComparativeGenomicsTests.cs` (confirmed: only a pointer comment remains).

- **Sourced, not code-echoed:** M1/M4/C2/S2 lock the *full* exact set against the paper's Example 1 or an
  independent brute force. P1 is a property test (every returned cluster contiguous in every genome) over a
  seeded random input — guards over-reporting generally.
- **Defects fixed this session (green-washing risk):** three tests used permissive `Contains.Item` /
  `Does.Not.Contain` where the exact set is known and computable:
  - **M3** (`RepeatedLabels_WindowStillMatches`): was `Contains.Item({a,b,c})` → rewritten to
    `Is.EquivalentTo({a,b,c},{a,b,c,x},{b,c,x})` (brute-forced).
  - **M5** (`ForeignGroupInsideWindow_BreaksCluster`): was `Does.Not.Contain({a,b,c})` → rewritten to
    `Is.Empty` (brute-forced: no size-≥3 common interval).
  - **S3** (`ConservedInTwoButSplitInThird`): was `Does.Not.Contain({1,2})`+`Contains.Item({2,3,4})` →
    rewritten to `Is.EquivalentTo({2,3},{3,4},{2,3,4})` (brute-forced).
  These now fail against an implementation that over- or under-reports, closing the green-wash gap.
- **Coverage:** the single public method and every Stage-A branch are exercised — golden vector (M1),
  one-genome-only negative (M2), paralog/repeated-label sequence path (M3), size filter (M4), foreign-group
  window break (M5), <2 genomes empty (S1), k=3 all-conserved (S2), k=3 split-in-one (S3), determinism (C1),
  trivial-whole-set-only (C2), null args (C3), property/INV-01 (P1).
- **Honest green:** FULL unfiltered suite = **6605 passed, 0 failed, 0 skipped-with-intent (1 benchmark
  skipped as designed)**; `dotnet build` 0 errors (4 pre-existing warnings, all in unrelated files —
  `ApproximateMatcher_EditDistance_Tests.cs` etc., untouched).

### Findings / defects (Stage B)

No implementation defect. Three test assertions were weak (could pass against an over-reporting
implementation); all three were rewritten to exact brute-forced sets and verified green. No code change
needed.

## Verdict & follow-ups

- **Stage A: PASS** — description matches primary literature verbatim; all numbers independently
  reproduced.
- **Stage B: PASS-WITH-NOTES** — implementation correct; three weak test assertions strengthened to exact
  sourced values (test-quality defect, fixed in-session).
- **End-state: ✅ CLEAN** — algorithm fully functional, tests lock sourced values, full suite green.
- **Test-quality gate: PASS** (after fixing M3/M5/S3).
