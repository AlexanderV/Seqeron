# Validation Report: PHYLO-COMP-001 — Tree Comparison (Robinson–Foulds, MRCA, Patristic)

- **Validated:** 2026-06-12   **Area:** Phylogenetics
- **Canonical method(s):** `PhylogeneticAnalyzer.RobinsonFouldsDistance(PhyloNode, PhyloNode)`,
  `PhylogeneticAnalyzer.FindMRCA(PhyloNode, string, string)`,
  `PhylogeneticAnalyzer.PatristicDistance(PhyloNode, string, string)`
- **Source:** `src/Seqeron/Algorithms/Seqeron.Genomics.Phylogenetics/PhylogeneticAnalyzer.cs`
- **Tests:** `tests/Seqeron/Seqeron.Genomics.Tests/PhylogeneticAnalyzer_TreeComparison_Tests.cs`
- **Stage A verdict:** PASS-WITH-NOTES
- **Stage B verdict:** PASS
- **End state:** CLEAN (no defect; full suite green)

---

## Stage A — Description

### Sources opened & what they confirm
- **Wikipedia, "Robinson–Foulds metric"** (https://en.wikipedia.org/wiki/Robinson%E2%80%93Foulds_metric).
  Confirms verbatim: RF "is defined as (A + B) where A is the number of partitions of data
  implied by the first tree but not the second tree and B is the number of partitions of data
  implied by the second tree but not the first tree." Partitions are obtained "by removing each
  branch" (edge → bipartition/split). Operation count = "the number of edges in T1 that are not
  in T2 plus the number of edges in T2 that are not in T1" — i.e. the symmetric difference,
  **both directions**. Rooted trees: "Rooted trees can be examined by attaching a dummy leaf to
  the root node."
- **Robinson & Foulds (1981)**, doi:10.1016/0025-5564(81)90043-2 — proved RF is a proper metric
  (identity, symmetry, triangle inequality).
- **Wikipedia, "Most recent common ancestor"** — MRCA = "the most recent individual from which
  all organisms of a set are inferred to have descended"; equivalently the deepest node that is
  an ancestor of all members (LCA/concestor).

### Formula check
- **Canonical (unrooted) RF** = number of nontrivial **bipartitions (splits)** present in one
  tree but not the other, summed both ways: `RF = |A\B| + |B\A|`. Trivial splits (a single leaf,
  or all leaves) are excluded. Identical trees → 0. Max for **unrooted** binary trees = `2(n−3)`.
- **This implementation deliberately uses the ROOTED CLADE (cluster) variant**, not unrooted
  bipartitions. Each internal node defines a clade = the set of taxa in its subtree;
  `RF_clade = |clades(T1) △ clades(T2)|`. Trivial clades (a single leaf, and the all-leaves
  root clade) are excluded. This is the rooted analogue of RF, equivalent (via the dummy-leaf
  construction in the Wikipedia source) to the unrooted RF of the trees with a dummy leaf added
  at the root. The evidence doc (`docs/Evidence/PHYLO-COMP-001-Evidence.md`, §7) and the source
  XML-doc state this scope explicitly. **The implementation choice is sourced and internally
  consistent — this is the NOTE: it is rooted clade distance, NOT unrooted bipartition RF; the
  two differ numerically and the spec/code correctly say which.**
- **Max RF (rooted variant)** = `2(n−2)`. A rooted binary tree on n leaves has n−1 internal
  nodes; excluding the root (all-leaves clade) leaves n−2 nontrivial clades. Dummy-leaf check:
  unrooted tree on n+1 leaves has max `2((n+1)−3) = 2(n−2)`. Consistent.
- **MRCA**: deepest node ancestral to both taxa; for a single taxon vs itself, the taxon node.
- **Patristic distance**: sum of branch lengths on the path through the MRCA =
  `d(x,MRCA) + d(y,MRCA)`.

### Edge-case semantics
- Identical trees → RF 0; completely different topologies → max (`2(n−2)` rooted). Defined.
- RF even: clade symmetric difference over the same leaf set is even (for binary trees). Defined.
- MRCA null root → null; non-existent taxon → null; MRCA(x,x) → x's leaf node. Defined & sourced.
- Patristic non-existent taxon → NaN; PD(x,x) = 0. Defined.

### Independent cross-check (hand computation)
Standard 4-taxa trees:
- Tree A `((A,B),(C,D))` → nontrivial clades {A,B},{C,D}.
- Tree B `((A,C),(B,D))` → nontrivial clades {A,C},{B,D}.
- `RF_clade = |{AB,CD}\{AC,BD}| + |{AC,BD}\{AB,CD}| = 2 + 2 = 4 = 2(n−2)` for n=4. ✓ (RF-M06)
- (For reference, the **unrooted** RF of the same two trees is 2 = 2(n−3) — single split each
  side — confirming the implementation is the rooted variant, as documented.)

3-taxa: `((A,B),C)` clade {A,B} vs `((A,C),B)` clade {A,C} → `RF = 1 + 1 = 2 = 2(3−2)`. ✓ (RF-M04/RF-S02)
Identical trees → 0. ✓ MRCA(A,B) in `((A,B),(C,D))` = node AB (2 taxa); MRCA(A,C) = root (4 taxa). ✓

### Findings / divergences
- **NOTE (not a defect):** the implementation computes **rooted clade/cluster distance**, not the
  classic unrooted bipartition RF. This is explicitly documented and the test expectations
  (RF=4 for the 4-taxa swap, max `2(n−2)`) are consistent with that choice. The spec correctly
  states which variant is used. A user expecting unrooted RF on these inputs would get a
  different number, so the divergence is documented here for clarity.

---

## Stage B — Implementation

### Code path reviewed
- `RobinsonFouldsDistance` (`PhylogeneticAnalyzer.cs:704`) → `GetClades`/`CollectClades`
  (`:713`–`:743`).
- `FindMRCA` (`:749`) → `FindMRCAInternal` (`:764`).
- `PatristicDistance` (`:785`) → `DistanceToTaxon` (`:796`).

### Formula realised correctly?
- **RF:** `GetClades` collects, for every internal node, the sorted set of subtree taxon names
  joined as a `"A|B|..."` key into a `HashSet<string>`. The nontrivial filter
  (`subtreeTaxa.Count > 1 && < totalLeaves`, `:737`) correctly excludes single-leaf and
  all-leaves (root) clades. `RobinsonFouldsDistance` returns
  `clades1.Except(clades2).Count() + clades2.Except(clades1).Count()` (`:709`) — the symmetric
  difference counted in **both directions**. ✓ The sorted-name key means two clades collide in
  the set iff they contain exactly the same taxon set (correct), and — unlike the previous
  "smaller-side canonical split" approach — it does **not** conflate complementary clades for
  n ≥ 5 (this defect was fixed per the spec's 2026-03-08 change log; verified by inspection).
- **MRCA:** standard post-order LCA. `FindMRCAInternal` returns the node where the two taxa are
  first found in different subtrees; the wrapper (`:758`) converts a leaf result with distinct
  taxa back to `null` (so one-missing / both-missing both yield null), while preserving the
  self-MRCA leaf for `taxon1 == taxon2`. ✓ Matches Stage-A semantics including MRCA-at-root and
  MRCA(leaf, itself).
- **Patristic:** locates MRCA, then sums branch lengths from MRCA down to each taxon
  (`DistanceToTaxon` adds the child branch length at each step); NaN when MRCA is null. ✓

### Cross-verification table recomputed vs code (via tests)
| Case | Expected (source/hand) | Code result |
|------|------------------------|-------------|
| RF identical 4-taxa | 0 | 0 ✓ |
| RF `((A,B),(C,D))` vs `((A,C),(B,D))` | 4 = 2(n−2) | 4 ✓ |
| RF 3-taxa different topology | 2 | 2 ✓ |
| RF even / non-negative (3 sizes) | 0,2,4 | match ✓ |
| MRCA(A,A) | leaf A | leaf A ✓ |
| MRCA(A,B) siblings | parent {A,B} | {A,B} ✓ |
| MRCA(A,C) cross-clade | root {A,B,C,D} | root ✓ |
| MRCA non-existent / null root | null | null ✓ |
| PD(A,B) 4-taxa | 1.0 | 1.0 ✓ |
| PD(A,C) 4-taxa | 5.0 | 5.0 ✓ |
| PD(C,D) 4-taxa | 2.0 | 2.0 ✓ |
| PD non-existent | NaN | NaN ✓ |

### Variant/delegate consistency
- `Bootstrap` reuses the same `GetClades` clade representation for support counting — consistent
  with the RF clade definition. ✓

### Test quality audit
- 34 tests in the unit; MUST/SHOULD tests assert **exact sourced values** (RF 0/2/4, PD 1.0/3.0/5.0/2.0),
  not just "≥ 0" or "no throw". Edge cases (null root, non-existent taxon one/both, self-MRCA,
  cross-clade root MRCA, evenness, max bound) are covered. Deterministic. ✓

### Findings / defects
- None. Off-by-direction (counting only one side), trivial-split inclusion, and
  complementary-clade collision — the three classic RF defects called out in the protocol — are
  all absent: both directions counted, trivial clades excluded, direct-clade keys avoid
  complement collisions.

---

## Verdict & follow-ups
- **Stage A: PASS-WITH-NOTES** — definitions and worked examples confirmed against Wikipedia and
  Robinson & Foulds (1981). The single note is that the implementation uses the **rooted clade
  (cluster) distance** variant (max `2(n−2)`), not unrooted bipartition RF (max `2(n−3)`); this
  is explicitly and correctly documented in spec, evidence, and code.
- **Stage B: PASS** — code faithfully realises the rooted-clade RF, MRCA, and patristic-distance
  definitions; all cross-check values reproduced; classic RF defects absent.
- **End state: CLEAN.** No code changes. Unit tests: 34 passed. Full suite:
  `Seqeron.Genomics.Tests` 4461 passed / 0 failed.
- **Follow-up (optional, non-blocking):** if an unrooted bipartition RF is ever required, add it
  as a separate method rather than changing the documented rooted-clade behaviour.
