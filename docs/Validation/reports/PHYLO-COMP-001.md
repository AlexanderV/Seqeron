# Validation Report: PHYLO-COMP-001 — Tree Comparison (Robinson–Foulds + MRCA)

- **Validated:** 2026-06-17   **Area:** Phylogenetics
- **Re-validation scope:** (1) unrooted-bipartition RF added beside rooted-clade RF
  (`CalculateUnrootedRobinsonFoulds` / `CalculateNormalizedUnrootedRobinsonFoulds`, commit 4846294);
  (2) N-ary (multifurcating) refactor of `PhyloNode` to a `Children` list, with RF + MRCA traversing
  N children (commit c4f0190).
- **Canonical method(s):** `RobinsonFouldsDistance` (rooted clade), `CalculateUnrootedRobinsonFoulds`
  + `CalculateNormalizedUnrootedRobinsonFoulds` (unrooted bipartition), `FindMRCA`, `PatristicDistance`.
- **Stage A verdict:** ✅ PASS
- **Stage B verdict:** ✅ PASS
- **End-state:** ✅ CLEAN (no defect; existing tests pass the hard mutation gate; 0 code/test change)

---

## Stage A — Description

### Sources opened (this session)
- **Robinson & Foulds (1981)** "Comparison of phylogenetic trees", *Math. Biosci.* 53:131–147 —
  RF = symmetric difference of the sets of **bipartitions/splits** (partitions induced by removing
  internal edges); trivial splits (terminal/leaf edges) carry no information and are excluded.
- **Wikipedia "Robinson–Foulds metric"** (WebFetched) — verbatim: RF "is defined as (A + B) where A is
  the number of partitions of data implied by the first tree but not the second tree and B is the
  number of partitions of data implied by the second tree but not the first tree"; "straightforward to
  normalize RF distances so they range between zero and one."
- **Rice CS comp571 tree-metrics worked example** (cited in the test file) — two 5-taxon trees differing
  by one internal split each → unrooted RF = 1 + 1 = 2.

### Formula check
- **Unrooted RF** = |splits(T1) △ splits(T2)| over non-trivial bipartitions. ✔ matches Wikipedia A+B.
- **Rooted-clade RF** = |clades(T1) △ clades(T2)|, clades = subtree taxon-sets, root + leaves excluded.
  This is the standard rooted variant (root-sensitive). ✔
- **Multifurcation:** contracting an internal edge removes exactly that edge's split from Σ(T) — an
  unresolved (polytomy) node contributes FEWER non-trivial splits/clades. ✔ (R&F 1981 edge-contraction).
- **Normalized unrooted RF denominator 2n−6 = 2(n−3):** a binary unrooted tree on n leaves has n−3
  internal edges (non-trivial splits), so the symmetric difference is at most 2(n−3). ✔ (independently
  derived; matches the code's denominator). For n=3 the denominator is 0 (one unrooted topology) → 0.
- **Parity:** RF is always even for two binary trees on the same taxa (each differing split in one tree
  forces a differing split in the other). ✔

### Independently hand-derived cases (required by the prompt)
All derived **from the bipartition/clade definition by hand**, not from code output:

| Case | Trees | Hand-derived RF |
|------|-------|-----------------|
| (a) identical | T1 vs T1 | rooted 0, unrooted 0 |
| (b) single NNI | `((B,C),(A,(D,E)))` vs `((A,B),(C,(D,E)))` — splits {BC\|ADE, DE\|ABC} vs {AB\|CDE, DE\|ABC}; shared DE\|ABC; BC\|ADE & AB\|CDE unique | **unrooted RF = 2** |
| (c) **root-invariance** | X=`((A,B),(C,(D,E)))`, Y=`(((A,B),C),(D,E))` — SAME unrooted tree, different root edge. Unrooted splits of BOTH = {AB\|CDE, ABC\|DE} → **unrooted RF = 0**. Rooted clades: X={AB,DE,CDE}, Y={AB,DE,ABC} → {CDE} only-X + {ABC} only-Y → **rooted RF = 2 ≠ 0** | **unrooted 0 / rooted 2** |
| (d) **collapse** (binary vs multifurcating, one edge contracted) | binary `(((A,B),C),(D,E))` vs collapsed `((A,B,C),(D,E))`. Rooted clades binary={AB, ABC, DE}, collapsed={ABC, DE} → lost {AB} → **rooted RF = 1**. Unrooted splits binary={AB\|CDE, ABC\|DE}, collapsed={ABC\|DE} → lost {AB\|CDE} → **unrooted RF = 1** | **rooted 1 / unrooted 1** (= splits lost) |
| (e) **MRCA over ≥3-child polytomy** | `((A,B,C),(D,E))`: MRCA(A,B)=MRCA(B,C)=MRCA(A,C)= the 3-child polytomy node; MRCA(A,D)=MRCA(C,E)= root; MRCA(A,Z)=null (Z missing) | as stated |
| (extra) max 5-taxon | `((A,B),(C,(D,E)))` vs `(((A,C),E),(B,D))` — disjoint internal splits → **unrooted RF = 4 = 2(n−3)**; normalized 4/4 = 1 | 4 / 1.0 |

### Edge-case semantics
- < 3 leaves → `ArgumentException` (no non-trivial bipartition possible). ✔
- Mismatched leaf sets → `ArgumentException` (RF defined only on a common taxon set). ✔
- null tree → `ArgumentNullException`. ✔
- n=3 normalized → denominator 0 → defined as 0. ✔

**Stage A verdict: ✅ PASS** — both metrics and the multifurcation handling are mathematically sound;
the normalization denominator and the root-invariance rationale are correct.

---

## Stage B — Implementation

### Code path reviewed
- `RobinsonFouldsDistance` + `GetClades`/`CollectClades` — `PhylogeneticAnalyzer.cs:845–861, 1011–1036`.
- `CalculateUnrootedRobinsonFoulds` + `GetBipartitions`/`CollectBipartitions`/`CanonicalSplitKey`
  — `PhylogeneticAnalyzer.cs:899–1009`.
- `CalculateNormalizedUnrootedRobinsonFoulds` — `:936–946`.
- `FindMRCA`/`FindMRCAInternal` — `:1042–1086`; `PatristicDistance`/`DistanceToTaxon` — `:1091–1119`.
- `ParseNewick` parses multifurcations faithfully (`:662–726`) and populates `node.Taxa`.

### Formula realised correctly? (evidence)
- Unrooted RF reads the rooted binary tree as unrooted by forming the subtree side at every internal
  node and **canonicalising the split to its smaller side** (`CanonicalSplitKey`), so the duplicate
  split contributed by the two root children dedups in the `HashSet` and the root itself
  (side = total) is excluded. Non-trivial filter `side ≥ 2 && total−side ≥ 2`. ✔ This is the correct
  root-invariant reduction.
- Multifurcation: `CollectBipartitions`/`CollectClades`/`FindMRCAInternal`/`DistanceToTaxon`/`GetLeaves`/
  `CalculateTreeLength` all iterate `node.Children` (N-ary), so a polytomy contributes one
  side/clade — exactly the edge-contraction semantics. ✔

### Cross-verification vs code
Ran the unit fixture (`PhylogeneticAnalyzer_TreeComparison_Tests`, 46 tests) — all green. Every value
above is asserted exactly: NNI=2, root-invariance unrooted 0 / rooted 2, collapse rooted 1 + unrooted 1,
max 4, normalized 0.5 / 1.0 / 0, MRCA-over-polytomy, polytomy stats (5 leaves, length 2.1, A→C PD 0.4).

### Test-quality audit — HARD mutation gate (all run this session)
1. **Root-invariance gate** — mutated `CalculateUnrootedRobinsonFoulds` to use `GetClades` (rooted-clade
   masquerade). Result: **5 tests fail**, including `URF-ROOT`
   (`UnrootedRobinsonFoulds_DifferByRootPositionOnly_IsZero_WhileRootedRfIsNonZero`), the single-NNI,
   maximally-different, and both normalized tests. ✔ A rooted impl of the unrooted method DOES fail.
2. **Multifurcation gate** — mutated every N-ary traversal (`CollectClades`, `CollectBipartitions`,
   `FindMRCAInternal`, `DistanceToTaxon`) to `Children.Take(2)` (ignore 3rd+ child). Result: **all 4
   multifurcation tests fail** — `RF-MULTI-ROOTED`, `RF-MULTI-UNROOTED`, `MRCA-POLYTOMY`,
   `STATS-POLYTOMY`. ✔ No test passes if multifurcation handling is wrong (no green-washing); the
   collapse-RF tests assert exactly the number of splits/clades lost.
   Source restored byte-for-byte after each mutation (0 `MUTATION` markers remain; clean `git diff`).

### Findings / defects
**None.** Code is correct; tests assert exact sourced integers/normalized values and survive the gate.

---

## Verdict & follow-ups
- **Stage A: ✅ PASS.  Stage B: ✅ PASS.  End-state: ✅ CLEAN.**
- No code change, no test change. Full unfiltered `dotnet test`: **6779 passed / 0 failed**, build 0 errors.
