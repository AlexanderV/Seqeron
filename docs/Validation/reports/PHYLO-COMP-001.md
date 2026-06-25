# Validation Report: PHYLO-COMP-001 — Tree Comparison (Robinson–Foulds + MRCA + Patristic)

- **Validated:** 2026-06-24   **Area:** Phylogenetics
- **Re-validation scope:** independent re-confirmation of (1) rooted-clade RF
  (`RobinsonFouldsDistance`), (2) unrooted-bipartition RF
  (`CalculateUnrootedRobinsonFoulds` / `CalculateNormalizedUnrootedRobinsonFoulds`, commit 4846294),
  (3) the N-ary (multifurcating) `PhyloNode.Children` refactor with RF/MRCA/patristic traversing
  N children (commit c4f0190).
- **Canonical method(s):** `RobinsonFouldsDistance` (rooted clade),
  `CalculateUnrootedRobinsonFoulds` + `CalculateNormalizedUnrootedRobinsonFoulds` (unrooted
  bipartition), `FindMRCA`, `PatristicDistance`.
- **Stage A verdict:** ✅ PASS
- **Stage B verdict:** ✅ PASS
- **End-state:** ✅ CLEAN (no defect; no code/test change this session)

---

## Stage A — Description

### Sources opened (this session)
- **Wikipedia "Robinson–Foulds metric"** (WebFetched 2026-06-24) — verbatim: RF "is defined as
  (A + B) where A is the number of partitions of data implied by the first tree but not the second
  tree and B is the number of partitions of data implied by the second tree but not the first
  tree"; "it is straightforward to normalize RF distances so they range between zero and one";
  "Rooted trees can be examined by attaching a dummy leaf to the root node." (The article carries a
  neutrality notice; the A+B definition and normalization claim are uncontroversial and match the
  primary paper.)
- **Robinson & Foulds (1981)** "Comparison of phylogenetic trees", *Math. Biosci.* 53:131–147 —
  RF = symmetric difference of the sets of bipartitions/splits induced by removing internal edges;
  trivial splits (terminal/leaf edges) carry no information and are excluded; contracting an
  internal edge removes exactly that edge's split from Σ(T).
- **Rice CS comp571 tree-metrics worked example** (cited in the test file) — two 5-taxon trees
  differing by one internal split each → unrooted RF = 1 + 1 = 2.

### Formula check
- **Unrooted RF** = |splits(T1) △ splits(T2)| over non-trivial bipartitions. ✔ matches Wikipedia A+B.
- **Rooted-clade RF** = |clades(T1) △ clades(T2)|, clades = subtree taxon-sets with root and leaves
  excluded. Standard root-sensitive rooted variant. ✔
- **Max unrooted RF = 2(n−3) = 2n−6**: a binary unrooted tree on n leaves has n−3 internal edges
  (non-trivial splits), so the symmetric difference is at most 2(n−3). Independently derived; matches
  the normalization denominator in code. For n=3 the denominator is 0 (one unrooted topology) → 0. ✔
- **Multifurcation:** contracting an internal edge removes exactly that edge's split/clade — a
  polytomy contributes FEWER non-trivial splits/clades (edge-contraction semantics, R&F 1981). ✔
- **Parity:** RF is always even for two binary trees on the same taxa. ✔

### Independently hand-derived cases (from the bipartition/clade definitions, NOT from code output)

| Case | Trees | Hand-derived RF |
|------|-------|-----------------|
| (a) identical | `((B,C),(A,(D,E)))` vs itself | rooted 0, unrooted 0 |
| (b) single NNI | `((B,C),(A,(D,E)))` vs `((A,B),(C,(D,E)))` — splits {BC\|ADE, DE\|ABC} vs {AB\|CDE, DE\|ABC}; shared DE\|ABC; BC\|ADE & AB\|CDE unique | **unrooted RF = 2** |
| (c) **root-invariance** | X=`((A,B),(C,(D,E)))`, Y=`(((A,B),C),(D,E))` — SAME unrooted tree, different root edge. Unrooted splits of BOTH = {AB\|CDE, ABC\|DE} → **unrooted RF = 0**. Rooted clades X={AB,DE,CDE}, Y={AB,DE,ABC} → △={CDE,ABC} → **rooted RF = 2 ≠ 0** | **unrooted 0 / rooted 2** |
| (d) **collapse** | binary `(((A,B),C),(D,E))` vs collapsed `((A,B,C),(D,E))`. Rooted clades binary={AB,ABC,DE}, collapsed={ABC,DE} → lost {AB} → **rooted RF = 1**. Unrooted splits binary={AB\|CDE, ABC\|DE}, collapsed={ABC\|DE} → lost {AB\|CDE} → **unrooted RF = 1** | **rooted 1 / unrooted 1** |
| (e) **max 5-taxon** | `((A,B),(C,(D,E)))` vs `(((A,C),E),(B,D))` — splits {AB\|CDE, DE\|ABC} vs {AC\|BDE, ACE\|BD}, disjoint → **unrooted RF = 4 = 2(n−3)**; normalized 4/4 = 1 | 4 / 1.0 |
| (f) **MRCA over polytomy** | `((A,B,C),(D,E))`: MRCA(A,B)=MRCA(B,C)=MRCA(A,C)= the 3-child node; MRCA(A,D)=MRCA(C,E)= root; MRCA(A,Z)=null (Z missing) | as stated |

### Edge-case semantics
- < 3 leaves → `ArgumentException` (no non-trivial bipartition possible). ✔
- Mismatched leaf sets → `ArgumentException` (RF defined only on a common taxon set). ✔
- null tree → `ArgumentNullException`. ✔
- n=3 normalized → denominator 2n−6 = 0 → defined as 0. ✔
- MRCA: null root → null; non-existent taxon → null. Patristic: same taxon → 0; non-existent → NaN. ✔

**Stage A verdict: ✅ PASS** — both RF metrics, the multifurcation handling, the normalization
denominator, and the root-invariance rationale are mathematically sound and externally sourced.

---

## Stage B — Implementation

### Code path reviewed (`src/Seqeron/Algorithms/Seqeron.Genomics.Phylogenetics/PhylogeneticAnalyzer.cs`)
- `RobinsonFouldsDistance` + `GetClades`/`CollectClades` — `:879–895, 1045–1070`.
  Clade key = sorted, `|`-joined subtree taxa; non-trivial filter `Count > 1 && Count < totalLeaves`
  (excludes leaves and root). Symmetric difference via `Except`. ✔
- `CalculateUnrootedRobinsonFoulds` + `GetBipartitions`/`CollectBipartitions`/`CanonicalSplitKey`
  — `:933–1043`. Null/leaf-count/leaf-set-equality guards present; split canonicalised to the
  smaller side (lexicographic tie-break) so {S, complement} dedups in the `HashSet` regardless of
  root placement; the root side (= total) is excluded by the `total − side ≥ 2` non-trivial filter. ✔
- `CalculateNormalizedUnrootedRobinsonFoulds` — `:970–980`: `rf / (2n−6)`, returns 0 when
  denominator ≤ 0 (n ≤ 3). ✔
- `FindMRCA`/`FindMRCAInternal` — `:1076–1120`: N-ary; ≥2 children with a hit ⇒ this node is MRCA;
  leaf-result with distinct taxa ⇒ one taxon missing ⇒ null. ✔
- `PatristicDistance`/`DistanceToTaxon` — `:1125–1153`: sum of branch lengths from MRCA down to each
  taxon over N children; NaN when MRCA null or taxon absent. ✔
- `GetLeaves`/`CalculateTreeLength`/`GetTreeDepth` — all iterate `node.Children` (N-ary). ✔

### Formula realised correctly? (evidence)
- Unrooted RF reads the binary tree as unrooted by canonicalising each internal node's subtree side
  to the smaller side, so the duplicate split contributed by the two root children dedups and the
  root itself is excluded — the correct root-invariant reduction.
- Multifurcation: every traversal is N-ary, so a polytomy contributes exactly one side/clade —
  precisely the edge-contraction semantics.

### Cross-verification table recomputed vs code (unit fixture, 46 tests, all green)
| Case | Method | Hand value | Test assertion |
|------|--------|-----------|----------------|
| identical | unrooted/rooted | 0 / 0 | URF-M01, RF-M01 = 0 ✔ |
| single NNI | unrooted | 2 | URF-M02 = 2 ✔ |
| root-invariance | unrooted / rooted | 0 / 2 | URF-ROOT = (0, 2) ✔ |
| collapse | rooted / unrooted | 1 / 1 | RF-MULTI-ROOTED=1, RF-MULTI-UNROOTED=1 ✔ |
| max 5-taxon | unrooted / normalized | 4 / 1.0 | URF-M04=4, URF-N01=1.0 ✔ |
| single NNI normalized | normalized | 0.5 | URF-N01 = 0.5 ✔ |
| n=3 normalized | normalized | 0 | URF-N02 = 0 ✔ |
| 3-taxa diff rooted | rooted | 2 | RF-M04/RF-S02 = 2 ✔ |
| 4-taxa max rooted | rooted | 4 | RF-M06 = 4 ✔ |
| MRCA over polytomy | FindMRCA | polytomy/root/null | MRCA-POLYTOMY ✔ |
| PD 3-taxa | patristic | 1.0,3.0,3.0 | PD-M05 ✔ |
| PD 4-taxa | patristic | 1.0,5.0,2.0 | PD-S01 ✔ |
| polytomy stats | leaves/length/depth/PD | 5 / 2.1 / 2 / 0.3,0.4 | STATS-POLYTOMY ✔ |

### Variant/delegate consistency
Normalized unrooted RF = raw unrooted RF / (2n−6); both share `GetBipartitions`. Rooted and unrooted
share `GetLeaves`. Symmetric (RF-M02, URF-M03, MRCA symmetry, PD-M03) confirmed. ✔

### Test-quality audit
Assertions check exact sourced integers/normalized doubles (RF=0/1/2/4, normalized 0/0.5/1.0), not
"no-throw" tautologies. The root-invariance test (`URF-ROOT`) asserts BOTH unrooted 0 AND rooted 2 on
the same pair — locks the distinguishing property. Collapse tests assert the exact count of splits/
clades lost (1). Multifurcation tests assert the polytomy node really has 3 children before checking
behaviour. Edge cases (null, <3 leaves, mismatched leaf sets, missing taxon) all asserted. Prior
session ran hard mutation gates (rooted-masquerade and `Children.Take(2)`) that broke the relevant
tests; logic is unchanged since, so those guarantees still hold.

### Findings / defects
**None.**

---

## Verdict & follow-ups
- **Stage A: ✅ PASS.  Stage B: ✅ PASS.  End-state: ✅ CLEAN.**
- No code change, no test change this session. Unit fixture
  (`PhylogeneticAnalyzer_TreeComparison_Tests`): **46 passed / 0 failed**; build 0 warnings, 0 errors.
- Full unfiltered suite not run (no code touched). Independent hand-derivations match every asserted
  value, including the load-bearing root-invariance case (unrooted RF 0 vs rooted RF 2) and the
  max RF = 2(n−3) bound.
