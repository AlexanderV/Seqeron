# Evidence Artifact: ONCO-PHYLO-001

**Test Unit ID:** ONCO-PHYLO-001
**Algorithm:** Tumor Phylogeny Reconstruction — clonal tree from CCF clusters (sum rule + lineage-precedence rule)
**Date Collected:** 2026-06-15

---

## Online Sources

### Popic V, Salari R, Hajirasouliha I, Kashef-Haghighi D, West RB, Batzoglou S (2015) — LICHeE

**URL:** https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0647-8 (DOI 10.1186/s13059-015-0647-8); full text PDF retrieved at https://arxiv.org/pdf/1412.8574
**Accessed:** 2026-06-15
**Authority rank:** 1 (peer-reviewed paper, *Genome Biology*)

**How retrieved:** WebSearch query "LICHeE Popic 2015 Genome Biology fast scalable inference multi-sample cancer lineages VAF perfect phylogeny constraints" → fetched the arXiv preprint PDF (arxiv.org/pdf/1412.8574); the PDF was parsed (pypdf) and the constraint passages quoted below were extracted from the rendered text.

**Key Extracted Points:**

1. **Perfect-phylogeny ordering constraints (verbatim):** "Firstly, (1) a mutation present in a given set of samples cannot be a successor of a mutation that is present in a smaller subset of these samples … (2) a given mutation cannot have a VAF higher than that of its predecessor mutation (except due to CNVs), since all cells containing this mutation will also contain the predecessor. Finally, (3) the sum of the VAFs of mutations disjointly present in distinct subclones cannot exceed the VAF of a common predecessor mutation present in these subclones, since the subclones with the descendent mutations must contain the parent mutations."
2. **Ancestor ≥ descendant edge rule (verbatim, Eq. 2):** an edge `(u → v)` is added only if for all samples i: "(1) `u.VAFi ≥ v.VAFi − ϵuv` and (2) if `u.VAFi = 0`, `v.VAFi = 0`", where `ϵuv` is the VAF noise error margin.
3. **Sum rule (verbatim, Eq. 5):** a valid lineage tree T "must be a spanning tree of the network that satisfies the following requirement `∀ nodes u∈T : ∀i∈ samples : Σ_{v s.t. (u→v)∈T} v.VAFi ≤ u.VAFi + ϵ`. That is, the sum of the VAF centroids of all the children must not exceed the centroid of the parent."
4. **Why inequality, not equality:** "we use inequality here since our method does not require all the true lineage branches to have been observed."

### Zheng L, Dang H, Niknafs N, et al. (2022) — PICTograph

**URL:** https://pmc.ncbi.nlm.nih.gov/articles/PMC9344857/ (DOI 10.1093/bioinformatics/btac367, *Bioinformatics* 38(15):3677–3683)
**Accessed:** 2026-06-15
**Authority rank:** 1 (peer-reviewed paper, *Bioinformatics*)

**How retrieved:** WebSearch query "CITUP tumor phylogeny CCF sum rule pigeonhole ancestor descendant cancer cell fraction clonal tree" → WebFetch of the PMC full-text article with a prompt extracting the sum condition and lineage-precedence rule.

**Key Extracted Points:**

1. **Sum Condition (CCF form, verbatim):** "the CCF of an ancestral clone must be greater than or equal to the sum of CCFs of its descendants." Trees where "the sum of the CCFs of the descendants of a parent node exceeded the parent node's CCF by more than ε₂ were excluded" (default ε₂ = 0.2).
2. **Lineage Precedence (CCF form, verbatim):** "the CCF of any mutation cannot exceed the CCF of its ancestor." Operationally the descendant cluster must be present in a subset of the samples in which the ancestral cluster is present, and the descendant CCF is at most ε₁ greater than the ancestral CCF in each sample (default ε₁ = 0.1).
3. **Sum Condition generalizes the pigeonhole principle:** the constraint that the summed children CCF cannot exceed the parent CCF is the cell-fraction analogue of the pigeonhole principle applied per node.

### LICHeE reference implementation (viq854/lichee, MIT) — opened 2026-09-28 (B24 F18)

**URL:** raw.githubusercontent.com/viq854/lichee/master/LICHeE/src/lineage/{PHYNetwork,PHYTree,PHYNode,PHYEdge,Parameters,LineageEngine,SNVDataStore,SNVGroup,AAFClusterer}.java; binary LICHeE/release/lichee.jar (run under OpenJDK 21 through a harness calling `PHYNetwork` → `getLineageTrees` → `evaluateLineageTrees`, the pipeline of `LineageEngine.buildLineage` steps 4–6, with `-cp` settings `VAF_MAX = 1.0`, `VAF_ERROR_MARGIN = ε`).

**Key extracted points:**
1. Constraint network (`PHYNetwork` constructor): nodes levelled by number of samples present; within-profile pairs, level → next non-empty lower level, orphan nodes → closest higher level (from level+2) else root; `checkAndAddEdge` keeps one direction (smaller one-sided VAF excess; ties → n2→n1). Complete network (`-c`, `ALL_EDGES`) is the fallback when no tree is found.
2. Tree search: Gabow & Myers (1978) enumeration (`grow`) with `PHYTree.checkConstraint` = sum rule `Σ_children > u + VAF_ERROR_MARGIN ⇒ reject`; caps `MAX_NUM_TREES = 100000`, `MAX_NUM_GROW_CALLS = 1e8`.
3. Ranking: `PHYTree.computeErrorScore` = √(Σ_nodes Σ_samples max(0, Σ_children − u)²); `Collections.sort` (stable); top tree = index 0. QP consistency check off by default (`NUM_TREES_FOR_CONSISTENCY_CHECK = 0`).
4. No valid tree ⇒ LICHeE reports none (after `fixNetwork` removal of non-robust clusters and the `ALL_EDGES` retry) — it never returns a sum-rule-violating tree. (`fixNetwork` ported in F43, see below.)

### LICHeE per-cluster error margins — opened 2026-10-10 (B24 F42)

**Source:** github.com/viq854/lichee master `26c2a7018aa3c5f777bf5c894e440017f7e50e2b`: `PHYNetwork.getAAFErrorMargin` (l. 237–263), `checkAndAddEdge` (l. 180–232), `PHYTree.checkConstraint`, `PHYNode.getStdDev`, `AAFClusterer.Cluster.recomputeCentroidAndStdDev` and `AAFClusterer.em` (cluster SD).

1. `getAAFErrorMargin(from, to, i)`: `se_parent = 1.96·from.getStdDev(i)/√|from.members|` (root: `VAF_ERROR_MARGIN`), same for the child; returns `se_parent + se_child` if it exceeds `VAF_ERROR_MARGIN`, else `VAF_ERROR_MARGIN` — i.e. `max(ε, se_u + se_v)`; ε is a floor, not added (except that the root's own se is ε). The `STATIC_ERROR_MARGIN` switch is commented out: LICHeE always uses it.
2. Used only in `checkAndAddEdge` (`comp_12 += n1.AAF(i) ≥ n2.AAF(i) − margin(n1,n2,i)`, and symmetric); the one-sided error totals that choose the orientation use no margin. `PHYTree.checkConstraint` (sum rule) uses the static `VAF_ERROR_MARGIN` (`errMargin += getAAFErrorMargin` is commented out).
3. SD = population SD (divisor n) of the member AAFs per sample (`Math.sqrt(Σ (x − mean)²/members.size())`, both in `recomputeCentroidAndStdDev` and in `em`); centroid = Σ/n. `PHYNode.getStdDev` returns 0 for a sample outside the presence profile; the root's SD is 0. With a clusters file (`--clustersFile`) LICHeE sets SD = 0 (`c.setStdDev(new double[numSamples])`), i.e. the static margin.
4. Harness (`F4243Harness.java`, scratch, compiled against `lichee.jar` + `LICHeE/lib/*.jar`, OpenJDK 21): builds one `SNVGroup` per presence profile (first-appearance order) from member CCF rows, one `Cluster` per input cluster via `recomputeCentroidAndStdDev`, robustness as in `SNVGroup.setSubPopulations`, then runs `LineageEngine.buildLineage` steps 4–6 verbatim (`VAF_MAX = 1`, `VAF_ERROR_MARGIN = ε`).
5. Cross-check: 9 000 random member-level inputs (1–3 samples, 2–6 clusters × 2–5 robust members, member noise SD 0.01–0.2, ε ∈ {0, 0.02, 0.05, 0.1}; seeds 1 and 2): 9 000 / 9 000 identical (feasibility, network mode, tree count, top tree, error score); 57 594 / 57 594 centroids/SDs bit-identical. Margins changed the outcome vs the static ε on 16 inputs.

| Fixture | ε | lichee.jar with margins | lichee.jar static (SD 0) |
|---------|---|-------------------------|--------------------------|
| p1: A = {0.44×4},{0.56×4} (mean 0.5, SD 0.060000000000000026); B = 2 × [0.46,0.46,0.46,0.57] | 0.05 | 0 valid trees | 1 tree root→B→A, error 0.06928203230275505 |
| t02150 (3 clusters, 2 samples) | 0 | complete network, 2 trees, root→{1,2,3}, error 0 | 1 tree, 2→1, root→{2,3} |
| t01061 (4 clusters, 3 samples) | 0.1 | 0 valid trees | 1 tree root→2→3→4→1, error 0.10016236818286589 |
| t00298 (6 clusters, 2 samples) | 0 | complete network, 28 trees, root→{1,3,4,5}, 4→2, 1→6, error 0 | 5 trees, root→{1,4}, 4→2, 2→3, 1→6, 6→5 |

(Member rows in `OncologyAnalyzer_ReconstructPhylogenyClusterSummaries_Tests.cs`.) A margin-admitted edge u→v with `v.CCF[i] − u.CCF[i] > ε` can never pass the static-ε sum rule, so the margins act through orientation (both directions pass → the smaller one-sided excess wins) and orphan attachment.

### LICHeE fixNetwork — opened 2026-10-10 (B24 F43)

**Source:** same commit: `LineageEngine.buildLineage` l. 110–124, `PHYNetwork.fixNetwork`, `SNVGroup.setSubPopulations` (cluster robustness), `SNVGroup.removeCluster`, `SNVGroup.equals` (no `hashCode`), `SNVDataStore` (`SNVEntry.isRobust`, clusters-file loader), `Parameters` (`MIN_ROBUST_CLUSTER_SUPPORT = 2`, `MIN_VAF_PRESENT = MAX_VAF_ABSENT = 0.005`).

1. Loop: if no tree, `do { numNodes = net.numNodes; net = net.fixNetwork(); trees = net.getLineageTrees(); delta = numNodes − net.numNodes; } while (delta ≠ 0 && no trees)`; if still none: `ALL_EDGES = true; net = new PHYNetwork(groups, …)` — built from the same `SNVGroup` objects, which `removeCluster` mutated, so dropped clusters stay dropped.
2. `fixNetwork`: iterates `nodesById.values()` (Integer keys ⇒ ascending node id), picks the first non-robust cluster, replacing it only when `members.size()` is strictly smaller ⇒ smallest non-robust cluster, ties → lowest node id. Removes one cluster per call and rebuilds from `new ArrayList<>(HashSet<SNVGroup>)`.
3. Robustness: cluster robust ⇔ ≥ 2 robust member SNVs (`setSubPopulations`); SNV robust unless some sample VAF ∈ [`MAX_VAF_ABSENT`, `MIN_VAF_PRESENT`) — an empty band at the defaults, so with SNV input every member is robust and robust ⇔ n ≥ 2. Clusters read from `--clustersFile` are never `setRobust()` ⇒ all removable.
4. Group order after a rebuild = `HashSet<SNVGroup>` iteration (identity hash codes); the port keeps first-appearance order = LICHeE with a constant identity hash (`-XX:+UnlockExperimentalVMOptions -XX:hashCode=2`).
5. Cross-check (harness as in F42; `MIN_ROBUST_CLUSTER_SUPPORT` rule applied to per-member robust flags): 4 000 random inputs (1–3 samples, 2–6 clusters × 1–3 members, ~10 % non-robust members, seed 3; 1 958 with a removal): 4 000 / 4 000 identical with `hashCode=2`; 3 844 / 4 000 with the default identity hash (the 156 differ between LICHeE's own two runs).

| Fixture (ε) | Clusters (member rows; robust members) | lichee.jar (both hash settings) |
|-------------|----------------------------------------|---------------------------------|
| f1 (0) | A [0.5,0.5]×1; B [0.55,0]×1 | removed [1]; 1 tree root→B, error 0 |
| f2 (0) | A [0.5,0.5]×3 (1 robust); B [0.55,0]×1 | removed [2]; root→A |
| f3 (0.02) | F18 complete-network clusters ×2 each; D [0.01,0.01,0.01]×1 | removed [4]; complete network; root→{1,2}, 2→3; error 0.003292532308117998 |
| f4 (0) | A [0.5,0.5]×2, B [0.55,0]×3, no robust member (clusters file) | removed [1]; root→B |
| f4r (0) | same, all members robust | no tree |
| f6 (0) | A [0.5,0.5]×2; B [0.55,0]×1; C [0,0.55]×1 | removed [2, 3]; root→A |
| f7 (0) | A [0.5,0.5]×2; B [0.55,0]×1; C [0,0.55]×2; D [0,0.3]×1 | removed [2, 4]; complete network; no tree |

### LICHeE default ε in `-cp` mode and the grow-call cap — opened 2026-10-10 (B24 F44)

**Source:** same commit: `Parameters.java` (`VAF_ERROR_MARGIN = 0.1`, `VAF_MAX = 0.5`, `MAX_NUM_GROW_CALLS = 100000000`), `LineageEngine.java` option parsing l. 369–376 (`-cp`: `CP = true; VAF_MAX = 1.0; MAX_ALLOWED_VAF = 1.0`; `-e` only when given: `VAF_ERROR_MARGIN = parse(-e)`; option help "VAF error margin (default: 0.1)"), `PHYNetwork.grow` l. 443–530.

1. **Default ε:** `-cp` does not touch `VAF_ERROR_MARGIN`, so LICHeE's CCF-input default is ε = 0.1. The library default `DefaultPhylogenyTolerance` is now 0.1 (was 0).
2. **Cap:** `grow` increments `numGrowCalls` on entry; after a constraint-passing edge it returns when `numGrowCalls ≥ MAX_NUM_GROW_CALLS`, *without* undoing that edge or re-adding its `ff` edges, and its caller resumes after its own `grow(t)` (pop/restore, remove edge, bridge test, next stack edge). Every later recursion is behind the same check and the counter never decreases, so no further `grow` call is entered and no tree is added: `spanningTrees` = the trees completed in calls 1..cap (a tree completed *at* call cap is kept). The corrupted `t`/`edges` state is never read afterwards (`evaluateLineageTrees` scores the cloned trees; `fixNetwork`/`ALL_EDGES` build new networks with a fresh counter). The port's immediate stop is therefore result-equivalent.
3. Harness `F4445Harness.java` (= F4243 + `MAX_NUM_GROW_CALLS` from `-Dcap`, a `ties` count and optional `HashMap` group order): input c1 (ε 0.05, 7 clusters × 2 identical robust members, profile 11: [0.9,0.8], [0.5,0.45], [0.4,0.3], [0.3,0.35], [0.2,0.25], [0.1,0.15], [0.15,0.1]):

| cap | lichee.jar | parents of 1..7 (−1 = root) |
|-----|------------|------------------------------|
| 40 | 0 trees (default and `ALL_EDGES` searches both capped) | — |
| 80 | 7 trees, error 0.07071067811865477 | −1, 1, 2, 1, 4, 7, −1 |
| 320 | 17 trees, error 0.050000000000000044 | −1, 1, 2, 1, 4, 5, −1 |
| 10⁸ | 176 trees (30 at the minimum), error 0 | −1, 1, 2, 1, 4, −1, 5 |

4. Re-derived unit fixtures at the new default (jar `-e 0.1`): S3 [0.5,0.5]/[0.55,0] → 1 tree root→A→B, error 0.050000000000000044 (ε 0: none); 4 private 0.05…0.053 → 24 trees (15 at ε 0), same top chain; t02150 margins → default network, 1 tree, root→2, 2→{1,3}, error 0.016000000000000014 (static with the Java centroids: same); t00298 margins → complete network, 44 trees (28 at ε 0), same top tree; static t00298 → 6 trees (5), same top tree; f1 → nothing removed, root→1→2, error 0.050000000000000044; f2/f4/f4r → root→1→2, error 0.050000000000000044; f6 → root→1→{2,3}, f7 → root→1→{2,3}, 3→4, both error 0.07071067811865482; trunk fixtures unchanged ([1,1],[1,1],[0.4,0] → 2→1→3; [0.97],[0.5] → root→1→2). Fixtures whose purpose is ε = 0 (S3, F42 t02150/t00298, F43 f1–f7, the strict trunk case) now pass `tolerance: 0` explicitly — their ε = 0 jar values are unchanged.

### LICHeE equal-score tie order — opened 2026-10-10 (B24 F45)

**Source:** same commit: `SNVDataStore.loadSNVFileWithClusters` (`tag2SNVs = new HashMap<String, ArrayList<SNVEntry>>()`; per cluster line `if (tag2SNVs.containsKey(profile)) … else tag2SNVs.put(profile, …)`; no removal on this path), `LineageEngine.buildLineage` step 2 (`for (String groupTag : snvsByTag.keySet()) groups.add(…)`), `PHYNetwork` constructor (node ids assigned in group order; `nodes` `HashMap<Integer,…>`, `nodesById` `HashMap<Integer,…>`, `edges` `HashMap<PHYNode,…>` with `PHYNode.hashCode() = nodeId`), `PHYTree.compareTo` + `Collections.sort` (stable), OpenJDK 21 `java.util.HashMap` / `String.hashCode`.

1. **Which maps decide the order:** the only map whose iteration order reaches the result is `tag2SNVs` (key = presence-profile `String`, e.g. "101"). Node ids follow its `keySet()` order; ids fix the `edges` insertion order (hence the `grow` stack and the enumeration order), the argument order of intra-group `checkAndAddEdge` (orientation when both directions tie), the `computeErrorScore` summation order (sorted by id) and the `fixNetwork` scan (`nodesById.values()`, Integer keys = ascending id). The other maps have Integer keys or `PHYNode` keys with `hashCode = nodeId` (deterministic) or are only read order-insensitively (bridge test). The single identity-hash structure is `fixNetwork`'s `HashSet<SNVGroup>` (F43, post-removal only).
2. **Determinism:** `String.hashCode` is a fixed 31-polynomial ⇒ the order is the same in every JVM run (g45_s2: 2 000 inputs on a default JVM, no hash flags, identical to the port). Exact emulation is therefore possible; nothing is BLOCKED.
3. **Emulation:** `OncologyAnalyzer.JavaStringHashMapOrder.KeyOrder` ports `putVal` (tail append; note `computeIfAbsent` would insert at the bin head — the F42/F43 harness's `LinkedHashMap.computeIfAbsent` was insertion-ordered and is replaced), `resize` (16 → ×2 when size > 0.75·n, lo/hi split keeping order), `treeifyBin` (bin ≥ 9 nodes: resize if n < 64, else `TreeNode.treeify`), `putTreeVal` (new node linked after its tree parent), `balanceInsertion`/rotations, `moveRootToFront`, `split` (re-treeify > 6, `untreeify` ≤ 6). Check vs `java.util.HashMap.keySet()` (HmOrder.java): 8 000 / 8 000 random key sets identical (1–20-char 0/1 keys, up to 202 keys, 1 404 sets with treeified bins, 2 732 treeify events).
4. **Harness `F4445Harness -Dhm=1`** (`containsKey` + `put` like `SNVDataStore`, `removalStep` per `fixNetwork` call): port identical on r1 3 000 + r2 6 000 (F42 sets), m1 4 000 (F43 set, `-XX:hashCode=2`, removal order per step), g45 static 3 000 + member 3 000 (1–5 samples, 2–8 clusters, ε ∈ {0, 0.02, 0.05, 0.1, 0.2}) + g45_s2 2 000 = 21 000 / 21 000. Insertion-order grouping differs from the jar on 5 / 3 000 (static) and 85 / 3 000 (member-level) g45 inputs.

| Fixture (ε) | Clusters | Profiles: first appearance → HashMap order | lichee.jar (HashMap order) | Insertion order would give |
|-------------|----------|---------------------------------------------|----------------------------|-----------------------------|
| t02042 (0.1) | [0.25,0,0.8], [0.2,0.68,0], [0.2,0,0], [0.73,0,0.9] | 101,110,100 → 110,100,101 | 3 trees tied at 0; root→{2,4}, 4→1, 4→3 | 2→3 |
| t02235 (0.2) | [0,0.74,0], [0.4,0,0.6], [0.3,0.21,0], [0.8,0,0.8], [0.3,0,0] | 010,101,110,100 → 110,100,101,010 | 3 trees tied at 0.10000000000000009; root→{1,3,4}, 4→2, 4→5 | 3→5 |
| t00750 (0) | [0,0,0.23,0.3,0], [0,0,0.44,0.34,0.9], [0.1,0.3,0,0,0.1], [0,0,0,0,0.07], [0,0,0.3,0.6,0] | 00110,00111,11001,00001 → 11001,00110,00111,00001 | 2 trees tied at 0; root→{2,3,5}, 5→1, 3→4 | 2→4 |
| t01561 (0) | [0.46,0,0.6], [0,0,0.4], [0,0,0.3], [0,0.45,0.3] | 101,001,011 → 011,001,101 | 2 trees tied at 0; root→{1,4}, 1→2, 2→3 | 4→3 |
| t01396 (0.05) | [0,0,0.5,0,0.1], [0,0.19,0.48,0,0.19], [0.41,0.42,0,0.7,0], [0.4,0.42,0.04,0.67,0] | 00101,01101,11010,11110 → 11010,00101,11110,01101 | 1 tree, error 0.03741657386773935 | error 0.03741657386773934 |

5. **F43 fixture f6 re-derived:** with HashMap order (11, 01, 10) C = [0,0.55] is node 2 and is removed first: jar `removalStep [3]`, `removalStep [2]` → `RemovedClusterIds` = [3, 2] (was [2, 3] under the insertion-order harness; final tree unchanged). Other F42/F43/F44 fixtures: unchanged.

### Werner B et al. (2017), *Sci Rep* 7:44991 — trunk definition (WebSearch snippet)

"alterations that are in the trunk of the tree must be present in all cells of the tumour" ⇒ truncal ⇔ CCF = 1 in every sample.

---

## Documented Corner Cases and Failure Modes

### From Popic et al. (2015) — LICHeE

1. **Branch with insufficient parent budget:** when three candidate sibling groups have CCFs whose sum exceeds the parent CCF, "no more than one of these groups can be a descendant of the parent group, without violating the VAF phylogenetic constraint"; conflicting branches must be removed/re-placed.
2. **Under-determined private mutations:** private (single-sample) mutation groups are "under-constrained (i.e. multiple tree nodes can serve as ancestors)"; their exact placement is ambiguous, so a deterministic tie-break is required.
3. **Noise margin:** equality of the sum rule is relaxed to `≤ u + ϵ` to tolerate VAF/CCF measurement noise.

### From Zheng et al. (2022) — PICTograph

1. **Presence-pattern constraint:** a descendant cluster must be present only in (a subset of) the samples where the ancestor is present; this is the CCF analogue of LICHeE constraint (1).

---

## Test Datasets

### Dataset: Linear (chain) clonal evolution — derived from the lineage-precedence + sum rules

**Source:** Popic et al. (2015) Eq. 2 (ancestor ≥ descendant) and Eq. 5 (sum rule); single sample.

| Cluster | CCF (sample 1) | Expected placement (deepest valid ancestor) |
|---------|----------------|---------------------------------------------|
| Normal (root) | 1.0 | — |
| A | 1.0 | child of Normal (A.CCF ≤ Normal.CCF) |
| B | 0.6 | child of A (0.6 ≤ 1.0; A budget 1.0 ≥ 0.6) |
| C | 0.3 | child of B (0.3 ≤ 0.6; B budget 0.6 ≥ 0.3) |

Two trees are valid (C under A or under B; B can only be A's child because root: 1.0+0.6 > 1). LICHeE (lichee.jar) ranks root→A→{B, C} first (2 valid trees, error 0). Trunk = {A} (CCF 1). Branches = {B, C}.

### Dataset: Branching clonal evolution (two samples) — derived from constraints (1),(2),(3)

**Source:** Popic et al. (2015) constraints (1) presence pattern, (2) ancestor ≥ descendant per sample, (3) sum rule per sample.

| Cluster | CCF sample 1 | CCF sample 2 | Expected placement |
|---------|--------------|--------------|--------------------|
| Normal (root) | 1.0 | 1.0 | — |
| A (trunk) | 1.0 | 1.0 | child of Normal |
| B | 0.6 | 0.0 | child of A (B private to s1; not a descendant of C because C.s1=0 < B.s1=0.6) |
| C | 0.0 | 0.7 | child of A (C private to s2; not a descendant of B because B.s2=0 < C.s2=0.7) |

Edges: Normal→A, A→B, A→C. A is the trunk (parent of two sibling branches B, C). Sum rule under A holds per sample: s1 = 0.6+0.0 = 0.6 ≤ 1.0; s2 = 0.0+0.7 = 0.7 ≤ 1.0. B and C cannot be ancestor/descendant of each other (each is absent in the other's sample → constraint (1)/(2) violated both directions).

### Dataset: Sum-rule violation forces branching instead of nesting

**Source:** Popic et al. (2015) Eq. 5 — children CCF sum may not exceed parent CCF.

| Cluster | CCF (sample 1) | Expected placement |
|---------|----------------|--------------------|
| Normal | 1.0 | — |
| A | 1.0 | child of Normal |
| B | 0.6 | child of A (A budget 1.0 ≥ 0.6) |
| C | 0.6 | child of A — NOT child of B, because B already has nothing but B.CCF=0.6 ≥ C.CCF=0.6 would be allowed by lineage rule; however attaching C under A keeps A's children sum 0.6+0.6=1.2 > 1.0 → violates sum rule |

Expected: with B and C both 0.6 they cannot both be children of the same parent (sum 1.2 > 1.0). LICHeE orients the equal pair C→B (`checkAndAddEdge` tie → later→earlier), so the only valid tree is Normal→A→C→B (lichee.jar: 1 valid tree). This shows the sum rule converting a would-be sibling pair into a chain.

### Dataset: No valid tree / FP-exact sum rule / complete network (lichee.jar outputs)

| ε | Clusters | LICHeE result |
|---|----------|---------------|
| 0 | A=[0.5,0.5], B=[0.55,0] | none (A→B violates Eq. 2; root 0.5+0.55 > 1) |
| 0.1 | same | A→B, error 0.050000000000000044 |
| 0 | A=[0.3,1.0], B=[0.2,0.2], C=[0,0.8] | A→{B,C}, 1 tree, error 0 (former greedy: root→B, Eq. 5 violated) |
| 0.02 | A=[0.482,0,0.443137], B=[0.519,0.796779,0.56], C=[0.19,0.418589,0.24] | complete network; root→{A,B}, B→C; 1 tree; error 0.003292532308117998 |
| 0 | 0.05, 0.051, 0.052, 0.053 (1 sample) | 15 trees; top = chain root→0.053→0.052→0.051→0.05 |

---

## Assumptions

1. **Tie-break among equal-score trees** — LICHeE returns the first tree of its Gabow–Myers enumeration among those with the minimal error score; the enumeration order depends on node order, which LICHeE takes from a Java `HashMap<String>` of presence profiles. Since B24 F45 the port reproduces that `HashMap` iteration order exactly (`JavaStringHashMapOrder`; within a profile: input order) — no longer an assumption. (Superseded 2026-09-28: the former "deepest valid ancestor" greedy, which could return sum-rule-violating trees — B24 F18.)
2. ~~**ASSUMPTION: Noise margin ε = 0**~~ — superseded 2026-10-10 (B24 F44): the default is LICHeE's `-e` default 0.1 (also in `-cp` mode); `tolerance: 0` gives the strict inequalities.

---

## Recommendations for Test Coverage

1. **MUST Test:** Linear chain from descending single-sample CCFs (Normal→A→B→C). — Evidence: Popic 2015 Eq. 2; Zheng 2022 lineage precedence.
2. **MUST Test:** Branching tree from two private single-sample clusters (A→B, A→C siblings). — Evidence: Popic 2015 constraints (1)(2)(3).
3. **MUST Test:** Sum rule rejects two equal-CCF siblings under one parent and forces a chain. — Evidence: Popic 2015 Eq. 5.
4. **MUST Test:** Trunk identification = clusters on the root path with CCF = 1 (≥ 1 − ε) in all samples; branch identification = all non-trunk clusters. — Evidence: Werner 2017; Popic 2015.
5. **MUST Test:** Ancestor CCF ≥ descendant CCF holds on every reconstructed edge (invariant). — Evidence: Popic 2015 Eq. 2.
6. **MUST Test:** Per-node sum rule holds on the reconstructed tree (invariant). — Evidence: Popic 2015 Eq. 5; Zheng 2022 sum condition.
7. **SHOULD Test:** Empty input → tree with only the root, no trunk/branch mutations. — Rationale: boundary.
8. **SHOULD Test:** Single cluster → child of root, is the trunk. — Rationale: boundary.
9. **SHOULD Test:** Null input / CCF out of [0,1] / NaN / inconsistent sample count → exceptions. — Rationale: documented validation.
10. **COULD Test:** Tolerance ε > 0 admits a near-violating edge that ε = 0 would reject. — Rationale: optional parameter.

---

## References

1. Popic V, Salari R, Hajirasouliha I, Kashef-Haghighi D, West RB, Batzoglou S. (2015). Fast and scalable inference of multi-sample cancer lineages. *Genome Biology* 16:91. https://doi.org/10.1186/s13059-015-0647-8 (full-text PDF: https://arxiv.org/pdf/1412.8574)
2. Zheng L, Dang H, Niknafs N, et al. (2022). Estimation of cancer cell fractions and clone trees from multi-region sequencing of tumors (PICTograph). *Bioinformatics* 38(15):3677–3683. https://doi.org/10.1093/bioinformatics/btac367 (open access: https://pmc.ncbi.nlm.nih.gov/articles/PMC9344857/)

---

## Change History

- **2026-06-15**: Initial documentation.
- **2026-09-28**: B24 F18/F19 — LICHeE reference code + jar cross-check; datasets corrected to LICHeE output; no-valid-tree case; CCF-based trunk (Werner 2017).
- **2026-10-10**: B24 F42 — LICHeE per-cluster `1.96·sd/√n` edge margins (`getAAFErrorMargin`), jar-locked.
- **2026-10-10**: B24 F43 — LICHeE `fixNetwork` (non-robust cluster removal, then `ALL_EDGES` on the reduced set), jar-locked.
- **2026-10-10**: B24 F44 — default ε = LICHeE 0.1 (`-cp` keeps `VAF_ERROR_MARGIN`); post-10⁸-cap equivalence proof + reduced-cap jar locks; fixtures re-derived at `-e 0.1`.
- **2026-10-10**: B24 F45 — group order = Java `HashMap<String>` iteration of the profile tags (exact emulation incl. treeified bins); 21 000 / 21 000 jar inputs, 8 000 / 8 000 key sets; f6 removal order [3, 2].
