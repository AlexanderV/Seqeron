# Validation Report: PHYLO-TREE-001 — Phylogenetic Tree Construction (UPGMA + Neighbor-Joining)

- **Validated:** 2026-06-24 (independent re-confirmation)   **Area:** Phylogenetics
- **Re-validation context:** Re-confirm the N-ary (multifurcating) refactor (commit `c4f0190`):
  `BuildUPGMA` (binary ultrametric), `BuildNeighborJoining` (final trifurcation),
  `PhyloNode.Children`. Done in a fresh context, formulas re-derived from primary sources,
  both Wikipedia worked examples re-run by hand, code probed via the test suite.
- **Canonical method(s):** `PhylogeneticAnalyzer.BuildTree` / `BuildTreeFromMatrix` →
  `BuildUPGMA` (`PhylogeneticAnalyzer.cs:306`), `BuildNeighborJoining` (`:427`),
  `ToNewick` (`:577`), `ParseNewick` (`:657`).
- **Stage A verdict:** ✅ PASS
- **Stage B verdict:** ✅ PASS
- **End-state:** ✅ CLEAN (no code defect; source byte-for-byte unchanged since `c4f0190`)

---

## Stage A — Description

### Sources opened & what they confirm

- **Saitou & Nei 1987**, "The Neighbor-joining Method…", *Mol. Biol. Evol.* 4(4):406–425 — NJ
  joins, at each step, the pair minimizing total tree length; branch lengths from the additive
  system; final 3 nodes → one central node with three branches (the unrooted-tree centre):
  `L_iu = (d_ij + d_ik − d_jk)/2` and symmetric forms. (DOI: 10.1093/oxfordjournals.molbev.a040454)
- **Studier & Keppler 1988** — O(n³) reformulation of NJ via the Q-matrix; same join criterion.
- **Wikipedia "Neighbor joining"** — `Q(i,j) = (n−2)·d(i,j) − Σ_k d(i,k) − Σ_k d(j,k)`;
  `δ(f,u) = ½ d(f,g) + 1/(2(n−2))·[Σ_k d(f,k) − Σ_k d(g,k)]`, `δ(g,u) = d(f,g) − δ(f,u)`;
  update `d(u,k) = ½[d(f,k) + d(g,k) − d(f,g)]`. 5-taxon additive worked example as cross-check.
- **Sokal & Michener 1958 / Wikipedia "UPGMA"** — proportional (size-weighted) averaging
  `d((A∪B),X) = (|A|·d(A,X)+|B|·d(B,X))/(|A|+|B|)`; new-cluster height `= d(i,j)/2`; incremental
  branch length `= height(new) − height(child)`; strictly binary, ultrametric. 5S-rRNA 5-bacteria
  example as cross-check.

### Formula check (code ↔ source)

| Quantity | Source | Code |
|----------|--------|------|
| NJ Q-matrix | `(n−2)d(i,j) − r_i − r_j` | `:470` `(m-2)*dist[i,j] - r[i] - r[j]` ✓ |
| NJ branch length δ(f,u) | `½d(f,g) + (r_f − r_g)/(2(n−2))` | `:493` `distIJ/2 + (r[i]-r[j])/(2*(m-2))` ✓ (r = full row sums, equivalent to Σ form) |
| NJ distance update | `½[d(f,k)+d(g,k)−d(f,g)]` | `:514` ✓ |
| NJ final trifurcation | `δ_i=(d_ij+d_ik−d_jk)/2` (+sym) | `:537–539` ✓; stop at `active.Count > 3` (`:445`), centre = 3 `Children` (`:544`) |
| UPGMA proportional avg | `(|A|d_AX+|B|d_BX)/(|A|+|B|)` | `:397` `(dIK*size_i+dJK*size_j)/newSize` ✓ |
| UPGMA height / branch | `d(i,j)/2`; `height(new)−height(child)` | `:373`, `:376–377` ✓ |

NJ negative branch lengths are preserved (not clamped) per Saitou-Nei; UPGMA branch lengths
are `Math.Max(0, …)` which is correct since UPGMA heights are monotone non-decreasing.

### Edge-case semantics

`<2` sequences / null / unequal lengths → `ArgumentException` (`:141`, `:149`). Identical
sequences → zero distances. Gap/ambiguous-only columns → `comparableSites==0 → distance 0`
(`:256`). 2-taxon degenerate NJ → single edge split `d/2` (`:552`). All defined & sourced.

### Independent cross-check (hand-recomputed in Python from formulas only, NOT read off code)

**NJ — Wikipedia 5-taxon additive matrix** (a..e: ab5 ac9 ad9 ae8 bc10 bd10 be9 cd8 ce7 de3):

| step | n | join | Q_min | branch lengths |
|------|---|------|-------|----------------|
| 1 | 5 | a+b → u1 | −50 | δ(a)=2, δ(b)=3 |
| 2 | 4 | u1+c → u2 | −28 | δ(c)=4, δ(u1)=3 |
| final | 3 | centre {u2,d,e} | — | δ(u2)=2, δ(d)=2, δ(e)=1 |

Centre is exactly 3 children. Patristic reconstruction recovers all 10 input distances exactly
(additive matrix → NJ recovers the true tree).

**UPGMA — Wikipedia 5S-rRNA matrix** (a..e):

merge (a,b)@17 (h=8.5) → e+(a,b)@22 (h=11) → (c,d)@28 (h=14) → final@33 (h=16.5).
Branches δ(a)=δ(b)=8.5, δ((a,b)→u)=2.5, δ(e)=11, δ(c)=δ(d)=14, δ(v,root)=5.5, δ(w,root)=2.5.
Ultrametric: all 5 tips at 16.5 from root; strictly binary (4 internal binary nodes).
Matches the published Wikipedia dendrogram.

### Findings / divergences

None. Description is mathematically correct and faithful to the primary sources.

---

## Stage B — Implementation

### Code path reviewed

`BuildUPGMA` (`:306–413`), `BuildNeighborJoining` (`:427–572`), `ToNewick` (`:577–620`),
`ParseNewick` (`:657–760`). **Source is byte-for-byte identical to the prior validation**
(`git log` shows no change to `PhylogeneticAnalyzer.cs` since `c4f0190`).

### Formula realised correctly?

The actual produced values were confirmed via the test suite (which pins hand-derived values):
- NJ root has `Children.Count==3`; Newick = `(((a:2.0000,b:3.0000):3.0000,c:4.0000):2.0000,d:2.0000,e:1.0000);`
  — top-level descendant list has THREE members; δ(d)=2, δ(e)=1, δ(((a,b),c))=2; round-trips
  byte-for-byte through ParseNewick preserving the trifurcation.
- UPGMA strictly binary (4 internal nodes, all `Children.Count==2`); Newick =
  `(((a:8.5000,b:8.5000):2.5000,e:11.0000):5.5000,(c:14.0000,d:14.0000):2.5000);`.

Both reproduce my Stage-A hand-derived numbers exactly.

### Cross-verification table recomputed vs code

| Example | My hand value | Test assertion (code-produced) | Match |
|---------|---------------|-------------------------------|-------|
| NJ join order / Q | a+b(−50), u1+c(−28) | S02b/S02d | ✓ |
| NJ δ(a),δ(b) | 2, 3 | S02c | ✓ |
| NJ trifurcation δ(u2),δ(d),δ(e) | 2, 2, 1 | S02d ExactTrifurcation | ✓ |
| NJ patristic = input | all 10 | S02 additivity | ✓ |
| UPGMA heights | 8.5,11,14,16.5 | S01b/S01c | ✓ |
| UPGMA Newick | as above | S01d ExactNewick | ✓ |

### Variant/delegate consistency

`BuildTree` (from sequences) and `BuildTreeFromMatrix` (direct matrix) route to the same
`BuildUPGMA`/`BuildNeighborJoining`. `Left`/`Right` are convenience accessors over `Children[0/1]`;
all traversals (`GetLeaves`, `PatristicDistance`, `FindMRCA`, clades, bipartitions) use `Children`,
so multifurcations are handled uniformly.

### Test quality audit

The exact-value tests (S01d, S02d, S02e) assert hand-derived Newick strings and branch lengths
at 1e-10 tolerance, with explicit comments stating values were derived from the algorithm, NOT
captured from code. `Children.Count` guards on both the NJ trifurcation (==3) and UPGMA binarity
(all ==2) lock the headline refactor behaviour. The prior session's mutation checks
(`>3`→`>2` collapse; proportional→naive averaging) confirmed these tests are real (each kills ≥3 tests).

### Findings / defects

None. Tie-break note (unchanged from prior): at NJ step 2 the code deterministically selects
(u1,c) at Q=−28; Wikipedia notes a tie Q(u1,c)=Q(d,e)=−28. Both are valid NJ topologies on an
additive matrix and both reproduce the input distances; tests assert the deterministic produced
topology.

---

## Verdict & follow-ups

- **Stage A:** ✅ PASS — formulas and both worked examples independently re-verified against
  Saitou-Nei 1987 / Studier-Keppler 1988 / Sokal-Michener; hand-computation matches published trees.
- **Stage B:** ✅ PASS — code (unchanged) reproduces the NJ trifurcation (δ=2,2,1), the binary
  UPGMA ultrametric tree, and the byte-for-byte Newick round trip; exact sourced tests pass.
- **End-state:** ✅ CLEAN. No code change. Phylogenetic tree tests: **63 passed / 0 failed**.
