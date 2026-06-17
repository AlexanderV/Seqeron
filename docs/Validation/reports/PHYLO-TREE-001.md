# Validation Report: PHYLO-TREE-001 — Phylogenetic Tree Construction (UPGMA + Neighbor-Joining)

- **Validated:** 2026-06-17   **Area:** Phylogenetics
- **Re-validation context:** Phase-3 ledger row #9 — N-ary (multifurcating) tree refactor (commit c4f0190).
  `PhyloNode` now holds `List<PhyloNode> Children`; **Neighbor-Joining now emits its natural
  unrooted trifurcation** at the final central node (the classic Saitou-Nei result), enabling a true
  NJ round-trip through Newick; UPGMA must remain strictly bifurcating.
- **Canonical method(s):** `PhylogeneticAnalyzer.BuildTree` / `BuildTreeFromMatrix` →
  `BuildUPGMA` (`PhylogeneticAnalyzer.cs:306`), `BuildNeighborJoining` (`:427`),
  `ToNewick` (`:566`), `ParseNewick` (`:632`).
- **Stage A verdict:** ✅ PASS
- **Stage B verdict:** ✅ PASS
- **End-state:** ✅ CLEAN (no code defect; 3 strict sourced tests added to close a test-quality gap)

---

## Stage A — Description

### Sources opened

- **Saitou & Nei 1987**, "The Neighbor-joining Method: A New Method for Reconstructing Phylogenetic
  Trees", *Mol. Biol. Evol.* 4(4):406–425. The NJ algorithm: at each step the pair minimizing the
  total tree length is joined; branch lengths from the additive system. Final step with three nodes
  remaining → a single central node with three branches (unrooted-tree centre), with
  `L_iu = (d_ij + d_ik − d_jk)/2` and symmetric forms for j, k.
  (https://doi.org/10.1093/oxfordjournals.molbev.a040454)
- **Wikipedia "Neighbor joining"** — Q-matrix `Q(i,j) = (n−2)·d(i,j) − Σ_k d(i,k) − Σ_k d(j,k)`;
  branch length of joined member f to new node u:
  `δ(f,u) = ½ d(f,g) + 1/(2(n−2)) · [Σ_k d(f,k) − Σ_k d(g,k)]`, the other `δ(g,u) = d(f,g) − δ(f,u)`;
  distance update `d(u,k) = ½ [d(f,k) + d(g,k) − d(f,g)]`. Worked **5-taxon example** with a known
  answer used as the cross-check. (https://en.wikipedia.org/wiki/Neighbor_joining)
- **Sokal & Michener 1958** / **Wikipedia "UPGMA"** — proportional (size-weighted) averaging
  `d((A∪B),X) = (|A|·d(A,X) + |B|·d(B,X)) / (|A|+|B|)`; new cluster height `= d(i,j)/2`; incremental
  branch length `= height(new) − height(child)`; strictly **binary, ultrametric**. Worked **5S-rRNA
  5-bacteria example** used as the cross-check. (https://en.wikipedia.org/wiki/UPGMA)

### Formula check

All formulas in the code match the sources exactly: Q-matrix (`:470`), NJ branch lengths (`:482–483`),
NJ distance update (`:503`), the final-step trifurcation system (`:526–528`), UPGMA proportional
averaging (`:397`), UPGMA height/incremental-branch (`:373–377`).

### Independent cross-check (hand-recomputed, NOT read off the code)

**NJ — Wikipedia 5-taxon additive matrix** (a..e):

```
   a  b  c  d  e
a  0  5  9  9  8
b  5  0 10 10  9
c  9 10  0  8  7
d  9 10  8  0  3
e  8  9  7  3  0
```

Hand-run (Python, formulas only):

| step | n | join | Q_min | branch lengths |
|------|---|------|-------|----------------|
| 1 | 5 | a+b → u1 | −50 | δ(a)=2, δ(b)=3 |
| 2 | 4 | u1+c → u2 | −28 | δ(c)=4, δ(u1→u2)=3 |
| final | 3 | centre {u2, d, e} | — | δ(u2)=2, δ(d)=2, δ(e)=1 |

Final-step trifurcation from `δ_i=(d_ij+d_ik−d_jk)/2` on `d(d,e)=3, d(d,u2)=4, d(e,u2)=3`
→ δ(d)=2, δ(e)=1, δ(u2)=2. The centre has **exactly 3 children**. Manual patristic reconstruction
on this tree reproduces **all 10 input distances exactly** (additive matrix → NJ recovers the true
tree, Saitou-Nei guarantee).

**UPGMA — Wikipedia 5S-rRNA matrix** (a..e):

```
   a  b  c  d  e
a  0 17 21 31 23
b 17  0 30 34 21
c 21 30  0 28 39
d 31 34 28  0 43
e 23 21 39 43  0
```

Hand-run: merge (a,b)@17 (h=8.5) → e@22 (h=11) → (c,d)@28 (h=14) → final@33 (h=16.5).
Branch lengths δ(a)=δ(b)=8.5, δ(u,v)=2.5, δ(e)=11, δ(c)=δ(d)=14, δ(v,r)=5.5, δ(w,r)=2.5.
Ultrametric root-to-tip = 16.5 for all 5 tips. Strictly binary (4 internal binary nodes for 5 leaves).
These match the published Wikipedia dendrogram exactly.

### Findings / divergences

None. Description is mathematically correct and faithful to the primary sources.

---

## Stage B — Implementation

### Code path reviewed

`BuildNeighborJoining` (`:427–561`) and `BuildUPGMA` (`:306–413`); `ToNewick`/`ToNewickRecursive`
(`:566–609`); `ParseNewick`/`ParseNewickRecursive` (`:632–726`).

### Formula realised correctly? (evidence — probe against the real compiled code)

A standalone probe (`/tmp/probe`) referencing the production assembly produced, for the NJ matrix:

```
Root children count: 3
  child ((a,b),c)  branch=2   (internal clade, leaves {a,b,c})
  child d          branch=2   (leaf)
  child e          branch=1   (leaf)
Newick: (((a:2.0000,b:3.0000):3.0000,c:4.0000):2.0000,d:2.0000,e:1.0000);
RoundTrip equal: True
patristic(a,b)=5 … (d,e)=3  → all 10 == input
```

For the UPGMA matrix:

```
Root children: 2
Newick: (((a:8.5000,b:8.5000):2.5000,e:11.0000):5.5000,(c:14.0000,d:14.0000):2.5000);
```

Both reproduce my **hand-derived** numbers exactly: NJ centre is a genuine 3-child trifurcation with
δ=(2,2,1); inner branches a=2,b=3,u1=3,c=4; UPGMA is strictly binary with the correct ultrametric
heights. NJ round-trips through Newick **byte-for-byte** with the trifurcation preserved on re-parse.

### Test quality audit (HARD gate applied)

Pre-existing tests already had a real trifurcation guard
(`BuildTree_NJ_FinalNodeIsTrifurcation_RoundTripsThroughNewick`, asserts `Children.Count==3` + round
trip) and a real UPGMA-binary guard (`…EveryInternalNodeHasExactlyTwoChildren`). **Gap found:** no test
locked the **exact trifurcation branch lengths** (δ(d)=2, δ(e)=1, δ(internal)=2) nor the **exact
Newick strings** for either method — these exact values are the headline of the refactor and a wrong
trifurcation formula or a binarised serialisation could pass the looser checks. Closed with **3 strict
sourced tests** (all expected values hand-derived in Stage A, not captured from code):

- `BuildTree_NJ_WikipediaExample_ExactTrifurcationBranchLengths` — `Children.Count==3`; pins
  δ(d)=2, δ(e)=1, δ(((a,b),c))=2; pins the internal child's leaf set = {a,b,c}.
- `BuildTree_NJ_WikipediaExample_ExactNewickRoundTrip` — pins the exact string
  `(((a:2.0000,b:3.0000):3.0000,c:4.0000):2.0000,d:2.0000,e:1.0000);` (top level = 3 members) and a
  byte-for-byte build→ToNewick→ParseNewick→ToNewick round trip with 3 re-parsed children.
- `BuildTree_UPGMA_WikipediaExample_ExactNewickAndStrictlyBinary` — pins the exact ultrametric string
  `(((a:8.5000,b:8.5000):2.5000,e:11.0000):5.5000,(c:14.0000,d:14.0000):2.5000);` and every internal
  node `Children.Count==2` (4 internal nodes).

**Mutation verification (anti-green-washing):**
- NJ `while (active.Count > 3)` → `> 2` (collapse trifurcation to a bifurcating join): kills the 2 new
  NJ tests **plus** the pre-existing trifurcation test (3 fail). A bifurcating NJ cannot reproduce the
  3-member top-level Newick — exactly the defect the gate targets.
- UPGMA proportional averaging → naive `(dIK+dJK)/2`: kills the new UPGMA exact-Newick test plus the
  two pre-existing branch-length/ultrametric tests (3 fail).

Both mutations reverted; production source is **byte-for-byte unchanged** (no code defect).

### Findings / defects

No code defect. One test-quality gap (exact trifurcation branch lengths + exact Newick strings not
locked) → fixed with 3 strict sourced tests. Tie-break note: at NJ step 2 the code deterministically
selects (u1,c) (Q=−28). Wikipedia notes a tie Q(u1,c)=Q(d,e)=−28; both are valid NJ topologies for an
additive matrix and both reproduce the input distances (additive guarantee). Tests assert the
deterministic produced topology and exact values, consistent with the protocol's "code obeys tests
derived from the algorithm" rule.

---

## Verdict & follow-ups

- **Stage A:** ✅ PASS — formulas and the two worked examples verified against Saitou-Nei 1987 and
  Sokal-Michener/UPGMA; hand-computation matches the published trees.
- **Stage B:** ✅ PASS — code reproduces the hand-derived NJ trifurcation (δ=2,2,1), the binary UPGMA
  tree, and the byte-for-byte Newick round trip; strict tests added; mutation-confirmed real.
- **End-state:** ✅ CLEAN. Full unfiltered suite **6779 passed / 0 failed**; build 0 errors.
