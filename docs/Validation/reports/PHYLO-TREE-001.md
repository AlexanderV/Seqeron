# Validation Report: PHYLO-TREE-001 — Phylogenetic Tree Construction (UPGMA & Neighbor-Joining)

- **Validated:** 2026-06-12   **Area:** Phylogenetics
- **Canonical method(s):** `PhylogeneticAnalyzer.BuildTree(IReadOnlyDictionary<string,string>, DistanceMethod, TreeMethod)`, with `BuildTreeFromMatrix(...)`, internal `BuildUPGMA(...)`, `BuildNeighborJoining(...)`
- **Source file:** `src/Seqeron/Algorithms/Seqeron.Genomics.Phylogenetics/PhylogeneticAnalyzer.cs`
- **Test file:** `tests/Seqeron/Seqeron.Genomics.Tests/PhylogeneticAnalyzer_TreeConstruction_Tests.cs`
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS
- **End state:** CLEAN

---

## Stage A — Description

### Sources opened
- Wikipedia "Neighbor joining" — Q-matrix, branch-length, and distance-update formulas; full 5-taxon worked example.
- Wikipedia "UPGMA" — average-linkage (size-weighted) cluster update, height = d/2, 5-taxon (5S rRNA) worked example.
- Cross-references: Saitou & Nei (1987) MBE 4(4):406-425 (NJ original), Sokal & Michener (1958) (UPGMA original) — cited in evidence doc, consistent with Wikipedia statements.

### Formula check (confirmed verbatim against Wikipedia)

**UPGMA**
- Join the two clusters with the SMALLEST distance. ✓
- New height = d(i,j)/2 (ultrametric). ✓
- Average-linkage update weighted by cluster sizes:
  d(A∪B, X) = (|A|·d(A,X) + |B|·d(B,X)) / (|A|+|B|). ✓
- Produces a ROOTED ultrametric tree. ✓

**Neighbor-Joining**
- Q(i,j) = (n−2)·d(i,j) − Σ_k d(i,k) − Σ_k d(j,k); join MINIMUM Q. ✓ (the (n−2) factor is present.)
- δ(f,u) = ½·d(f,g) + 1/(2(n−2))·(Σd(f,k) − Σd(g,k)); δ(g,u) = d(f,g) − δ(f,u). ✓
- Update: d(u,k) = ½(d(f,k) + d(g,k) − d(f,g)). ✓
- Produces an UNROOTED additive tree (rooted by convention; here midpoint d/2 on the final edge, which preserves patristic distances). ✓

### Independent hand computation

**NJ (matrix a-b=5, a-c=9, a-d=9, a-e=8, b-c=10, b-d=10, b-e=9, c-d=8, c-e=7, d-e=3; n=5)**
Row sums: r(a)=31, r(b)=34, r(c)=34, r(d)=30, r(e)=27.
Q(a,b) = 3·5 − 31 − 34 = −50 (unique minimum); Q(d,e) = 3·3 − 30 − 27 = −48.
→ a,b joined first.
δ(a,u) = ½·5 + (31−34)/(2·3) = 2.5 − 0.5 = **2.0**; δ(b,u) = 5 − 2 = **3.0**.
Subsequent joins reproduce the additive matrix; all 10 patristic distances equal the input. **Matches Wikipedia exactly.**

**UPGMA (matrix a-b=17, a-c=21, a-d=31, a-e=23, b-c=30, b-d=34, b-e=21, c-d=28, c-e=39, d-e=43)**
1. min = d(a,b)=17 → u at height 8.5; δ(a,u)=δ(b,u)=**8.5**.
2. d((ab),e)=(23+21)/2=22, d((ab),c)=(21+30)/2=25.5, d((ab),d)=(31+34)/2=32.5; min = 22 → v at 11; δ(e,v)=**11**, δ(u,v)=11−8.5=**2.5**.
3. d(c,d)=28 → w at 14; δ(c,w)=δ(d,w)=**14**.
4. Final join at 33 → root r at 16.5; δ(v,r)=16.5−11=**5.5**, δ(w,r)=16.5−14=**2.5**.
All tips equidistant from root: **16.5**. **Matches Wikipedia exactly.**

### Edge-case semantics (sourced)
2 taxa → trivial binary tree; <2 / null / empty → ArgumentException; unequal lengths → ArgumentException; identical sequences → distance 0; ties (NJ step 2 here is Q=−28 tied between (u,c) and (d,e)) → either choice is a valid NJ result per Wikipedia. All defined and standard.

**Stage A verdict: PASS** — every formula and worked-example value matches the authoritative sources; no defects in the description.

---

## Stage B — Implementation

### Code path reviewed
`PhylogeneticAnalyzer.cs`: `BuildUPGMA` (L248-355), `BuildNeighborJoining` (L360-466), distance functions (L141-240).

### Formula realised correctly?
- **UPGMA height** `newHeight = minDist / 2` (L315). ✓
- **UPGMA incremental branch** `Math.Max(0, newHeight − clusterHeights[child])` (L318-319) — ultrametric, non-negative. ✓
- **UPGMA size-weighted update** `(dIK·clusterSizes[minI] + dJK·clusterSizes[minJ]) / newSize` (L339), using full cluster sizes (sizes updated only after, L349). This is genuine average-linkage UPGMA, NOT simple averaging. ✓
- **NJ Q-matrix** `(m-2)*dist[i,j] - r[i] - r[j]` (L401) — (n−2) factor present, minimum tracked. ✓
- **NJ branch length** `(distIJ/2) + (r[minI]-r[minJ])/(2*(m-2))` (L413), other side `distIJ - branchI` (L414). ✓
- **NJ distance update** `(dist[minI,k] + dist[minJ,k] - dist[minI,minJ]) / 2` (L434). ✓
- **NJ unrooted/negative branches**: negative lengths preserved (no clamp, L426-427); final two-node join split d/2 (L449), midpoint rooting preserving patristic distances. ✓

### Cross-verification table recomputed vs code (via tests)

| Quantity | External value | Code (test) |
|----------|---------------|-------------|
| UPGMA δ(a,u),δ(b,u) | 8.5 | 8.5 (S01b) ✓ |
| UPGMA δ(e,v) | 11 | 11 ✓ |
| UPGMA δ(u,v) | 2.5 | 2.5 ✓ |
| UPGMA δ(c,w),δ(d,w) | 14 | 14 ✓ |
| UPGMA δ(v,r),δ(w,r) | 5.5 / 2.5 | 5.5 / 2.5 ✓ |
| UPGMA tip→root | 16.5 (all) | 16.5 (S01c) ✓ |
| NJ first join | (a,b), Q=−50 | (a,b) (S02b) ✓ |
| NJ δ(a,u),δ(b,u) | 2 / 3 | 2 / 3 (S02c) ✓ |
| NJ patristic = input | all 10 match | match (S02) ✓ |

### Variant/delegate consistency
`BuildTree` and `BuildTreeFromMatrix` share the same `BuildUPGMA`/`BuildNeighborJoining` core. Distance methods (Hamming/PDistance/JC/K2P) all produce valid trees (C03).

### Numerical robustness
Worked-example arithmetic is integer/half-integer, exact in IEEE-754; tests assert at 1e-10. JC saturation returns +∞ (defined). No div-by-zero on stated ranges (n≥2; NJ loop only runs while active>2 so (m−2)≥1).

### Test quality audit
32 test runs assert exact externally-sourced values (not tautologies), deterministic, covering MUST/SHOULD/COULD plus edge cases (2/3 taxa, identical seqs, gap-only columns, validation throws, NJ tie handling, ultrametric, additive patristic). Strong.

**Stage B verdict: PASS** — code faithfully realises every validated formula; all cross-check values reproduce the authoritative worked examples.

---

## Verdict & follow-ups

- **STATE: CLEAN.** No defects found. Common scientific traps checked and absent: UPGMA uses size-weighted (not simple) average linkage; NJ Q-matrix includes the (n−2) factor; NJ branch-length and distance-update formulas exact; NJ is treated as unrooted with patristic-preserving midpoint root and negatives unclamped.
- **Tests:** PHYLO-TREE-001 filter (`~TreeConstruction`) = 32 passed, 0 failed. Full `Seqeron.Genomics.Tests` = 4461 passed, 0 failed (matches baseline).
- **Code changed:** none.
