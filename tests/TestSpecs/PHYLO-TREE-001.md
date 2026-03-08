# PHYLO-TREE-001: Tree Construction

**Test Unit ID:** PHYLO-TREE-001
**Category:** Phylogenetics
**Created:** 2026-02-01
**Status:** ☑ Complete
**Evidence:** [PHYLO-TREE-001-Evidence.md](../docs/Evidence/PHYLO-TREE-001-Evidence.md)

---

## 1. Scope

### Canonical Method
- `PhylogeneticAnalyzer.BuildTree(IReadOnlyDictionary<string, string> sequences, DistanceMethod, TreeMethod)`

### Supported Tree Methods
- UPGMA (Unweighted Pair Group Method with Arithmetic Mean)
- Neighbor-Joining

### Returns
- `PhylogeneticTree` record with: Root, Taxa, DistanceMatrix, Method

---

## 2. Test Classification

### MUST Tests (Core Functionality)

| ID | Test | Rationale | Source |
|----|------|-----------|--------|
| M01 | BuildTree returns valid PhylogeneticTree | Basic contract | Definition |
| M02 | Tree contains all input taxa as leaves | INV-01 | Wikipedia |
| M03 | UPGMA method produces tree with Method="UPGMA" | Contract | Implementation |
| M04 | NeighborJoining method produces tree with Method="NeighborJoining" | Contract | Implementation |
| M05 | Two sequences create binary tree with both as leaves | Minimum case | Algorithm |
| M06 | Three sequences create valid binary tree structure | Basic case | Algorithm |
| M07 | Four sequences: structure, MRCA of identical pair, distance=0 | Standard case + structural | Wikipedia example |
| M08 | Throws on single sequence | Validation | Definition |
| M09 | Throws on unequal sequence lengths | Alignment required | Definition |
| M10 | All branch lengths are non-negative | INV-05 | UPGMA property |
| M11 | UPGMA produces rooted tree | INV-U02 | Wikipedia UPGMA |
| M12 | Identical sequences produce zero-distance subtree | Algorithm | UPGMA/NJ |
| M13 | Case-insensitive sequence handling | Robustness | Implementation |

### SHOULD Tests (Enhanced Coverage)

| ID | Test | Rationale | Source |
|----|------|-----------|--------|
| S01 | Wikipedia UPGMA example: correct clustering order | Validation | Wikipedia UPGMA working example |
| S01b | Wikipedia UPGMA example: leaf + internal branch lengths (δ(a,u)=8.5, δ(u,v)=2.5, δ(v,root)=5.5, δ(c,w)=14, δ(w,root)=2.5, δ(e,v)=11) | Validation | Wikipedia UPGMA working example |
| S01c | Wikipedia UPGMA example: ultrametric (all tips at 16.5 from root) | INV-U01 | Wikipedia UPGMA dendrogram |
| S02 | Wikipedia NJ example: patristic distances match input matrix (INV-N01 topology guarantee) | INV-N01 | Wikipedia NJ conclusion |
| S02b | Wikipedia NJ example: a,b joined first (Q=-50 unique minimum) | Validation | Wikipedia NJ first joining |
| S02c | Wikipedia NJ example: δ(a,u)=2, δ(b,u)=3 | Validation | Wikipedia NJ first branch estimation |
| S03 | UPGMA ultrametric property on general input | INV-U01 | Wikipedia UPGMA |
| S04 | DistanceMatrix in result is symmetric | Consistency | Matrix definition |
| S05 | Tree leaves match exact input taxon names (EquivalentTo, no extras) | Traceability | Definition |
| ~~S06~~ | ~~NJ additive matrix topology guarantee~~ | Removed | Duplicate of S02 |

### COULD Tests (Extended Coverage)

| ID | Test | Rationale | Source |
|----|------|-----------|--------|
| C01 | Large input (50+ sequences) completes in reasonable time | Performance | O(n³) complexity |
| C02 | Gap-only columns are handled correctly | Edge case | Implementation |
| C03 | Different distance methods produce valid trees | Flexibility | Implementation |
| C04 | Tree total length is sum of all branch lengths | Consistency | Definition |
| C05 | Both tree methods produce trees with same leaf set | Consistency | Definition |

---

## 3. Edge Cases

| Case | Input | Expected | Source |
|------|-------|----------|--------|
| Minimum input | 2 sequences | Binary tree | Definition |
| Single sequence | 1 sequence | ArgumentException | Validation |
| Unequal lengths | Mixed lengths | ArgumentException | Alignment requirement |
| Empty dictionary | 0 sequences | ArgumentException | Validation |
| Null input | null | ArgumentException | Validation |
| Identical sequences | All same | Valid tree, zero distances | Algorithm |
| Case variation | "acgt", "ACGT" | Treated as identical | Implementation |
| Gap-only columns | "--" in all seqs | Skipped; comparableSites=0 → distance 0 | Implementation |

---

## 4. Test Files

| File | Content |
|------|---------|
| PhylogeneticAnalyzer_TreeConstruction_Tests.cs | All PHYLO-TREE-001 tests (M01–M13, S01–S05, C01–C05) — 32 test runs |
| PhylogeneticSnapshotTests.cs | Snapshot tests for UPGMA and NJ tree output |
| PhylogeneticProperties.cs | Property-based tests for tree invariants |

---

## 5. Decisions

| Decision | Choice | Rationale |
|----------|--------|-----------|
| Null handling | Throw ArgumentException | Fail fast |
| Empty input | Throw ArgumentException | Minimum 2 sequences required |
| Default method | UPGMA | Simpler, widely used |
| Default distance | JukesCantor | Standard for DNA |

---

## 6. Deviations and Assumptions

### Deviations from Sources

None. The implementation strictly follows the algorithms as described in the authoritative sources:

| Aspect | Source Requirement | Implementation |
|--------|-------------------|----------------|
| UPGMA branch lengths | Incremental: height(new) − height(child) (Wikipedia UPGMA) | Tracks cluster heights; computes incremental branch lengths |
| UPGMA ultrametric property | All tips equidistant from root (Wikipedia UPGMA) | Verified by test S01c using Wikipedia example (16.5 for all tips) |
| NJ branch length formula | δ(f,u) = d(f,g)/2 + (Σd(f,k) − Σd(g,k))/(2(n−2)) (Wikipedia NJ) | Exact formula implemented |
| NJ negative branch lengths | "often assigns negative lengths to some branches" (Wikipedia NJ) | Negative values preserved, not clamped |
| NJ distance update formula | d(u,k) = (d(f,k) + d(g,k) − d(f,g))/2 (Wikipedia NJ eq. 3) | Exact formula implemented |
| NJ additive matrix guarantee | Correct topology for additive matrices (Wikipedia NJ) | Verified by test S02 using Wikipedia example |
| NJ final two-node join | Produces unrooted tree; rooted by convention | Last edge split d/2 each (midpoint rooting preserves patristic distances) |

### Assumptions

None. All algorithmic behavior is defined by the external sources listed in Section 1.

---

## 7. Coverage Classification (2026-03-08)

Applied systematic coverage classification:

| Test | Status | Action |
|------|--------|--------|
| M01–M06 | ✅ Covered | No changes |
| M07 | ⚠ Weak → ✅ | Strengthened: added MRCA sibling check for identical pair + distance=0 assertion |
| M08–M13 | ✅ Covered | No changes |
| S01 | ✅ Covered | No changes |
| S01b | ⚠ Weak → ✅ | Tightened tolerance 0.01→1e-10 (exact integer arithmetic); added internal branch lengths (δ(u,v)=2.5, δ(w,root)=2.5, δ(v,root)=5.5) |
| S01c | ⚠ Weak → ✅ | Tightened tolerance 0.01→1e-10 |
| S02 | ⚠ Weak → ✅ | Tightened tolerance 0.01→1e-10; absorbed INV-N01 label from removed S06 |
| S02b | ✅ Covered | No changes |
| S02c | ⚠ Weak → ✅ | Tightened tolerance 0.01→1e-10 |
| S03 | ✅ Covered | No changes |
| S04 | ✅ Covered | No changes |
| S05 | ⚠ Weak → ✅ | Replaced individual Does.Contain with Is.EquivalentTo (guards against extra leaves) |
| S06 | 🔁 Duplicate → Removed | Identical to S02: same data (WikipediaNJMatrix), same assertion (patristic=input) |
| C01 | ❌ Missing → ✅ | Implemented: 50 sequences, 100bp, asserts completion + correct taxa count |
| C02 | ❌ Missing → ✅ | Implemented: gap-only columns skipped, distance=0 for identical non-gap sites |
| C03 | ✅ Covered | Relabeled from previous C01 (was mislabeled) |
| C04 | ✅ Covered | Relabeled from previous C02 (extra test beyond spec) |
| C05 | ✅ Covered | Relabeled from previous C03 (extra test beyond spec) |

**Summary: 0 missing, 0 weak, 0 duplicate. Total: 32 test runs (was 31: −1 duplicate, +2 new).**

### Theory Verification

All test assertion values are derived from authoritative sources (Wikipedia UPGMA/NJ worked examples), not from current implementation output. Manual arithmetic trace confirms:

- **UPGMA:** Heights 8.5→11→14→16.5 from Wikipedia example (integer arithmetic, exact in IEEE 754)
- **NJ:** δ(a,u)=2, δ(b,u)=3 from Wikipedia first step; all 10 patristic distances = input matrix (additive)
- **Tolerances:** 1e-10 is appropriate — all intermediate values are representable exactly as doubles
