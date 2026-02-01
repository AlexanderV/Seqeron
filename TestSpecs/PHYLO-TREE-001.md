# PHYLO-TREE-001: Tree Construction

**Test Unit ID:** PHYLO-TREE-001  
**Category:** Phylogenetics  
**Created:** 2026-02-01  
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
| M07 | Four sequences create valid tree with all taxa | Standard case | Wikipedia example |
| M08 | Throws on single sequence | Validation | Definition |
| M09 | Throws on unequal sequence lengths | Alignment required | Definition |
| M10 | All branch lengths are non-negative | INV-05 | UPGMA property |
| M11 | UPGMA produces rooted tree | INV-U02 | Wikipedia UPGMA |
| M12 | Identical sequences produce zero-distance subtree | Algorithm | UPGMA/NJ |
| M13 | Case-insensitive sequence handling | Robustness | Implementation |

### SHOULD Tests (Enhanced Coverage)

| ID | Test | Rationale | Source |
|----|------|-----------|--------|
| S01 | Wikipedia UPGMA example produces expected clustering | Validation | Wikipedia |
| S02 | Wikipedia NJ example produces expected branch lengths | Validation | Wikipedia |
| S03 | Tree depth is logarithmic in number of taxa | Property | Binary tree |
| S04 | DistanceMatrix in result matches input calculations | Consistency | Implementation |
| S05 | Tree leaves match exact input taxon names | Traceability | Definition |

### COULD Tests (Extended Coverage)

| ID | Test | Rationale | Source |
|----|------|-----------|--------|
| C01 | Large input (50+ sequences) completes in reasonable time | Performance | O(n³) complexity |
| C02 | Gap-only columns are handled correctly | Edge case | Implementation |
| C03 | Different distance methods produce valid trees | Flexibility | Implementation |

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

---

## 4. Test Consolidation Plan

### Current State
- PhylogeneticAnalyzerTests.cs contains mixed tests for PHYLO-TREE-001, PHYLO-NEWICK-001, PHYLO-COMP-001

### Target State
- Create PhylogeneticAnalyzer_TreeConstruction_Tests.cs for PHYLO-TREE-001
- Keep PHYLO-NEWICK-001 and PHYLO-COMP-001 tests in PhylogeneticAnalyzerTests.cs temporarily
- Remove duplicate BuildTree tests from original file

### Tests to Extract
From PhylogeneticAnalyzerTests.cs:
- BuildTree_UPGMA_ReturnsValidTree
- BuildTree_NeighborJoining_ReturnsValidTree
- BuildTree_ContainsAllTaxa
- BuildTree_TwoSequences_CreatesBinaryTree
- BuildTree_ThrowsOnSingleSequence
- BuildTree_ThrowsOnUnequalLengths
- BuildTree_CaseInsensitive
- GetLeaves_ReturnsAllLeafNodes
- CalculateTreeLength_SumsAllBranches
- GetTreeDepth_ReturnsCorrectDepth

---

## 5. Decisions

| Decision | Choice | Rationale |
|----------|--------|-----------|
| Null handling | Throw ArgumentException | Fail fast |
| Empty input | Throw ArgumentException | Minimum 2 sequences required |
| Default method | UPGMA | Simpler, widely used |
| Default distance | JukesCantor | Standard for DNA |

---

## 6. Open Questions

None - algorithm behavior is well-defined in sources.
