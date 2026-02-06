# Test Specification: PHYLO-COMP-001

**Test Unit ID:** PHYLO-COMP-001
**Title:** Tree Comparison
**Algorithm Group:** Phylogenetics
**Status:** ☑ Complete
**Last Updated:** 2026-02-01

---

## Scope

### Canonical Methods (Deep Testing Required)
| Method | Class | Complexity |
|--------|-------|------------|
| `RobinsonFouldsDistance(PhyloNode, PhyloNode)` | PhylogeneticAnalyzer | O(n) |
| `FindMRCA(PhyloNode, string, string)` | PhylogeneticAnalyzer | O(n) |
| `PatristicDistance(PhyloNode, string, string)` | PhylogeneticAnalyzer | O(n) |

### Helper Methods (Smoke Testing)
| Method | Class | Notes |
|--------|-------|-------|
| `GetLeaves(PhyloNode)` | PhylogeneticAnalyzer | Tree traversal helper |
| `CalculateTreeLength(PhyloNode)` | PhylogeneticAnalyzer | Sum of branch lengths |
| `GetTreeDepth(PhyloNode)` | PhylogeneticAnalyzer | Tree height |

---

## Test Cases

### MUST Tests (Required - Evidence-Based)

#### Robinson-Foulds Distance

| ID | Test Case | Evidence Source | Expected |
|----|-----------|-----------------|----------|
| RF-M01 | Identical trees have RF distance 0 | Wikipedia: "RF distance of 0 indicates identical trees" | RF = 0 |
| RF-M02 | RF distance is symmetric | Robinson & Foulds (1981): metric property | RF(T1,T2) = RF(T2,T1) |
| RF-M03 | RF distance is non-negative | Robinson & Foulds (1981): metric property | RF ≥ 0 |
| RF-M04 | Different topologies have positive RF | Wikipedia: symmetric difference | RF > 0 |
| RF-M05 | RF distance is even | Wikipedia: symmetric difference of two sets | RF % 2 = 0 |

#### FindMRCA

| ID | Test Case | Evidence Source | Expected |
|----|-----------|-----------------|----------|
| MRCA-M01 | Same taxon returns itself | Wikipedia MRCA: "most recent common ancestor" | Node name = taxon |
| MRCA-M02 | Sibling taxa return parent | Wikipedia: MRCA is deepest common ancestor | Parent node |
| MRCA-M03 | Distant taxa return deeper ancestor | Wikipedia: MRCA definition | Ancestor node containing both |
| MRCA-M04 | MRCA contains both taxa | Wikipedia: MRCA properties | Both taxa in subtree |
| MRCA-M05 | Null root returns null | Edge case | null |

#### PatristicDistance

| ID | Test Case | Evidence Source | Expected |
|----|-----------|-----------------|----------|
| PD-M01 | Same taxon has distance 0 | Definition: no path to traverse | Distance = 0.0 |
| PD-M02 | Sibling distance = sum of branch lengths | Definition | Sum of both branches |
| PD-M03 | Distance is symmetric | Metric property | PD(x,y) = PD(y,x) |
| PD-M04 | Non-existent taxon returns NaN | Edge case | double.NaN |
| PD-M05 | Distance is non-negative | Metric property | PD ≥ 0 |

### SHOULD Tests (Recommended)

| ID | Test Case | Rationale |
|----|-----------|-----------|
| RF-S01 | Three taxa tree comparison | Simple non-trivial case |
| RF-S02 | Four taxa with different groupings | Classic RF example |
| MRCA-S01 | MRCA of all taxa is root | Property verification |
| PD-S01 | Verify path through specific MRCA | Algorithm correctness |

### COULD Tests (Optional)

| ID | Test Case | Rationale |
|----|-----------|-----------|
| RF-C01 | Large tree performance | Complexity verification |
| MRCA-C01 | Very deep tree MRCA | Stack depth testing |
| PD-C01 | Complex tree path verification | Additional confidence |

---

## Invariants to Verify

### Robinson-Foulds
1. ∀T: RF(T, T) = 0
2. ∀T1, T2: RF(T1, T2) = RF(T2, T1)
3. ∀T1, T2: RF(T1, T2) ≥ 0
4. ∀T1, T2: RF(T1, T2) is even

### MRCA
1. ∀x: MRCA(x, x) represents x
2. ∀x, y: MRCA(x, y) = MRCA(y, x)
3. ∀x, y: MRCA(x, y) is an ancestor of both x and y

### Patristic Distance
1. ∀x: PD(x, x) = 0
2. ∀x, y: PD(x, y) = PD(y, x)
3. ∀x, y: PD(x, y) ≥ 0
4. ∀x, y: PD(x, y) = d(x, MRCA) + d(y, MRCA)

---

## Audit Notes

### Existing Test Analysis

**File:** `PhylogeneticAnalyzerTests.cs`

| Test | Classification | Action |
|------|----------------|--------|
| GetLeaves_ReturnsAllLeafNodes | Covered (helper) | Keep as smoke test |
| CalculateTreeLength_SumsAllBranches | Covered (helper) | Keep as smoke test |
| GetTreeDepth_ReturnsCorrectDepth | Covered (helper) | Keep as smoke test |
| FindMRCA_FindsCommonAncestor | Covered | Enhance with explicit assertions |
| FindMRCA_SameTaxon_ReturnsTaxonItself | Covered | Keep |
| PatristicDistance_CalculatesTreePathDistance | Weak | Strengthen with explicit values |
| PatristicDistance_SameTaxon_ReturnsZero | Covered | Keep |
| RobinsonFouldsDistance_IdenticalTrees_ReturnsZero | Covered | Keep |
| RobinsonFouldsDistance_DifferentTrees_ReturnsPositive | Weak | Strengthen |
| Bootstrap tests | Different scope | Keep separate |
| GetLeaves_NullRoot_ReturnsEmpty | Covered (edge) | Keep |
| CalculateTreeLength_NullRoot_ReturnsZero | Covered (edge) | Keep |

### Consolidation Plan
- Rename file to `PhylogeneticAnalyzer_TreeComparison_Tests.cs`
- Add missing MUST tests
- Strengthen weak tests with explicit expected values
- Keep existing good tests
- Remove duplicate coverage

---

## Test Data

### Standard Test Trees

```
Tree A (4 taxa):
    root
   /    \
  AB     CD
 / \    / \
A   B  C   D

Tree B (4 taxa - different topology):
    root
   /    \
  AC     BD
 / \    / \
A   C  B   D
```

### Branch Length Test Tree
```
    root (BL=0)
   /          \
  X (BL=1.5)   Y (BL=2.0)
 / \          / \
A   B        C   D
(0.5)(0.5)  (1.0)(1.0)

Patristic distances:
- A to B: 0.5 + 0.5 = 1.0
- A to C: 0.5 + 1.5 + 2.0 + 1.0 = 5.0
- C to D: 1.0 + 1.0 = 2.0
```

---

## Open Questions / Decisions

1. **Q:** Should RF distance normalize by maximum possible?
   **A:** ASSUMPTION: Current implementation returns raw count, not normalized.

2. **Q:** Behavior with multifurcating trees?
   **A:** ASSUMPTION: Implementation assumes binary trees.

---

## Definition of Done

- [x] Evidence document created
- [x] Algorithm documentation created
- [x] TestSpec created
- [ ] Tests implemented
- [ ] All MUST tests pass
- [ ] Zero warnings
- [ ] Checklist updated
