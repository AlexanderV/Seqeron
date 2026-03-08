# Test Specification: PHYLO-COMP-001

**Test Unit ID:** PHYLO-COMP-001
**Title:** Tree Comparison
**Algorithm Group:** Phylogenetics
**Status:** ☑ Complete
**Last Updated:** 2026-03-08

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
| RF-M03 | RF distance is non-negative — exact values on diverse tree sizes | Robinson & Foulds (1981): metric property | RF = 0, 0, 2 |
| RF-M04 | Different topologies have positive RF — exact value on 3-taxa trees | Wikipedia: symmetric difference | RF = 2 |
| RF-M05 | RF distance is even — exact values on multiple tree sizes | Wikipedia: symmetric difference of two sets | 4%2=0, 2%2=0, 0%2=0 |
| RF-M06 | Four-taxa max-difference RF = 4 | Rooted RF: 2(n-2) = 4 for n=4 | RF = 4 |

#### FindMRCA

| ID | Test Case | Evidence Source | Expected |
|----|-----------|-----------------|----------|
| MRCA-M01 | Same taxon returns itself | Wikipedia MRCA: "most recent common ancestor" | Node name = taxon |
| MRCA-M02 | Sibling taxa return parent | Wikipedia: MRCA is deepest common ancestor | Parent node |
| MRCA-M03 | Distant taxa return deeper ancestor | Wikipedia: MRCA definition | Ancestor node containing both |
| MRCA-M04 | MRCA contains both taxa | Wikipedia: MRCA properties | Both taxa in subtree |
| MRCA-M05 | Null root returns null | Edge case | null |
| MRCA-M06 | Non-existent taxon returns null | Edge case: Evidence doc | null |

#### PatristicDistance

| ID | Test Case | Evidence Source | Expected |
|----|-----------|-----------------|----------|
| PD-M01 | Same taxon has distance 0 | Definition: no path to traverse | Distance = 0.0 |
| PD-M02 | Sibling distance = sum of branch lengths | Definition | Sum of both branches |
| PD-M03 | Distance is symmetric | Metric property | PD(x,y) = PD(y,x) |
| PD-M04 | Non-existent taxon returns NaN | Edge case | double.NaN |
| PD-M05 | Distance is non-negative — exact values on three-taxa tree | Metric property | PD(A,B)=1.0, PD(A,C)=3.0, PD(B,C)=3.0 |

### SHOULD Tests (Recommended)

| ID | Test Case | Rationale |
|----|-----------|-----------|
| RF-S01 | Three taxa tree comparison (same topology) | RF = 0 |
| RF-S02 | Three taxa different topology | RF = 2 (exact, max for n=3 rooted) |
| RF-S03 | RF bounded by 2(n-2) | Max bound verification |  
| MRCA-S01 | MRCA of all taxa is root | Property verification |
| PD-S01 | All documented PD values for four-taxa tree | A-B=1.0, A-C=5.0, C-D=2.0 |
| PD-S02 | PD triangle inequality | Metric property |

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
4. ∀T1, T2: RF(T1, T2) is even (for binary trees with same leaf set)
5. ∀T1, T2: RF(T1, T2) ≤ 2(n-2) for rooted binary trees with n taxa

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

**File:** `PhylogeneticAnalyzer_TreeComparison_Tests.cs`

All MUST, SHOULD, and exact-value tests implemented and passing.
Helper methods retained as smoke tests. Bootstrap tests kept in same file for coverage.

### Changes (2026-03-08)

**Implementation fixes:**
- `FindMRCA`: fixed to return `null` when either taxon does not exist (was returning the found node)
- `RobinsonFouldsDistance`: replaced canonical-split approach with direct clade comparison, matching the documented rooted RF definition and preventing complementary-clade collisions for n ≥ 5
- Bootstrap method updated to use clade terminology consistently

**Tests strengthened (⚠ Weak → exact values):**
- RF-M03: `Is.GreaterThanOrEqualTo(0)` → exact values on 3 diverse tree pairs
- RF-M04: same trees as RF-M06 with `Is.GreaterThan(0)` → different trees (3-taxa) with exact RF = 2
- RF-M05: single pair evenness → exact values + evenness on 3 tree sizes
- PD-M05: `Is.GreaterThanOrEqualTo(0)` → exact values on 3-taxa tree (1.0, 3.0, 3.0)

**Tests added (❌ Missing):**
- MRCA-S01: cross-clade taxa pairs all return root (4 pairs verified)
- MRCA-M06: non-existent taxon returns null (one-missing and both-missing cases)

**Labels fixed (🔁 Duplicate/mislabeled):**
- Sequence-built RF test: "RF-S02" → "RF integration"
- MRCA symmetry test: "MRCA-S01" → "MRCA invariant" (was mislabeled; real MRCA-S01 now implemented)

### Previous Consolidation
- Renamed file to `PhylogeneticAnalyzer_TreeComparison_Tests.cs`
- All MUST tests added
- Exact-value tests for RF and PD added
- Triangle inequality test for PD added
- RF max bound test added
- Three-taxa different-topology RF test added

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

## Definition of Done

- [x] Evidence document created
- [x] Algorithm documentation created
- [x] TestSpec created
- [x] Tests implemented
- [x] All MUST tests pass
- [x] Zero warnings
- [x] Checklist updated
