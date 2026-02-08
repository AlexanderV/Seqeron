# PHYLO-NEWICK-001: Newick I/O

## Test Unit Summary

| Field | Value |
|-------|-------|
| **ID** | PHYLO-NEWICK-001 |
| **Area** | Phylogenetic |
| **Canonical Methods** | `ToNewick(PhyloNode, bool)`, `ParseNewick(string)` |
| **Complexity** | O(n) where n = tree size |
| **Status** | ☑ Complete |

---

## Methods Under Test

| Method | Class | Type | Test Depth |
|--------|-------|------|------------|
| `ToNewick(PhyloNode, bool)` | PhylogeneticAnalyzer | Canonical | Deep |
| `ParseNewick(string)` | PhylogeneticAnalyzer | Canonical | Deep |

---

## Invariants (from Evidence)

| ID | Invariant | Source |
|----|-----------|--------|
| N1 | Newick string MUST end with semicolon | Wikipedia Grammar |
| N2 | Leaf count in parsed tree MUST match leaf names in string | Semantic |
| N3 | Round-trip (ToNewick → ParseNewick) MUST preserve topology | Semantic |
| N4 | Round-trip MUST preserve leaf names | Semantic |
| N5 | Branch lengths MUST be non-negative when parsed | PHYLIP implicit |
| N6 | ParseNewick on empty/null MUST throw ArgumentException | Implementation |
| N7 | ToNewick on null SHOULD return empty string | Implementation |

---

## Test Specification

### MUST Tests (Required for DoD)

| Test ID | Test Name | Invariant | Source |
|---------|-----------|-----------|--------|
| M01 | ToNewick_SimpleTree_EndsWithSemicolon | N1 | Wikipedia Grammar |
| M02 | ToNewick_WithBranchLengths_IncludesColons | - | Wikipedia Examples |
| M03 | ToNewick_WithoutBranchLengths_NoColons | - | Wikipedia Examples |
| M04 | ToNewick_ContainsAllLeafNames | N2 | Semantic |
| M05 | ParseNewick_SimpleBinaryTree_ParsesCorrectly | - | Wikipedia Example `(A,B);` |
| M06 | ParseNewick_WithBranchLengths_ExtractsValues | N5 | Wikipedia Example `(A:0.1,B:0.2);` |
| M07 | ParseNewick_NestedBinaryTree_ParsesRecursively | - | Wikipedia Example `((A,B),(C,D));` |
| M08 | ParseNewick_LeafCountMatchesInput | N2 | PHYLIP Examples |
| M09 | RoundTrip_PreservesLeafNames | N4 | Semantic invariant |
| M10 | RoundTrip_PreservesTopology | N3 | Semantic invariant |
| M11 | ParseNewick_EmptyString_ThrowsArgumentException | N6 | Edge case |
| M12 | ParseNewick_NullString_ThrowsArgumentException | N6 | Edge case |
| M13 | ToNewick_NullNode_ReturnsEmpty | N7 | Edge case |

### SHOULD Tests (Recommended)

| Test ID | Test Name | Rationale | Source |
|---------|-----------|-----------|--------|
| S01 | ParseNewick_InternalNodeNames_ExtractsName | Full format support | Wikipedia Example `(A,B)E;` |
| S02 | ParseNewick_FullFormat_ParsesAllFields | Complete parsing | Wikipedia Example full format |
| S03 | RoundTrip_PreservesBranchLengths_WithTolerance | Precision | Semantic |
| S04 | ToNewick_LargeTree_ProducesValidFormat | Scalability | ASSUMPTION |
| S05 | ParseNewick_WhitespaceAtEnd_HandlesGracefully | Robustness | Wikipedia Notes |

### COULD Tests (Optional)

| Test ID | Test Name | Rationale | Source |
|---------|-----------|-----------|--------|
| C01 | ParseNewick_SingleTaxon_ParsesCorrectly | Rare format | PHYLIP Example `A;` |
| C02 | ParseNewick_MultifurcatingTree_ParsesChildren | Beyond binary | Wikipedia Grammar |

---

## Edge Cases (from Evidence)

| Case | Input | Expected | Covered By |
|------|-------|----------|------------|
| Empty string | `""` | Throw ArgumentException | M11 |
| Null string | `null` | Throw ArgumentException | M12 |
| Null node | `ToNewick(null)` | Return empty string | M13 |
| Whitespace only | `"   "` | Throw ArgumentException | S05 |
| Missing semicolon | `"(A,B)"` | Parse should handle | ASSUMPTION (impl adds) |
| Single leaf | `"A;"` | Single node tree | C01 |

---

## Audit of Existing Tests

### Current State (PhylogeneticAnalyzerTests.cs)

| Test | Status | Action |
|------|--------|--------|
| `ToNewick_SimpleTree_ProducesValidFormat` | Keep | Rename to match pattern |
| `ToNewick_WithBranchLengths_IncludesColons` | Keep | Good coverage |
| `ToNewick_WithoutBranchLengths_NoColons` | Keep | Good coverage |
| `ParseNewick_SimpleTree_ParsesCorrectly` | Keep | Good basic test |
| `ParseNewick_WithBranchLengths_ExtractsValues` | Keep | Good coverage |
| `ParseNewick_NestedTree_ParsesRecursively` | Keep | Good coverage |
| `ParseNewick_RoundTrip_PreservesStructure` | Keep | Rename for clarity |
| `ParseNewick_EmptyString_Throws` | Keep | Edge case coverage |

### Consolidation Plan

1. **Extract PHYLO-NEWICK-001 tests** from `PhylogeneticAnalyzerTests.cs` into new canonical file
2. ~~**Add missing tests**: M04, M08, M10, M12, M13, S01-S05~~ ✅ All added
3. **Rename tests** to follow `Method_Scenario_ExpectedResult` pattern
4. **Add Assert.Multiple** for invariant groups
5. **Remove PHYLO-COMP-001 tests** from the original file (leave for separate processing)

---

## Test File Structure

**Canonical File:** `PhylogeneticAnalyzer_NewickIO_Tests.cs`

```
PhylogeneticAnalyzer_NewickIO_Tests
├── ToNewick Tests
│   ├── ToNewick_SimpleTree_EndsWithSemicolon
│   ├── ToNewick_WithBranchLengths_IncludesColons
│   ├── ToNewick_WithoutBranchLengths_NoColons
│   ├── ToNewick_ContainsAllLeafNames
│   ├── ToNewick_NullNode_ReturnsEmptyString
│   └── ToNewick_LargeTree_ProducesValidFormat
├── ParseNewick Tests
│   ├── ParseNewick_SimpleBinaryTree_ParsesCorrectly
│   ├── ParseNewick_WithBranchLengths_ExtractsValues
│   ├── ParseNewick_NestedBinaryTree_ParsesRecursively
│   ├── ParseNewick_LeafCountMatchesInput
│   ├── ParseNewick_InternalNodeNames_ExtractsName
│   ├── ParseNewick_FullFormat_ParsesAllFields
│   └── ParseNewick_SingleTaxon_ParsesSingleNode
├── Round-Trip Tests
│   ├── RoundTrip_PreservesLeafNames
│   ├── RoundTrip_PreservesTopology
│   └── RoundTrip_PreservesBranchLengths
└── Edge Cases
    ├── ParseNewick_EmptyString_ThrowsArgumentException
    ├── ParseNewick_NullString_ThrowsArgumentException
    └── ParseNewick_WhitespaceOnly_ThrowsArgumentException
```

---

## Open Questions / Decisions

| Question | Decision | Rationale |
|----------|----------|-----------|
| Support multifurcating trees? | No (implementation is binary) | Implementation constraint |
| Support quoted names? | No (not implemented) | Document as limitation |
| Support comments in `[]`? | No | Not implemented |

---

## References

- Evidence: [PHYLO-NEWICK-001-Evidence.md](../docs/Evidence/PHYLO-NEWICK-001-Evidence.md)
- Algorithm Doc: [docs/algorithms/Phylogenetics/Newick_Format.md](../docs/algorithms/Phylogenetics/Newick_Format.md)
