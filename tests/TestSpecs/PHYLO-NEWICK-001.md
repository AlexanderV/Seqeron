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
| N1 | Newick string MUST end with semicolon | Wikipedia Grammar: `Tree → Subtree ";"` |
| N2 | Leaf count in parsed tree MUST match leaf names in string | Wikipedia Grammar: leaves are `Name` productions |
| N3 | Round-trip (ToNewick → ParseNewick) MUST preserve topology | Wikipedia Grammar: unambiguous binary structure |
| N4 | Round-trip MUST preserve leaf names | Wikipedia Grammar: `Leaf → Name` |
| N5 | Branch lengths are real numbers after `:` | Wikipedia Grammar: `Length → empty \| ":" number` |
| N6 | ParseNewick on empty/null MUST throw ArgumentException | Wikipedia Grammar: at minimum `Subtree ";"` required |
| N7 | ToNewick on null SHOULD return empty string | Defensive programming |
| N8 | Internal node names MUST appear after `)` in output | Wikipedia Grammar: `Internal → "(" BranchSet ")" Name` |
| N9 | Number formatting MUST use `.` as decimal separator | Olsen: `branch_length ==> signed_number \| unsigned_number` (no locale) |

---

## Test Specification

### MUST Tests (Required for DoD)

| Test ID | Test Name | Invariant | Source |
|---------|-----------|-----------|--------|
| M01 | ToNewick_SimpleTree_EndsWithSemicolon | N1 | Wikipedia Grammar: `Tree → Subtree ";"` |
| M02 | ToNewick_WithBranchLengths_IncludesColons | N5 | Wikipedia Examples: `(A:0.1,B:0.2);` |
| M03 | ToNewick_WithoutBranchLengths_NoColons | N5 | Wikipedia Grammar: `Length → empty` |
| M04 | ToNewick_ContainsAllLeafNames | N2 | Wikipedia Grammar: `Leaf → Name` |
| M05 | ParseNewick_SimpleBinaryTree_ParsesCorrectly | - | Wikipedia Example: `(A,B);` |
| M06 | ParseNewick_WithBranchLengths_ExtractsValues | N5 | Wikipedia Example: `(A:0.1,B:0.2);` |
| M07 | ParseNewick_NestedBinaryTree_ParsesRecursively | - | Wikipedia Example: `((A,B),(C,D));` |
| M08 | ParseNewick_LeafCountMatchesInput | N2 | PHYLIP Examples |
| M09 | RoundTrip_PreservesLeafNames | N4 | Wikipedia Grammar: round-trip invariant |
| M10 | RoundTrip_PreservesTopology | N3 | Wikipedia Grammar: round-trip invariant |
| M11 | ParseNewick_EmptyString_ThrowsArgumentException | N6 | Wikipedia Grammar: minimum structure required |
| M12 | ParseNewick_NullString_ThrowsArgumentException | N6 | Wikipedia Grammar: minimum structure required |
| M13 | ToNewick_NullNode_ReturnsEmpty | N7 | Defensive edge case |
| M14 | ParseNewick_RootBranchLength_ExtractsValue | N5 | Olsen Grammar: `tree ==> descendant_list [root_label] [:branch_length] ;` |
| M15 | ParseNewick_RootNameAndBranchLength_ParsesCorrectly | N5, N8 | Olsen Grammar (same as M14) |
| M16 | RoundTrip_FullFormat_PreservesInternalNames | N8 | Wikipedia Grammar: `Internal → "(" BranchSet ")" Name` |
| M17 | ToNewick_WithValidInternalNames_EmitsNames | N8 | Wikipedia Grammar: `Internal → "(" BranchSet ")" Name` |
| M18 | ToNewick_WithMetacharacterNames_OmitsInvalidNames | N9 | Olsen: unquoted labels prohibit metacharacters |

### SHOULD Tests (Recommended)

| Test ID | Test Name | Rationale | Source |
|---------|-----------|-----------|--------|
| S01 | ParseNewick_InternalNodeNames_ExtractsName | Full format support | Wikipedia Example: `(A,B)E;` all names |
| S02 | ParseNewick_FullFormat_ParsesAllFields | Complete parsing | Wikipedia Example: full format |
| S03 | RoundTrip_PreservesBranchLengths_WithTolerance | Precision | Wikipedia Grammar: `Length → ":" number` |
| S04 | ToNewick_LargeTree_ProducesValidFormat | Grammar validity at scale | Wikipedia Grammar: balanced `()`, trailing `;`, leaf `Name` |
| S05 | ParseNewick_WhitespaceAtEnd_HandlesGracefully | Robustness | Wikipedia Notes: "Whitespace may appear anywhere except within unquoted string or Length" |
| S06 | ParseNewick_MissingSemicolon_ParsesCorrectly | Lenient parsing | Wikipedia Grammar: `Tree → Subtree ";"` requires `;`; parser is lenient |

### COULD Tests (Optional)

| Test ID | Test Name | Rationale | Source |
|---------|-----------|-----------|--------|
| C01 | ParseNewick_SingleTaxon_ParsesCorrectly | Rare format | PHYLIP Example: `A;` |

---

## Edge Cases (from Evidence)

| Case | Input | Expected | Source | Covered By |
|------|-------|----------|--------|------------|
| Empty string | `""` | Throw ArgumentException | Wikipedia Grammar: min structure required | M11 |
| Null string | `null` | Throw ArgumentException | Wikipedia Grammar: min structure required | M12 |
| Null node | `ToNewick(null)` | Return empty string | Defensive | M13 |
| Whitespace only | `"   "` | Throw ArgumentException | Wikipedia Notes: whitespace is not a tree | S05 |
| Missing semicolon | `"(A,B)"` | Parse leniently (strip `;` if present) | Wikipedia Grammar requires `;`; lenient parser | S06 |
| Single leaf | `"A;"` | Single node tree | PHYLIP Examples | C01 |
| Root branch length | `"(A,B):0.0;"` | Parse root branch length | Olsen Grammar: `tree ==> ... [:branch_length] ;` | M14 |
| Root name + length | `"(A,B)Root:0.5;"` | Parse root name and length | Olsen Grammar | M15 |
| Invalid internal name | `ToNewick` of UPGMA tree | Omit metacharacter names | Olsen: unquoted label restrictions | M18 |

---

## Coverage Classification

| Test ID | Test Name | Status | Notes |
|---------|-----------|--------|-------|
| M01 | ToNewick_SimpleTree_EndsWithSemicolon | ✅ Covered | N1: semicolon termination |
| M02 | ToNewick_WithBranchLengths_IncludesColons | ✅ Covered | N5: `Name:number` format verified via regex |
| M03 | ToNewick_WithoutBranchLengths_NoColons | ✅ Covered | N5: absence of colons when disabled |
| M04 | ToNewick_ContainsAllLeafNames | ✅ Covered | N2: all 4 leaf names present |
| M05 | ParseNewick_SimpleBinaryTree_ParsesCorrectly | ✅ Covered | Structure + leaf names (A, B) |
| M06 | ParseNewick_WithBranchLengths_ExtractsValues | ✅ Covered | N5: exact values 0.1, 0.2 |
| M07 | ParseNewick_NestedBinaryTree_ParsesRecursively | ✅ Covered | Full structure: leaf names + positions |
| M08 | ParseNewick_LeafCountMatchesInput | ✅ Covered | N2: 2/3/4/5 taxa test cases |
| M09 | RoundTrip_PreservesLeafNames | ✅ Covered | N4: sorted name lists match |
| M10 | RoundTrip_PreservesTopology | ✅ Covered | N3: sibling relationships (A,B) and (C,D) |
| M11 | ParseNewick_EmptyString_ThrowsArgumentException | ✅ Covered | N6: ArgumentException with message |
| M12 | ParseNewick_NullString_ThrowsArgumentException | ✅ Covered | N6: ArgumentException |
| M13 | ToNewick_NullNode_ReturnsEmptyString | ✅ Covered | N7: returns empty |
| M14 | ParseNewick_RootBranchLength_ExtractsValue | ✅ Covered | N5: root BL=0.0, children BL=0.1/0.2 |
| M15 | ParseNewick_RootNameAndBranchLength_ParsesCorrectly | ✅ Covered | N5+N8: name="Root", BL=0.5 |
| M16 | RoundTrip_FullFormat_PreservesInternalNames | ✅ Covered | N8: Root, AB, CD names survive |
| M17 | ToNewick_WithValidInternalNames_EmitsNames | ✅ Covered | N8: ``)AB:`, ``)CD:`, ``Root;`` patterns |
| M18 | ToNewick_WithMetacharacterNames_OmitsInvalidNames | ✅ Covered | N9: metacharacter names suppressed |
| S01 | ParseNewick_InternalNodeNames_ExtractsName | ✅ Covered | Exact name "Root" |
| S02 | ParseNewick_FullFormat_ParsesAllFields | ✅ Covered | All names + all 6 branch lengths |
| S03 | RoundTrip_PreservesBranchLengths | ✅ Covered | 6 individual BLs within ±0.0001 |
| S04 | ToNewick_LargeTree_ProducesValidFormat | ✅ Covered | Balanced parens, 8 taxa, semicolon |
| S05a | ParseNewick_WhitespaceOnly_ThrowsArgumentException | ✅ Covered | ArgumentException |
| S05b | ParseNewick_TrailingWhitespace_ParsesCorrectly | ✅ Covered | Parses with trailing spaces |
| S06 | ParseNewick_MissingSemicolon_ParsesCorrectly | ✅ Covered | Leaf names A, B verified |
| C01 | ParseNewick_SingleTaxon_ParsesSingleNode | ✅ Covered | Name="A", IsLeaf=true |

**Summary:** 26/26 tests ✅ Covered. No missing, no duplicates.

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
│   ├── ToNewick_LargeTree_ProducesValidFormat
│   ├── ToNewick_WithValidInternalNames_EmitsNames
│   └── ToNewick_WithMetacharacterNames_OmitsInvalidNames
├── ParseNewick Tests
│   ├── ParseNewick_SimpleBinaryTree_ParsesCorrectly
│   ├── ParseNewick_WithBranchLengths_ExtractsValues
│   ├── ParseNewick_NestedBinaryTree_ParsesRecursively
│   ├── ParseNewick_LeafCountMatchesInput
│   ├── ParseNewick_InternalNodeNames_ExtractsName
│   ├── ParseNewick_FullFormat_ParsesAllFields
│   ├── ParseNewick_SingleTaxon_ParsesSingleNode
│   ├── ParseNewick_RootBranchLength_ExtractsValue
│   └── ParseNewick_RootNameAndBranchLength_ParsesCorrectly
├── Round-Trip Tests
│   ├── RoundTrip_PreservesLeafNames
│   ├── RoundTrip_PreservesTopology
│   ├── RoundTrip_PreservesBranchLengths
│   └── RoundTrip_FullFormat_PreservesInternalNames
└── Edge Cases
    ├── ParseNewick_EmptyString_ThrowsArgumentException
    ├── ParseNewick_NullString_ThrowsArgumentException
    ├── ParseNewick_WhitespaceOnly_ThrowsArgumentException
    ├── ParseNewick_TrailingWhitespace_ParsesCorrectly
    └── ParseNewick_MissingSemicolon_ParsesCorrectly
```

---

## Documented Limitations (per Scope)

These are intentional scope boundaries, documented with their spec basis.

| Limitation | Spec Reference | Rationale |
|------------|---------------|-----------|
| Binary trees only | Wikipedia Grammar: `BranchSet → Branch \| Branch "," BranchSet` supports N children; implementation restricted to 2 | UPGMA/NJ produce bifurcating trees; multifurcation not needed for current algorithms |
| No quoted names | Olsen: `quoted_label ==> ' string '`; not implemented | Out of scope; all current taxa use simple alphanumeric names |
| No `[]` comments | Wikipedia Notes: "Comments are enclosed in square brackets"; not implemented | Out of scope; no use case for comment parsing |
| No underscore→blank | PHYLIP: "underscore stands for a blank"; not implemented | Out of scope; no taxa use underscores |
| F4 branch length precision | Wikipedia Grammar: `Length → ":" number`; spec imposes no precision limit; F4 chosen for readability | Adequate for UPGMA/NJ round-trip (±0.00005); matches typical bioinformatics output |

---

## References

- Evidence: [PHYLO-NEWICK-001-Evidence.md](../docs/Evidence/PHYLO-NEWICK-001-Evidence.md)
- Algorithm Doc: [docs/algorithms/Phylogenetics/Newick_Format.md](../docs/algorithms/Phylogenetics/Newick_Format.md)
- Wikipedia: [Newick format](https://en.wikipedia.org/wiki/Newick_format) — Grammar, examples, notes
- PHYLIP: [The Newick tree format](https://phylipweb.github.io/phylip/newicktree.html) — Felsenstein, original specification
- Olsen: [Interpretation of Newick's 8:45 Tree Format](https://phylipweb.github.io/phylip/newick_doc.html) — Formal grammar (1990)
