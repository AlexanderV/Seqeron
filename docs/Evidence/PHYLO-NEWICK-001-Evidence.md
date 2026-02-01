# PHYLO-NEWICK-001: Newick I/O - Evidence Document

## Test Unit Information
- **ID:** PHYLO-NEWICK-001
- **Area:** Phylogenetic
- **Canonical Methods:** `PhylogeneticAnalyzer.ToNewick()`, `PhylogeneticAnalyzer.ParseNewick()`
- **Date:** 2026-02-01

---

## Authoritative Sources

### Primary Sources

1. **Wikipedia - Newick format**
   - URL: https://en.wikipedia.org/wiki/Newick_format
   - Description: Comprehensive documentation of Newick tree format specification
   - Key content: Grammar rules, examples, edge cases, escaping rules

2. **PHYLIP - The Newick tree format (Felsenstein)**
   - URL: https://phylipweb.github.io/phylip/newicktree.html
   - Description: Original specification by Joe Felsenstein (one of the format authors)
   - Key content: Origin, examples, uniqueness considerations

3. **Gary Olsen - Interpretation of Newick's 8:45 Tree Format**
   - URL: https://phylipweb.github.io/phylip/newick_doc.html
   - Referenced via: Wikipedia citation
   - Key content: Formal grammar description

---

## Extracted Test Data from Sources

### Wikipedia Examples

| Format | Newick String | Description |
|--------|---------------|-------------|
| No names | `(,,(,));` | All nodes unnamed |
| Leaf names only | `(A,B,(C,D));` | Only leaf nodes named |
| All names | `(A,B,(C,D)E)F;` | All nodes (internal + leaf) named |
| Branch lengths only | `(:0.1,:0.2,(:0.3,:0.4):0.5);` | All but root have distances |
| All with root distance | `(:0.1,:0.2,(:0.3,:0.4):0.5):0.0;` | All nodes have distances |
| Popular format | `(A:0.1,B:0.2,(C:0.3,D:0.4):0.5);` | Distances and leaf names |
| Full format | `(A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5)F;` | Distances and all names |
| Rooted on leaf | `((B:0.2,(C:0.3,D:0.4)E:0.5)F:0.1)A;` | Rare: rooted on leaf |

### PHYLIP Examples

| Newick String | Description |
|---------------|-------------|
| `(B,(A,C,E),D);` | Simple tree with multifurcation |
| `(B:6.0,(A:5.0,C:3.0,E:4.0):5.0,D:11.0);` | Tree with branch lengths |
| `(B:6.0,(A:5.0,C:3.0,E:4.0)Ancestor1:5.0,D:11.0);` | Internal node with name |
| `A;` | Single leaf tree |
| `((A,B),(C,D));` | Binary symmetric tree |
| `(Alpha,Beta,Gamma,Delta,,Epsilon,,,);` | Tree with empty names |

### Complex Real-World Examples (PHYLIP)

```
((raccoon:19.19959,bear:6.80041):0.84600,((sea_lion:11.99700,seal:12.00300):7.52973,((monkey:100.85930,cat:47.14069):20.59201,weasel:18.87953):2.09460):3.87382,dog:25.46154);
```

---

## Grammar Specification (Wikipedia/Olsen)

```
Tree → Subtree ";"
Subtree → Leaf | Internal
Leaf → Name
Internal → "(" BranchSet ")" Name
BranchSet → Branch | Branch "," BranchSet
Branch → Subtree Length
Name → empty | string
Length → empty | ":" number
```

---

## Documented Corner Cases (from Sources)

### From Wikipedia

| Case | Description | Source |
|------|-------------|--------|
| Empty string | Invalid - must contain tree structure | Wikipedia Grammar |
| Missing semicolon | Format ends with semicolon | Wikipedia Grammar |
| Whitespace | Allowed anywhere except in names/numbers | Wikipedia Notes |
| Underscore → blank | Underscore in unquoted string becomes blank | Wikipedia Notes |
| Single quotes | String can be quoted with single quotes | Wikipedia Notes |
| Nested comments | Comments in `[]` allowed in some dialects | Wikipedia Notes |
| Unnamed nodes | `(,,(,));` - all nodes can be empty | Wikipedia Examples |
| Single taxon | `A;` - valid single leaf tree | PHYLIP Examples |

### From PHYLIP

| Case | Description | Source |
|------|-------------|--------|
| Multifurcation | `(A,B,C,D);` - more than 2 children | PHYLIP spec |
| Non-uniqueness | Multiple representations for same tree | PHYLIP Non-Uniqueness |
| Left-right order | Order of descendants is arbitrary | PHYLIP Non-Uniqueness |
| Rooted on leaf | `((B,C),A);` - rare but valid | PHYLIP Rooted trees |

---

## Implementation-Specific Constraints

Based on current implementation (`PhylogeneticAnalyzer.cs`):

1. **Binary trees only**: Implementation builds binary trees (UPGMA/NJ produce bifurcating trees)
2. **Branch lengths**: Supported in both parsing and export
3. **Internal node names**: Supported in parsing, auto-generated in export
4. **No quoted names**: Current implementation does not handle quoted names
5. **No comments**: Comments in `[]` not supported
6. **Standard format**: Follows popular format `(A:0.1,B:0.2);`

---

## Invariants (Derived from Sources)

| ID | Invariant | Source |
|----|-----------|--------|
| N1 | Newick string MUST end with semicolon | Wikipedia Grammar |
| N2 | Leaf count in parsed tree MUST match leaf names in string | Semantic |
| N3 | Round-trip (ToNewick → ParseNewick) MUST preserve topology | Semantic |
| N4 | Round-trip MUST preserve leaf names | Semantic |
| N5 | Branch lengths MUST be non-negative | PHYLIP implicit |
| N6 | ParseNewick on empty string MUST throw | Implementation |
| N7 | ToNewick on null node SHOULD return empty or throw | Implementation |

---

## Test Categories (Evidence-Based)

### MUST Tests (Core Functionality)

1. **Semicolon termination** - ToNewick output ends with `;` (N1)
2. **Parse simple binary tree** - `(A,B);` parses correctly
3. **Parse with branch lengths** - `(A:0.1,B:0.2);` extracts values
4. **Parse nested tree** - `((A,B),(C,D));` creates correct structure
5. **Leaf count preservation** - Parsed tree has correct number of leaves
6. **Round-trip topology** - ToNewick → ParseNewick preserves structure (N3, N4)
7. **Empty string throws** - ParseNewick("") throws exception (N6)

### SHOULD Tests (Extended Functionality)

1. **Internal node names** - `(A,B)Root;` parses internal name
2. **Mixed format** - `(A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5)F;` full format
3. **Whitespace handling** - `( A , B );` with whitespace
4. **Branch length precision** - Values preserved with adequate precision

### COULD Tests (Optional/Dialect-Specific)

1. **Single taxon** - `A;` single leaf tree
2. **Empty names** - `(,);` unnamed leaves (not fully supported)
3. **Quoted names** - `'A B':0.1` quoted with spaces (not implemented)

---

## Quality Criteria

- [x] Sources are authoritative (Wikipedia, PHYLIP)
- [x] Test cases traceable to documented examples
- [x] Corner cases from grammar specification
- [x] Implementation constraints documented
- [x] Invariants clearly defined

---

## References

1. Wikipedia contributors. "Newick format." Wikipedia, The Free Encyclopedia.
   https://en.wikipedia.org/wiki/Newick_format

2. Felsenstein, J. "The Newick tree format." PHYLIP documentation.
   https://phylipweb.github.io/phylip/newicktree.html

3. Olsen, G. "Interpretation of Newick's 8:45 Tree Format." 30 August 1990.
   https://phylipweb.github.io/phylip/newick_doc.html
