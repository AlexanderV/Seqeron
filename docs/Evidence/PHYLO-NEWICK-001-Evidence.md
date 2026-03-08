# PHYLO-NEWICK-001: Newick I/O - Evidence Document

## Test Unit Information
- **ID:** PHYLO-NEWICK-001
- **Area:** Phylogenetic
- **Canonical Methods:** `PhylogeneticAnalyzer.ToNewick()`, `PhylogeneticAnalyzer.ParseNewick()`
- **Date:** 2026-03-08

---

## Authoritative Sources

### Primary Sources

1. **Wikipedia - Newick format**
   - URL: https://en.wikipedia.org/wiki/Newick_format
   - Description: Comprehensive documentation of Newick tree format specification
   - Key content: Grammar rules, examples, edge cases, escaping rules, label restrictions

2. **PHYLIP - The Newick tree format (Felsenstein)**
   - URL: https://phylipweb.github.io/phylip/newicktree.html
   - Description: Original specification by Joe Felsenstein (one of the format authors)
   - Key content: Origin, examples, uniqueness considerations, underscore convention

3. **Gary Olsen - Interpretation of Newick's 8:45 Tree Format (1990)**
   - URL: https://phylipweb.github.io/phylip/newick_doc.html
   - Description: Formal grammar description
   - Key content: Formal BNF grammar, quoting rules, label restrictions, root branch length

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

### Olsen Example

```
(((One:0.2,Two:0.3):0.3,(Three:0.5,Four:0.3):0.2):0.3,Five:0.7):0.0;
```
Note: root has branch length `:0.0` — this is the TreeAlign convention documented by Olsen.

---

## Grammar Specification

### Wikipedia Grammar

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

### Olsen Grammar

```
tree ==> descendant_list [ root_label ] [ : branch_length ] ;
descendant_list ==> ( subtree { , subtree } )
subtree ==> descendant_list [internal_node_label] [: branch_length]
        ==> leaf_label [: branch_length]
label ==> unquoted_label | quoted_label
unquoted_label ==> string_of_printing_characters
quoted_label ==> ' string_of_printing_characters '
branch_length ==> signed_number | unsigned_number
```

Key difference: Olsen grammar explicitly allows root to have `[:branch_length]` and `[root_label]`.

### Label Restrictions (Olsen, Wikipedia Notes)

Unquoted labels may NOT contain: blanks, parentheses `()`, square brackets `[]`, single quotes `'`, colons `:`, semicolons `;`, or commas `,`.

Underscore `_` in unquoted labels is converted to blank (PHYLIP convention).

---

## Documented Corner Cases (from Sources)

### From Wikipedia

| Case | Description | Source |
|------|-------------|--------|
| Whitespace | Allowed anywhere except in unquoted names/numbers | Wikipedia Notes |
| Underscore → blank | Underscore in unquoted string becomes blank | Wikipedia Notes |
| Single quotes | String can be quoted with single quotes | Wikipedia Notes |
| Nested comments | Comments in `[]` allowed in some dialects | Wikipedia Notes |
| Unnamed nodes | `(,,(,));` - all nodes can be empty | Wikipedia Examples |
| Semicolon required | `Tree → Subtree ";"` — semicolon terminates | Wikipedia Grammar |

### From PHYLIP

| Case | Description | Source |
|------|-------------|--------|
| Multifurcation | `(A,B,C,D);` - more than 2 children | PHYLIP spec |
| Non-uniqueness | Multiple representations for same tree | PHYLIP Non-Uniqueness |
| Left-right order | Order of descendants is arbitrary | PHYLIP Non-Uniqueness |
| Single taxon | `A;` - valid single leaf tree | PHYLIP Examples |

### From Olsen

| Case | Description | Source |
|------|-------------|--------|
| Root branch length | `tree ==> ... [:branch_length] ;` — root can have length | Olsen Grammar |
| Root label | `tree ==> ... [root_label] ...` — root can have name | Olsen Grammar |

---

## Implementation Compliance

Based on current implementation (`PhylogeneticAnalyzer.cs`):

| Spec Requirement | Status | Notes |
|-----------------|--------|-------|
| Semicolon termination | ✅ Compliant | `ToNewick` appends `;` |
| Internal node names in output | ✅ Compliant | Emits valid unquoted labels per grammar |
| Invalid label suppression | ✅ Compliant | UPGMA/NJ names with metacharacters omitted per Olsen label rules |
| Branch length format | ✅ Compliant | Uses `.` decimal separator via InvariantCulture |
| Root branch length parsing | ✅ Compliant | Handles Olsen `[:branch_length]` after root subtree |
| Parsing: balanced parens | ✅ Compliant | Recursive descent matches `(` `)` |
| Parsing: InvariantCulture | ✅ Compliant | `double.TryParse` with InvariantCulture |
| Parsing: scientific notation | ✅ Compliant | Handles `e`, `E`, `+`, `-` in numbers |

---

## Invariants (Derived from Sources)

| ID | Invariant | Source |
|----|-----------|--------|
| N1 | Newick string MUST end with semicolon | Wikipedia Grammar: `Tree → Subtree ";"` |
| N2 | Leaf count in parsed tree MUST match leaf names in string | Wikipedia Grammar: leaves are `Name` productions |
| N3 | Round-trip (ToNewick → ParseNewick) MUST preserve topology | Wikipedia Grammar: unambiguous structure |
| N4 | Round-trip MUST preserve leaf names | Wikipedia Grammar: `Leaf → Name` |
| N5 | Branch lengths are real numbers after `:` | Wikipedia Grammar: `Length → ":" number` |
| N6 | ParseNewick on empty string MUST throw | Wikipedia Grammar: minimum structure required |
| N7 | ToNewick on null node SHOULD return empty | Defensive programming |
| N8 | Internal node names after `)` in output | Wikipedia Grammar: `Internal → "(" BranchSet ")" Name` |
| N9 | Number formatting uses `.` decimal separator | Olsen grammar: numbers are locale-independent |

---

## Documented Limitations (per Scope)

| Limitation | Spec Reference | Rationale |
|------------|---------------|-----------|
| Binary trees only | Wikipedia: `BranchSet → Branch \| Branch "," BranchSet` supports N children | UPGMA/NJ produce bifurcating trees |
| No quoted names | Olsen: `quoted_label ==> ' string '` | Out of scope |
| No `[]` comments | Wikipedia Notes: comments in square brackets | Out of scope |
| No underscore→blank | PHYLIP: underscore convention | Out of scope |
| F4 precision | Spec imposes no precision limit | Adequate for UPGMA/NJ; ±0.00005 |

---

## Quality Criteria

- [x] Sources are authoritative (Wikipedia, PHYLIP, Olsen)
- [x] Test cases traceable to documented examples
- [x] Corner cases from grammar specification
- [x] Implementation constraints documented as limitations with spec references
- [x] Invariants clearly defined with source citations
- [x] No assumptions — all decisions backed by external sources
- [x] InvariantCulture enforced for locale-independent number formatting

---

## References

1. Wikipedia contributors. "Newick format." Wikipedia, The Free Encyclopedia.
   https://en.wikipedia.org/wiki/Newick_format

2. Felsenstein, J. "The Newick tree format." PHYLIP documentation.
   https://phylipweb.github.io/phylip/newicktree.html

3. Olsen, G. "Interpretation of Newick's 8:45 Tree Format." 30 August 1990.
   https://phylipweb.github.io/phylip/newick_doc.html
