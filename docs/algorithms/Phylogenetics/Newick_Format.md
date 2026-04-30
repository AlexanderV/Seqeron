# Newick Format

| Field | Value |
|-------|-------|
| Algorithm Group | Phylogenetics |
| Test Unit ID | PHYLO-NEWICK-001 |
| Related Projects | N/A |
| Implementation Status | Simplified |
| Last Reviewed | 2026-04-30 |

## 1. Overview

The Newick format is a compact text representation for phylogenetic trees. In this repository it appears as both a serializer (`ToNewick`) and a parser (`ParseNewick`) for the `PhyloNode` tree type. The implementation supports the core parenthesized tree grammar, internal names when they are valid unquoted labels, optional branch lengths, optional root branch lengths during parsing, and invariant-culture numeric formatting. It is simplified relative to the full ecosystem of Newick dialects because it only handles binary trees and omits quoted labels and comments.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

Newick, also called the New Hampshire format, is a standard text serialization for rooted or unrooted phylogenetic trees. Its governing logic is grammatical rather than biological: a tree is represented as recursively nested subtrees, optionally annotated with labels and branch lengths.

### 2.2 Core Model

The repository follows the basic Olsen-style grammar summarized in the original document:

```text
Tree      -> Subtree ";"
Subtree   -> Leaf | Internal
Leaf      -> Name
Internal  -> "(" BranchSet ")" Name
BranchSet -> Branch | Branch "," BranchSet
Branch    -> Subtree Length
Name      -> empty | string
Length    -> empty | ":" number
```

Serializer and parser both rely on the same core constructs: parentheses for internal nodes, commas for sibling separation, optional labels after a subtree, and branch lengths introduced by `:`.

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | `ToNewick` output ends with `;` | The serializer appends a semicolon after recursive traversal |
| INV-02 | Branch lengths are rendered with `.` as decimal separator | The serializer uses `CultureInfo.InvariantCulture` |
| INV-03 | Internal node names are emitted only when they are valid unquoted Newick labels | The serializer suppresses names containing Newick metacharacters |
| INV-04 | `ParseNewick` accepts an optional root branch length after the main subtree | The parser checks for `:` after parsing the root subtree |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `[ToNewick] node` | `PhyloNode` | required | Root node of the tree to serialize | `null` returns an empty string |
| `[ToNewick] includeBranchLengths` | `bool` | `true` | Whether child branch lengths are emitted | Applies to non-root edges |
| `[ParseNewick] newick` | `string` | required | Newick-formatted tree string | Null, empty, or whitespace-only input throws `ArgumentException` |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `[ToNewick] return value` | `string` | Serialized Newick string ending with `;` |
| `[ParseNewick] return value` | `PhyloNode` | Parsed binary tree rooted at a `PhyloNode` |

### 3.3 Preconditions and Validation

The parser trims the input string and strips a trailing semicolon if present. It then performs recursive descent over parentheses, commas, labels, and branch lengths. Only binary branch sets are supported by the current tree model because `PhyloNode` exposes only `Left` and `Right` children. Serializer-side internal-node label emission is restricted to valid unquoted labels; internal labels containing spaces or Newick metacharacters are omitted rather than quoted. `ParseNewick(...)` itself is more permissive and reads labels as raw runs up to the next structural delimiter rather than enforcing Olsen-style unquoted-label restrictions.

## 4. Algorithm

### 4.1 High-Level Steps

1. For serialization, return an empty string when the root node is `null`.
2. Traverse the tree recursively in depth-first order.
3. Emit leaf names directly.
4. Emit internal nodes as `(left,right)` with optional child branch lengths.
5. Append an internal node name only when it is a valid unquoted label.
6. Append the final semicolon.
7. For parsing, recursively descend through `(`, `,`, `)`, labels, and numeric branch lengths.
8. Optionally parse a root branch length after the main subtree.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

The serializer formats branch lengths with `ToString("F4", CultureInfo.InvariantCulture)`. The parser accepts digits, decimal points, signs, and scientific notation markers (`e`, `E`) inside numeric fields. Valid unquoted labels exclude blanks, parentheses, square brackets, single quotes, colons, semicolons, and commas.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `ToNewick` | `O(n)` | `O(h)` | `n` = number of nodes, `h` = recursion depth |
| `ParseNewick` | `O(n)` | `O(h)` | `n` = input length, `h` = parse-tree depth |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [PhylogeneticAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Phylogenetics/PhylogeneticAnalyzer.cs)

- `PhylogeneticAnalyzer.ToNewick(PhyloNode, bool)`: Serializes a `PhyloNode` tree to Newick text.
- `PhylogeneticAnalyzer.ParseNewick(string)`: Parses a Newick string into a `PhyloNode` tree.

### 5.2 Current Behavior

The current serializer always emits binary trees because the `PhyloNode` type contains only `Left` and `Right` child slots. Internal node names containing Newick metacharacters are silently suppressed rather than quoted, but leaf labels are emitted verbatim because the serializer's unquoted-label check is applied only to internal-node names. The parser accepts optional root branch lengths and scientific-notation numbers, and it parses labels as raw character runs until a Newick delimiter is encountered.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Parenthesized recursive tree structure with commas between siblings.
- Optional branch lengths using `:`.
- Optional internal node names when those names are valid unquoted labels.

**Intentionally simplified:**

- Only binary trees are supported; **consequence:** multifurcating `BranchSet` forms from the broader Newick grammar are out of scope for this parser and serializer.
- Labels are never quoted; **consequence:** internal labels requiring quoting are omitted on serialization, while parsing remains permissive and does not enforce the same unquoted-label restrictions.

**Not implemented:**

- Square-bracket comments and underscore-to-blank translation; **users should rely on:** no current alternative in this class.

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| `ToNewick(null)` | Returns an empty string | Explicit null branch in the serializer |
| Empty or whitespace-only parse input | Throws `ArgumentException` | Parser rejects missing tree text |
| Root branch length after the subtree | Parsed into `root.BranchLength` | Dedicated post-root `:` handling in `ParseNewick` |
| Scientific notation in branch lengths | Accepted | `ParseNumber` allows `e`, `E`, `+`, and `-` |

### 6.2 Limitations

This implementation omits several Newick dialect features: multifurcations, quoted labels, comments, and underscore-to-blank conversion. It is also asymmetric about labels: internal names are filtered to a conservative unquoted subset during serialization, while parsing accepts permissive raw label tokens up to the next delimiter. It is appropriate for the binary trees produced by the repository's UPGMA and Neighbor-Joining builders, but not for general-purpose Newick interoperability across all tools.

## 7. Examples and Related Material

- [PHYLO-NEWICK-001](../../../tests/TestSpecs/PHYLO-NEWICK-001.md) documents the repository's Newick-format test specification.

## 8. References

1. Wikipedia contributors. Newick format. Wikipedia. https://en.wikipedia.org/wiki/Newick_format
2. Felsenstein, J. The Newick tree format. https://phylipweb.github.io/phylip/newicktree.html
3. Olsen, G. 1990. Interpretation of Newick's 8:45 Tree Format.

