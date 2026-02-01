# Newick Format

## Overview

The Newick format (also known as New Hampshire format) is a standard notation for representing phylogenetic trees in text form. Named after Newick's restaurant in Dover, New Hampshire, where the format was finalized in 1986.

## Specification

### Grammar (from Gary Olsen)

```
Tree     → Subtree ";"
Subtree  → Leaf | Internal
Leaf     → Name
Internal → "(" BranchSet ")" Name
BranchSet → Branch | Branch "," BranchSet
Branch   → Subtree Length
Name     → empty | string
Length   → empty | ":" number
```

### Format Examples

| Pattern | Example | Description |
|---------|---------|-------------|
| Unnamed | `(,,(,));` | All nodes unnamed |
| Leaf names | `(A,B,(C,D));` | Only leaves named |
| All names | `(A,B,(C,D)E)F;` | Internal nodes also named |
| With lengths | `(A:0.1,B:0.2);` | Branch lengths included |
| Full | `(A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5)F;` | Names + lengths |

## Implementation

### Location

- **Class:** `PhylogeneticAnalyzer`
- **Namespace:** `Seqeron.Genomics`
- **File:** `PhylogeneticAnalyzer.cs`

### Methods

#### ToNewick

```csharp
public static string ToNewick(PhyloNode node, bool includeBranchLengths = true)
```

Converts a `PhyloNode` tree to Newick format string.

**Parameters:**
- `node`: Root node of the tree
- `includeBranchLengths`: Whether to include `:length` notation

**Returns:** Newick-formatted string ending with semicolon

**Algorithm:**
1. Recursive depth-first traversal
2. Leaf nodes → emit name
3. Internal nodes → emit `(left,right)` with optional lengths
4. Append semicolon at end

**Complexity:** O(n) where n = number of nodes

#### ParseNewick

```csharp
public static PhyloNode ParseNewick(string newick)
```

Parses a Newick format string into a tree structure.

**Parameters:**
- `newick`: Newick-formatted string

**Returns:** Root `PhyloNode` of the parsed tree

**Throws:** `ArgumentException` if string is empty or null

**Algorithm:**
1. Strip trailing semicolon
2. Recursive descent parsing
3. `(` starts internal node
4. `,` separates siblings
5. `:` precedes branch length
6. `)` closes internal node
7. Non-delimiter characters form names

**Complexity:** O(n) where n = string length

### PhyloNode Structure

```csharp
public class PhyloNode
{
    public string Name { get; set; }
    public double BranchLength { get; set; }
    public PhyloNode? Left { get; set; }
    public PhyloNode? Right { get; set; }
    public bool IsLeaf => Left == null && Right == null;
    public List<string> Taxa { get; set; }
}
```

## Invariants

| ID | Invariant | Verification |
|----|-----------|--------------|
| N1 | Output ends with semicolon | `newick.EndsWith(";")` |
| N2 | Leaf count preserved | Count leaves after parsing |
| N3 | Round-trip preserves topology | Compare leaf sets |
| N4 | Round-trip preserves names | Compare sorted name lists |
| N5 | Branch lengths ≥ 0 | Non-negative values |

## Limitations

The current implementation has these limitations relative to the full Newick specification:

1. **Binary trees only**: Does not parse multifurcating trees (3+ children)
2. **No quoted names**: Names with spaces or special characters not supported
3. **No comments**: `[comment]` syntax not parsed
4. **No underscore conversion**: Underscores not converted to spaces

## Test Coverage

See [PHYLO-NEWICK-001.md](../../TestSpecs/PHYLO-NEWICK-001.md) for test specification.

## References

1. Wikipedia. "Newick format." https://en.wikipedia.org/wiki/Newick_format
2. Felsenstein, J. "The Newick tree format." https://phylipweb.github.io/phylip/newicktree.html
3. Olsen, G. "Interpretation of Newick's 8:45 Tree Format." (1990)
