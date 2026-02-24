# Suffix Tree (Ukkonen)

**Algorithm Group:** Pattern Matching / String Indexing
**Projects:** `SuffixTree.Core`, `SuffixTree`, `SuffixTree.Persistent`

---

## 1. Overview

A suffix tree is a compressed trie of all suffixes of a string. It supports
O(m) substring search and serves as the core index structure for pattern
matching, repeat detection, and sequence comparison in this repository.

The implementation follows Ukkonen's online construction algorithm (1995).
Two interchangeable backends share a single `ISuffixTree` interface:

| Backend | Project | Target | Backing |
|---------|---------|--------|---------|
| **In-memory** | `SuffixTree` | .NET 8 | Heap-allocated node objects |
| **Persistent** | `SuffixTree.Persistent` | .NET 9 | Memory-mapped files (MMF) |

Shared algorithms are written once against `ISuffixTreeNavigator<TNode>`
with a `struct` constraint, enabling JIT specialization for each backend.

---

## 2. Complexity

| Operation | Time | Notes |
|-----------|------|-------|
| Build | O(n) amortized | Ukkonen's online algorithm |
| Contains | O(m) in-memory; O(m log d) persistent | Child lookup is hash/inline vs binary search |
| FindAllOccurrences | O(m + k) in-memory; O(m log d + k) persistent | k = number of returned positions |
| CountOccurrences | O(m) in-memory; O(m log d) persistent | Leaf count is precomputed; no subtree DFS |
| LongestRepeatedSubstring | First call: O(\|LRS\|) in-memory; persistent O(h + \|LRS\|) typical, O(n + \|LRS\|) fallback; then O(1) cached | Persistent can fall back to deepest-node DFS if header metadata is unavailable |
| LongestCommonSubstring | O(m + h) in-memory; O(m log d + h) persistent | Streaming match + one leaf-position recovery |
| FindAllLongestCommonSubstrings | O(m + Σ subtree(best matches)); worst case O(n·m) | Collects leaves for each maximal match candidate |
| FindExactMatchAnchors | O(m + a·h); worst case O(n·m) | a = anchors emitted, each needs leaf-position recovery |
| EnumerateSuffixes | O(n²) total | Lazy DFS, O(n) per suffix |
| GetAllSuffixes | O(n²) | Materialized sorted list |

Where: n = text length, m = pattern/query length, k = result count, d = max branching factor on visited nodes, h = max depth walked to recover a leaf position, a = anchor count.

---

## 3. Construction (Ukkonen's Algorithm)

Both backends use the same logic:

1. Process characters left to right, extending the tree with each character.
2. Maintain an **active point** (active node, active edge, active length)
   and a **remainder** counter.
3. Three extension rules per phase:
   - **Rule 1:** Implicit leaf extension (open-ended leaf edges grow automatically).
   - **Rule 2:** New leaf creation + optional edge split.
   - **Rule 3:** Character already present — showstopper.
4. **Suffix links** connect node for "xα" to node for "α", ensuring amortized
   O(1) jumps during construction.
5. A special **terminator** (integer key `−1`) is appended to force all suffixes
   to be explicit leaves.

Post-construction passes (both backends):
- **Bottom-up leaf counting** — iterative post-order traversal assigns
  `LeafCount` to every node. `CountOccurrences` reads this in O(1).
- **Deepest internal node** — identified during the same pass. Its depth
  equals the LRS length (O(1) query). Persistent format stores deepest-node
  offset and precomputed LRS depth in the v6 header (bytes 72 and 80).
- **Suffix link validation** — in DEBUG builds, the in-memory backend checks
  that internal non-root nodes have valid suffix links.

The persistent builder additionally handles **compact→large transitions**
mid-build: when the next allocation would exceed `0xFFFFFFFE`, it switches
from 24-byte (compact) to 32-byte (large) nodes. A jump table bridges
cross-zone suffix links and child arrays (see
[Persistent README](../../../src/SuffixTree/Algorithms/SuffixTree.Persistent/README.md)).

---

## 4. Algorithms

### 4.1 Contains — O(m)

Walk from root, matching pattern characters against edge labels.
The in-memory implementation uses a hybrid strategy:
`SequenceEqual` (SIMD) for edges ≥ 8 characters, scalar loop otherwise.
Zero allocations via `ReadOnlySpan<char>`.

### 4.2 FindAllOccurrences — O(m + k)

Match pattern to find a terminal node, then iterative stack-based DFS
collects all leaf positions. Each leaf's position =
`textLength − node.DepthFromRoot − edgeLength`.

### 4.3 CountOccurrences — O(m)

Match pattern, then return `node.LeafCount` — precomputed during
construction by the bottom-up pass.

### 4.4 LongestRepeatedSubstring — first call non-constant, then O(1) cached

The deepest internal node (by `DepthFromRoot + edgeLength`) is found
during construction. Its total depth is the LRS length.
- In-memory: cached in private fields.
- Persistent: deepest-node offset + LRS depth are stored in the v6 header (72/80), then cached;
  if metadata is unavailable, a fallback DFS is used once.

### 4.5 LongestCommonSubstring — O(m + h) in-memory; O(m log d + h) persistent

Uses `SuffixTreeAlgorithms.FindAllLcs<TNode, TNav>`:

Walk the query string character-by-character against the tree. On mismatch,
follow suffix links and rescan to maintain the longest match. Track the
best match length and positions throughout.

Variants:
- `LongestCommonSubstring(other)` → string only
- `LongestCommonSubstringInfo(other)` → `(string, posInText, posInOther)`;
  returns `("", -1, -1)` if none
- `FindAllLongestCommonSubstrings(other)` → canonical LCS string + all positions in both strings

### 4.6 FindExactMatchAnchors — O(m + a·h), worst case O(n·m)

Uses `SuffixTreeAlgorithms.FindExactMatchAnchors<TNode>`:

Same suffix-link streaming as LCS, but with **peak tracking**. Emits
right-maximal matches when the match length drops below a minimum after
being above it. Produces non-overlapping anchors suitable for alignment
chaining (MUMmer/LAGAN-style MEMs).

### 4.7 Suffix Enumeration

DFS traversal in sorted child-key order (ascending). Concatenates edge
labels to produce suffixes.
- `EnumerateSuffixes()` — lazy `IEnumerable<string>`, avoids O(n²) peak memory.
- `GetAllSuffixes()` — materialized `IReadOnlyList<string>`.

### 4.8 Traverse (Visitor Pattern)

`Traverse(ISuffixTreeVisitor)` — deterministic DFS in sorted key order.
Calls `VisitNode`, `EnterBranch`, `ExitBranch` for each node/edge.
Used by `SuffixTreeSerializer` for structural hashing.

---

## 5. Interface Hierarchy

### ISuffixTree

```
ITextSource Text
int NodeCount
int LeafCount
int MaxDepth
bool IsEmpty

bool Contains(string / ReadOnlySpan<char>)
IReadOnlyList<int> FindAllOccurrences(string / ReadOnlySpan<char>)
int CountOccurrences(string / ReadOnlySpan<char>)
string LongestRepeatedSubstring()
ReadOnlyMemory<char> LongestRepeatedSubstringMemory()
IEnumerable<string> EnumerateSuffixes()
IReadOnlyList<string> GetAllSuffixes()
string LongestCommonSubstring(string / ReadOnlySpan<char>)
(string, int, int) LongestCommonSubstringInfo(string)
(string, IReadOnlyList<int>, IReadOnlyList<int>) FindAllLongestCommonSubstrings(string)
string PrintTree()
void Traverse(ISuffixTreeVisitor)
IReadOnlyList<(int, int, int)> FindExactMatchAnchors(string, int)
```

### ISuffixTreeNavigator\<TNode\>

```
ITextSource Text
TNode Root
TNode NullNode
bool IsNull(TNode)
bool IsRoot(TNode)
int GetEdgeSymbol(TNode, int)
int LengthOf(TNode)
TNode GetSuffixLink(TNode)
bool TryGetChild(TNode, int, out TNode)
void CollectLeaves(TNode, int, List<int>)
int FindAnyLeafPosition(TNode, int)
```

### ITextSource

```
int Length
char this[int]
string Substring(int, int)
ReadOnlySpan<char> Slice(int, int)
```

Implementations: `StringTextSource` (wraps `string`),
`MemoryMappedTextSource` and `AsciiMemoryMappedTextSource` (MMF-backed).

---

## 6. API Reference

### In-Memory Tree

```csharp
// Build
var tree = SuffixTree.Build("banana");             // from string
var tree = SuffixTree.Build(textSource);            // from ITextSource
var tree = SuffixTree.Build(readOnlyMemory);        // from ReadOnlyMemory<char>
var tree = SuffixTree.Build(readOnlySpan);          // from ReadOnlySpan<char>
bool ok  = SuffixTree.TryBuild("text", out var t);  // non-throwing
var empty = SuffixTree.Empty;                       // singleton empty tree

// Search
bool found        = tree.Contains("ana");
var positions     = tree.FindAllOccurrences("ana"); // IReadOnlyList<int>
int count         = tree.CountOccurrences("ana");

// Algorithms
string lrs        = tree.LongestRepeatedSubstring();
string lcs        = tree.LongestCommonSubstring("bandana");
var (s, p1, p2)   = tree.LongestCommonSubstringInfo("bandana");
var allLcs        = tree.FindAllLongestCommonSubstrings("bandana");
var anchors       = tree.FindExactMatchAnchors("bandana", minLength: 3);

// Enumeration
var suffixes      = tree.GetAllSuffixes();          // IReadOnlyList<string>
var lazySuffixes  = tree.EnumerateSuffixes();        // IEnumerable<string>

// Diagnostics
string viz        = tree.PrintTree();
tree.Traverse(visitor);
```

### Persistent (MMF-Backed) Tree

```csharp
using SuffixTree.Persistent;

// Build into MMF file (hybrid v6)
using var tree = (IDisposable)PersistentSuffixTreeFactory.Create(
    new StringTextSource("banana"), "tree.dat");
var st = (ISuffixTree)tree;

// Build in heap memory (no file)
using var tree = (IDisposable)PersistentSuffixTreeFactory.Create(
    new StringTextSource("banana"));

// Load existing file (read-only, v6 with compact/large base auto-detection)
using var loaded = (IDisposable)PersistentSuffixTreeFactory.Load("tree.dat");

// All ISuffixTree methods work identically
st.Contains("ana");
st.FindExactMatchAnchors("bandana", 3);
```

### Serialization

```csharp
// Deterministic structural hash (SHA256)
byte[] hash = SuffixTreeSerializer.CalculateLogicalHash(tree);

// Stream export/import (format v2)
SuffixTreeSerializer.Export(tree, stream);
var imported = SuffixTreeSerializer.Import(stream, new HeapStorageProvider());

// File-based export/import (MMF)
using var saved  = (IDisposable)SuffixTreeSerializer.SaveToFile(tree, "saved.tree");
using var loaded = (IDisposable)SuffixTreeSerializer.LoadFromFile("saved.tree");
```

---

## 7. Project Structure

```
src/SuffixTree/Algorithms/
├── SuffixTree.Core/              ← Interfaces + shared algorithms (.NET 8)
│   ├── ISuffixTree.cs
│   ├── ISuffixTreeNavigator.cs
│   ├── ITextSource.cs
│   ├── StringTextSource.cs
│   └── SuffixTreeAlgorithms.cs   ← LCS + Anchors (generic, JIT-specialized)
├── SuffixTree/                   ← In-memory implementation (.NET 8)
│   ├── SuffixTree.cs             ← Build, factory methods
│   ├── SuffixTree.Construction.cs ← Ukkonen's algorithm
│   ├── SuffixTree.Search.cs      ← Contains, FindAll, Count
│   ├── SuffixTree.Algorithms.cs  ← LRS, LCS, Anchors
│   ├── SuffixTree.Navigator.cs   ← ISuffixTreeNavigator<SuffixTreeNode>
│   ├── SuffixTree.Diagnostics.cs ← PrintTree, Traverse, ComputeStatistics
│   └── SuffixTreeNode.cs         ← Hybrid children storage (inline ≤ 4 → Dictionary)
└── SuffixTree.Persistent/       ← Disk-backed implementation (.NET 9)
    ├── PersistentSuffixTree.cs
    ├── PersistentSuffixTreeBuilder.cs
    ├── PersistentSuffixTreeFactory.cs
    ├── PersistentSuffixTreeNavigator.cs
    ├── PersistentSuffixTreeNode.cs
    ├── IStorageProvider.cs
    ├── HeapStorageProvider.cs
    ├── MappedFileStorageProvider.cs
    ├── AsciiMemoryMappedTextSource.cs
    ├── MemoryMappedTextSource.cs
    ├── NodeLayout.cs
    ├── HybridLayout.cs
    ├── StorageFormat.cs
    ├── PersistentConstants.cs
    └── SuffixTreeSerializer.cs
```

### Dependencies

```
SuffixTree.Core  ←──  SuffixTree (in-memory)
       ↑
       └──────────  SuffixTree.Persistent
```

---

## 8. MCP Integration

The `SuffixTree.Mcp.Core` project exposes 12 tools via Model Context Protocol:

| Tool | Method |
|------|--------|
| `suffix_tree_contains` | Pattern existence check |
| `suffix_tree_count` | Occurrence count |
| `suffix_tree_find_all` | All positions |
| `suffix_tree_lrs` | Longest repeated substring |
| `suffix_tree_lcs` | Longest common substring |
| `suffix_tree_stats` | Tree statistics |
| `find_longest_repeat` | DNA longest tandem repeat |
| `find_longest_common_region` | DNA common region |
| `calculate_similarity` | K-mer Jaccard similarity |
| `hamming_distance` | Hamming distance |
| `edit_distance` | Levenshtein edit distance |
| `count_approximate_occurrences` | Approximate pattern matching |

Each tool validates input and returns a structured result record.
See [MCP docs](../../mcp/README.md) for connection and usage guides.

---

## 9. Tests

As of 2026-02-24 (`dotnet test` on the three SuffixTree test projects):

- `SuffixTree.Tests`: 353 passed
- `SuffixTree.Persistent.Tests`: 471 passed
- `SuffixTree.Mcp.Core.Tests`: 59 passed
- **Total: 883 passed**

For detailed, maintained coverage matrix by suite/scenario, see
[`tests/SuffixTree/SUFFIX_TREE_TEST_MATRIX.md`](../../../tests/SuffixTree/SUFFIX_TREE_TEST_MATRIX.md).

---

## 10. Benchmarks

The `SuffixTree.Benchmarks` project (BenchmarkDotNet) includes scenarios for:

| Category | Benchmarks |
|----------|------------|
| **Build** | DNA 10K, DNA 100K, random 10K, random 100K |
| **Contains** | Found, not found, long pattern, very long pattern, single char |
| **FindAll** | Single, multiple, many, not found |
| **Count** | Single, multiple, many |
| **LRS** | Short, medium, long, repetitive |
| **LCS** | Short, medium, long, no match |
| **Hairpin** | Build, search, full pipeline |

Run with:

```bash
cd apps/SuffixTree.Benchmarks
dotnet run -c Release
```

The `SuffixTree.Console` app provides an exhaustive stress harness:
three phases testing small strings (all substrings), large strings
(random + all suffixes), and 13 edge cases.

---

## 11. References

- Ukkonen, E. (1995). *On-line construction of suffix trees.* Algorithmica, 14(3), 249–260.
- Gusfield, D. (1997). *Algorithms on Strings, Trees, and Sequences.* Cambridge University Press.
- Delcher, A. et al. (1999). *Alignment of whole genomes.* Nucleic Acids Research (MUMmer — suffix tree anchor approach).
- https://visualgo.net/en/suffixtree

---

## 12. Related Documentation

- Persistent format details: [SuffixTree.Persistent/README.md](../../../src/SuffixTree/Algorithms/SuffixTree.Persistent/README.md)
- Exact pattern search: [Exact_Pattern_Search.md](Exact_Pattern_Search.md)
- MCP tool docs: [docs/mcp/](../../mcp/README.md)
