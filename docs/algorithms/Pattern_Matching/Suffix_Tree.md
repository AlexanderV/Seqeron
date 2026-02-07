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
| Build | O(n) | Ukkonen's online algorithm |
| Contains | O(m) | Walk edges; SIMD for ≥ 8-char edges |
| FindAllOccurrences | O(m + k) | k = number of matches |
| CountOccurrences | O(m) | Precomputed leaf counts |
| LongestRepeatedSubstring | O(1) | Precomputed during construction |
| LongestCommonSubstring | O(n + m) | Suffix-link streaming |
| FindAllLongestCommonSubstrings | O(m + k) | All occurrences of LCS |
| FindExactMatchAnchors | O(n + m) | Suffix-link streaming with peak tracking |
| EnumerateSuffixes | O(n²) total | Lazy DFS, O(n) per suffix |
| GetAllSuffixes | O(n²) | Materialized sorted list |

Where: n = text length, m = pattern/query length, k = result count.

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
  equals the LRS length (O(1) query). Persistent format stores the offset
  in the v5 header at byte 72.
- **Suffix link validation** — in DEBUG builds, every internal non-root node
  is checked for a valid suffix link.

The persistent builder additionally handles **compact→large transitions**
mid-build: when the next allocation would exceed `0xFFFFFFFE`, it switches
from 28-byte (compact) to 40-byte (large) nodes. A jump table bridges
cross-zone suffix links and child arrays (see
[Persistent README](../../src/SuffixTree/Algorithms/SuffixTree.Persistent/README.md)).

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

### 4.4 LongestRepeatedSubstring — O(1)

The deepest internal node (by `DepthFromRoot + edgeLength`) is found
during construction. Its total depth is the LRS length.
- In-memory: cached in private fields.
- Persistent: offset stored in v5 header byte 72; read on first access, then cached.

### 4.5 LongestCommonSubstring — O(n + m)

Uses `SuffixTreeAlgorithms.FindLongestCommonSubstring<TNode>`:

Walk the query string character-by-character against the tree. On mismatch,
follow suffix links and rescan to maintain the longest match. Track the
best match length and positions throughout.

Variants:
- `LongestCommonSubstring(other)` → string only
- `LongestCommonSubstringInfo(other)` → `(string, posInText, posInOther)`;
  returns `("", -1, -1)` if none
- `FindAllLongestCommonSubstrings(other)` → all occurrences of the LCS

### 4.6 FindExactMatchAnchors — O(n + m)

Uses `SuffixTreeAlgorithms.FindExactMatchAnchors<TNode>`:

Same suffix-link streaming as LCS, but with **peak tracking**. Emits
right-maximal matches when the match length drops below a minimum after
being above it. Produces non-overlapping anchors suitable for alignment
chaining (MUMmer/LAGAN-style MEMs).

### 4.7 Suffix Enumeration

DFS traversal in sorted child-key order (ascending). Concatenates edge
labels to produce suffixes.
- `EnumerateSuffixes()` — lazy `IEnumerable<string>`, avoids O(n²) peak memory.
- `GetAllSuffixes()` — materialized `List<string>`.

### 4.8 Traverse (Visitor Pattern)

`Traverse(ITreeVisitor)` — deterministic recursive DFS in sorted key order.
Calls `OnEnterNode`, `OnBeforeBranch`, `OnAfterBranch` for each node/edge.
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
IEnumerable<int> FindAllOccurrences(string / ReadOnlySpan<char>)
int CountOccurrences(string / ReadOnlySpan<char>)
string LongestRepeatedSubstring()
IEnumerable<string> EnumerateSuffixes()
List<string> GetAllSuffixes()
string LongestCommonSubstring(string / ReadOnlySpan<char>)
(string, int, int) LongestCommonSubstringInfo(string)
IEnumerable<(int, int, int)> FindAllLongestCommonSubstrings(string)
string PrintTree()
void Traverse(ITreeVisitor)
IEnumerable<(int, int, int)> FindExactMatchAnchors(string, int)
```

### ISuffixTreeNavigator\<TNode\>

```
ITextSource TextSource
TNode Root
TNode NullNode
bool IsNull(TNode)
bool IsRoot(TNode)
int GetEdgeChar(TNode, int)
int GetEdgeLength(TNode)
int GetDepth(TNode)
int GetDepthBeforeEdge(TNode)
TNode GetSuffixLink(TNode)
TNode GetChild(TNode, int)
IEnumerable<int> GetLeafPositions(TNode)
int GetAnyLeafPosition(TNode)
```

### ITextSource

```
int Length
char this[int]
string Substring(int, int)
ReadOnlySpan<char> AsSpan(int, int)
```

Implementations: `StringTextSource` (wraps `string`), `MemoryMappedTextSource`
(zero-copy pointer into MMF).

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
var positions     = tree.FindAllOccurrences("ana"); // IEnumerable<int>
int count         = tree.CountOccurrences("ana");

// Algorithms
string lrs        = tree.LongestRepeatedSubstring();
string lcs        = tree.LongestCommonSubstring("bandana");
var (s, p1, p2)   = tree.LongestCommonSubstringInfo("bandana");
var allLcs        = tree.FindAllLongestCommonSubstrings("bandana");
var anchors       = tree.FindExactMatchAnchors("bandana", minLength: 3);

// Enumeration
var suffixes      = tree.GetAllSuffixes();          // List<string>
var lazySuffixes  = tree.EnumerateSuffixes();        // IEnumerable<string>

// Diagnostics
string viz        = tree.PrintTree();
tree.Traverse(visitor);
```

### Persistent (MMF-Backed) Tree

```csharp
using SuffixTree.Persistent;

// Build into MMF file (hybrid v5)
using var tree = (IDisposable)PersistentSuffixTreeFactory.Create(
    new StringTextSource("banana"), "tree.dat");
var st = (ISuffixTree)tree;

// Build in heap memory (no file)
using var tree = (IDisposable)PersistentSuffixTreeFactory.Create(
    new StringTextSource("banana"));

// Load existing file (read-only, auto-detects v3/v4/v5)
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
See [MCP docs](../../docs/mcp/README.md) for connection and usage guides.

---

## 9. Tests

**599 tests total** across three test projects:

### In-Memory Tests — 236 tests

```
SuffixTree.Tests/
├── Core/                  (52 tests)
│   ├── BuildTests.cs          — 21 factory methods, input validation, TryBuild
│   ├── DiagnosticsTests.cs    — 5 PrintTree + ToString
│   ├── InvariantTests.cs      — 23 structural invariants (suffix existence, count parity)
│   └── StatisticsTests.cs     — 3 node/leaf count correctness
├── Search/                (44 tests)
│   ├── ContainsTests.cs       — 12 substring presence, span parity, case sensitivity
│   ├── CountOccurrencesTests.cs — 13 count correctness, overlapping, span parity
│   └── FindAllOccurrencesTests.cs — 19 position correctness, lazy evaluation, span parity
├── Algorithms/            (74 tests)
│   ├── BruteForceVerificationTests.cs — 15 LCS/LRS vs O(n²m)/O(n³) brute force
│   ├── LiteratureExamplesTests.cs     — 6 Ukkonen 1995, Fibonacci strings
│   ├── LongestCommonSubstringTests.cs — 28 LCS, FindAll, Info, span parity
│   ├── LongestRepeatedSubstringTests.cs — 13 LRS correctness + occurrence proof
│   └── SuffixEnumerationTests.cs      — 12 sorted enumeration, lazy, uniqueness
├── Regression/            (7 tests)
│   └── RegressionTests.cs     — Edge-offset reset, suffix-link traversal, stale offset, null char
├── Robustness/            (23 tests)
│   ├── EdgeCaseTests.cs       — 13 long patterns, single char, whitespace, palindromes
│   ├── StressTests.cs         — 5 100K-char builds, 10K queries, LCS stress
│   └── ThreadSafetyTests.cs   — 5 parallel reads, builds, buffer safety
└── Compatibility/         (36 tests)
    ├── BinaryDataTests.cs     — 20 null bytes, control chars, all 256 values
    └── UnicodeTests.cs        — 16 Cyrillic, Chinese, Greek, emoji, surrogates, ZWJ
```

### Persistent Tests — 339 tests

```
SuffixTree.Persistent.Tests/
├── Core/                  (37 tests)
│   ├── SuffixTreeTestBase.cs      — 8 base tests (inherited by 3 backends)
│   ├── HeapSuffixTreeTests.cs     — 8 heap backend
│   ├── LargeFormatSuffixTreeTests.cs — 8 large (v3) backend
│   ├── MmfSuffixTreeTests.cs      — 13 MMF backend + persistence, stress, concurrency
│   └── BuilderGuardTests.cs       — 8 double-build, deepest node, memory release
├── Format/                (109 tests)
│   ├── StorageFormatTests.cs          — 20 compact/large parity, auto-detection, constants
│   └── HybridTransitionZoneTests.cs   — 89 exhaustive sweep ×8 texts, cross-zone links
├── Parity/                (17 tests)
│   ├── ParityTests.cs            — 11 full API parity (10 strings + 10 random)
│   └── TopologyParityTests.cs    — 6 node-by-node topology comparison
├── Safety/                (78 tests)
│   ├── DisposeBehaviorTests.cs       — 12 ODE, idempotent, read-only, exception safety
│   ├── ConcurrencyTests.cs          — 5 parallel reads/dispose (1000–5000 threads)
│   ├── ProviderSafetyTests.cs       — 28 ODE, bounds checks, data preservation
│   ├── NodeApiSafetyTests.cs        — 14 TryGetChild, jumped nodes, API visibility
│   ├── NullGuardTests.cs            — 3 null guards
│   ├── MemoryMappedTextSourceTests.cs — 9 char access, dispose safety, TOCTOU fix
│   └── OverflowAndValidationTests.cs — 7 overflow, slice, chunked text
├── Serialization/         (14 tests)
│   ├── LogicalPersistenceTests.cs — 6 hash parity, round-trip, file save/load
│   ├── SerializerHashTests.cs    — 2 hash stability, byte-order regression
│   ├── ImportTruncationTests.cs  — 2 truncation detection
│   └── SetSizeAndExportTests.cs  — 4 dispose guard, negative size, chunked export
└── Validation/            (49 tests)
    ├── LoadValidationTests.cs        — 17 header validation (magic, version, bounds, SIZE)
    ├── MathematicalInvariantsTests.cs — 10 LeafCount = text.Length, NodeCount bounds
    ├── TreeContractTests.cs          — 12 lex order, depth semantics, visitor balance
    ├── EmptyPatternContractTests.cs  — 4 empty pattern contracts
    └── BoundaryAndEdgeCaseTests.cs   — 6 resizing, edge invariant, pathological patterns
```

### MCP Tests — 24 tests

```
SuffixTree.Mcp.Core.Tests/
├── SuffixTreeContainsTests.cs          — 2
├── SuffixTreeCountTests.cs             — 2
├── SuffixTreeFindAllTests.cs           — 2
├── SuffixTreeLcsTests.cs               — 2
├── SuffixTreeLrsTests.cs               — 2
├── SuffixTreeStatsTests.cs             — 2
├── FindLongestRepeatTests.cs           — 2
├── FindLongestCommonRegionTests.cs     — 2
├── CalculateSimilarityTests.cs         — 2
├── HammingDistanceTests.cs             — 2
├── EditDistanceTests.cs                — 2
└── CountApproximateOccurrencesTests.cs — 2
```

Each MCP tool has one input-validation test and one functional test.

---

## 10. Benchmarks

The `SuffixTree.Benchmarks` project (BenchmarkDotNet) measures 22 scenarios:

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

- Persistent format details: [SuffixTree.Persistent/README.md](../../src/SuffixTree/Algorithms/SuffixTree.Persistent/README.md)
- Exact pattern search: [Exact_Pattern_Search.md](Exact_Pattern_Search.md)
- MCP tool docs: [docs/mcp/](../../docs/mcp/README.md)
