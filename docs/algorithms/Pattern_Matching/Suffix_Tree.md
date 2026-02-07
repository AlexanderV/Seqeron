# Suffix Tree (Ukkonen)

**Algorithm Group:** Pattern Matching / String Indexing
**Implementation:** `SuffixTree` project (`SuffixTree/`)

---

## 1. Overview

A suffix tree is a compressed trie of all suffixes of a string. It supports fast substring search and is a core index structure used by several pattern matching tasks in this repository.

This implementation follows Ukkonen's online construction algorithm and is optimized for high performance in C#.

## 2. Complexity

| Operation | Complexity | Notes |
|----------|------------|-------|
| Build | O(n) | Ukkonen's algorithm |
| Contains | O(m) | Match along edges |
| Find all occurrences | O(m + k) | k = number of matches |
| Count occurrences | O(m) | Uses precomputed leaf counts |
| Longest repeated substring | O(1) | Precomputed during construction |
| Longest common substring | O(m) | Uses suffix links |
| Find exact match anchors | O(n + m) | Suffix-link streaming (MUMmer-style) |
| Find all longest common substrings | O(m + k) | All positions of LCS |

## 3. Implementation Highlights

- Online construction with active point and remainder tracking.
- Suffix links for amortized O(1) jumps between internal nodes.
- Suffix-link streaming for O(n+m) exact match anchor finding (MUMmer/LAGAN-style MEMs).
- Edge compression (tree, not trie) for memory efficiency.
- Internal terminator key to avoid conflicts with any input character.
- Read-only queries are thread-safe after build.
- Shared algorithms via `ISuffixTreeNavigator<TNode>` (generic specialization, zero overhead).
- Two implementations behind a single `ISuffixTree` interface:
  - **In-memory** (`SuffixTree.Build`) — heap-allocated, fastest for moderate data.
  - **Persistent** (`PersistentSuffixTreeFactory.Create`) — MMF-backed, scales to multi-GB text.
- **Adaptive storage format**: persistent trees automatically select a compact 32-bit offset layout
  for texts up to ~50 M characters (28-byte nodes, ~30% smaller files) and switch to 64-bit offsets
  for larger texts. Format is detected on load — no user configuration needed.
- Serialization via `SuffixTreeSerializer`: export (text + SHA256 hash), import (rebuild via Ukkonen).

## 4. API Sketch

### In-Memory Tree

```csharp
var tree = SuffixTree.Build("banana");
bool exists = tree.Contains("ana");
var positions = tree.FindAllOccurrences("ana");
int count = tree.CountOccurrences("ana");

string lrs = tree.LongestRepeatedSubstring();
string lcs = tree.LongestCommonSubstring("bandana");
var (sub, posText, posOther) = tree.LongestCommonSubstringInfo("bandana");

// O(n+m) exact match anchors via suffix-link streaming
var anchors = tree.FindExactMatchAnchors("bandana", minLength: 3);
foreach (var (posInText, posInQuery, length) in anchors)
    Console.WriteLine($"Match: text[{posInText}..] = query[{posInQuery}..] len={length}");
```

### Persistent (MMF-Backed) Tree

```csharp
using SuffixTree.Persistent;

// Build directly into a memory-mapped file
// Format is selected automatically: compact (32-bit) for texts ≤ 50M chars,
// large (64-bit) for bigger texts. Detected on load — no config needed.
using var tree = (IDisposable)PersistentSuffixTreeFactory.Create(
    new StringTextSource("banana"), "tree.dat");
var st = (ISuffixTree)tree;
st.Contains("ana"); // O(m), data stays on disk
st.FindExactMatchAnchors("bandana", 3); // suffix links preserved

// Load existing file (read-only, instant startup, auto-detects format)
using var loaded = (IDisposable)PersistentSuffixTreeFactory.Load("tree.dat");
```

### Serialization (Portable Format)

```csharp
using SuffixTree.Persistent;

// Export any ISuffixTree to a stream (text + SHA256 hash)
using var ms = new MemoryStream();
SuffixTreeSerializer.Export(tree, ms);

// Import: rebuilds via Ukkonen — full suffix link support
ms.Position = 0;
var imported = SuffixTreeSerializer.Import(ms, new HeapStorageProvider());

// Convenience: save to / load from MMF file
var saved = SuffixTreeSerializer.SaveToFile(tree, "portable.tree");
var loaded = SuffixTreeSerializer.LoadFromFile("portable.tree");
```

Key entry points:
- `SuffixTree.Build(...)`, `SuffixTree.TryBuild(...)`, `SuffixTree.Empty`
- `PersistentSuffixTreeFactory.Create(...)`, `PersistentSuffixTreeFactory.Load(...)`
- `SuffixTreeSerializer.Export`, `Import`, `SaveToFile`, `LoadFromFile`
- `ISuffixTree` — unified interface for both implementations
- Search: `Contains`, `FindAllOccurrences`, `CountOccurrences`
- Algorithms: `LongestRepeatedSubstring`, `LongestCommonSubstring`, `LongestCommonSubstringInfo`, `FindAllLongestCommonSubstrings`, `FindExactMatchAnchors`
- Diagnostics: `PrintTree`, `Traverse`

## 5. Related Documentation

- Exact matching with suffix tree: `docs/algorithms/Pattern_Matching/Exact_Pattern_Search.md`
- Persistent suffix tree (MMF): `src/SuffixTree/Algorithms/SuffixTree.Persistent/README.md`
- MCP tool docs: `docs/mcp/tools/core/suffix_tree_*`

## 6. Benchmarks and Tests

- Performance benchmarks: `SuffixTree.Benchmarks/`
- Stress harness: `SuffixTree.Console/`
- In-memory tests: `SuffixTree.Tests/` (245 tests)
- Persistent tests: `SuffixTree.Persistent.Tests/` (157 tests — parity, serialization, MMF, format parity, auto-detection)

## 7. References

- Ukkonen, E. (1995). On-line construction of suffix trees. Algorithmica.
- Gusfield, D. (1997). Algorithms on Strings, Trees, and Sequences.
- Delcher, A. et al. (1999). Alignment of whole genomes (MUMmer — suffix tree anchor approach).
- https://visualgo.net/en/suffixtree
