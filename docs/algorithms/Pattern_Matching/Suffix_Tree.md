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

## 3. Implementation Highlights

- Online construction with active point and remainder tracking.
- Suffix links for amortized O(1) jumps between internal nodes.
- Edge compression (tree, not trie) for memory efficiency.
- Internal terminator key to avoid conflicts with any input character.
- Read-only queries are thread-safe after build.

## 4. API Sketch

```csharp
var tree = SuffixTree.Build("banana");
bool exists = tree.Contains("ana");
var positions = tree.FindAllOccurrences("ana");
int count = tree.CountOccurrences("ana");

string lrs = tree.LongestRepeatedSubstring();
string lcs = tree.LongestCommonSubstring("bandana");

foreach (var suffix in tree.EnumerateSuffixes())
{
    Console.WriteLine(suffix);
}
```

Key entry points:
- `SuffixTree.Build(...)`, `SuffixTree.TryBuild(...)`, `SuffixTree.Empty`
- `ISuffixTree` for abstraction
- Search: `Contains`, `FindAllOccurrences`, `CountOccurrences`
- Algorithms: `LongestRepeatedSubstring`, `LongestCommonSubstring`, `FindAllLongestCommonSubstrings`
- Diagnostics: `PrintTree`

## 5. Related Documentation

- Exact matching with suffix tree: `docs/algorithms/Pattern_Matching/Exact_Pattern_Search.md`
- MCP tool docs: `docs/mcp/tools/core/suffix_tree_*`

## 6. Benchmarks and Tests

- Performance benchmarks: `SuffixTree.Benchmarks/`
- Stress harness: `SuffixTree.Console/`
- Unit tests: `SuffixTree.Tests/`

## 7. References

- Ukkonen, E. (1995). On-line construction of suffix trees. Algorithmica.
- Gusfield, D. (1997). Algorithms on Strings, Trees, and Sequences.
- https://visualgo.net/en/suffixtree
