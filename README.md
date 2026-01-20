# SuffixTree

A high-performance implementation of **Ukkonen's suffix tree algorithm** in C#.

[![.NET](https://img.shields.io/badge/.NET-8.0-512BD4)](https://dotnet.microsoft.com/)
[![Tests](https://img.shields.io/badge/tests-201%20passed-brightgreen)]()
[![License](https://img.shields.io/badge/license-MIT-blue)](LICENSE)

## Overview

A **suffix tree** is a compressed trie containing all suffixes of a given string. It is one of the most powerful data structures in stringology, enabling efficient solutions to numerous string processing problems.

This implementation uses **Ukkonen's online algorithm** (1995), which constructs the tree in **O(n)** time and allows **O(m)** substring search, where *n* is the length of the text and *m* is the length of the pattern.

### Why Suffix Trees?

Suffix trees provide optimal time complexity for many string operations:

| Operation | Time Complexity | Notes |
|-----------|---------------|-------|
| Build tree | O(n) | Ukkonen's algorithm |
| Substring search | O(m) | Character-by-character matching |
| Find all occurrences | O(m + k) | m = pattern, k = occurrences |
| Count occurrences | O(m) | Uses precomputed leaf counts |
| Longest repeated substring | O(1)* | *O(n) on first call, cached |
| Longest common substring | O(m) | Uses suffix links for O(m) traversal |

*where n = text length, m = pattern length, k = number of occurrences*

## Features

- **O(n) construction** — linear time tree building using Ukkonen's algorithm
- **O(m) substring search** — efficient pattern matching
- **Full Unicode support** — works with any characters (Cyrillic, Chinese, emoji, etc.)
- **Binary data support** — handles null characters and control characters
- **Terminator key** — uses internal key (-1) to avoid conflicts with any input character
- **Suffix links** — optimized traversal during construction and LCS queries
- **Thread-safe reads** — built tree can be queried from multiple threads
- **Memory efficient** — uses edge compression (not a trie)
- **.NET 8.0** — modern SDK-style project

## Installation

Clone the repository and build:

```bash
git clone https://github.com/your-username/SuffixTree.git
cd SuffixTree
dotnet build
```

Or add the project reference to your solution.

## Usage

```csharp
// Build a suffix tree
var tree = SuffixTree.Build("banana");

// Access original text
string text = tree.Text;  // "banana"

// Tree statistics
int nodes = tree.NodeCount;    // Total nodes in tree
int leaves = tree.LeafCount;   // Number of leaf nodes
int depth = tree.MaxDepth;     // Maximum depth in characters

// Check if substring exists
bool found = tree.Contains("nan");  // true
bool notFound = tree.Contains("xyz");  // false

// Zero-allocation search with Span
bool exists = tree.Contains("ana".AsSpan());

// Find all occurrences
var positions = tree.FindAllOccurrences("ana");  // [1, 3]
int count = tree.CountOccurrences("ana");  // 2

// Longest repeated substring
string lrs = tree.LongestRepeatedSubstring();  // "ana"

// Longest common substring with positions
var (substring, posInText, posInOther) = tree.LongestCommonSubstringInfo("bandana");

// Find ALL positions of longest common substring
var (lcs, positionsInText, positionsInOther) = tree.FindAllLongestCommonSubstrings("xabcyabcz");

// Get all suffixes (sorted)
var suffixes = tree.GetAllSuffixes();

// Lazy enumeration for large strings
foreach (var suffix in tree.EnumerateSuffixes())
{
    Console.WriteLine(suffix);
}

// Debug: print tree structure
Console.WriteLine(tree.PrintTree());
```

## Project Structure

```
SuffixTree.sln
├── SuffixTree/              # Core library
│   ├── SuffixTree.cs        # Ukkonen's algorithm implementation
│   ├── SuffixTreeNode.cs    # Internal node class
│   └── ISuffixTree.cs       # Public interface
├── SuffixTree.Console/      # Stress test console app
│   └── Program.cs           # Exhaustive verification tests
├── SuffixTree.Tests/        # Unit tests (NUnit) - 201 tests
│   ├── Algorithms/          # LCS, LRS, suffix enumeration (62 tests)
│   ├── Compatibility/       # Unicode, binary data (36 tests)
│   ├── Core/                # Build, statistics, diagnostics (21 tests)
│   ├── Regression/          # Bug prevention tests (7 tests)
│   ├── Robustness/          # Edge cases, stress, threading (26 tests)
│   └── Search/              # Contains, Count, FindAll (45 tests)
└── SuffixTree.Benchmarks/   # Performance benchmarks
    └── SuffixTreeBenchmarks.cs
```

## Benchmarks

Performance measured on Intel Core i7-1185G7 @ 3.00GHz, .NET 8.0:

### Build Performance

| Text Size | Time | Allocated |
|-----------|------|-----------|
| 100 chars | 7 μs | 19 KB |
| 10K chars | 1.36 ms | 2 MB |
| 100K chars | 52 ms | 20 MB |
| 50K DNA | 20 ms | 8 MB |

### Search Performance

| Operation | 100 chars | 10K chars | 100K chars | 50K DNA |
|-----------|-----------|-----------|------------|---------|
| Contains | 18 ns | 25 ns | 33 ns | 37 ns |
| Count | — | 16 ns | 27 ns | 21 ns |
| FindAll | 50 ns | 1.4 μs | 25 μs | 3 μs |
| LRS | 8 ns | 877 ns | 161 μs | 9 ns |
| LCS | 183 ns | 476 ns | 676 ns | 980 ns |

*LRS = Longest Repeated Substring, LCS = Longest Common Substring*

### Key Observations

- **Build time scales linearly** with text size (O(n) confirmed)
- **Contains is extremely fast** — 18-37 ns regardless of tree size (hybrid SIMD optimization)
- **Count uses precomputed leaf counts** — O(m) instead of O(m+k), ~1,800x faster on large texts
- **LRS uses depth tracking** — O(1) after first call, 59x-26,000x faster than tree traversal
- **LCS uses suffix links** — O(m) traversal, 100x+ faster than naive approach
- **Zero allocations** for Contains and Count operations
- **DNA sequences** (small alphabet) build faster due to better tree compression

Run benchmarks yourself:
```bash
cd SuffixTree.Benchmarks
dotnet run -c Release
```

## Algorithm Details

The implementation follows Ukkonen's canonical algorithm with:

1. **Active Point** — tracks current position in tree (`activeNode`, `activeEdge`, `activeLength`)
2. **Suffix Links** — enable O(1) jumps between internal nodes
3. **Remainder** — counts pending suffixes to be inserted
4. **Rule 1** — no matching edge exists → create new leaf
5. **Rule 2** — mismatch in middle of edge → split edge, create internal node
6. **Rule 3** — character matches → suffix exists implicitly (showstopper)
7. **Walk-Down** — traverse edges when `activeLength >= edgeLength`
8. **Terminator** — internal key (`-1`) ensures all suffixes end at leaves without conflicting with any input character

### Understanding Suffix Trees

A suffix tree for string S is a tree where:
- Each **edge** is labeled with a non-empty substring of S
- Each **internal node** (except root) has at least 2 children
- No two edges from the same node start with the same character
- Each **leaf** represents a unique suffix of S
- Concatenating edge labels from root to leaf gives the suffix

**Example:** For string "banana", the suffixes are:
```
banana  (position 0)
anana   (position 1)
nana    (position 2)
ana     (position 3)
na      (position 4)
a       (position 5)
```

### Naive vs Ukkonen's Algorithm

**Naive approach** (O(n²)):
```
for each suffix s[i..n]:
    insert suffix into trie
```
This is slow because we process each suffix independently.

**Ukkonen's approach** (O(n)):
```
for each character s[i]:
    extend ALL current suffixes simultaneously
```
The key insight: when we add character `c` to the tree, ALL existing suffixes need to be extended with `c`. Ukkonen's algorithm does this in amortized O(1) time per suffix using clever tricks.

### The Three Extension Rules

When extending suffix tree with new character `c`:

**Rule 1: Create new leaf**
```
Current position is at a node, no edge starts with 'c'
→ Create new leaf edge labeled 'c'

     [node]                    [node]
        |           →            |  \
      (existing)              (existing) c
```

**Rule 2: Split edge**
```
Current position is in middle of edge, next char ≠ 'c'
→ Split edge, create new internal node, add leaf for 'c'

     [node]                    [node]
        |                         |
      "abc"          →          "ab"
        |                      /    \
      [leaf]               "c"      "x..."
                          [new]     [old]
```

**Rule 3: Do nothing (showstopper)**
```
Character 'c' already exists at current position
→ Suffix is already implicitly in tree, stop phase

This is the key optimization! When we find 'c' already exists,
ALL remaining suffixes will also match (they're shorter).
```

### The Active Point

The **active point** tracks where we are in the tree:

```csharp
activeNode   // Which node we're at (or descended from)
activeEdge   // Which edge we're on (first character)
activeLength // How far along that edge
```

**Why is this needed?**

When Rule 3 fires, we don't insert anything - we just note that the suffix exists implicitly. But we need to remember WHERE it exists for the next phase. The active point tracks this location.

**Example:** Building "aab"
```
After 'a':   activeNode=root, activeEdge='a', activeLength=1
             (we're 1 char into the 'a' edge)
             
After 'aa':  activeNode=root, activeEdge='a', activeLength=2
             (we're 2 chars into the 'a' edge - it grew!)
             
After 'aab': Split! The 'aa' edge becomes 'a' with children 'a' and 'b'
```

### Suffix Links

**Problem:** After inserting suffix s[j..i], we need to insert s[j+1..i].

**Naive solution:** Walk from root each time → O(n) per suffix → O(n²) total

**Suffix links solution:** Jump directly!

```
A suffix link connects internal node for "xα" to node for "α"
(where x is a single character and α is a string)

Example: node for "an" links to node for "n"
         node for "ana" links to node for "na"
```

**In code:**
```csharp
class Node {
    Node SuffixLink;  // Points to node representing shorter suffix
}
```

After inserting at node X, follow suffix link to node Y, then continue inserting. This gives us O(1) amortized time per suffix!

### The Remainder Variable

`remainder` counts how many suffixes still need explicit insertion.

**Why?** When Rule 3 fires, we stop but DON'T decrement remainder. These "pending" suffixes will be handled in future phases.

```
Building "aab":

Phase 1 ('a'): Insert 'a', remainder=0
Phase 2 ('a'): Rule 3! 'a' edge exists. remainder=1
Phase 3 ('b'): 
  - Insert "ab" (split 'aa' edge), remainder=1
  - Follow suffix link to root
  - Insert "b", remainder=0
```

### Walk-Down Operation

When `activeLength >= edgeLength`, we must "walk down" to the next node:

```csharp
if (activeLength >= edgeLength) {
    activeNode = edge;           // Move to end of edge
    activeLength -= edgeLength;  // Remaining length
    activeEdge = next character; // Start of next edge
}
```

**Example:** If activeLength=5 but edge "abc" has length 3:
```
Before: activeNode=root, activeEdge='a', activeLength=5
        (conceptually at position 5 along "abc...")
        
After:  activeNode=node_after_abc, activeEdge=?, activeLength=2
        (now 2 chars into the next edge)
```

### Terminator Character

**Problem:** Some suffixes may remain implicit after processing.

**Example:** In "aa", suffix "a" is implicit (it's a prefix of "aa").

**Solution:** Add unique terminator '$' (we use '\0'):
```
"aa" → "aa$"
Now suffix "a$" is different from "aa$", both are explicit leaves.
```

This ensures every suffix ends at a unique leaf, giving us a proper suffix tree (not a suffix trie).

### Complete Algorithm Pseudocode

```
function BuildSuffixTree(s):
    s = s + '$'  // Add terminator
    create root node
    activeNode = root
    activeEdge = none
    activeLength = 0
    remainder = 0
    
    for each character c in s:
        remainder++
        lastCreatedNode = null
        
        while remainder > 0:
            if activeLength == 0:
                activeEdge = c
            
            if no edge starting with activeEdge:
                // Rule 1: create leaf
                create leaf edge from activeNode
                setSuffixLink(lastCreatedNode, activeNode)
            else:
                edge = getEdge(activeNode, activeEdge)
                if walkDownNeeded(edge):
                    continue  // Restart with updated active point
                
                if characterAtPosition == c:
                    // Rule 3: showstopper
                    activeLength++
                    setSuffixLink(lastCreatedNode, activeNode)
                    break
                
                // Rule 2: split edge
                splitNode = splitEdge(edge, activeLength)
                create leaf from splitNode for c
                setSuffixLink(lastCreatedNode, splitNode)
                lastCreatedNode = splitNode
            
            remainder--
            
            if activeNode == root and activeLength > 0:
                // Rule 1 variant: stay at root, adjust active point
                activeLength--
                activeEdge = s[currentPosition - remainder + 1]
            else:
                // Follow suffix link (or go to root)
                activeNode = activeNode.suffixLink ?? root
```

### How It Works

Ukkonen's algorithm builds the suffix tree incrementally, processing one character at a time:

```
Phase i: Extend tree with character s[i]
  - For each suffix s[j..i] where j = 0..i:
    - If suffix already exists implicitly → do nothing (Rule 3)
    - If at a node with no matching edge → create new leaf
    - If in middle of edge with mismatch → split edge, create leaf
  - Use suffix links to jump between suffixes efficiently
```

### Key Invariants

- **Remainder** — number of suffixes to be inserted
- **Active Point** — where the next insertion will happen
- **Suffix Links** — connect internal nodes representing s[i..j] to nodes representing s[i+1..j]

### Example: Building tree for "banana"

```
Step 1: "b"      → Create leaf for 'b'
Step 2: "ba"     → Extend 'b' edge
Step 3: "ban"    → Extend 'b' edge
Step 4: "bana"   → Extend 'b' edge, implicit 'a' suffix
Step 5: "banan"  → Split edge, create suffix links
Step 6: "banana" → Continue splitting and linking
Step 7: "$"      → Add terminator, finalize all suffixes
```

### Tree Structure for "banana$"

```
root
├── a
│   ├── na
│   │   ├── na$ (suffix: "anana")
│   │   └── $ (suffix: "ana")
│   └── $ (suffix: "a")
├── banana$ (suffix: "banana")
├── na
│   ├── na$ (suffix: "nana")
│   └── $ (suffix: "na")
└── $ (suffix: "")
```

## Performance

Benchmarks on random strings (Intel i7, .NET 8.0):

| String Length | Build Time | Memory |
|---------------|------------|--------|
| 1,000 | ~0.2 ms | ~100 KB |
| 10,000 | ~2 ms | ~2.2 MB |
| 100,000 | ~65 ms | ~23 MB |
| 1,000,000 | ~700 ms | ~230 MB |

Search time is always O(m) regardless of tree size. Count and LRS operations use O(1) precomputed data after tree construction.

## Building

```bash
dotnet build
```

## Testing

```bash
# Run all 201 unit tests
dotnet test

# Run stress tests (1M+ substring verifications)
dotnet run --project SuffixTree.Console
```

### Test Coverage

The test suite is organized into 6 categories:

| Category | Tests | Description |
|----------|-------|-------------|
| **Algorithms** | 62 | LCS, LRS, suffix enumeration |
| **Search** | 45 | Contains, Count, FindAll |
| **Compatibility** | 36 | Unicode, binary data, null chars |
| **Robustness** | 26 | Edge cases, stress, thread safety |
| **Core** | 21 | Build, statistics, diagnostics |
| **Regression** | 7 | Previously fixed bugs |

### Verified Correctness

- **201 unit tests** — all passing (NUnit 4.2.2)
- **Full Unicode support** — Cyrillic, Chinese, Greek, emoji, surrogate pairs
- **Binary data support** — null characters, control characters
- **Thread safety** — concurrent read operations verified
- **Classic test cases** — "mississippi", "banana", "abracadabra"

## API Reference

### `SuffixTree.Build(string value)`

Creates a new suffix tree from the given string.

```csharp
var tree = SuffixTree.Build("hello world");
```

**Parameters:**
- `value` — the string to build the tree from (cannot be null)

**Returns:** A new `SuffixTree` instance

**Throws:** `ArgumentNullException` if value is null

---

### `Contains(string value)`

Checks if the given string is a substring of the tree content.

```csharp
bool found = tree.Contains("world");  // true
```

**Parameters:**
- `value` — the substring to search for (cannot be null)

**Returns:** `true` if substring exists, `false` otherwise

**Throws:** `ArgumentNullException` if value is null

**Time Complexity:** O(m) where m is the length of the search string

---

### `Contains(ReadOnlySpan<char> value)`

Zero-allocation overload for performance-critical scenarios.

```csharp
bool found = tree.Contains("world".AsSpan());  // true
```

---

### `CountOccurrences(string pattern)` / `CountOccurrences(ReadOnlySpan<char> pattern)`

Counts how many times a pattern appears in the text.

```csharp
var tree = SuffixTree.Build("banana");
int count = tree.CountOccurrences("ana");  // 2
```

**Time Complexity:** O(m) where m is pattern length (uses precomputed leaf counts)

---

### `FindAllOccurrences(string pattern)` / `FindAllOccurrences(ReadOnlySpan<char> pattern)`

Finds all starting positions where the pattern occurs.

```csharp
var tree = SuffixTree.Build("banana");
var positions = tree.FindAllOccurrences("ana");  // [1, 3]
```

**Time Complexity:** O(m + k) where m is pattern length, k is number of occurrences

---

### `LongestRepeatedSubstring()`

Finds the longest substring that appears at least twice.

```csharp
var tree = SuffixTree.Build("banana");
string lrs = tree.LongestRepeatedSubstring();  // "ana"
```

**Time Complexity:** O(1) after first call (O(n) initial computation, result is cached)

---

### `LongestCommonSubstring(string other)`

Finds the longest common substring between the tree's text and another string.

```csharp
var tree = SuffixTree.Build("abcdefgh");
string lcs = tree.LongestCommonSubstring("xxcdefxx");  // "cdef"
```

**Time Complexity:** O(m) where m is length of other string

---

### `LongestCommonSubstringInfo(string other)`

Finds the longest common substring with position information.

```csharp
var tree = SuffixTree.Build("abcdefgh");
var (substring, posInText, posInOther) = tree.LongestCommonSubstringInfo("xxcdefxx");
// substring="cdef", posInText=2, posInOther=2
```

**Time Complexity:** O(m) where m is length of other string

---

### `FindAllLongestCommonSubstrings(string other)`

Finds all positions where the longest common substring occurs.

```csharp
var tree = SuffixTree.Build("abcabc");
var (substring, positionsInText, positionsInOther) = tree.FindAllLongestCommonSubstrings("xabcyabcz");
// substring="abc", positionsInOther=[1, 5]
```

**Time Complexity:** O(m) where m is length of other string

---

### `NodeCount`, `LeafCount`, `MaxDepth`

Tree statistics properties.

```csharp
var tree = SuffixTree.Build("banana");
int nodes = tree.NodeCount;    // Total nodes including root
int leaves = tree.LeafCount;   // Number of leaf nodes (suffixes + 1)
int depth = tree.MaxDepth;     // Maximum depth in characters
```

**Time Complexity:** O(n) for each property

---

### `GetAllSuffixes()`

Returns all suffixes in lexicographic order (useful for debugging).

```csharp
var tree = SuffixTree.Build("banana");
var suffixes = tree.GetAllSuffixes();
// ["a", "ana", "anana", "banana", "na", "nana"]
```

---

### `PrintTree()`

Returns a detailed string representation of the tree structure (for debugging).

```csharp
Console.WriteLine(tree.PrintTree());
```

---

### `ToString()`

Returns a summary of the tree content.

```csharp
Console.WriteLine(tree);  // "SuffixTree (length: 12, content: "hello world")"
```

## Applications

Suffix trees are used in:

- **Text editors** — find/replace, autocomplete
- **Bioinformatics** — DNA/protein sequence analysis
- **Data compression** — LZ77, LZW algorithms
- **Plagiarism detection** — finding similar text passages
- **Search engines** — substring indexing
- **Spell checkers** — approximate string matching

## References

- Ukkonen, E. (1995). *On-line construction of suffix trees*. Algorithmica, 14(3), 249-260.
- Gusfield, D. (1997). *Algorithms on Strings, Trees, and Sequences*. Cambridge University Press.
- [Suffix Tree Visualization](https://visualgo.net/en/suffixtree)

## Requirements

- .NET 8.0 SDK

## License

See [LICENSE](LICENSE) file.

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.
