# SuffixTree

A high-performance implementation of **Ukkonen's suffix tree algorithm** in C#.

[![.NET](https://img.shields.io/badge/.NET-8.0-512BD4)](https://dotnet.microsoft.com/)
[![Tests](https://img.shields.io/badge/tests-143%20passed-brightgreen)]()
[![License](https://img.shields.io/badge/license-MIT-blue)](LICENSE)

## Overview

A **suffix tree** is a compressed trie containing all suffixes of a given string. It is one of the most powerful data structures in stringology, enabling efficient solutions to numerous string processing problems.

This implementation uses **Ukkonen's online algorithm** (1995), which constructs the tree in **O(n)** time and allows **O(m)** substring search, where *n* is the length of the text and *m* is the length of the pattern.

### Why Suffix Trees?

Suffix trees provide optimal time complexity for many string operations:

| Operation | Time Complexity |
|-----------|-----------------|
| Build tree | O(n) |
| Substring search | O(m) |
| Find all occurrences | O(m + k) |
| Longest repeated substring | O(n) |
| Longest common substring | O(n + m) |

*where n = text length, m = pattern length, k = number of occurrences*

## Features

- **O(n) construction** — linear time tree building using Ukkonen's algorithm
- **O(m) substring search** — efficient pattern matching
- **Full Unicode support** — works with any characters (Cyrillic, Chinese, emoji, etc.)
- **Terminator character** — ensures all suffixes are explicit (proper suffix tree)
- **Suffix links** — optimized traversal during construction
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

// Check if substring exists
bool found = tree.Contains("nan");  // true
bool notFound = tree.Contains("xyz");  // false

// All suffixes are searchable
tree.Contains("banana");  // true
tree.Contains("anana");   // true
tree.Contains("nana");    // true
tree.Contains("ana");     // true
tree.Contains("na");      // true
tree.Contains("a");       // true

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
├── SuffixTree.Tests/        # Unit tests (NUnit)
│   └── UnitTests.cs         # 143 comprehensive tests
└── SuffixTree.Benchmarks/   # Performance benchmarks
    └── SuffixTreeBenchmarks.cs
```

## Benchmarks

Performance measured on Intel Core i7-1185G7 @ 3.00GHz, .NET 8.0:

### Build Performance

| Text Size | Time | Allocated |
|-----------|------|-----------|
| 100 chars | 7.1 μs | 23 KB |
| 10K chars | 1.84 ms | 3.5 MB |
| 100K chars | 67.5 ms | 36 MB |
| 50K DNA | 50.6 ms | 15 MB |

### Search Performance

| Operation | 100 chars | 10K chars | 100K chars | 50K DNA |
|-----------|-----------|-----------|------------|---------|
| Contains | 26 ns | 56 ns | 58 ns | 118 ns |
| Count | — | 906 ns | 33.7 μs | 2.1 μs |
| FindAll | 49 ns | 1.07 μs | 33.1 μs | 2.65 μs |
| LRS | 1.0 μs | 494 μs | 20.3 ms | 7.2 ms |

*LRS = Longest Repeated Substring*

### Key Observations

- **Build time scales linearly** with text size (O(n) confirmed)
- **Contains is extremely fast** — 26-118 ns regardless of tree size
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
3. **Rule 1** — handle extension at root node
4. **Rule 3** — implicit extension (showstopper)
5. **Walk-Down** — traverse edges when `activeLength >= edgeLength`
6. **Terminator** — `'\0'` character ensures all suffixes end at leaves

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
| 1,000 | < 1 ms | ~50 KB |
| 10,000 | ~2 ms | ~500 KB |
| 100,000 | ~20 ms | ~5 MB |
| 1,000,000 | ~200 ms | ~50 MB |

Search time is always O(m) regardless of tree size.

## Building

```bash
dotnet build
```

## Testing

```bash
# Run unit tests (88 tests)
dotnet test

# Run stress tests (1M+ substring verifications)
dotnet run --project SuffixTree.Console
```

### Test Coverage

The test suite includes:

- **Basic operations** — build, contains, edge cases
- **Algorithm verification** — suffix links, walk-down, edge splitting
- **Stress tests** — random strings up to 5000 characters
- **Unicode tests** — Cyrillic, Chinese, emoji, mixed content
- **Performance tests** — 10,000+ character strings
- **Regression tests** — all previously fixed bugs

### Verified Correctness

- **88 unit tests** — all passing
- **1,137,909 substrings** — verified in stress tests
- **313 test strings** — including edge cases like "mississippi", "abracadabra"

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

**Time Complexity:** O(m + k) where m is pattern length, k is number of occurrences

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

**Time Complexity:** O(n)

---

### `LongestCommonSubstring(string other)`

Finds the longest common substring between the tree's text and another string.

```csharp
var tree = SuffixTree.Build("abcdefgh");
string lcs = tree.LongestCommonSubstring("xxcdefxx");  // "cdef"
```

**Time Complexity:** O(m) where m is length of other string

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
