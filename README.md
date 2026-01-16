# SuffixTree

A high-performance implementation of **Ukkonen's suffix tree algorithm** in C#.

## Overview

Suffix tree is a compressed trie of all suffixes of a given string. This implementation uses Ukkonen's online algorithm, which constructs the tree in **O(n)** time and allows **O(m)** substring search, where *n* is the length of the text and *m* is the length of the pattern.

## Features

- **O(n) construction** — linear time tree building using Ukkonen's algorithm
- **O(m) substring search** — efficient pattern matching
- **Full Unicode support** — works with any characters (Cyrillic, Chinese, emoji, etc.)
- **Terminator character** — ensures all suffixes are explicit (proper suffix tree)
- **Suffix links** — optimized traversal during construction
- **.NET 8.0** — modern SDK-style project

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
│   └── SuffixTree.cs        # Ukkonen's algorithm implementation
├── SuffixTree.Console/      # Stress test console app
│   └── Program.cs           # Exhaustive verification tests
└── SuffixTree.Tests/        # Unit tests (NUnit)
    └── UnitTests.cs         # 88 comprehensive tests
```

## Algorithm Details

The implementation follows Ukkonen's canonical algorithm with:

1. **Active Point** — tracks current position in tree (`activeNode`, `activeEdge`, `activeLength`)
2. **Suffix Links** — enable O(1) jumps between internal nodes
3. **Rule 1** — handle extension at root node
4. **Rule 3** — implicit extension (showstopper)
5. **Walk-Down** — traverse edges when `activeLength >= edgeLength`
6. **Terminator** — `'\0'` character ensures all suffixes end at leaves

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

## Requirements

- .NET 8.0 SDK

## License

See [LICENSE](LICENSE) file.
