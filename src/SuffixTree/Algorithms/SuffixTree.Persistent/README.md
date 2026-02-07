# SuffixTree.Persistent

A high-performance, disk-backed persistent suffix tree implementation for .NET.

Based on Ukkonen's algorithm, this library provides an efficient way to index and search large volumes of text without loading the entire tree structure into heap memory. By leveraging Memory-Mapped Files (MMF), it allows for O(m) substring searches on trees that exceed available RAM.

## Features

- **Disk-Backed Storage**: Uses Memory-Mapped Files specifically optimized for suffix tree structures.
- **Persistence**: Build once, save to disk, and reload instantly in subsequent sessions.
- **Hybrid Storage Format**: Construction starts with the compact layout and automatically promotes to the large layout only if the storage actually exceeds the uint32 address space:
  - **Compact (v4)**: 32-bit offsets, 28-byte nodes, 8-byte child entries. ~30% smaller files, better CPU cache locality. Used for all trees that fit within ~4 GB of storage.
  - **Large (v3)**: 64-bit offsets, 40-byte nodes, 12-byte child entries. No practical size limit. Activated automatically when compact overflows.
  - Format is detected automatically on load — no user configuration needed.
- **Scalability**: Capable of handling gigabytes of text by offloading the node graph to disk.
- **Thread-Safe Read Access**: Supports multiple concurrent readers for the same tree file.
- **Suffix Links Preserved**: Ukkonen's algorithm writes suffix links natively, enabling `FindExactMatchAnchors` on persistent trees.
- **Logical Identity**: Deterministic checksumming (SHA256) ensures distinct trees can be uniquely identified regardless of their memory layout or storage format.
- **Portable Serialization (v2)**: Export stores text + SHA256 hash; import rebuilds via Ukkonen's, guaranteeing 100% functionality including suffix links.

## Usage Examples

### Creating a Persistent Tree

```csharp
using SuffixTree.Persistent;

string text = "abracadabra";
string filePath = "my_tree.suffix";

// This builds the tree directly into a memory-mapped file
using (var tree = PersistentSuffixTreeFactory.Create(text, filePath) as IDisposable)
{
    var st = (ISuffixTree)tree;
    bool found = st.Contains("bra"); // O(m)
    int count = st.CountOccurrences("a"); // 5
}
```

### Loading an Existing Tree

```csharp
// Load in read-only mode for lightning-fast startup and memory efficiency
using (var tree = PersistentSuffixTreeFactory.Load(filePath) as IDisposable)
{
    var st = (ISuffixTree)tree;
    var positions = st.FindAllOccurrences("abra");
}
```

## Internal Architecture

The library follows Clean Architecture principles:

- **[IStorageProvider](file:///d:/Prototype/src/SuffixTree/Algorithms/SuffixTree.Persistent/IStorageProvider.cs)**: An abstraction layer for memory access.
  - `HeapStorageProvider`: Uses a byte array for in-memory binary operations.
  - `MappedFileStorageProvider`: Uses `MemoryMappedFile` for disk-based persistence.
- **[PersistentSuffixTreeNode](file:///d:/Prototype/src/SuffixTree/Algorithms/SuffixTree.Persistent/PersistentSuffixTreeNode.cs)**: A binary handle (struct) to a fixed-size node in storage.
- **[PersistentSuffixTreeBuilder](file:///d:/Prototype/src/SuffixTree/Algorithms/SuffixTree.Persistent/PersistentSuffixTreeBuilder.cs)**: Implements Ukkonen's algorithm for stream-based construction.
- **[PersistentSuffixTree](file:///d:/Prototype/src/SuffixTree/Algorithms/SuffixTree.Persistent/PersistentSuffixTree.cs)**: The search engine implementing `ISuffixTree`.

### Binary Layout (Node)

Nodes are fixed-size blocks whose layout depends on the automatically selected storage format.
The format is recorded in the file header (version field) and detected on load.

#### Compact Layout (v4, default) — 28-byte nodes

| Offset | Size | Field | Description |
| :--- | :--- | :--- | :--- |
| 0 | 4 | Start | Edge start position in text |
| 4 | 4 | End | Edge end position |
| 8 | **4** | SuffixLink | Offset to suffix link node (uint32) |
| 12 | 4 | Depth | Node depth from root |
| 16 | 4 | LeafCount | Number of leaves in subtree |
| 20 | **4** | ChildrenHead | Offset to first child entry (uint32) |
| 24 | 4 | ChildCount | Number of children |

Child entry: Key (4 B) + ChildNodeOffset (**4 B**) = **8 bytes**.

#### Large Layout (v3) — 40-byte nodes

| Offset | Size | Field | Description |
| :--- | :--- | :--- | :--- |
| 0 | 4 | Start | Edge start position in text |
| 4 | 4 | End | Edge end position |
| 8 | **8** | SuffixLink | Offset to suffix link node (int64) |
| 16 | 4 | Depth | Node depth from root |
| 20 | 4 | LeafCount | Number of leaves in subtree |
| 24 | **8** | ChildrenHead | Offset to first child entry (int64) |
| 32 | 4 | ChildCount | Number of children |
| 36 | 4 | (padding) | Alignment padding |

Child entry: Key (4 B) + ChildNodeOffset (**8 B**) = **12 bytes**.

#### Performance Notes

- The Compact format saves ~30% disk space and improves CPU cache hit rates because more nodes fit in L1/L2/L3 cache lines.
- The `NodeLayout` class uses `[AggressiveInlining]` for offset read/write helpers; the `if (OffsetIs64Bit)` branch is perfectly predicted by the CPU (same direction every time), adding zero measurable overhead.
- Format promotion is transparent: the builder starts with Compact and only retries with Large if an allocation would exceed the uint32 limit. The cost of the retry (at most 2× build time) is amortized over the lifetime of the tree.

### Algorithmic Stability
- **Iterative Algorithms**: Leaf counting and tree traversal are implemented using stacks rather than recursion. This prevents `StackOverflowException` when processing extremely deep or repetitive trees.

## Serialization & Portability

The `SuffixTreeSerializer` utility provides logical identity, portable export/import, and convenient file-based persistence.

**Format v2**: Export stores only the source text and a SHA256 structural hash. Import rebuilds the tree from text via `PersistentSuffixTreeBuilder` (Ukkonen's algorithm), guaranteeing 100% functionality including suffix links for `FindExactMatchAnchors`.

### 1. Calculating Logical Hash
Get a unique fingerprint of the tree structure that is identical across different memory layouts.

```csharp
byte[] hash = SuffixTreeSerializer.CalculateLogicalHash(tree);
string hashString = BitConverter.ToString(hash).Replace("-", "");
Console.WriteLine($"Tree Fingerprint: {hashString}");
```

### 2. Stream-Based Export / Import
Export any `ISuffixTree` to a portable stream; import into any `IStorageProvider`.

```csharp
// Export (compact: text + hash only)
using (var outStream = File.Create("tree.bin"))
    SuffixTreeSerializer.Export(tree, outStream);

// Import — rebuilds via Ukkonen, suffix links created natively
using (var inStream = File.OpenRead("tree.bin"))
{
    var storage = new HeapStorageProvider();
    var imported = SuffixTreeSerializer.Import(inStream, storage);
    imported.FindExactMatchAnchors("query", 3); // works — suffix links intact
}
```

### 3. File-Based Save / Load (MMF)
Convenience methods that create a native persistent tree file.

```csharp
// Save any ISuffixTree to a memory-mapped file
using var saved = SuffixTreeSerializer.SaveToFile(tree, "saved.tree") as IDisposable;

// Load from file — read-only, instant startup
using var loaded = SuffixTreeSerializer.LoadFromFile("saved.tree") as IDisposable;
var st = (ISuffixTree)loaded;
st.FindExactMatchAnchors("query", 3); // full functionality
```

## Quality Assurance

The library is verified by 157 tests:
- **Differential Parity**: Cross-validation against the reference in-memory implementation for all `ISuffixTree` methods.
- **Format Parity**: Compact and Large formats produce identical query results and logical hashes.
- **Anchor Parity**: `FindExactMatchAnchors` results match between in-memory and persistent trees.
- **Logical Parity**: SHA256 checksums match across Heap, MMF, imported, Compact, and Large storage backends.
- **Serialization Round-Trip**: Export → Import and SaveToFile → LoadFromFile preserve full functionality.
- **Auto-Detection**: `Load()` correctly identifies v3 (Large) and v4 (Compact) files from the header.
- **Stress Testing**: Validates stability with large datasets (1MB+) and concurrent access.

## License

This project is part of the Prototype SuffixTree suite.
