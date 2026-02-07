# SuffixTree.Persistent

A high-performance, disk-backed persistent suffix tree implementation for .NET 9.

Built on Ukkonen's algorithm, the library indexes and searches large volumes of text
without loading the tree into managed heap memory. Memory-mapped files (MMF) keep the
node graph on disk while providing O(m) substring lookups.

---

## Features

| Feature | Details |
|---------|---------|
| **Disk-backed storage** | Memory-mapped files optimized for suffix tree node access. |
| **Persistence** | Build once → save to disk → reload instantly in subsequent sessions. |
| **Hybrid storage format (v5)** | Starts in compact (32-bit) mode; promotes to large (64-bit) mid-build only if storage exceeds ~4 GB. A jump table bridges cross-zone suffix links and child arrays. Format detected automatically on load. |
| **Scalable** | Handles gigabytes of text by offloading the node graph to disk. |
| **Thread-safe reads** | Multiple concurrent readers on the same tree file. |
| **Suffix links preserved** | Written natively by Ukkonen's algorithm, enabling `FindExactMatchAnchors`. |
| **Deterministic identity** | SHA256 structural hash — identical for same content regardless of layout. |
| **Portable serialization (v2)** | Export stores text + SHA256 hash; import rebuilds via Ukkonen's, guaranteeing full functionality. |

---

## Quick Start

### Build a Persistent Tree

```csharp
using SuffixTree.Persistent;

string text = "abracadabra";
string filePath = "my_tree.suffix";

// Build directly into a memory-mapped file (hybrid v5 format)
using var tree = (IDisposable)PersistentSuffixTreeFactory.Create(
    new StringTextSource(text), filePath);

var st = (ISuffixTree)tree;
st.Contains("bra");            // O(m)
st.CountOccurrences("a");      // 5
st.LongestRepeatedSubstring(); // "abra"
```

### Load an Existing Tree

```csharp
// Read-only mode — instant startup, minimal memory footprint
using var tree = (IDisposable)PersistentSuffixTreeFactory.Load("my_tree.suffix");
var st = (ISuffixTree)tree;

var positions = st.FindAllOccurrences("abra"); // [0, 7]
st.FindExactMatchAnchors("bandana", 3);        // suffix links preserved
```

### Export / Import (Portable)

```csharp
// Export any ISuffixTree to a stream (text + SHA256 hash)
using var ms = new MemoryStream();
SuffixTreeSerializer.Export(tree, ms);

// Import — rebuilds via Ukkonen, full suffix link support
ms.Position = 0;
var imported = SuffixTreeSerializer.Import(ms, new HeapStorageProvider());

// File-based convenience
using var saved  = (IDisposable)SuffixTreeSerializer.SaveToFile(tree, "portable.tree");
using var loaded = (IDisposable)SuffixTreeSerializer.LoadFromFile("portable.tree");
```

---

## Architecture

```
┌──────────────────────────────────────────────────────────┐
│                     ISuffixTree                          │
│  (Contains, FindAll, Count, LRS, LCS, Anchors, ...)     │
└──────────────────┬───────────────────────────────────────┘
                   │
        ┌──────────┴──────────┐
        │                     │
  SuffixTree             PersistentSuffixTree
  (in-memory)            (disk-backed)
        │                     │
  SuffixTreeNode         PersistentSuffixTreeNode
  (heap objects)         (struct handles → byte offsets)
                              │
                     ┌────────┴────────┐
                     │                 │
              HeapStorage       MappedFileStorage
              (byte[])          (MemoryMappedFile)
```

### Key Components

| Component | Responsibility |
|-----------|---------------|
| `IStorageProvider` | Abstraction for byte-level read/write. Two implementations: heap array and MMF. |
| `PersistentSuffixTreeBuilder` | Ukkonen's algorithm writing nodes directly to storage. Handles compact→large transitions. |
| `PersistentSuffixTreeFactory` | Static factory: `Create` (build), `Load` (read-only open). |
| `PersistentSuffixTree` | Query engine implementing `ISuffixTree`. Delegates shared algorithms via `ISuffixTreeNavigator<TNode>`. |
| `PersistentSuffixTreeNode` | Struct handle — reads/writes fixed-size fields at a byte offset in storage. |
| `PersistentSuffixTreeNavigator` | Struct navigator for generic algorithm dispatch. All methods `[AggressiveInlining]`. |
| `NodeLayout` | Immutable descriptor of field sizes/offsets. Two singletons: `Compact` (v4) and `Large` (v3). |
| `HybridLayout` | Zone detection + jump table resolution for v5 hybrid trees. |
| `MemoryMappedTextSource` | `ITextSource` backed by a raw pointer into MMF. Zero-copy character access. |
| `SuffixTreeSerializer` | SHA256 hash, export/import (format v2), file save/load. |

---

## Binary Format

### File Header

All format versions share bytes 0–47. Hybrid v5 extends to byte 79.

#### Base Header (48 bytes) — v3 / v4 / v5

| Offset | Size | Field | Description |
|--------|------|-------|-------------|
| 0–7 | 8 | `MAGIC` | `0x5452454558494646` — ASCII `"SUFFIXTR"` |
| 8–11 | 4 | `VERSION` | `3` = Large, `4` = Compact, `5` = Hybrid |
| 12–15 | 4 | `TEXT_LEN` | Character count |
| 16–23 | 8 | `ROOT` | int64 byte offset of root node |
| 24–31 | 8 | `TEXT_OFF` | int64 byte offset of stored text |
| 32–35 | 4 | `NODE_COUNT` | Total node count |
| 36–39 | 4 | _(reserved)_ | — |
| 40–47 | 8 | `SIZE` | Expected file size (validated on load) |

#### v5 Extended Header (bytes 48–79)

| Offset | Size | Field | Description |
|--------|------|-------|-------------|
| 48–55 | 8 | `TRANSITION` | Byte offset where compact→large transition occurred (`-1` if pure compact) |
| 56–63 | 8 | `JUMP_START` | Start of the jump table region |
| 64–71 | 8 | `JUMP_END` | End of the jump table region |
| 72–79 | 8 | `DEEPEST_NODE` | Byte offset of the deepest internal node (for O(1) LRS) |

### Node Layouts

#### Compact (v4) — 28-byte Nodes

| Offset | Size | Field | Type |
|--------|------|-------|------|
| 0–3 | 4 | Start | uint32 |
| 4–7 | 4 | End | uint32 |
| 8–11 | 4 | SuffixLink | uint32 |
| 12–15 | 4 | DepthFromRoot | uint32 |
| 16–19 | 4 | LeafCount | uint32 |
| 20–23 | 4 | ChildrenHead | uint32 |
| 24–27 | 4 | ChildCount | int32 |

**Child entry:** Key (4 B) + Offset (4 B) = **8 bytes**.
**NULL sentinel:** `0xFFFFFFFF`. Max valid offset: `0xFFFFFFFE` (~4.29 GB).

#### Large (v3) — 40-byte Nodes

| Offset | Size | Field | Type |
|--------|------|-------|------|
| 0–3 | 4 | Start | uint32 |
| 4–7 | 4 | End | uint32 |
| 8–15 | 8 | SuffixLink | int64 |
| 16–19 | 4 | DepthFromRoot | uint32 |
| 20–23 | 4 | LeafCount | uint32 |
| 24–31 | 8 | ChildrenHead | int64 |
| 32–35 | 4 | ChildCount | int32 |
| 36–39 | 4 | _(padding)_ | — |

**Child entry:** Key (4 B) + Offset (8 B) = **12 bytes**.
**NULL sentinel:** `-1L`.

### Hybrid v5 Format

Construction always starts in compact mode. When the next allocation would
exceed `uint.MaxValue - 1`, the builder transitions to large mode:

1. **Transition offset** — recorded in the header at byte 48.
2. **Jump table** — a contiguous block of `int64` entries materialized at finalization:
   - *Suffix-link jumps:* Compact nodes whose suffix links point into the large zone.
     The compact `uint32` field points to the jump entry; the entry stores the `int64` target.
   - *Child-array jumps:* Compact parents with any large-zone children.
     `ChildrenHead` points to the jump entry; child entries use 12-byte (large) format.
     The high bit of `ChildCount` is set as a "jumped" flag.
3. **Zone detection** — `HybridLayout.GetLayoutForOffset(offset)` returns Compact if
   `offset < transitionOffset`, else Large.
4. **Jump resolution** — `HybridLayout.ResolveSuffixLink(offset)` dereferences `int64`
   targets for offsets falling within the jump table range.

### Serialization Format v2

Used by `SuffixTreeSerializer.Export` / `Import`:

| Field | Type | Description |
|-------|------|-------------|
| Magic | int64 | `0x53544C4F47494332` — ASCII `"STLOGIC2"` |
| Version | int32 | `2` |
| Text length | 7-bit encoded int | Character count |
| Text chars | char[] | Raw characters |
| Node count | int32 | Expected after rebuild |
| Hash length | int32 | Always `32` (SHA256) |
| Hash | byte[32] | SHA256 of text + tree structure |

Import rebuilds the tree from scratch via Ukkonen's algorithm, guaranteeing
correct suffix links. Both node count and structural hash are validated after rebuild.

---

## Performance

### Format Savings

The compact layout saves ~30% disk space vs. large and improves CPU cache hit rates
(more nodes per cache line). The `NodeLayout` accessor methods use `[AggressiveInlining]`;
the `if (OffsetIs64Bit)` branch is perfectly predicted by the CPU (same direction for
every node in single-zone trees), adding zero measurable overhead.

### Iterative Algorithms

All tree traversals (leaf counting, child serialization, depth computation, `FindAll`,
`PrintTree`) use explicit stacks instead of recursion — no `StackOverflowException`
even for extremely deep or repetitive trees.

### O(1) Longest Repeated Substring

The deepest internal node is identified during the builder's post-order pass and
its offset is stored in the v5 header (byte 72). `LongestRepeatedSubstring()` reads
this single node — no traversal on query.

### Chunked Text Storage

Text is written to the storage in 4096-character chunks, avoiding large contiguous
allocations during build. On load, `MemoryMappedTextSource` provides zero-copy character
access via a raw pointer into the MMF.

### Construction Memory

After `FinalizeTree`, the builder clears all in-memory collections (suffix link map,
child lists, deferred jump entries) to release GC pressure before the tree is returned.

---

## Load Validation

`PersistentSuffixTree.Load` performs comprehensive validation before returning a tree:

| Check | Description |
|-------|-------------|
| Storage ≥ header size | Detects truncated files |
| Magic bytes | Must be `"SUFFIXTR"` |
| Version ∈ {3, 4, 5} | Rejects unknown formats |
| Root offset in bounds | Root must be within storage |
| Root offset ≥ 0 | Rejects negative offsets |
| Text offset in bounds | Text region must be within storage |
| Text length ≥ 0 | Rejects negative lengths |
| Text region fits | `textOff + textLen × 2 ≤ Size` |
| v5: transition in bounds | Transition offset must be within storage |
| v5: jump end ≥ jump start | Jump table must be non-negative length |
| v5: deepest node in bounds | Offset must be within storage (if ≥ 0) |
| v5: deepest node ≥ −1 | Rejects negative offsets other than sentinel |
| Header SIZE = storage size | Detects mismatch / corruption |

---

## Quality Assurance

**339 tests** organized into 6 categories:

### Core/ — 37 tests
Base functionality across all three storage backends (Heap, MMF, Large-format)
plus builder guard rails.

| File | Tests | Focus |
|------|-------|-------|
| `SuffixTreeTestBase.cs` | 8 (×3) | Contains, FindAll, Count, LRS, LCS, empty, Unicode, suffixes |
| `HeapSuffixTreeTests.cs` | 8 | Heap backend |
| `LargeFormatSuffixTreeTests.cs` | 8 | Large (v3) backend |
| `MmfSuffixTreeTests.cs` | 13 | MMF backend + persistence, 50 KB / 1 MB stress, concurrency, corruption |
| `BuilderGuardTests.cs` | 8 | Double-build guard, deepest node tracking, O(1) LRS, memory release |

### Format/ — 109 tests
Storage format correctness, compact vs. large parity, hybrid transition zone sweep.

| File | Tests | Focus |
|------|-------|-------|
| `StorageFormatTests.cs` | 20 | Compact/Large build+query, format detection, size savings, layout constants |
| `HybridTransitionZoneTests.cs` | 89 | Exhaustive transition-point sweep (×8 texts), cross-zone suffix links, stress |

### Parity/ — 17 tests
Differential testing between persistent and in-memory reference implementations.

| File | Tests | Focus |
|------|-------|-------|
| `ParityTests.cs` | 11 | Full API parity (10 strings + 10 random) |
| `TopologyParityTests.cs` | 6 | Node-by-node topology comparison (6 texts) |

### Safety/ — 78 tests
Dispose, concurrency, providers, node API, null guards, text source, overflow.

| File | Tests | Focus |
|------|-------|-------|
| `DisposeBehaviorTests.cs` | 12 | ODE for all 17 methods, idempotent dispose, read-only, PrintTree, exception safety |
| `ConcurrencyTests.cs` | 5 | 1000 parallel MMF reads, 5000 heap reads, concurrent dispose trials |
| `ProviderSafetyTests.cs` | 28 | ODE, negative allocate, data preservation, bounds checks (both providers) |
| `NodeApiSafetyTests.cs` | 14 | TryGetChild, jumped handling, offset overflow, internal API visibility |
| `NullGuardTests.cs` | 3 | Null guards for factory methods |
| `MemoryMappedTextSourceTests.cs` | 9 | Char access, slicing, pointer cleanup, dispose safety, TOCTOU fix |
| `OverflowAndValidationTests.cs` | 7 | Integer overflow, slice overflow, negative length, zero-cap, chunked text |

### Serialization/ — 14 tests
Export/import round-trips, hash consistency, truncation detection.

| File | Tests | Focus |
|------|-------|-------|
| `LogicalPersistenceTests.cs` | 6 | Hash parity (Heap/MMF/Reference), byte-identical exports, round-trip |
| `SerializerHashTests.cs` | 2 | Same-tree hash stability, little-endian regression |
| `ImportTruncationTests.cs` | 2 | Truncated stream, zero-chars-after-header |
| `SetSizeAndExportTests.cs` | 4 | Post-dispose guard, negative size, valid update, 10 K chunked export |

### Validation/ — 49 tests
Load validation, mathematical invariants, tree contracts, boundary cases.

| File | Tests | Focus |
|------|-------|-------|
| `LoadValidationTests.cs` | 17 | Truncated storage, bad magic, unknown version, bounds, SIZE mismatch |
| `MathematicalInvariantsTests.cs` | 10 | LeafCount = text.Length, NodeCount bounds, terminator isolation |
| `TreeContractTests.cs` | 12 | Lex suffix order, null guard, LCS positions, depth, visitor balance |
| `EmptyPatternContractTests.cs` | 4 | Empty → all positions, count = text.Length |
| `BoundaryAndEdgeCaseTests.cs` | 6 | Frequent resizing, edge count invariant, pathological patterns |

---

## License

Part of the Seqeron bioinformatics suite. See [LICENSE](../../../../LICENSE).
