# SuffixTree.Persistent

A high-performance, disk-backed persistent suffix tree implementation for .NET 9.

Built on Ukkonen's algorithm, the library indexes and searches large volumes of text
without loading the tree into managed heap memory. Memory-mapped files (MMF) keep the
node graph on disk while providing O(m log d) substring lookups (`d` = branching factor).

---

## Features

| Feature | Details |
|---------|---------|
| **Disk-backed storage** | Memory-mapped files optimized for suffix tree node access. |
| **Persistence** | Build once → save to disk → reload instantly in subsequent sessions. |
| **Hybrid storage format (v6)** | Starts in compact (32-bit) mode; promotes to large (64-bit) mid-build only if storage exceeds ~4 GB. A jump table bridges cross-zone suffix links and child arrays. Format detected automatically on load. |
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

// Build directly into a memory-mapped file (hybrid v6 format)
using var tree = (IDisposable)PersistentSuffixTreeFactory.Create(
    new StringTextSource(text), filePath);

var st = (ISuffixTree)tree;
st.Contains("bra");            // O(m log d)
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
| `NodeLayout` | Immutable descriptor of field sizes/offsets. Two singletons: `Compact` (v6, 24-byte node) and `Large` (v6, 32-byte node). |
| `HybridLayout` | Zone detection + jump table resolution for hybrid v6 trees. |
| `MemoryMappedTextSource` | `ITextSource` backed by a raw pointer into MMF. Zero-copy character access. |
| `SuffixTreeSerializer` | SHA256 hash, export/import (format v2), file save/load. |

---

## Binary Format

### File Header

Current runtime format is **v6**. Header size is **88 bytes**.

#### v6 Header (bytes 0–87)

| Offset | Size | Field | Description |
|--------|------|-------|-------------|
| 0–7 | 8 | `MAGIC` | `0x5452454558494646` — ASCII `"SUFFIXTR"` |
| 8–11 | 4 | `VERSION` | `6` |
| 12–15 | 4 | `TEXT_LEN` | Character count |
| 16–23 | 8 | `ROOT` | int64 byte offset of root node |
| 24–31 | 8 | `TEXT_OFF` | int64 byte offset of stored text |
| 32–35 | 4 | `NODE_COUNT` | Total node count |
| 36–39 | 4 | `FLAGS` | Bit field (`FLAG_TEXT_ASCII` for 1-byte text storage) |
| 40–47 | 8 | `SIZE` | Expected file size (validated on load) |
| 48–55 | 8 | `TRANSITION` | Byte offset where compact→large transition occurred (`-1` if pure compact) |
| 56–63 | 8 | `JUMP_START` | Start of the jump table region |
| 64–71 | 8 | `JUMP_END` | End of the jump table region |
| 72–79 | 8 | `DEEPEST_NODE` | Byte offset of deepest internal node (`-1` if unavailable) |
| 80–83 | 4 | `LRS_DEPTH` | Precomputed LRS depth (`DepthFromRoot + edgeLength`) |
| 84–87 | 4 | `BASE_NODE_SIZE` | Base node size (`24`=Compact, `32`=Large) |

### Node Layouts

#### Compact (v6) — 24-byte Nodes

| Offset | Size | Field | Type |
|--------|------|-------|------|
| 0–3 | 4 | Start | uint32 |
| 4–7 | 4 | End | uint32 |
| 8–11 | 4 | SuffixLink | uint32 |
| 12–15 | 4 | LeafCount | uint32 |
| 16–19 | 4 | ChildrenHead | uint32 |
| 20–23 | 4 | ChildCount | int32 |

**Child entry:** Key (4 B) + Offset (4 B) = **8 bytes**.
**NULL sentinel:** `0xFFFFFFFF`. Max valid offset: `0xFFFFFFFE` (~4.29 GB).

#### Large (v6) — 32-byte Nodes

| Offset | Size | Field | Type |
|--------|------|-------|------|
| 0–3 | 4 | Start | uint32 |
| 4–7 | 4 | End | uint32 |
| 8–15 | 8 | SuffixLink | int64 |
| 16–19 | 4 | LeafCount | uint32 |
| 20–27 | 8 | ChildrenHead | int64 |
| 28–31 | 4 | ChildCount | int32 |

**Child entry:** Key (4 B) + Offset (8 B) = **12 bytes**.
**NULL sentinel:** `-1L`.

### Hybrid v6 Format

Construction always starts in compact mode. When the next allocation would
exceed `uint.MaxValue - 1`, the builder transitions to large mode:

1. **Transition offset** — recorded in the header at byte 48.
2. **Jump table** — a contiguous block of `int64` entries materialized at finalization:
   - *Suffix-link jumps:* Compact nodes whose suffix links point into the large zone.
     The compact `uint32` field points to the jump entry; the entry stores the `int64` target.
   - *Child-array jumps:* Compact parents with any large-zone children.
     `ChildrenHead` points to the jump entry; child entries use 12-byte (large) format.
     The high bit of `ChildCount` is set as a "jumped" flag.
3. **Zone detection** — `HybridLayout.LayoutForOffset(offset)` returns Compact if
   `offset < transitionOffset`, else Large.
4. **Jump resolution** — `HybridLayout.ResolveJump(offset)` dereferences `int64`
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

The compact layout saves ~25% disk space vs. large and improves CPU cache hit rates
(more nodes per cache line). The `NodeLayout` accessor methods use `[AggressiveInlining]`;
the `if (OffsetIs64Bit)` branch is perfectly predicted by the CPU (same direction for
every node in single-zone trees), adding zero measurable overhead.

### Iterative Algorithms

All tree traversals (leaf counting, child serialization, depth computation, `FindAll`,
`PrintTree`) use explicit stacks instead of recursion — no `StackOverflowException`
even for extremely deep or repetitive trees.

### Fast-Path Longest Repeated Substring

The deepest internal node is identified during the builder's post-order pass and
its offset/depth are stored in the v6 header (bytes 72 and 80). `LongestRepeatedSubstring()`
uses this metadata on the fast path; if metadata is unavailable, it falls back to a one-time DFS.

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
| Version = 6 | Rejects unknown formats |
| Root offset in bounds | Root must be within `[HEADER_SIZE_V6, Size)` |
| Text offset in bounds | Text region must be within storage |
| Text length ≥ 0 | Rejects negative lengths |
| Text region fits | `textOff + textLen × bytesPerChar ≤ Size` (`bytesPerChar` from flags) |
| v6: transition in bounds | Transition offset must be within storage (when present) |
| v6: jump end ≥ jump start | Jump table must be non-negative length |
| v6: deepest node in bounds | Offset must be within storage (if not sentinel) |
| Header SIZE = storage size | Detects mismatch / corruption |

---

## Quality Assurance

As of 2026-02-24 (`dotnet test tests/SuffixTree/SuffixTree.Persistent.Tests/SuffixTree.Persistent.Tests.csproj -c Release`):

- **471 tests passed**, 0 failed, 0 skipped.

Test suites cover core API behavior, format/hybrid transitions, parity with in-memory tree,
safety/dispose/concurrency, serialization invariants, and load/contract validation.
See the maintained matrix:
[`tests/SuffixTree/SUFFIX_TREE_TEST_MATRIX.md`](../../../../tests/SuffixTree/SUFFIX_TREE_TEST_MATRIX.md).

---

## License

Part of the Seqeron bioinformatics suite. See [LICENSE](../../../../LICENSE).
