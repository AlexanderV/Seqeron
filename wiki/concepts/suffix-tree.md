---
type: concept
title: "Suffix tree (Ukkonen construction + dual in-memory / persistent backends)"
tags: [pattern-matching, algorithm, sequence-comparison]
mcp_tools:
  - suffix_tree_contains
  - suffix_tree_count
  - suffix_tree_find_all
  - suffix_tree_lrs
  - suffix_tree_lcs
  - suffix_tree_stats
sources:
  - docs/algorithms/Pattern_Matching/Suffix_Tree.md
source_commit: fdb2411237de02e2dfdb4ef7ed2c0d6a77cb52c6
created: 2026-07-15
updated: 2026-07-15
graph:
  relationships:
    - predicate: relates_to
      object: concept:exact-pattern-search
      source: suffix-tree
      evidence: "Suffix_Tree.md §4.1–4.3 — Contains/FindAllOccurrences/CountOccurrences (the exact-pattern-search PAT-EXACT-001 primitives) are the search algorithms provided by this suffix-tree data structure via ISuffixTree."
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:longest-repeated-substring
      source: suffix-tree
      evidence: "Suffix_Tree.md §2/§4.4 — LongestRepeatedSubstring is the deepest-internal-node query on this tree, computed during construction (deepest node identified in the bottom-up pass, O(1) cached); the single-string repeat family drives this tree."
      confidence: high
      status: current
---

# Suffix tree (Ukkonen construction + dual in-memory / persistent backends)

A **suffix tree** is a compressed trie of *every* suffix of a string `T`, with a
terminal sentinel so each suffix ends in an explicit leaf. It is the shared **index
data structure** beneath Seqeron's pattern matching, repeat detection, and
sequence-comparison work: once built in `O(n)`, it answers substring queries in `O(m)`
(pattern length), independent of `n`. This page is the **data structure itself**
(construction, storage, interface). The *query applications* that ride on it are
separate units:

- exact search (`Contains` / `CountOccurrences` / `FindAllOccurrences`) → [[exact-pattern-search]] (PAT-EXACT-001);
- single-string repeats (deepest node with ≥ 2 leaves) → [[longest-repeated-substring]] (GENOMIC-REPEAT-001);
- two-string common substring (generalized tree) → [[longest-common-substring]] (GENOMIC-COMMON-001);
- word/anchor engine for dot-plots → [[dot-plot-word-match]].

The implementation follows **Ukkonen's online construction (1995)**, `O(n)` amortized.
The validation/evidence pattern for the algorithms it hosts is
[[algorithm-validation-evidence]].

## Two interchangeable backends, one interface

Both backends implement a single `ISuffixTree` interface, so callers are backend-blind:

| Backend | Project | Target | Backing store |
|---------|---------|--------|---------------|
| **In-memory** | `SuffixTree` | .NET 8 | heap-allocated node objects; hybrid child storage (inline ≤ 4 children → `Dictionary`) |
| **Persistent** | `SuffixTree.Persistent` | .NET 9 | memory-mapped file (MMF), hybrid **v6** format; also builds heap-only with no file |

Shared query algorithms (LCS, exact-match anchors) are written **once** against a
generic `ISuffixTreeNavigator<TNode>` with a **`struct` constraint**, so the JIT
**specializes** the same code for each backend's node type with no virtual-dispatch
cost. `SuffixTree.Core` holds the interfaces + shared algorithms; the in-memory and
persistent projects both depend up onto it.

## Ukkonen construction (both backends, same logic)

1. Process characters **left to right**, extending the tree by one character per phase.
2. Maintain an **active point** — `(active node, active edge, active length)` — plus a
   **remainder** counter of suffixes still owed.
3. Three **extension rules** per phase:
   - **Rule 1** — implicit leaf extension: open-ended leaf edges grow automatically (a
     leaf stays a leaf), so no per-character work is needed for them.
   - **Rule 2** — create a new leaf, splitting an edge if the next character diverges
     mid-edge (introduces a new internal node).
   - **Rule 3** — the character is already present on the path → **showstopper**: end the
     phase, increment remainder, advance active length.
4. **Suffix links** connect the node for `xα` to the node for `α`, giving amortized
   `O(1)` jumps between successive suffix insertions — the trick that makes the whole
   build linear.
5. A special **terminator** (integer key `−1`) is appended so every suffix becomes an
   **explicit leaf** (no suffix is a prefix of another), which is what makes leaf-based
   counting exact.

### Post-construction passes

Run once, right after the build, and cached:

- **Bottom-up leaf counting** — iterative post-order traversal writes a `LeafCount` to
  every node. `CountOccurrences` then reads it in `O(1)` (no subtree DFS) — this is why
  [[exact-pattern-search]]'s count primitive is `O(m)`.
- **Deepest internal node** — identified in the same pass. Its string depth
  (`DepthFromRoot + edgeLength`) *is* the [[longest-repeated-substring]] length, so LRS
  is an `O(1)` cached lookup after the first call.
- **Suffix-link validation** — DEBUG-only: the in-memory backend asserts that internal
  non-root nodes have valid suffix links.

The persistent v6 header **precomputes and stores** the deepest-node offset and LRS depth
(header bytes 72 and 80), so a loaded file answers LRS without re-walking; if that
metadata is missing it falls back to a one-time deepest-node DFS. The persistent builder
also handles a **compact → large node transition mid-build**: when the next allocation
would exceed `0xFFFFFFFE`, it switches from **24-byte compact** to **32-byte large**
nodes, bridging cross-zone suffix links and child arrays with a jump table.

## Complexity

`n` = text length, `m` = query length, `k` = result count, `d` = branching factor at
visited nodes, `h` = depth walked to recover a leaf position, `a` = anchors emitted.
Persistent child lookup is **binary search** (`log d`) vs the in-memory **hash/inline**
(`O(1)`), which is the only asymptotic gap between the backends:

| Operation | In-memory | Persistent |
|-----------|-----------|------------|
| Build | `O(n)` amortized | `O(n)` amortized |
| `Contains` | `O(m)` | `O(m log d)` |
| `FindAllOccurrences` | `O(m + k)` | `O(m log d + k)` |
| `CountOccurrences` | `O(m)` (precomputed `LeafCount`) | `O(m log d)` |
| `LongestRepeatedSubstring` | `O(\|LRS\|)` first call, then `O(1)` | `O(h + \|LRS\|)` typ., `O(n + \|LRS\|)` fallback, then `O(1)` |
| `LongestCommonSubstring` | `O(m + h)` | `O(m log d + h)` |
| `EnumerateSuffixes` / `GetAllSuffixes` | `O(n²)` total | `O(n²)` total |

`FindAllLongestCommonSubstrings` and `FindExactMatchAnchors` are `O(n·m)` worst case
(they collect leaves per maximal-match candidate).

## Algorithm surface hosted on the tree

Beyond the exact-search trio (own unit [[exact-pattern-search]]) and LRS
([[longest-repeated-substring]]), the tree exposes:

- **`LongestCommonSubstring` / `LongestCommonSubstringInfo` / `FindAllLongestCommonSubstrings`**
  — walk the query char-by-char; on mismatch, follow suffix links and rescan to keep the
  longest running match (`SuffixTreeAlgorithms.FindAllLcs`). Surfaces as
  [[longest-common-substring]].
- **`FindExactMatchAnchors(query, minLength)`** — the same suffix-link streaming plus
  **peak tracking**: emit right-maximal, non-overlapping matches when the match length
  falls back below `minLength` after exceeding it. These are **MEMs** for
  alignment-chaining (MUMmer / LAGAN style).
- **Suffix enumeration** — DFS in ascending child-key order concatenating edge labels;
  `EnumerateSuffixes()` is lazy (avoids the `O(n²)` peak), `GetAllSuffixes()` materializes
  a sorted list.
- **`Traverse(ISuffixTreeVisitor)`** — deterministic DFS in sorted key order
  (`VisitNode` / `EnterBranch` / `ExitBranch`), used by the serializer for structural
  hashing.

`Contains` uses a **hybrid comparison**: `SequenceEqual` (SIMD) for edge fragments ≥ 8
characters, scalar loop otherwise, zero-allocation via `ReadOnlySpan<char>`.
`ITextSource` abstracts the backing text — `StringTextSource` (wraps a `string`) or the
MMF-backed `MemoryMappedTextSource` / `AsciiMemoryMappedTextSource`.

## Serialization

- **Logical hash** — `SuffixTreeSerializer.CalculateLogicalHash` produces a deterministic
  **SHA-256** structural hash (via the visitor traversal), so two trees can be compared
  for structural equality regardless of backend.
- **Stream export/import** — format **v2**, `Export`/`Import` against any
  `IStorageProvider` (e.g. `HeapStorageProvider`).
- **File export/import** — `SaveToFile` / `LoadFromFile` produce/consume MMF-backed trees;
  `PersistentSuffixTreeFactory.Load` auto-detects compact vs large base in the v6 header.

## API entry points

```csharp
var tree = SuffixTree.Build("banana");              // string / ITextSource / ReadOnlyMemory / ReadOnlySpan
bool ok  = SuffixTree.TryBuild("text", out var t);  // non-throwing
var empty = SuffixTree.Empty;                        // singleton

// persistent (MMF file or heap-only), same ISuffixTree surface
using var p = (IDisposable)PersistentSuffixTreeFactory.Create(new StringTextSource("banana"), "tree.dat");
```

Source files: `SuffixTree.Construction.cs` (Ukkonen), `SuffixTree.Search.cs` (Contains /
FindAll / Count), `SuffixTree.Algorithms.cs` + `SuffixTreeAlgorithms.cs` (LRS / LCS /
anchors), `SuffixTree.Persistent/` (MMF backend). MCP: `SuffixTree.Mcp.Core` exposes 12
tools — `suffix_tree_contains` / `_count` / `_find_all` / `_lrs` / `_lcs` / `_stats`
directly wrap this structure; the remaining DNA-flavored tools (`find_longest_repeat`,
`find_longest_common_region`, `calculate_similarity`, `hamming_distance`, `edit_distance`,
`count_approximate_occurrences`) belong to their own concepts.

## Sources and deviations

Ukkonen (1995) *On-line construction of suffix trees*; Gusfield (1997) *Algorithms on
Strings, Trees, and Sequences* (ch. 5–7 application family); Delcher et al. (1999) MUMmer
(the suffix-tree anchor approach behind `FindExactMatchAnchors`). **No algorithm
deviations**; the design choices are engineering ones — the dual heap/MMF backends behind
one interface, the `struct`-constrained navigator for JIT specialization, the hybrid
inline/`Dictionary` child storage, the SIMD edge comparison, and the v6 compact→large
node transition — none of which change the classical construction or query semantics.

The engine was later hardened by a contract-freeze refactoring campaign — see
[[suffix-tree-controlled-refactoring]] for the frozen invariants (the persistent v6 binary
format, the `ISuffixTree` public surface, and the six `suffix_tree_*` MCP tool names) and the
per-cycle characterization gates that keep each refactor phase behaviour-preserving.
