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
| Build | O(n) amortized | Ukkonen's online algorithm |
| Contains | O(m) in-memory; O(m log d) persistent | Child lookup is hash/inline vs binary search |
| FindAllOccurrences | O(m + k) in-memory; O(m log d + k) persistent | k = number of returned positions |
| CountOccurrences | O(m) in-memory; O(m log d) persistent | Leaf count is precomputed; no subtree DFS |
| LongestRepeatedSubstring | First call: O(\|LRS\|) in-memory; persistent O(h + \|LRS\|) typical, O(n + \|LRS\|) fallback; then O(1) cached | Persistent can fall back to deepest-node DFS if header metadata is unavailable |
| FindAllLongestRepeatedSubstrings | O(nodes · log σ + Σ (L + occ)) | One `Traverse`, then one `FindAllOccurrences` per deepest internal node (L = LRS length) |
| LongestCommonSubstring | O(m + h) in-memory; O(m log d + h) persistent | Streaming match + one leaf-position recovery |
| FindAllLongestCommonSubstrings | O(m + Σ subtree(best matches)); worst case O(n·m) | Collects leaves for each maximal match candidate |
| FindExactMatchAnchors | O(m + a·h); worst case O(n·m) | a = anchors emitted, each needs leaf-position recovery |
| FindMaximalExactMatches | O(n + m + R) | R = right-maximal matches ≥ minLength (leaves visited; ≥ output size) — MUMmer 3 bound |
| FindMaximalUniqueMatches | O(n + m + k log k) | k = MUM candidates (≤ m) |
| EnumerateSuffixes | O(n²) total | Lazy DFS, O(n) per suffix |
| GetAllSuffixes | O(n²) | Materialized sorted list |

Where: n = text length, m = pattern/query length, k = result count, d = max branching factor on visited nodes, h = max depth walked to recover a leaf position, a = anchor count.

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
  equals the LRS length (O(1) query). Persistent format stores deepest-node
  offset and precomputed LRS depth in the v6 header (bytes 72 and 80).
- **Suffix link validation** — in DEBUG builds, the in-memory backend checks
  that internal non-root nodes have valid suffix links.

The persistent builder additionally handles **compact→large transitions**
mid-build: when the next allocation would exceed `0xFFFFFFFE`, it switches
from 24-byte (compact) to 32-byte (large) nodes. A jump table bridges
cross-zone suffix links and child arrays (see
[Persistent README](../../../src/SuffixTree/Algorithms/SuffixTree.Persistent/README.md)).

Text persistence is lossless at the code-unit level: the persistent file (non-ASCII
texts), the `SuffixTreeSerializer` v2 payload and the structural hash store raw
UTF-16LE code units, so any `string` accepted by the in-memory tree — including
lone surrogates and surrogate pairs straddling the 4096-char write chunks — reloads
and re-imports unchanged (review 2026-09 F17). Empty texts round-trip through MMF files (F16).

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

### 4.4 LongestRepeatedSubstring — first call non-constant, then O(1) cached

The deepest internal node (by `DepthFromRoot + edgeLength`) is found
during construction. Its total depth is the LRS length.
- In-memory: cached in private fields.
- Persistent: deepest-node offset + LRS depth are stored in the v6 header (72/80), then cached;
  if metadata is unavailable, a fallback DFS is used once.
- Ties (several distinct repeats of maximal length): the representative is unspecified and
  implementation-specific (in-memory and persistent may return different, equally long repeats,
  e.g. `baab` → `a` vs `b`); the length is always the maximum (brute-force verified, 2026-09 review).
  `LongestRepeatedSubstring` itself is unchanged; for **all ties** use §4.4.1.

#### 4.4.1 FindAllLongestRepeatedSubstrings — every tie, every position

`FindAllLongestRepeatedSubstrings()` (on `ISuffixTree` as a default interface member, overridden
by both `SuffixTree` and `PersistentSuffixTree`; implemented once in
`SuffixTreeAlgorithms.FindAllLongestRepeatedSubstrings(ISuffixTree)`) returns every distinct
longest repeated substring with all its 0-based start positions, ascending (occurrences may
overlap), ordered by first occurrence; empty when no character repeats. Method: a substring
occurs ≥ 2 times iff its locus is at or above an internal node (Gusfield 1997 §7.1), so the
longest repeats are exactly the path labels of the internal nodes of maximal string depth — one
`Traverse` collects them (distinct nodes ⇒ distinct labels), and each label's leaves are its
occurrences (`FindAllOccurrences`, sorted). The output is identical for every tree implementation
(heap / hybrid / memory-mapped / reloaded persistent parity tests). Examples: `abcxbcaxcab` →
(`ab`, [0, 9]), (`bc`, [1, 4]), (`ca`, [5, 8]); `xyzzyxxz` → (`x`, [0, 5, 6]), (`y`, [1, 4]),
(`z`, [2, 3, 7]); `banana` → (`ana`, [1, 3]). Verified against an O(n²) brute force on 400 random
texts (ASCII and non-ASCII, n ≤ 120) plus Python-checked literals.

### 4.5 LongestCommonSubstring — O(m + h) in-memory; O(m log d + h) persistent

Uses `SuffixTreeAlgorithms.FindAllLcs<TNode, TNav>`:

Walk the query string character-by-character against the tree. On mismatch,
follow suffix links and rescan to maintain the longest match. Track the
best match length and positions throughout.

Variants:
- `LongestCommonSubstring(other)` → string only
- `LongestCommonSubstringInfo(other)` → `(string, posInText, posInOther)`;
  returns `("", -1, -1)` if none. Tie-break: the substring whose occurrence in `other` comes
  first; `posInOther` is that first occurrence, `posInText` is *one* occurrence in the text
  (any-leaf walk — not necessarily leftmost, may differ between in-memory and persistent).
- `FindAllLongestCommonSubstrings(other)` → canonical LCS string (same tie-break) + all its
  positions in both strings, each list ascending and duplicate-free. Other substrings of the
  same maximal length are not part of this result — use the next variant.
- `FindAllDistinctLongestCommonSubstrings(other)` → **every** distinct LCS string (all length
  ties), each with all positions in the text and in `other` (ascending, duplicate-free), ordered by
  first occurrence in `other` (entry 0 = `FindAllLongestCommonSubstrings`). With L = max ms(i),
  every occurrence in `other` of a length-L common substring ends exactly at an i with ms(i) = L,
  so grouping those ends by substring yields all ties; text positions are the leaves below each
  substring's locus. Default interface body for external implementers: O(n·m) DP over `Text`
  (identical output). Validated against a Python brute force (1200 random pairs, 382 with ties)
  and in-memory = persistent.

### 4.6 FindExactMatchAnchors — O(m + a·h), worst case O(n·m)

Uses `SuffixTreeAlgorithms.FindExactMatchAnchors<TNode>`:

Same suffix-link streaming as LCS (matching statistics `ms(i)`), but with **peak
tracking**: for every maximal run of query end positions with `ms(i) ≥ minLength`, one
anchor — the first peak of the run — is emitted when the run ends. Each anchor is a
maximal exact match (MEM; Kurtz et al. 2004), with one arbitrary text occurrence.

Declared differences from MUMmer (verified by brute force, 2026-09 review):
- a **subset** of the MEMs ≥ `minLength` (one per run), not all MEMs and not MUMs
  (no uniqueness test, Delcher et al. 1999);
- anchors of adjacent runs can **overlap** in the query by < `minLength` characters:
  text `aba`, query `ababa`, `minLength = 3` → `(0,0,3)`, `(0,2,3)` (ms = 1,2,3,2,3);
- `minLength ≤ 0` returns an empty list.

For the **complete** MEM set (every MEM, every text occurrence) or MUMs use §4.9
(`FindMaximalExactMatches` / `FindMaximalUniqueMatches`). `FindExactMatchAnchors` is kept
unchanged as the lightweight anchor seed of `AnchorBasedAligner`; its anchors are always a
subset of `FindMaximalExactMatches` (test `FindExactMatchAnchors_IsSubsetOfFullMemSet`).

### 4.7 Suffix Enumeration

DFS traversal in sorted child-key order (ascending). Concatenates edge
labels to produce suffixes.
- `EnumerateSuffixes()` — lazy `IEnumerable<string>`, avoids O(n²) peak memory.
- `GetAllSuffixes()` — materialized `IReadOnlyList<string>`.

### 4.8 Traverse (Visitor Pattern)

`Traverse(ISuffixTreeVisitor)` — deterministic DFS in sorted key order.
Calls `VisitNode`, `EnterBranch`, `ExitBranch` for each node/edge.
Used by `SuffixTreeSerializer` for structural hashing.

### 4.9 Maximal exact / unique matches (MUMmer 3) — O(n + m + R)

`FindMaximalExactMatches(query, minLength)` and
`FindMaximalUniqueMatches(query, minLength, MumUniqueness)` (shared code in
`SuffixTreeAlgorithms`, identical output from `SuffixTree` and `PersistentSuffixTree`).
The tree's text is the **reference**; all coordinates are 0-based; forward strand only.

**Definitions** (Kurtz et al. 2004; Delcher et al. 1999):

| Set | Triple (r, q, len), len ≥ minLength, text[r..r+len) = query[q..q+len), and … | MUMmer 3 |
|-----|------|------|
| MEM | left-maximal (r = 0, q = 0 or text[r−1] ≠ query[q−1]) and right-maximal (a string ends or text[r+len] ≠ query[q+len]); **every** reference occurrence | `mummer -maxmatch -l L` |
| MUM, `Reference` | a MEM whose string occurs exactly once in the reference (may repeat in the query) | `mummer -mumreference -l L` (MUMmer default; `-mumcand`; "MAM" in MUMmer 4) |
| MUM, `Both` (default) | a MEM whose string occurs exactly once in the reference **and** exactly once in the query | `mummer -mum -l L` |

**Algorithm.** The query is streamed through the tree with suffix links (matching statistics,
Chang & Lawler 1994), keeping for each query start q the locus of the longest prefix of
query[q..] found in the text (length ms(q)) and — as MUMmer 3 `findmaxmat.c` — the locus of its
prefix of length `minLength`.
- *MEM*: the subtree below the `minLength` locus holds exactly the text suffixes matching ≥
  `minLength` characters. Walking down the matching path, a leaf in a sibling subtree that
  leaves the path at string depth d matches exactly d characters, and a leaf below the ms(q)
  locus matches ms(q) — so each leaf is a right-maximal match with known length; a constant-time
  character test keeps the left-maximal ones. Cost O(n + m + R), R = number of right-maximal
  matches ≥ `minLength` (the same bound as MUMmer 3; every path node visited has a sibling leaf
  that is counted in R).
- *MUM `Reference`* (MUMmer 3 `findmumcand.c`, `checkiflocationisMUMcand`): the ms(q) locus lies
  inside a leaf edge (string unique in the reference), ms(q) ≥ `minLength`, left-maximal.
- *MUM `Both`* (MUMmer 3 `cleanMUMcand.c`, `mumuniqueinquery`): the `Reference` candidates sorted by
  reference start (longer first) are swept; a candidate whose reference interval ends at or before
  the rightmost end seen so far lies inside another candidate (its string recurs in the query)
  and is dropped, equal candidates are dropped together. The sweep starts at −1; MUMmer 3 starts
  it at 0, which additionally drops a length-1 MUM at reference position 0 (only reachable with
  `-l 1`) — that artefact is not reproduced.

**Output order:** ascending by query position, then reference position. (MUMmer prints `-maxmatch`
and `-mumreference` by query position, `-mum` by reference position; the sets are identical.)

**Guards:** `query == null` → `ArgumentNullException`; `minLength < 1` or an undefined
`MumUniqueness` → `ArgumentOutOfRangeException`; empty text/query or `minLength` above either
length → empty list.

**Validation (2026-09-30).** MUMmer 3.23 (Ubuntu `mummer 3.23+dfsg-8`), commands
`mummer -maxmatch|-mum|-mumreference -l L ref.fa qry.fa` (no `-b`/`-r`, so forward strand):
66 cases (6 classic + 60 seeded random pairs, lengths 50–2000, alphabets {AC, ACG, ACGT},
queries partly built from mutated reference fragments, L ∈ {2…20}) — 270 366 MEMs, 5 695
`-mumreference` and 2 573 `-mum` matches, all identical; plus 300 short random pairs with
L ∈ {1,2,3} identical except the documented length-1 `-mum` sweep artefact (2 cases). Every
set also equals an independent brute-force implementation of the definitions (324 pairs).
MUMmer-produced outputs are locked in `MaximalMatchTests` / `MaximalMatchParityTests`.

### 4.10 Maximal repeated pairs (MUMmer `repeat-match -f`) — O(n + z) enumeration

`FindMaximalRepeatedPairs(minLength, isUniqueSymbol = null)` (default member of `ISuffixTree`,
shared code `SuffixTreeAlgorithms.FindMaximalRepeatedPairs` over `Traverse`, so `SuffixTree` and
`PersistentSuffixTree` return identical lists).

**Definition** (Gusfield 1997 §7.12): a triple (i, j, L), 0-based, i < j, L ≥ `minLength`, with
text[i..i+L) = text[j..j+L), right-maximal (j + L = n or text[i+L] ≠ text[j+L]) and left-maximal
(i = 0 or text[i−1] ≠ text[j−1]). Copies may overlap (tandem repeats: `acgtacgtacgt`, L 3 →
(0, 4, 8), (0, 8, 4)). A pair (i, j) has exactly one maximal length. This is MUMmer 3
`repeat-match -f -n minLength` (forward strand; repeat-match prints 1-based positions in tree order).

**Unique symbols.** Characters for which `isUniqueSymbol` returns true match nothing, not even
themselves (separators, `N`): a match stops before them (right-maximal) and a preceding unique
character makes a pair left-maximal — the definition applied to the text with every such
occurrence replaced by a fresh symbol (as `RepeatFinder.FindDirectRepeats`, non-ACGT unique, and
Vmatch/REPuter separators). Without a predicate every character matches itself, like
repeat-match (which also pairs `N` with `N`).

**Algorithm** (Gusfield 1997 §7.12.3). Bottom-up over the tree, each node keeps its leaves in one
linked list per left character (plus a "unique" class for position 0 / a unique left neighbour).
Merging a child into its parent v of string depth d ≥ `minLength` pairs every leaf of the child
with every leaf already merged into v whose left class differs (or is unique): their LCP is exactly
d and the pair is left-maximal; then the lists are concatenated in O(1) per class (smaller class
map into the larger). With unique symbols, a subtree whose edge contains the first unique
character at string depth u (u(p) = distance from p to the next unique character) is "exploded":
all of its leaves pair with each other at length u and nothing deeper is emitted. Every maximal
pair is produced exactly once. O(n log σ′ + z) enumeration for z pairs plus O(z log z) for the sort.

**Output order:** ascending FirstPosition, then SecondPosition.

**Validation (2026-10-01).** MUMmer 3.23 `repeat-match -f -n L` (Ubuntu `mummer 3.23+dfsg-8`; source
`src/tigr/repeat-match.cc` from the Ubuntu orig tarball, `List_Maximal_Matches` / `List_Matches`):
40 random genomes 2 kb–200 kb ({AC, ACG, ACGT}, planted repeats, L 8–22) → 512 083 pairs identical;
239 small random texts (L 1–6) → 32 283 pairs identical. B04 `RepeatFinder.FindDirectRepeats(seq, L,
int.MaxValue, int.MinValue)` (independent suffix-array/LCP-interval enumeration, non-ACGT unique):
20 sequences 1–50 kb with 2 % N → 381 139 pairs identical, plus 265 small cases. O(n³) brute force
(with and without unique symbols): 1200 random texts, 124 565 pairs identical. In-memory =
persistent (heap, hybrid, MMF, reload) on every run.

### 4.11 Longest common substring of k strings / k-common substring — O(N log N)

`SuffixTree.FindLongestCommonSubstrings(texts, minSupport = k)` and
`SuffixTree.LongestCommonSubstringLengthsBySupport(texts)` (shared code
`SuffixTreeAlgorithms.FindLongestCommonSubstrings(texts, minSupport, buildTree)` /
`LongestCommonSubstringLengthsBySupport(texts, buildTree)`, so any tree implementation can be used,
e.g. `s => PersistentSuffixTreeFactory.CreatePersistent(new StringTextSource(s))`; the built tree is
disposed after use).

Generalized suffix tree (Gusfield 1997 §7.6): the texts are concatenated as t₁$₁…t_k$_k with k
distinct separator characters absent from every text (`BuildGeneralizedText`: Private Use Area
U+E000 upward, then any other unused code unit). Each separator occurs once, so no internal node's
path label contains one. C(v), the number of distinct texts with a leaf below v, is computed with
Hui's (1992) colour-set-size method in one depth-first pass: coloured leaves minus one per pair of
consecutive same-text leaves, charged to their lowest common ancestor (found by binary search on the
DFS stack by entry time). l(q) = max string depth of an internal node with C(v) ≥ q; the answer for
`minSupport` q ≥ 2 is the set of path labels of internal nodes of depth l(q) with C(v) ≥ q (a
length-l(q) substring with support ≥ q ending inside an edge would make the node below it deeper with
the same support). `minSupport` = 1 → the longest text(s). Results distinct, sorted ordinally;
empty when nothing reaches the support (e.g. an empty text with q = k).

**Validation.** Rosalind LCSM sample (GATTACA, TAGACCA, ATACA → sample answer `AC`; all longest:
AC, CA, TA; q = 2 → TACA; l = [7, 4, 2]); Python brute force over all substrings: 800 random sets
(k 1–6, lengths 0–25), every q; in-memory = persistent.

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
IReadOnlyList<int> FindAllOccurrences(string / ReadOnlySpan<char>)
int CountOccurrences(string / ReadOnlySpan<char>)
string LongestRepeatedSubstring()
ReadOnlyMemory<char> LongestRepeatedSubstringMemory()
IReadOnlyList<(string Substring, IReadOnlyList<int> Positions)> FindAllLongestRepeatedSubstrings()   // default member
IEnumerable<string> EnumerateSuffixes()
IReadOnlyList<string> GetAllSuffixes()
string LongestCommonSubstring(string / ReadOnlySpan<char>)
(string, int, int) LongestCommonSubstringInfo(string)
(string, IReadOnlyList<int>, IReadOnlyList<int>) FindAllLongestCommonSubstrings(string)
IReadOnlyList<(string, IReadOnlyList<int>, IReadOnlyList<int>)> FindAllDistinctLongestCommonSubstrings(string)   // default member
IReadOnlyList<(int, int, int)> FindMaximalRepeatedPairs(int, Func<char, bool>? = null)   // default member
string PrintTree()
void Traverse(ISuffixTreeVisitor)
IReadOnlyList<(int, int, int)> FindExactMatchAnchors(string, int)
IReadOnlyList<(int, int, int)> FindMaximalExactMatches(string, int)
IReadOnlyList<(int, int, int)> FindMaximalUniqueMatches(string, int, MumUniqueness = Both)
```

### ISuffixTreeNavigator\<TNode\>

```
ITextSource Text
TNode Root
TNode NullNode
bool IsNull(TNode)
bool IsRoot(TNode)
int GetEdgeSymbol(TNode, int)
int LengthOf(TNode)
TNode GetSuffixLink(TNode)
bool TryGetChild(TNode, int, out TNode)
void CollectLeaves(TNode, int, List<int>)
int FindAnyLeafPosition(TNode, int)
void GetChildren(TNode, List<TNode>)
```

### ITextSource

```
int Length
char this[int]
string Substring(int, int)
ReadOnlySpan<char> Slice(int, int)
```

Implementations: `StringTextSource` (wraps `string`),
`MemoryMappedTextSource` and `AsciiMemoryMappedTextSource` (MMF-backed).

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
var positions     = tree.FindAllOccurrences("ana"); // IReadOnlyList<int>
int count         = tree.CountOccurrences("ana");

// Algorithms
string lrs        = tree.LongestRepeatedSubstring();
string lcs        = tree.LongestCommonSubstring("bandana");
var (s, p1, p2)   = tree.LongestCommonSubstringInfo("bandana");
var allLcs        = tree.FindAllLongestCommonSubstrings("bandana");
var anchors       = tree.FindExactMatchAnchors("bandana", minLength: 3);
var mems          = tree.FindMaximalExactMatches("bandana", minLength: 3);   // mummer -maxmatch
var mums          = tree.FindMaximalUniqueMatches("bandana", 3, MumUniqueness.Both); // mummer -mum
var ties          = tree.FindAllDistinctLongestCommonSubstrings("bandana"); // every tied LCS
var pairs         = tree.FindMaximalRepeatedPairs(3);                      // repeat-match -f -n 3
var pairsN        = tree.FindMaximalRepeatedPairs(3, c => c == 'N');       // N never matches
var lcsK          = SuffixTree.FindLongestCommonSubstrings(new[] { "GATTACA", "TAGACCA", "ATACA" }); // AC, CA, TA
var lBySupport    = SuffixTree.LongestCommonSubstringLengthsBySupport(texts); // Gusfield l(q)

// Enumeration
var suffixes      = tree.GetAllSuffixes();          // IReadOnlyList<string>
var lazySuffixes  = tree.EnumerateSuffixes();        // IEnumerable<string>

// Diagnostics
string viz        = tree.PrintTree();
tree.Traverse(visitor);
```

### Persistent (MMF-Backed) Tree

```csharp
using SuffixTree.Persistent;

// Build into MMF file (hybrid v6)
using var tree = (IDisposable)PersistentSuffixTreeFactory.Create(
    new StringTextSource("banana"), "tree.dat");
var st = (ISuffixTree)tree;

// Build in heap memory (no file)
using var tree = (IDisposable)PersistentSuffixTreeFactory.Create(
    new StringTextSource("banana"));

// Load existing file (read-only, v6 with compact/large base auto-detection)
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
│   ├── MumUniqueness.cs          ← MUM mode (-mum / -mumreference)
│   └── SuffixTreeAlgorithms.cs   ← LCS, Anchors, MEM/MUM (generic, JIT-specialized)
├── SuffixTree/                   ← In-memory implementation (.NET 8)
│   ├── SuffixTree.cs             ← Build, factory methods
│   ├── SuffixTree.Construction.cs ← Ukkonen's algorithm
│   ├── SuffixTree.Search.cs      ← Contains, FindAll, Count
│   ├── SuffixTree.Algorithms.cs  ← LRS, LCS, Anchors, MEM/MUM
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
    ├── AsciiMemoryMappedTextSource.cs
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

The `SuffixTree.Mcp.Core` project exposes 18 tools via Model Context Protocol:

| Tool | Method |
|------|--------|
| `suffix_tree_contains` | Pattern existence check |
| `suffix_tree_count` | Occurrence count |
| `suffix_tree_find_all` | All positions |
| `suffix_tree_lrs` | Longest repeated substring |
| `suffix_tree_lcs` | Longest common substring |
| `suffix_tree_stats` | Tree statistics |
| `suffix_tree_all_lrs` | All tied longest repeated substrings with positions |
| `suffix_tree_find_mems` | Maximal exact matches (MUMmer `-maxmatch`) |
| `suffix_tree_find_mums` | Maximal unique matches (MUMmer `-mum` / `-mumreference`) |
| `suffix_tree_maximal_repeats` | Maximal repeated pairs (MUMmer `repeat-match -f`; optional unique symbols) |
| `suffix_tree_all_lcs` | All tied longest common substrings with positions in both texts |
| `suffix_tree_k_common_substrings` | Longest common substring of k strings / ≥ minSupport (Rosalind LCSM) |
| `find_longest_repeat` | DNA longest tandem repeat |
| `find_longest_common_region` | DNA common region |
| `calculate_similarity` | K-mer Jaccard similarity |
| `hamming_distance` | Hamming distance |
| `edit_distance` | Levenshtein edit distance |
| `count_approximate_occurrences` | Approximate pattern matching |

Each tool validates input and returns a structured result record.
See [MCP docs](../../mcp/README.md) for connection and usage guides.

---

## 9. Tests

As of 2026-02-24 (`dotnet test` on the three SuffixTree test projects):

- `SuffixTree.Tests`: 353 passed
- `SuffixTree.Persistent.Tests`: 471 passed
- `SuffixTree.Mcp.Core.Tests`: 59 passed
- **Total: 883 passed**

For detailed, maintained coverage matrix by suite/scenario, see
[`tests/SuffixTree/SUFFIX_TREE_TEST_MATRIX.md`](../../../tests/SuffixTree/SUFFIX_TREE_TEST_MATRIX.md).

---

## 10. Benchmarks

The `SuffixTree.Benchmarks` project (BenchmarkDotNet) includes scenarios for:

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
- Kurtz, S. et al. (2004). *Versatile and open software for comparing large genomes.* Genome Biology 5:R12 (MUMmer 3 — MEMs, MUMs, MUM-candidates). Source consulted: MUMmer 3.23 `src/kurtz/mm3src/findmaxmat.c`, `findmumcand.c`, `libbasedir/cleanMUMcand.c` (Ubuntu source package `mummer 3.23+dfsg`); MUMmer 4 `include/mummer/sparseSA.hpp` (`findMAM_each`, `findMUM_each`, `collectMEMs_each`).
- Khan, Z., Bloom, J.S., Kruglyak, L., Singh, M. (2009). *A practical algorithm for finding maximal exact matches in large sequence datasets using sparse suffix arrays.* Bioinformatics 25:1609–1616 (sparseMEM; MEM/MAM/MUM definitions reused by MUMmer 4).
- Gusfield, D. (1997), §7.6 (k-common substring problem, generalized suffix tree) and §7.12 / §7.12.3 (maximal repeated pairs, left-character lists, O(n + z)).
- Hui, L.C.K. (1992). *Color set size problem with applications to string matching.* CPM 1992, LNCS 644:230–243 (C(v) by LCA of consecutive same-colour leaves; as described in Gusfield 1997 §7.6 / §9.7 — the paper itself was not opened).
- MUMmer 3.23 `src/tigr/repeat-match.cc` (A. Delcher; `List_Maximal_Matches`, `List_Matches`, `-f`, `-n`) — Ubuntu orig tarball `mummer_3.23+dfsg.orig.tar.xz`.
- Rosalind LCSM "Finding a Shared Motif" (sample dataset GATTACA / TAGACCA / ATACA → `AC`; rosalind.info blocked, confirmed via WebSearch of solution repositories).
- Chang, W.I., Lawler, E.L. (1994). *Sublinear approximate string matching and biological applications.* Algorithmica 12:327–344 (matching statistics).
- https://visualgo.net/en/suffixtree

---

## 12. Related Documentation

- Persistent format details: [SuffixTree.Persistent/README.md](../../../src/SuffixTree/Algorithms/SuffixTree.Persistent/README.md)
- Exact pattern search: [Exact_Pattern_Search.md](Exact_Pattern_Search.md)
- MCP tool docs: [docs/mcp/](../../mcp/README.md)
