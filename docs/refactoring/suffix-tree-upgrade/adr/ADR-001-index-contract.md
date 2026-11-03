# ADR-001: Index contract and how the ESA plugs into the existing generic algorithms

Status: **proposed** (work package S1-0, 2026-11-03). Authorization: [PROCESS.md](../PROCESS.md). Plan: [PLAN.md](../PLAN.md) stage S1. Companions: [ADR-002](ADR-002-sa-construction.md) (construction), [ADR-003](ADR-003-esa-container-format.md) (container).
Paths are repo-relative; `file:line` refers to the tree at commit `b0922125` (branch `claude/stoic-maxwell-0olr2z`). Anything not checked against the code or by a measurement is marked **UNVERIFIED**.

## 1. Decisions already taken (not reopened here)

Managed C# only (no native binary, no P/Invoke, no copied libsais code); positions are `long` in the new contract from day one; contract first (a stable index contract, replaceable engines: ESA now, compressed FM / r-index and online Ukkonen later); the Ukkonen tree and the persistent v6 code stay untouched and are the test oracle; the ESA container is a separate format kind; S1a gate on 10 Mbp random DNA, 1 thread: SA+LCP build <= 3.3 s, index <= 16 B/char, peak RSS during build <= 24 B/char (ceiling goal 1.0 s, aspiration 0.5 s; libsais 0.328 s) ([PROGRESS-S1.md](../PROGRESS-S1.md)).

## 2. What the existing code constrains

| Fact | Evidence |
|---|---|
| `ISuffixTree` = `ISuffixTreeSearch` + `ISuffixTreeAnalysis` + `ISuffixTreeDiagnostics`; five members have default implementations that run on `Traverse` (`CountDistinctSubstrings*`, `FindAllLongestRepeatedSubstrings`, `FindMaximalRepeatedPairs`) | `SuffixTree.Core/ISuffixTree.cs:44-80` |
| Every position, count and length is `int`; `ITextSource` has `int Length` and `char this[int]` | `ISuffixTreeSearch.cs:19-29,55,71`; `ITextSource.cs:15,22` |
| The generic algorithms return `int` triples / `List<int>` | `SuffixTreeAlgorithms.cs:36,101,180,281` |
| Terminated tree: the text is followed by a unique terminator, symbol `-1`; its key sorts first; a suffix that is a prefix of another is a leaf under key `-1`; the empty suffix is a leaf under the root | `SuffixTree/SuffixTree.cs:45`; `PersistentConstants.cs:29`; `SuffixTree.Diagnostics.cs:189-191` (keys sorted ascending); visitors treat `p >= n` as the "empty-suffix leaf" (`SuffixTreeAlgorithms.Repeats.cs:131-134`, `SuffixTreeAlgorithms.KCommon.cs:255-256`) |
| Reported statistics: `NodeCount` counts every node of the terminated tree (including the empty-suffix leaf), `LeafCount = leaves - 1 = n`, `MaxDepth = deepest leaf - 1 = n`; empty text gives `(1, 0, 0)` | `SuffixTree.Diagnostics.cs:9-76` (`:11`, `:69-75`) |
| `CountOccurrences` of the empty pattern returns `n`; `FindAllOccurrences("")` returns `[0..n-1]`; non-empty counts are `LeafCount` of the locus (terminator leaves included) | `PersistentSuffixTree.Search.cs:35-36,58-60` |
| Two independent plug points for a backend: (a) `ISuffixTreeNavigator<TNode>` consumed by `SuffixTreeAlgorithms.*<TNode,TNav>` (LCS family, anchors, MEM, MUM); (b) `ISuffixTreeDiagnostics.Traverse(ISuffixTreeVisitor)` consumed by LRS family, distinct-substring counts, maximal repeated pairs, k-common substrings | `ISuffixTreeNavigator.cs:14-70`; `SuffixTreeAlgorithms.cs:567,579,677`; `SuffixTreeAlgorithms.Repeats.cs:62`; `SuffixTreeAlgorithms.KCommon.cs:192` |
| Both existing navigators are `internal` structs; the algorithms and the interface are `public`, so a new assembly can call them with its own navigator | `SuffixTree.Navigator.cs:13`; `PersistentSuffixTreeNavigator.cs:11`; `SuffixTreeAlgorithms.cs:22` |

## 3. The main design risk: exactly what the generic algorithms call

Complete list of navigator calls in `SuffixTreeAlgorithms.cs` and `SuffixTreeAlgorithms.Repeats.cs` (counted with `grep "nav\."`; nothing else is called):

| `ISuffixTreeNavigator` member | Call sites | What the ESA must give |
|---|---|---|
| `Text` | `:107,187,290`, `Repeats.cs:384` (`Length`, and `text[r-1]` for left-maximality at `:382`) | `ITextSource` view of the text |
| `Root`, `NullNode`, `IsNull`, `IsRoot` | `:112-113,197-201,295`; `IsRoot` only at `:475` | trivial |
| `GetEdgeSymbol(node, offset)` | `:223`, `:303` (leaf test: last symbol of the edge `< 0`), `:432` | symbol `text[SA[lb] + parentDepth + offset]`, `-1` when the index is `>= n` |
| `LengthOf(node)` | 10 sites, e.g. `:213-214,232,417,437,453,495` | `depth - parentDepth`: **the handle must know the parent's string depth** |
| `GetSuffixLink(node)` | **one** site, `:477`, inside `FollowSuffixLinkAndRescan` (`:464-516`), reached from `:405` (matching statistics), `:246-248` (MEM, two loci), `:311` (MUM) | see section 5 |
| `TryGetChild(node, key)` | `:229` (MEM), `:447` (extend), `:493` (rescan) | child by first edge symbol (`-1` = terminator child) |
| `CollectLeaves(node, depth, list)` | `:70,225,236`, `Repeats.cs:302` | `SA[lb..rb]` minus the empty suffix; `depth` ignored |
| `FindAnyLeafPosition(node, depth)` | `:63,305,524` | `SA[lb]`; `depth` ignored |
| `GetChildren(node, list)` | `:220` (MEM only) | children in any order |

Findings:
1. **`parent` is never called.** Node depth is never read from the node either: the algorithms keep `currentNodeDepth` themselves (`:15-19`, `:478`) and rely on the invariant "depth(suffixLink(v)) = depth(v) - 1". `DepthFromRoot` arguments exist only because v6 Slim does not store depth (`ISuffixTreeNavigator.cs:57-60`); an ESA ignores them.
2. **The suffix link is needed at exactly one point, and only for internal nodes.** `currentNode` becomes `currentEdge` only when an edge is fully consumed (`:437-442`); a leaf edge ends with the terminator `-1`, which a query character never equals, so a leaf is never `currentNode`. Therefore `GetSuffixLink` is called on the root (skipped at `:475`) or on an internal node; its target is always internal or root.
3. **Number of suffix-link calls is bounded by the query**: matching statistics grow by at most 1 per query character and drop by 1 per link, so at most `|query|` links (`AdvanceMatchingStatistics`, `:392-407`); MEM follows two loci per query start (`:245-248`), so at most `2|query|`; MUM one per start (`:310-311`).
4. **Rescan needs the link target's edge length** (`LengthOf` of the link target and of its children, `:495`; `MatchLocus` `:413-418`, MEM `:213`). An ESA node handle must therefore carry the parent's string depth, computed when the link target is created (`max(LCP[lb], LCP[rb+1])`, the standard parent-depth formula).
5. **Operations that need no suffix link**: `Contains`, `CountOccurrences`, `FindAllOccurrences`, `LongestRepeatedSubstring`, `FindAllLongestRepeatedSubstrings`, `CountDistinctSubstrings[ByLength]`, `FindMaximalRepeatedPairs`, k-common substrings, `GetAllSuffixes`. **Operations that need it** (matching-statistics family): `LongestCommonSubstring*`, `FindAllLongestCommonSubstrings`, `FindAllDistinctLongestCommonSubstrings`, `FindExactMatchAnchors`, `FindMaximalExactMatches`, `FindMaximalUniqueMatches`.
6. **Reuse is therefore cheap**: the ESA needs one `EsaNavigator : ISuffixTreeNavigator<EsaNode>` (struct, internal) and a `Traverse`. No algorithm is rewritten; all are instantiated for `<EsaNode, EsaNavigator>` and specialised by the JIT like the existing two.

## 4. The ESA as a terminated suffix tree (model, so that differential tests can be node-for-node)

* Physical arrays are built over `text$` with `n + 1` entries: `SA[0] = n` (the empty suffix), `SA[1..n]` = the libsais-convention suffix array of `text`; `LCP[i] = lcp(SA[i-1], SA[i])`, `LCP[0] = LCP[1] = 0`. SA-IS produces this layout natively (virtual sentinel, [ADR-002](ADR-002-sa-construction.md) section 2), so no shift is needed; libsais equality compares `SA[1..n]`.
* Internal nodes = lcp-intervals `[lb..rb]` with value `l > 0`; the root = `[0..n]`, `l = 0`; leaves = singletons `[i..i]` with string depth `n - SA[i] + 1` (terminator counted). Because of `$`, the root interval always has `l = 0`, so no virtual-root special case exists, and the tree is isomorphic to the Ukkonen tree of the same text: `NodeCount = 1 + I + (n+1)` (`I` = number of lcp-intervals with `l > 0`; `n = 0` returns 1, `SuffixTree.Diagnostics.cs:11`), `LeafCount = n`, `MaxDepth = n`, `Count(locus) = rb - lb + 1`.
* A suffix equal to a node's path label is the first singleton of that node and acts as its terminator child (edge `$`, length 1, key `-1`): `GetEdgeSymbol(leaf, len-1) = text[SA[i] + depth(leaf) - 1] = text[n] -> -1` (index `>= n`), which is what the MUM leaf test at `:303` checks.
* Node handle (internal, per `int`/`long` instantiation): `{ Lb, Rb, Depth, ParentDepth, Pos }` with `Pos = SA[Lb]` cached at creation (one SA read, then `GetEdgeSymbol` and `FindAnyLeafPosition` are free). `NullNode` = `Lb < 0`. 20 bytes for the int32 instantiation.
* `Traverse` is an iterative top-down walk of the lcp-interval tree in SA order (children in ascending key order, terminator first), emitting the events of `SuffixTree.Diagnostics.cs:182-229`: `VisitNode(start, end, leafCount = rb-lb+1, childCount, parentDepth)`, `EnterBranch(firstSymbol)`, `ExitBranch()`, with `start = occ + parentDepth`, `end = occ + depth` for internal nodes and `end = -1` (BOUNDLESS, `SuffixTreeNode.cs:20`) for leaves, `occ` = the text position of the occurrence the Ukkonen builder recorded for the edge.
  **Requirement R-HASH.** `SuffixTreeSerializer.CalculateLogicalHash` hashes the text and, for every node, `startIndex, endIndex, leafCount, childCount`, every `EnterBranch` key and an exit sentinel (`SuffixTreeSerializer.cs:25-47,219-241`), and `Import` validates a rebuilt Ukkonen tree against that stored hash (`SuffixTree.Persistent/README.md`, "Both node count and structural hash are validated after rebuild"). An ESA exported through the serializer is therefore importable only if its event stream is identical to Ukkonen's, **including the edge positions of internal nodes**. The consumers of `Traverse` in the library only need a valid occurrence for internal nodes (`SuffixTreeAlgorithms.cs:683-686,717-724` slices the text at `end`; `Repeats.cs:131-132` and `KCommon.cs:255` use exact positions only for leaves, where `p = start - depth` is the suffix start), so identity is a serializer requirement, not an algorithm requirement.
  **Hypothesis H-OCC (UNVERIFIED, first experiment of WP 1b.1):** Ukkonen inserts suffixes in increasing start order and a split keeps the upper part's start, so every edge records the *leftmost* occurrence of its subtree: `occ = min(SA[lb..rb])`. ESA computes `occ` for internal nodes in a bottom-up pre-pass (one transient `int` per internal node; `Traverse` is a diagnostic path). If H-OCC holds, `CalculateLogicalHash(esa) == CalculateLogicalHash(ukkonen)` is the strongest single oracle for the whole tree structure. If it does not, derive the actual rule from the experiment; if none exists, use `SA[lb]`, document that serializer `Export`/`Import` is not supported for ESA trees, and tell the owner.

Cost of the navigator members on the ESA (random-access reads; "r" = one potential cache miss):

| member | ESA realisation | cost |
|---|---|---|
| `TryGetChild` | enumerate child boundaries with the child table (`up/down/nextlIndex`), compare the symbol at depth `l`: `text[SA[b] + l]` (2 r per candidate) or, with the optional `NEXTSYM` byte array (`text[SA[i]+LCP[i]]`, +1 B/char, [ADR-003](ADR-003-esa-container-format.md)), 1 r | `<= sigma'` candidates |
| `GetChildren` | same enumeration, handles created with `ParentDepth = node.Depth` | `sigma'` |
| `CollectLeaves` | contiguous copy of `SA[lb..rb]`, dropping the empty suffix; this removes the subtree walk that costs 7.6 M vs 71 M positions/s today (`BENCHMARK.md:190`) | sequential |
| `GetSuffixLink` | section 5 | 4-8 r |
| everything else | O(1) from the handle | none |

## 5. Suffix links on an ESA (the question to resolve)

A plain ESA stores no suffix links. Fact **F1** (Gusfield 1997 / Abouelhoda, Kurtz, Ohlebusch 2004): for an internal node `v = [l..r]` with label `aY`, `SL(v)` is the node with label `Y`, and
`SL(v) = LCA(Psi(l), Psi(r))`, `Psi(i) = ISA[SA[i]+1]`.
Proof sketch: `l` and `r` are the first and last leaf of `v`, in different children, so `lcp(SA[l], SA[r]) = depth(v)` exactly; hence the suffixes `SA[l]+1` and `SA[r]+1` have lcp exactly `depth(v) - 1`, `Psi` is monotone on `[l..r]` (a common first symbol does not change order), and the interval of `Y` is the lcp-interval of value `depth(v) - 1` that contains `[Psi(l), Psi(r)]`. (The paper's section number for the linear-time suffix-link construction is UNVERIFIED; the identity above is derived here and is itself a test target.)

Options:

| # | Option | Per-link cost | Extra memory | Worst case | Verdict |
|---|---|---|---|---|---|
| O1 | Stored suffix-link table | 1-3 r | 4-8 B per internal node (needs both bounds of the target: about 5 B/char on random DNA with 0.62 internal nodes/char from `BENCHMARK.md:189`, up to 8 B/char) | O(1) | too fat for the 16 B budget; build pass costs about one lazy link per node |
| O2 | **Lazy via F1**: read `Psi(l)`, `Psi(r)` through `ISA` (optional section, 4 B/char, derivable from SA in O(n)), then locate the `(depth-1)`-interval around them | 3 r for the endpoints + the interval search | ISA 4 B/char + search aid | depends on the interval search (below) | **recommended generator** |
| O3 | Re-descend from the root along `text[Pos+1 .. Pos+depth)` with skip/count | node-depth(SL) child lookups | none | O(n) per link on `a^n`-like texts, quadratic overall | fallback only, flagged `Degraded` |
| O4 | ESA-only locus-based matching statistics (step `Psi` on the interval, no explicit nodes, no rescan) | like O2 | like O2 | like O2 | new algorithm, not a reuse; keep as the S3 `ms[]` candidate |
| O5 | Weiner links / backward search | rank structure (BWT) | n..2n bits | O(1) | S7 (FM) territory, not S1 |

Interval search for O2 (the only open sub-question): given `a = Psi(l) < b = Psi(r)` and the known value `t = depth(v) - 1`:
* (a) **LCP outward scan** from `[a, b]` while `LCP >= t`: trivial, O(width of the target).
* (b) **Parent climb** from leaf `a` with the child table until depth `<= t`: exact (F1 guarantees it stops at depth `t`), O(node-depth of leaf `a` - node-depth of the target); quadratic overall on `a^n`.
* (c) **RMQ + PSV/NSV helper over LCP** (a block-sparse structure, about 0.3-0.5 B/char, **UNVERIFIED** estimate): `p = RMQ(a+1..b)` (value `t`), `L = PSV_strict(p)`, `R+1 = NSV_strict(p)`: O(1) with a sparse table over blocks, O(block) otherwise. The same structure provides the public `Lcp(i, j)` and `Lca(i, j)` of PLAN.md:74-76, so its cost is shared.

**Recommendation (R-SL):**
1. Treat the suffix link as an **optional "navigation pack"** (sections `ISA` + `RMQ`, [ADR-003](ADR-003-esa-container-format.md)), not as part of the core index. Core ESA (text + SA + LCP + child + buckets) answers count/locate/contains/LRS/repeats/Traverse without it. The matching-statistics family needs the pack; if it is absent the engine builds it lazily on first use (O(n), thread-safe once-only, heap-resident) and reports `SuffixLink = Degraded` until then (capability flag, section 6).
2. Implement `GetSuffixLink` as **O2 with search (c)**; use (a) in the first working version (WP 1b.4) because it is the simplest and is the oracle for (c). (b) is a prototype for comparison only.
3. Complexity: per link 3 random reads + O(1) interval work; total cost of a streaming operation `O(|query| * (links + extend steps))`, the same asymptotics as the suffix tree. Memory: `+4.0 (ISA) + ~0.4 (RMQ)` B/char, i.e. Core about 10.1 B/char (SA 4, LCP about 1.1, child 4, text 1) and Full about 14.5-15 B/char, both under 16 (budget table in [ADR-003](ADR-003-esa-container-format.md) section 5). Optional `Psi` array (4 B/char, one independent read per endpoint instead of two dependent ones) may replace ISA for this purpose; decided by P1.
4. A stored table (O1) is **not** adopted; it is the fallback only if P1 shows the lazy link dominating MEM time (decision rule below).

**Prototypes that decide the open points (work package 1b.4a, before the navigator is finalised):**

| id | Question | Inputs | Decision rule |
|---|---|---|---|
| P1 | link kernel: (a) vs (b) vs (c), ISA vs Psi | `dna_10m`, `rep_10m`, `prot_1m`, `a^n` (1 M), Fibonacci (1 M), Thue-Morse (1 M) | ns/link (median, p99, max) and steps; choose (c) unless (a) or (b) is within 1.2x of (c) on all inputs **and** bounded on `a^n`/Fibonacci; if (c) > 600 ns median on `dna_10m`, evaluate O1 |
| P2 | end-to-end: MEM, 100 kbp query with 1 % substitutions, minLength 20 (`BENCHMARK.md:172-176`) | `dna_10m`, `rep_10m` | proposal for the S1b acceptance: <= 0.20 s on `dna_10m` (v6 persistent: 0.163 s, in-memory Ukkonen: 0.0623 s first call), stretch <= 0.10 s; MUM/anchors proportional |
| P3 | child-table packing (Abouelhoda single-array scheme, packing details UNVERIFIED) vs a plain 3-array reference; and whether the RMQ/PSV/NSV structure alone could replace the child table (it supports parent and children in O(1) each; saves about 3 B/char). Not reopening the child-table decision: only recorded as an optimisation candidate | `dna_10m`, `prot_1m` | child table stays unless the RMQ-only variant is >= 1.5x faster on count and smaller |
| P4 | **dropped 2026-11-03 (1a.5, review MAJ-5): optimized kernels are `int`-only, only the safe reference is generic; see [ADR-002](ADR-002-sa-construction.md) section 2 item 6.** Was: generic-struct index width (`int` vs `long`) code quality, see [ADR-002](ADR-002-sa-construction.md) section 2 | SA-IS inner loops | `int` instantiation within 3 % of a hand-written `int` clone |

## 6. The contract (C# sketch; final names fixed in WP 1b.0)

New files, additive, in `SuffixTree.Core` namespace `SuffixTree.Indexing` (no existing file is edited). The engine lives in a new project `SuffixTree.Esa` (`IsAotCompatible`, `IsTrimmable`, `AllowUnsafeBlocks`; same settings as `SuffixTree.Persistent.csproj`).

```csharp
namespace SuffixTree.Indexing;

public enum IndexOperation
{
    Count, ExtendRight, ExtendLeft, LazyLocate, TextAccess,
    Parent, Child, SuffixLink, StringDepth,
    SuffixArrayAccess, InverseSuffixArrayAccess, LcpAccess, Lcp /*lcp(i,j)*/, Lca,
}

public enum IndexSupport : byte
{
    Unsupported = 0, // the engine cannot do it; calling throws NotSupportedException (never fakes it)
    Degraded = 1,    // correct, but slower than the engine class norm, or needs a lazily built auxiliary structure
    Native = 2,
}

/// <summary>A place in the virtual suffix tree (explicit node or inside an edge).</summary>
public readonly struct IndexLocus : IEquatable<IndexLocus>
{
    public long Lo { get; }      // engine-defined; SA-interval engines: first rank
    public long Hi { get; }      // engine-defined; SA-interval engines: last rank
    public long Length { get; }  // string depth = number of matched characters
    public long Aux { get; }     // engine-defined (ESA: depth of the lowest explicit node at or below the locus)
    // equality, ToString, constructor omitted
}

public interface IIndexText                       // 64-bit counterpart of ITextSource
{
    long Length { get; }
    int SymbolAt(long position);                  // char code; -1 for position >= Length (terminator)
    int CopyTo(long start, Span<char> destination);
}

public interface ISuffixIndex
{
    string Engine { get; }                        // "esa", later "fm", "r-index", "ukkonen"
    long TextLength { get; }
    IIndexText Text { get; }
    IndexLocus Root { get; }                      // empty string; Count(Root) == TextLength
    IndexSupport Support(IndexOperation operation);
    string? SupportNote(IndexOperation operation);   // e.g. "O(width), builds ISA on first use"

    long Count(in IndexLocus locus);              // number of start positions of the locus string
    bool TryExtendRight(in IndexLocus locus, int symbol, out IndexLocus result);
    bool TryExtendRight(in IndexLocus locus, ReadOnlySpan<char> pattern, out IndexLocus result);
    bool TryExtendLeft(in IndexLocus locus, int symbol, out IndexLocus result);
    // lazy / paged locate: writes up to positions.Length start positions of the locus string, skipping 'skip'
    // in engine order (ESA: SA order, not ascending); returns the number written.
    int Locate(in IndexLocus locus, long skip, Span<long> positions);
}

public interface ISuffixTreeNavigation : ISuffixIndex
{
    bool IsExplicitNode(in IndexLocus locus);
    bool TryParent(in IndexLocus locus, out IndexLocus parent);      // nearest explicit ancestor above the locus
    bool TryChild(in IndexLocus node, int firstSymbol, out IndexLocus child);   // -1 = terminator child
    int GetChildren(in IndexLocus node, Span<IndexLocus> children);  // returns the child count (may exceed children.Length: then nothing is written)
    bool TrySuffixLink(in IndexLocus node, out IndexLocus link);
    int EdgeSymbol(in IndexLocus node, long offset);                 // symbol on the node's incoming edge
}

public interface ISuffixArrayAccess : ISuffixIndex                  // engines that really have SA/LCP/ISA
{
    long SuffixAt(long rank);          // SA[rank]
    long RankOf(long position);        // ISA[position]
    long LcpAt(long rank);             // LCP[rank]
    long Lcp(long rankA, long rankB);  // range minimum
}
```

Rules: (1) symbols are UTF-16 code units `0..65535`, `-1` = terminator (same as `ISuffixTreeNavigator.cs:33-35,48`); integer alphabets are added in S2 as a second overload set, not by changing these signatures. (2) `Support(op) == Unsupported` implies the method throws `NotSupportedException` naming the engine and the operation; no method returns a wrong or empty answer for an unsupported operation. (3) `Degraded` always comes with a `SupportNote` that states the cost. (4) Handles are values; engines must be safe for concurrent readers (the existing trees are, `README.md` of `SuffixTree.Persistent`, "Thread-safe reads"), lazily built auxiliary structures are published once with `Interlocked`/`Lazy`. (5) `Locate` order is engine order; callers sort if they need ascending positions (as for `ISuffixTreeSearch.cs:49-54`).

Capability matrix, now and planned (`N` native, `D` degraded, `-` unsupported):

| operation | ESA (S1) | FM / r-index (S7) | online Ukkonen (future engine) |
|---|---|---|---|
| Count, ExtendRight, LazyLocate, TextAccess | N | N (locate via samples) | N |
| ExtendLeft | **D** (interval positions `p-1` through ISA, min/max rank, O(width), needs the pack) | N | - |
| Child, Parent, StringDepth | N | D | N |
| SuffixLink | **D** until the pack exists, then N | D | N |
| SuffixArrayAccess / InverseSuffixArrayAccess / LcpAccess | N / N (pack) / N | D / D / D | - |
| Lcp(i,j), Lca | D without the RMQ helper, N with it | D | - |

## 7. Relation to the existing interfaces

`EsaSuffixTree : ISuffixTree, ISuffixTreeNavigation, ISuffixArrayAccess, IDisposable` (concrete class in `SuffixTree.Esa`; disposable like `PersistentSuffixTree`, and the k-common code already disposes trees that are `IDisposable`, `SuffixTreeAlgorithms.KCommon.cs:195-197`). The `int` surface is a facade over the `long` core:

| Existing member | Realisation | Notes |
|---|---|---|
| `Text`, `NodeCount`, `LeafCount`, `MaxDepth`, `IsEmpty` | header values (`NODE_COUNT`, `n`, `n`) | formulas in section 4; `n > int.MaxValue - 1` throws `NotSupportedException("use the long API")` at the facade, never overflows silently (cf. K0 gaps on unchecked int arithmetic) |
| `Contains`, `CountOccurrences`, `FindAllOccurrences` (string and span) | descent, then `rb-lb+1` or an SA slice | empty-pattern behaviour as `PersistentSuffixTree.Search.cs:10,18,35-36,58`; result order = SA order, allowed by `ISuffixTreeSearch.cs:49-54` |
| `LongestRepeatedSubstring` | `LRS_DEPTH`/`LRS_RANK` stored in the header (as v6 stores the deepest node, `PersistentConstants.cs:87-93`) | tie representative is implementation-specific by contract (`ISuffixTreeAnalysis.cs:13`) |
| `GetAllSuffixes`, `EnumerateSuffixes` | SA order, ranks `1..n` | same order as the Ukkonen tree (terminator-first = shorter prefix first) |
| `LongestCommonSubstring*`, `FindAllLongestCommonSubstrings`, `FindAllDistinctLongestCommonSubstrings`, `FindExactMatchAnchors`, `FindMaximalExactMatches`, `FindMaximalUniqueMatches` | `SuffixTreeAlgorithms.*<EsaNode, EsaNavigator>` (as `PersistentSuffixTree.Analysis.cs:99-190`) | needs the navigation pack (section 5) |
| `FindAllLongestRepeatedSubstrings`, `CountDistinctSubstrings`, `CountDistinctSubstringsByLength`, `FindMaximalRepeatedPairs` | the default interface methods already work through `Traverse`; `CountDistinctSubstrings[ByLength]` get native O(n) overrides from the LCP array (`sum(n-SA[i]) - sum(LCP)`; histogram of LCP for the per-length counts) | the override is verified against the `Traverse` default **and** Ukkonen |
| `PrintTree`, `Traverse` | iterative lcp-interval walk (section 4) | |

Results the contract itself leaves implementation-specific, so differential tests check properties instead of equality: the representative of `LongestRepeatedSubstring` on ties (`ISuffixTreeAnalysis.cs:13`); `PositionInText` of `LongestCommonSubstringInfo` (`:46-48`) and of `FindExactMatchAnchors` (`:77-84`, `SuffixTreeAlgorithms.cs:93-96`) (check: it is a valid occurrence and the match is maximal); order of `FindAllOccurrences`. Everything else (MEM, MUM, all-LCS lists, repeated pairs, LRS lists, counts, `NodeCount`, the `Traverse` event stream) must be identical.

## 8. Positions: `long` contract, `int` bridge

The new contract is `long` end to end. The generic algorithms are `int`-only (section 2), so in S1 the ESA reuses them only for `n <= int.MaxValue - 1`, and the SA/ISA/child arrays use `int` storage ([ADR-002](ADR-002-sa-construction.md): `int[]` limit 2,147,483,591 elements). Beyond that the facade throws; the `long` contract still works with the 64-bit SA width of the container ([ADR-003](ADR-003-esa-container-format.md)). Making MEM/MUM/LCS generic in the position type is S2/S3 work (a copy over `long`, not a change of the existing `int` code).

## 9. Consequences and risks

* Positive: no generic algorithm is rewritten; LRS/repeat families come from `Traverse` plus O(n) native overrides; count/locate become a binary-search-free interval + SA slice.
* Risk R1 (high): link cost dominates MEM/anchors (P1, P2). Mitigation: O1 table fallback, S3 locus-based `ms[]`.
* Risk R2 (medium): the node-handle size (20 B) and the `ITextSource` virtual call at `:382` make the reused algorithms slower than hand-inlined code; measure in P2.
* Risk R3 (medium): conventions (terminator leaf, empty-suffix leaf, `leafCount`, `End = -1`, internal-edge occurrence, H-OCC) are subtle; the full `Traverse` event-stream comparison and the logical-hash comparison with Ukkonen on random and adversarial strings are the guard (WP 1b.1).
* Risk R4 (low): `Aux`/`Lo`/`Hi` semantics are engine-defined; the conformance suite (WP 1b.0) only uses the public methods.
* UNVERIFIED: the ISA/RMQ memory figures, the 600 ns threshold, the MEM targets (no ESA exists yet), `NEXTSYM` benefit, and that `Traverse` leaf `leafCount` equals 1 for every consumer (the consumers shown ignore it).

## 10. Open questions for the owner

1. Which profile counts for the "index <= 16 B/char" gate: Core only (about 10-12 B/char) or Full with the navigation pack (about 15)? Proposal: both are reported; the gate is on Full.
2. `PLAN.md:82` states the 3.3 s build gate for "SA + LCP + child table"; `PROGRESS-S1.md:4` states it for "SA+LCP" (S1a). Proposal: S1a measures SA+LCP <= 3.3 s; S1b measures the whole default pipeline (SA + LCP + child) <= 3.3 s again.
