using System;
using System.Collections.Generic;

namespace SuffixTree;
/// <summary>
/// Shared suffix tree algorithms that operate on any node representation
/// via <see cref="ISuffixTreeNavigator{TNode}"/>.
/// <para>
/// Both in-memory and persistent suffix trees delegate LCS and
/// FindExactMatchAnchors to these methods, eliminating code duplication.
/// The <c>struct</c> constraint on the navigator ensures
/// the JIT specializes each call site — zero overhead vs hand-inlined code.
/// </para>
/// <para>
/// <b>v6 Slim compatibility:</b> these algorithms track node depth on-the-fly
/// using <c>currentNodeDepth</c> state instead of reading DepthFromRoot from
/// storage. After following a suffix link, depth decreases by exactly 1
/// (suffix link invariant). During rescan, depth accumulates edge lengths.
/// This eliminates the need for stored DepthFromRoot while preserving O(n+m).
/// </para>
/// </summary>
public static class SuffixTreeAlgorithms
{
    /// <summary>
    /// Finds the longest common substring between the tree's text and <paramref name="other"/>
    /// using O(n+m) suffix-link streaming traversal (matching statistics, Chang &amp; Lawler 1990;
    /// Gusfield 1997 §7.4/§7.8).
    /// </summary>
    /// <remarks>
    /// Tie-break: when several distinct substrings share the maximal length, the one whose
    /// occurrence in <paramref name="other"/> ends first (equivalently starts first) is returned,
    /// and only its positions are reported. With <paramref name="firstOnly"/> the text position is
    /// one occurrence reached by an arbitrary leaf walk (not necessarily the leftmost; it can differ
    /// between tree implementations). Without it, both position lists are ascending and duplicate-free.
    /// </remarks>
    public static (string Substring, List<int> PositionsInText, List<int> PositionsInOther)
        FindAllLcs<TNode, TNav>(ref TNav nav, string other, bool firstOnly)
        where TNav : struct, ISuffixTreeNavigator<TNode>
    {
        ArgumentNullException.ThrowIfNull(other);
        ReadOnlySpan<char> otherSpan = other.AsSpan();
        if (otherSpan.Length == 0 || nav.Text.Length == 0)
            return (string.Empty, new List<int>(), new List<int>());

        int maxLen = 0;
        // (Node, MatchEndInOther, DepthFromRoot of Node)
        var bestMatches = new List<(TNode Node, int MatchEndInOther, int DepthFromRoot)>();

        TNode currentNode = nav.Root;
        TNode currentEdge = nav.NullNode;
        int edgeOffset = 0;
        int currentMatchLen = 0;

        // Depth to END of currentNode's edge (= GetNodeDepth(currentNode)).
        // For root this is 0. Updated on suffix link follow (-1), edge
        // completion (+edgeLen), and rescan (+edgeLen per full edge).
        int currentNodeDepth = 0;

        for (int i = 0; i < otherSpan.Length; i++)
        {
            int c = otherSpan[i];

            while (true)
            {
                if (TryConsumeSymbol(ref nav, c, ref currentNode, ref currentEdge, ref edgeOffset, ref currentMatchLen, ref currentNodeDepth))
                    break;

                // Cannot extend — follow suffix link
                if (currentMatchLen == 0) break;
                FollowSuffixLinkAndRescan(ref nav, otherSpan, i, ref currentNode, ref currentEdge, ref edgeOffset, ref currentMatchLen, ref currentNodeDepth);
            }

            // Track best matches
            if (currentMatchLen > maxLen)
            {
                maxLen = currentMatchLen;
                bestMatches.Clear();
                TNode matchNode = nav.IsNull(currentEdge) ? currentNode : currentEdge;
                // DepthFromRoot of matchNode:
                // - currentNode: currentNodeDepth - LengthOf(currentNode)
                // - currentEdge: currentNodeDepth (= depth to END of parent = depth to START of child)
                int matchDFR = nav.IsNull(currentEdge)
                    ? currentNodeDepth - nav.LengthOf(currentNode)
                    : currentNodeDepth;
                bestMatches.Add((matchNode, i, matchDFR));
            }
            else if (currentMatchLen == maxLen && maxLen > 0 && !firstOnly)
            {
                TNode matchNode = nav.IsNull(currentEdge) ? currentNode : currentEdge;
                int matchDFR = nav.IsNull(currentEdge)
                    ? currentNodeDepth - nav.LengthOf(currentNode)
                    : currentNodeDepth;
                bestMatches.Add((matchNode, i, matchDFR));
            }
        }

        if (maxLen == 0)
            return (string.Empty, new List<int>(), new List<int>());

        var positionsInText = new List<int>();
        var positionsInOther = new List<int>();
        string substring = other.Substring(bestMatches[0].MatchEndInOther - maxLen + 1, maxLen);

        foreach (var match in bestMatches)
        {
            int positionInOther = match.MatchEndInOther - maxLen + 1;
            // When multiple distinct substrings share the same maximal length,
            // report positions only for the canonical returned substring.
            if (!otherSpan.Slice(positionInOther, maxLen).SequenceEqual(substring.AsSpan()))
                continue;

            positionsInOther.Add(positionInOther);

            if (firstOnly)
            {
                int leafPos = nav.FindAnyLeafPosition(match.Node, match.DepthFromRoot);
                if (leafPos >= 0)
                    positionsInText.Add(leafPos);
                break;
            }
            else
            {
                nav.CollectLeaves(match.Node, match.DepthFromRoot, positionsInText);
            }
        }

        if (!firstOnly)
        {
            // Leaves arrive in child-storage order, which differs between tree implementations;
            // report ascending, duplicate-free positions so every ISuffixTree returns the same lists.
            SortAndDeduplicateInPlace(positionsInText);
            SortAndDeduplicateInPlace(positionsInOther);
        }

        return (substring, positionsInText, positionsInOther);
    }

    /// <summary>
    /// Finds exact-match anchors between the tree's text and <paramref name="query"/>
    /// using O(n+m) suffix-link streaming with peak tracking.
    /// </summary>
    /// <remarks>
    /// Let ms(i) be the matching statistic at query end position i (length of the longest suffix of
    /// query[0..i] occurring in the text). For every maximal run of consecutive i with
    /// ms(i) ≥ <paramref name="minLength"/>, exactly one anchor is emitted: the first peak of ms in
    /// the run. Each anchor is a maximal exact match (MEM: not extendable left or right for any text
    /// occurrence), but the result is a subset of all MEMs ≥ minLength (one per run, one text
    /// occurrence per anchor — arbitrary, implementation-specific), and anchors from adjacent runs
    /// may overlap in the query by fewer than minLength characters (text "aba", query "ababa",
    /// minLength 3 → (0,0,3), (0,2,3)). minLength ≤ 0 returns an empty list.
    /// </remarks>
    public static IReadOnlyList<(int PositionInText, int PositionInQuery, int Length)>
        FindExactMatchAnchors<TNode, TNav>(ref TNav nav, string query, int minLength)
        where TNav : struct, ISuffixTreeNavigator<TNode>
    {
        ArgumentNullException.ThrowIfNull(query);
        ReadOnlySpan<char> querySpan = query.AsSpan();
        if (querySpan.Length == 0 || nav.Text.Length == 0 || minLength <= 0)
            return Array.Empty<(int, int, int)>();

        var results = new List<(int PositionInText, int PositionInQuery, int Length)>();

        TNode currentNode = nav.Root;
        TNode currentEdge = nav.NullNode;
        int edgeOffset = 0;
        int currentMatchLen = 0;
        int currentNodeDepth = 0;

        // Peak tracking
        int peakLen = 0;
        int peakEndInQuery = -1;
        TNode peakNode = nav.NullNode;
        int peakDepthFromRoot = 0;

        for (int i = 0; i < querySpan.Length; i++)
        {
            int c = querySpan[i];

            while (true)
            {
                if (TryConsumeSymbol(ref nav, c, ref currentNode, ref currentEdge, ref edgeOffset, ref currentMatchLen, ref currentNodeDepth))
                    break;

                // Cannot extend — follow suffix link
                if (currentMatchLen == 0) break;
                FollowSuffixLinkAndRescan(ref nav, querySpan, i, ref currentNode, ref currentEdge, ref edgeOffset, ref currentMatchLen, ref currentNodeDepth);
            }

            // Update peak tracking
            if (currentMatchLen >= minLength)
            {
                if (currentMatchLen > peakLen)
                {
                    peakLen = currentMatchLen;
                    peakEndInQuery = i;
                    peakNode = nav.IsNull(currentEdge) ? currentNode : currentEdge;
                    peakDepthFromRoot = nav.IsNull(currentEdge)
                        ? currentNodeDepth - nav.LengthOf(currentNode)
                        : currentNodeDepth;
                }
            }
            else if (peakLen >= minLength)
            {
                EmitAnchor(ref nav, results, peakNode, peakEndInQuery, peakLen, peakDepthFromRoot);
                peakLen = 0;
                peakEndInQuery = -1;
                peakNode = nav.NullNode;
            }
        }

        // Emit final run
        if (peakLen >= minLength && !nav.IsNull(peakNode))
        {
            EmitAnchor(ref nav, results, peakNode, peakEndInQuery, peakLen, peakDepthFromRoot);
        }

        return results;
    }

    private static bool TryConsumeSymbol<TNode, TNav>(
        ref TNav nav,
        int symbol,
        ref TNode currentNode,
        ref TNode currentEdge,
        ref int edgeOffset,
        ref int currentMatchLen,
        ref int currentNodeDepth)
        where TNav : struct, ISuffixTreeNavigator<TNode>
    {
        if (!nav.IsNull(currentEdge))
        {
            if (nav.GetEdgeSymbol(currentEdge, edgeOffset) != symbol)
                return false;

            edgeOffset++;
            currentMatchLen++;
            if (edgeOffset >= nav.LengthOf(currentEdge))
            {
                currentNodeDepth += nav.LengthOf(currentEdge);
                currentNode = currentEdge;
                currentEdge = nav.NullNode;
                edgeOffset = 0;
            }
            return true;
        }

        if (!nav.TryGetChild(currentNode, symbol, out var nextChild) || nav.IsNull(nextChild))
            return false;

        currentEdge = nextChild;
        edgeOffset = 1;
        currentMatchLen++;
        if (edgeOffset >= nav.LengthOf(currentEdge))
        {
            currentNodeDepth += nav.LengthOf(currentEdge);
            currentNode = currentEdge;
            currentEdge = nav.NullNode;
            edgeOffset = 0;
        }

        return true;
    }

    private static void FollowSuffixLinkAndRescan<TNode, TNav>(
        ref TNav nav,
        ReadOnlySpan<char> source,
        int sourceIndex,
        ref TNode currentNode,
        ref TNode currentEdge,
        ref int edgeOffset,
        ref int currentMatchLen,
        ref int currentNodeDepth)
        where TNav : struct, ISuffixTreeNavigator<TNode>
    {
        if (!nav.IsRoot(currentNode))
        {
            currentNode = nav.GetSuffixLink(currentNode);
            currentNodeDepth--;
        }
        currentMatchLen--;

        int nodeDepth = currentNodeDepth;
        int remaining = currentMatchLen - nodeDepth;

        if (remaining > 0)
        {
            int pos = sourceIndex - remaining;
            currentEdge = nav.NullNode;
            edgeOffset = 0;

            while (remaining > 0)
            {
                if (!nav.TryGetChild(currentNode, source[pos], out var nc) || nav.IsNull(nc))
                    break;
                int edgeLen = nav.LengthOf(nc);
                if (edgeLen <= remaining)
                {
                    pos += edgeLen;
                    remaining -= edgeLen;
                    currentNodeDepth += edgeLen;
                    currentNode = nc;
                }
                else
                {
                    currentEdge = nc;
                    edgeOffset = remaining;
                    remaining = 0;
                }
            }
        }
        else
        {
            currentEdge = nav.NullNode;
            edgeOffset = 0;
        }
    }

    private static void EmitAnchor<TNode, TNav>(
        ref TNav nav,
        List<(int PositionInText, int PositionInQuery, int Length)> results,
        TNode node, int endInQuery, int length, int depthFromRoot)
        where TNav : struct, ISuffixTreeNavigator<TNode>
    {
        int refPos = nav.FindAnyLeafPosition(node, depthFromRoot);
        if (refPos >= 0)
        {
            results.Add((refPos, endInQuery - length + 1, length));
        }
    }

    private static void SortAndDeduplicateInPlace(List<int> values)
    {
        if (values.Count < 2) return;
        values.Sort();
        int write = 1;
        for (int read = 1; read < values.Count; read++)
        {
            if (values[read] != values[write - 1])
                values[write++] = values[read];
        }

        if (write < values.Count)
            values.RemoveRange(write, values.Count - write);
    }

    /// <summary>
    /// Counts distinct substrings of the tree's text by length: element <c>i</c> of the result is
    /// the number of distinct substrings of length <c>i</c> for <c>i = 1..min(maxLength, n)</c>;
    /// element 0 is 1 (the empty substring). The terminator is not part of any substring.
    /// </summary>
    /// <remarks>
    /// Every distinct substring of length <c>i</c> ends at exactly one point at string depth
    /// <c>i</c> on exactly one edge, so the count equals the number of edges whose depth range
    /// covers <c>i</c> (Gusfield 1997 §7; Troyanskaya et al. 2002, Bioinformatics 18:679).
    /// Each edge adds +1 over its depth range in a difference array: O(nodes + maxLength).
    /// </remarks>
    /// <param name="tree">Any suffix tree implementation (traversed via <see cref="ISuffixTreeDiagnostics.Traverse"/>).</param>
    /// <param name="maxLength">Largest substring length to count (≥ 1).</param>
    public static long[] CountDistinctSubstringsByLength(ISuffixTree tree, int maxLength)
    {
        ArgumentNullException.ThrowIfNull(tree);
        ArgumentOutOfRangeException.ThrowIfLessThan(maxLength, 1);

        int n = tree.Text.Length;
        int m = Math.Min(maxLength, n);
        var visitor = new DistinctSubstringCounter(n, m);
        tree.Traverse(visitor);
        return visitor.ToCounts();
    }

    /// <summary>
    /// Total number of distinct non-empty substrings of the tree's text
    /// (sum of edge label lengths, terminator excluded). O(nodes).
    /// </summary>
    public static long CountDistinctSubstrings(ISuffixTree tree)
    {
        ArgumentNullException.ThrowIfNull(tree);
        var visitor = new DistinctSubstringCounter(tree.Text.Length, 0);
        tree.Traverse(visitor);
        return visitor.TotalEdgeLength;
    }

    /// <summary>
    /// Shared per-edge accumulator for distinct-substring counting. As a visitor it works with any
    /// tree via <see cref="ISuffixTreeDiagnostics.Traverse"/>; implementations with a cheaper node
    /// walk call <see cref="AddEdge"/> directly.
    /// </summary>
    public sealed class DistinctSubstringCounter : ISuffixTreeVisitor
    {
        private readonly int _n;
        private readonly int _m;
        private readonly long[] _diff;
        private bool _rootSeen;

        /// <param name="n">Text length (terminator excluded).</param>
        /// <param name="m">Largest substring length to count per length (0 = totals only).</param>
        public DistinctSubstringCounter(int n, int m)
        {
            ArgumentOutOfRangeException.ThrowIfNegative(n);
            ArgumentOutOfRangeException.ThrowIfNegative(m);
            _n = n;
            _m = m;
            _diff = new long[m + 2];
        }

        /// <summary>Total number of distinct non-empty substrings seen so far.</summary>
        public long TotalEdgeLength { get; private set; }

        /// <summary>
        /// Adds one edge covering text[start, end) whose label begins at string depth
        /// <paramref name="depth"/>. A leaf edge (end &lt; 0 or past the text) stops at the text end.
        /// </summary>
        public void AddEdge(int start, int end, int depth)
        {
            int realEnd = end < 0 || end > _n ? _n : end;
            int length = realEnd - start;
            if (length <= 0) return;

            TotalEdgeLength += length;

            int from = depth + 1;
            if (from > _m) return;
            int to = Math.Min(depth + length, _m);
            _diff[from]++;
            _diff[to + 1]--;
        }

        /// <summary>Per-length counts; element 0 is 1 (the empty substring).</summary>
        public long[] ToCounts()
        {
            var counts = new long[_m + 1];
            counts[0] = 1;
            long running = 0;
            for (int i = 1; i <= _m; i++)
            {
                running += _diff[i];
                counts[i] = running;
            }
            return counts;
        }

        /// <inheritdoc />
        public void VisitNode(int startIndex, int endIndex, int leafCount, int childCount, int depth)
        {
            if (!_rootSeen)
            {
                _rootSeen = true; // the root has no incoming edge
                return;
            }
            AddEdge(startIndex, endIndex, depth);
        }

        /// <inheritdoc />
        public void EnterBranch(int key) { }

        /// <inheritdoc />
        public void ExitBranch() { }
    }
}
