using System;
using System.Collections.Generic;

namespace SuffixTree;
/// <summary>
/// Shared suffix tree algorithms that operate on any node representation
/// via <see cref="ISuffixTreeNavigator{TNode}"/>.
/// <para>
/// Both in-memory and persistent suffix trees delegate LCS, FindExactMatchAnchors and
/// the maximal exact/unique match enumeration to these methods, eliminating code duplication.
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
    /// The complete MEM / MUM sets are <see cref="FindMaximalExactMatches{TNode, TNav}"/> and
    /// <see cref="FindMaximalUniqueMatches{TNode, TNav}"/>.
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

    /// <summary>
    /// Finds all maximal exact matches (MEMs) of length ≥ <paramref name="minLength"/> between the
    /// tree's text (the reference) and <paramref name="query"/> — the set reported by MUMmer 3
    /// <c>mummer -maxmatch -l minLength</c> on the forward strand (Kurtz et al. 2004, Genome Biol 5:R12).
    /// </summary>
    /// <remarks>
    /// A MEM is a triple (r, q, len) with text[r..r+len) = query[q..q+len) that is left-maximal
    /// (r = 0, q = 0 or text[r-1] ≠ query[q-1]) and right-maximal (either string ends or
    /// text[r+len] ≠ query[q+len]). Every reference occurrence is reported.
    /// <para>
    /// Algorithm (MUMmer 3 <c>findmaxmat.c</c>): for each query start q the suffix-link walk
    /// (matching statistics, Chang &amp; Lawler 1994) keeps two loci — the locus of the longest
    /// prefix of query[q..] occurring in the text (length ms(q)) and the locus of its prefix of
    /// length minLength. The subtree below the minLength locus holds exactly the text suffixes that
    /// match ≥ minLength characters; a leaf below the child that leaves the matching path at string
    /// depth d matches exactly d characters (right-maximal), a leaf below the ms(q) locus matches
    /// ms(q). Each leaf is then tested for left-maximality.
    /// Time O(|text| + |query| + R), R = number of right-maximal matches of length ≥ minLength
    /// (the leaves visited); space O(output).
    /// </para>
    /// Output order: ascending by query position, then by text position.
    /// </remarks>
    /// <exception cref="ArgumentNullException"><paramref name="query"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="minLength"/> &lt; 1.</exception>
    public static IReadOnlyList<(int PositionInText, int PositionInQuery, int Length)>
        FindMaximalExactMatches<TNode, TNav>(ref TNav nav, string query, int minLength)
        where TNav : struct, ISuffixTreeNavigator<TNode>
    {
        ArgumentNullException.ThrowIfNull(query);
        ArgumentOutOfRangeException.ThrowIfLessThan(minLength, 1);
        ReadOnlySpan<char> q = query.AsSpan();
        ITextSource text = nav.Text;
        if (q.Length < minLength || text.Length < minLength)
            return Array.Empty<(int, int, int)>();

        var results = new List<(int PositionInText, int PositionInQuery, int Length)>();
        var children = new List<TNode>();
        var leaves = new List<int>();
        var hits = new List<(int PositionInText, int Length)>();

        // Locus of the longest match query[s..s+mLen) (matching statistics ms(s) = mLen).
        TNode mNode = nav.Root, mEdge = nav.NullNode;
        int mOffset = 0, mLen = 0, mDepth = 0;
        // Locus of query[s..s+pLen), pLen = min(ms(s), minLength) — MUMmer's "ploc".
        TNode pNode = nav.Root, pEdge = nav.NullNode;
        int pOffset = 0, pLen = 0, pDepth = 0;

        for (int s = 0; s < q.Length; s++)
        {
            ExtendMatch(ref nav, q, s, q.Length - s, ref mNode, ref mEdge, ref mOffset, ref mLen, ref mDepth);
            ExtendMatch(ref nav, q, s, Math.Min(mLen, minLength), ref pNode, ref pEdge, ref pOffset, ref pLen, ref pDepth);

            if (mLen >= minLength)
            {
                // Node y whose edge holds the minLength locus; its edge spans string depths [yStart, yEnd).
                bool inEdge = !nav.IsNull(pEdge);
                TNode y = inEdge ? pEdge : pNode;
                int yStart = inEdge ? pDepth : pDepth - nav.LengthOf(pNode);
                int yEnd = yStart + nav.LengthOf(y);
                hits.Clear();

                while (mLen > yEnd)
                {
                    int pathSymbol = q[s + yEnd];
                    nav.GetChildren(y, children);
                    foreach (TNode child in children)
                    {
                        if (nav.GetEdgeSymbol(child, 0) == pathSymbol) continue;
                        leaves.Clear();
                        nav.CollectLeaves(child, yEnd, leaves);
                        AddLeftMaximal(text, q, s, leaves, yEnd, hits);
                    }

                    nav.TryGetChild(y, pathSymbol, out TNode next);
                    y = next;
                    yStart = yEnd;
                    yEnd = yStart + nav.LengthOf(y);
                }

                leaves.Clear();
                nav.CollectLeaves(y, yStart, leaves);
                AddLeftMaximal(text, q, s, leaves, mLen, hits);

                hits.Sort(static (a, b) => a.PositionInText.CompareTo(b.PositionInText));
                foreach (var (r, len) in hits)
                    results.Add((r, s, len));
            }

            // Advance both loci from query[s..] to query[s+1..] via suffix links.
            if (mLen > 0)
                FollowSuffixLinkAndRescan(ref nav, q, s + mLen, ref mNode, ref mEdge, ref mOffset, ref mLen, ref mDepth);
            if (pLen > 0)
                FollowSuffixLinkAndRescan(ref nav, q, s + pLen, ref pNode, ref pEdge, ref pOffset, ref pLen, ref pDepth);
        }

        return results;
    }

    /// <summary>
    /// Finds maximal unique matches (MUMs) of length ≥ <paramref name="minLength"/> between the tree's
    /// text (the reference) and <paramref name="query"/>, as MUMmer 3 on the forward strand:
    /// <see cref="MumUniqueness.Both"/> = <c>mummer -mum</c>, <see cref="MumUniqueness.Reference"/> =
    /// <c>mummer -mumreference</c> (Kurtz et al. 2004, Genome Biol 5:R12).
    /// </summary>
    /// <remarks>
    /// <b>Reference</b> (MUMmer 3 <c>findmumcand.c</c>, <c>checkiflocationisMUMcand</c>): for each query
    /// start q, the longest prefix of query[q..] occurring in the text is reported when its length
    /// is ≥ minLength, its locus lies on a leaf edge (the string occurs exactly once in the text) and
    /// the match is left-maximal. Equivalently: every MEM whose string occurs exactly once in the text.
    /// <para>
    /// <b>Both</b> (MUMmer 3 <c>cleanMUMcand.c</c>, <c>mumuniqueinquery</c>): the Reference candidates
    /// are sorted by text position (longer first on ties) and a candidate is dropped when its text
    /// interval ends at or before the rightmost end seen so far (it lies inside another candidate's
    /// text interval, so its string recurs in the query); equal candidates are dropped together.
    /// The result is every MEM whose string occurs exactly once in the text and exactly once in the
    /// query. (MUMmer 3 starts the sweep with a right end of 0 instead of -1, which also drops a
    /// length-1 candidate at text position 0; that artefact is not reproduced.)
    /// </para>
    /// Time O(|text| + |query| + k log k) for k candidates. Output order: ascending by query
    /// position, then by text position.
    /// </remarks>
    /// <exception cref="ArgumentNullException"><paramref name="query"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException">
    /// <paramref name="minLength"/> &lt; 1, or <paramref name="uniqueness"/> is not a defined value.
    /// </exception>
    public static IReadOnlyList<(int PositionInText, int PositionInQuery, int Length)>
        FindMaximalUniqueMatches<TNode, TNav>(ref TNav nav, string query, int minLength, MumUniqueness uniqueness)
        where TNav : struct, ISuffixTreeNavigator<TNode>
    {
        ArgumentNullException.ThrowIfNull(query);
        ArgumentOutOfRangeException.ThrowIfLessThan(minLength, 1);
        if (uniqueness is not (MumUniqueness.Both or MumUniqueness.Reference))
            throw new ArgumentOutOfRangeException(nameof(uniqueness), uniqueness, "Undefined MumUniqueness value.");
        ReadOnlySpan<char> q = query.AsSpan();
        ITextSource text = nav.Text;
        if (q.Length < minLength || text.Length < minLength)
            return Array.Empty<(int, int, int)>();

        var candidates = new List<(int PositionInText, int PositionInQuery, int Length)>();
        TNode node = nav.Root, edge = nav.NullNode;
        int offset = 0, len = 0, depth = 0;

        for (int s = 0; s < q.Length; s++)
        {
            ExtendMatch(ref nav, q, s, q.Length - s, ref node, ref edge, ref offset, ref len, ref depth);

            // Inside a leaf edge (a leaf edge ends with the terminator): the string occurs once in the text.
            if (len >= minLength && !nav.IsNull(edge) && nav.GetEdgeSymbol(edge, nav.LengthOf(edge) - 1) < 0)
            {
                int r = nav.FindAnyLeafPosition(edge, depth);
                if (r >= 0 && (s == 0 || r == 0 || text[r - 1] != q[s - 1]))
                    candidates.Add((r, s, len));
            }

            if (len > 0)
                FollowSuffixLinkAndRescan(ref nav, q, s + len, ref node, ref edge, ref offset, ref len, ref depth);
        }

        if (uniqueness == MumUniqueness.Reference || candidates.Count == 0)
            return candidates; // already ascending by query position (one candidate per start)

        var mums = KeepUniqueInQuery(candidates);
        mums.Sort(static (a, b) => a.PositionInQuery != b.PositionInQuery
            ? a.PositionInQuery.CompareTo(b.PositionInQuery)
            : a.PositionInText.CompareTo(b.PositionInText));
        return mums;
    }

    /// <summary>MUMmer 3 <c>mumuniqueinquery</c> (cleanMUMcand.c) sweep over reference-unique candidates.</summary>
    private static List<(int PositionInText, int PositionInQuery, int Length)> KeepUniqueInQuery(
        List<(int PositionInText, int PositionInQuery, int Length)> c)
    {
        c.Sort(static (a, b) => a.PositionInText != b.PositionInText
            ? a.PositionInText.CompareTo(b.PositionInText)
            : b.Length.CompareTo(a.Length));

        var kept = new List<(int PositionInText, int PositionInQuery, int Length)>();
        long rightmost = -1;
        bool ignorePrevious = false;
        for (int i = 0; i < c.Count; i++)
        {
            bool ignoreCurrent = false;
            long currentRight = (long)c[i].PositionInText + c[i].Length - 1;
            if (rightmost > currentRight)
            {
                ignoreCurrent = true;
            }
            else if (rightmost == currentRight)
            {
                ignoreCurrent = true;
                if (!ignorePrevious && c[i - 1].PositionInText == c[i].PositionInText)
                    ignorePrevious = true;
            }
            else
            {
                rightmost = currentRight;
            }

            if (i > 0 && !ignorePrevious)
                kept.Add(c[i - 1]);
            ignorePrevious = ignoreCurrent;
        }

        if (!ignorePrevious)
            kept.Add(c[^1]);
        return kept;
    }

    /// <summary>Extends the match query[start..start+len) along the tree while len &lt; maxLength.</summary>
    private static void ExtendMatch<TNode, TNav>(
        ref TNav nav, ReadOnlySpan<char> source, int start, int maxLength,
        ref TNode node, ref TNode edge, ref int edgeOffset, ref int length, ref int nodeDepth)
        where TNav : struct, ISuffixTreeNavigator<TNode>
    {
        while (length < maxLength)
        {
            if (!TryConsumeSymbol(ref nav, source[start + length], ref node, ref edge, ref edgeOffset, ref length, ref nodeDepth))
                return;
        }
    }

    private static void AddLeftMaximal(ITextSource text, ReadOnlySpan<char> query, int queryStart,
        List<int> textStarts, int length, List<(int PositionInText, int Length)> hits)
    {
        foreach (int r in textStarts)
        {
            if (queryStart == 0 || r == 0 || text[r - 1] != query[queryStart - 1])
                hits.Add((r, length));
        }
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
