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
    /// using O(n+m) suffix-link streaming traversal.
    /// </summary>
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
            DeduplicateInPlace(positionsInText);
            DeduplicateInPlace(positionsInOther);
        }

        return (substring, positionsInText, positionsInOther);
    }

    /// <summary>
    /// Finds exact-match anchors between the tree's text and <paramref name="query"/>
    /// using O(n+m) suffix-link streaming with peak tracking.
    /// </summary>
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

    private static void DeduplicateInPlace(List<int> values)
    {
        var seen = new HashSet<int>();
        int write = 0;
        for (int read = 0; read < values.Count; read++)
        {
            int value = values[read];
            if (!seen.Add(value))
                continue;

            values[write++] = value;
        }

        if (write < values.Count)
            values.RemoveRange(write, values.Count - write);
    }
}
