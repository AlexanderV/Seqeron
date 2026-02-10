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
        if (other.Length == 0 || nav.Text.Length == 0)
            return (string.Empty, new List<int>(), new List<int>());

        int maxLen = 0;
        var bestMatches = new List<(TNode Node, int MatchEndInOther)>();

        TNode currentNode = nav.Root;
        TNode currentEdge = nav.NullNode;
        int edgeOffset = 0;
        int currentMatchLen = 0;

        for (int i = 0; i < other.Length; i++)
        {
            int c = other[i];

            while (true)
            {
                if (!nav.IsNull(currentEdge))
                {
                    if (nav.GetEdgeSymbol(currentEdge, edgeOffset) == c)
                    {
                        edgeOffset++;
                        currentMatchLen++;
                        if (edgeOffset >= nav.LengthOf(currentEdge))
                        {
                            currentNode = currentEdge;
                            currentEdge = nav.NullNode;
                            edgeOffset = 0;
                        }
                        break;
                    }
                }
                else
                {
                    if (nav.TryGetChild(currentNode, c, out var nextChild) && !nav.IsNull(nextChild))
                    {
                        currentEdge = nextChild;
                        edgeOffset = 1;
                        currentMatchLen++;
                        if (edgeOffset >= nav.LengthOf(currentEdge))
                        {
                            currentNode = currentEdge;
                            currentEdge = nav.NullNode;
                            edgeOffset = 0;
                        }
                        break;
                    }
                }

                // Cannot extend — follow suffix link
                if (currentMatchLen == 0) break;

                if (!nav.IsRoot(currentNode))
                {
                    currentNode = nav.GetSuffixLink(currentNode);
                }
                currentMatchLen--;

                // Rescan from currentNode to restore edge position
                int nodeDepth = nav.GetNodeDepth(currentNode);
                int remaining = currentMatchLen - nodeDepth;

                if (remaining > 0)
                {
                    int pos = i - remaining;
                    currentEdge = nav.NullNode;
                    edgeOffset = 0;

                    while (remaining > 0)
                    {
                        if (!nav.TryGetChild(currentNode, other[pos], out var nc) || nav.IsNull(nc))
                            break;
                        int edgeLen = nav.LengthOf(nc);
                        if (edgeLen <= remaining)
                        {
                            pos += edgeLen;
                            remaining -= edgeLen;
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

            // Track best matches
            if (currentMatchLen > maxLen)
            {
                maxLen = currentMatchLen;
                bestMatches.Clear();
                bestMatches.Add((nav.IsNull(currentEdge) ? currentNode : currentEdge, i));
            }
            else if (currentMatchLen == maxLen && maxLen > 0 && !firstOnly)
            {
                bestMatches.Add((nav.IsNull(currentEdge) ? currentNode : currentEdge, i));
            }
        }

        if (maxLen == 0)
            return (string.Empty, new List<int>(), new List<int>());

        var positionsInText = new List<int>();
        var positionsInOther = new List<int>();

        foreach (var match in bestMatches)
        {
            positionsInOther.Add(match.MatchEndInOther - maxLen + 1);

            if (firstOnly)
            {
                int leafPos = nav.FindAnyLeafPosition(match.Node);
                if (leafPos >= 0)
                    positionsInText.Add(leafPos);
                break;
            }
            else
            {
                nav.CollectLeaves(match.Node, nav.GetDepthFromRoot(match.Node), positionsInText);
            }
        }

        string substring = other.Substring(bestMatches[0].MatchEndInOther - maxLen + 1, maxLen);
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
        if (query.Length == 0 || nav.Text.Length == 0 || minLength <= 0)
            return Array.Empty<(int, int, int)>();

        var results = new List<(int PositionInText, int PositionInQuery, int Length)>();

        TNode currentNode = nav.Root;
        TNode currentEdge = nav.NullNode;
        int edgeOffset = 0;
        int currentMatchLen = 0;

        // Peak tracking
        int peakLen = 0;
        int peakEndInQuery = -1;
        TNode peakNode = nav.NullNode;

        for (int i = 0; i < query.Length; i++)
        {
            int c = query[i];

            while (true)
            {
                if (!nav.IsNull(currentEdge))
                {
                    if (nav.GetEdgeSymbol(currentEdge, edgeOffset) == c)
                    {
                        edgeOffset++;
                        currentMatchLen++;
                        if (edgeOffset >= nav.LengthOf(currentEdge))
                        {
                            currentNode = currentEdge;
                            currentEdge = nav.NullNode;
                            edgeOffset = 0;
                        }
                        break;
                    }
                }
                else
                {
                    if (nav.TryGetChild(currentNode, c, out var nextChild) && !nav.IsNull(nextChild))
                    {
                        currentEdge = nextChild;
                        edgeOffset = 1;
                        currentMatchLen++;
                        if (edgeOffset >= nav.LengthOf(currentEdge))
                        {
                            currentNode = currentEdge;
                            currentEdge = nav.NullNode;
                            edgeOffset = 0;
                        }
                        break;
                    }
                }

                // Cannot extend — follow suffix link
                if (currentMatchLen == 0) break;

                if (!nav.IsRoot(currentNode))
                {
                    currentNode = nav.GetSuffixLink(currentNode);
                }
                currentMatchLen--;

                int nodeDepth = nav.GetNodeDepth(currentNode);
                int remaining = currentMatchLen - nodeDepth;

                if (remaining > 0)
                {
                    int pos = i - remaining;
                    currentEdge = nav.NullNode;
                    edgeOffset = 0;

                    while (remaining > 0)
                    {
                        if (!nav.TryGetChild(currentNode, query[pos], out var nc) || nav.IsNull(nc))
                            break;
                        int edgeLen = nav.LengthOf(nc);
                        if (edgeLen <= remaining)
                        {
                            pos += edgeLen;
                            remaining -= edgeLen;
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

            // Update peak tracking
            if (currentMatchLen >= minLength)
            {
                if (currentMatchLen > peakLen)
                {
                    peakLen = currentMatchLen;
                    peakEndInQuery = i;
                    peakNode = nav.IsNull(currentEdge) ? currentNode : currentEdge;
                }
            }
            else if (peakLen >= minLength)
            {
                EmitAnchor(ref nav, results, peakNode, peakEndInQuery, peakLen);
                peakLen = 0;
                peakEndInQuery = -1;
                peakNode = nav.NullNode;
            }
        }

        // Emit final run
        if (peakLen >= minLength && !nav.IsNull(peakNode))
        {
            EmitAnchor(ref nav, results, peakNode, peakEndInQuery, peakLen);
        }

        return results;
    }

    private static void EmitAnchor<TNode, TNav>(
        ref TNav nav,
        List<(int PositionInText, int PositionInQuery, int Length)> results,
        TNode node, int endInQuery, int length)
        where TNav : struct, ISuffixTreeNavigator<TNode>
    {
        int refPos = nav.FindAnyLeafPosition(node);
        if (refPos >= 0)
        {
            results.Add((refPos, endInQuery - length + 1, length));
        }
    }
}
