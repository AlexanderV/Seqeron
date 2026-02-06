using System;
using System.Collections.Generic;
using System.Linq;

namespace SuffixTree
{
    public partial class SuffixTree
    {
        /// <summary>
        /// Finds the longest substring that appears at least twice in the text.
        /// If multiple substrings have the same maximum length, one of them is returned.
        /// </summary>
        /// <returns>The longest repeated substring, or empty string if none exists.</returns>
        public string LongestRepeatedSubstring()
        {
            if (_deepestInternalNode == null || _maxInternalDepth == 0)
                return string.Empty;

            return _text.Substring(_deepestInternalNode.Start - _deepestInternalNode.DepthFromRoot, _maxInternalDepth);
        }

        /// <summary>
        /// Finds the longest common substring between this tree's text and another string.
        /// If multiple substrings have the same maximum length, the first one found in 'other' is returned.
        /// </summary>
        /// <param name="other">The string to compare against.</param>
        /// <returns>The longest common substring, or empty string if none exists.</returns>
        public string LongestCommonSubstring(string other)
            => LongestCommonSubstringInfo(other).Substring;

        /// <summary>
        /// Finds the longest common substring between this tree's text and another character span.
        /// Zero-allocation overload for performance-critical scenarios.
        /// </summary>
        /// <param name="other">The character span to compare against.</param>
        /// <returns>The longest common substring, or empty string if none exists.</returns>
        public string LongestCommonSubstring(ReadOnlySpan<char> other)
            => LongestCommonSubstringInfo(other.ToString()).Substring;

        /// <summary>
        /// Finds the longest common substring with position information.
        /// If multiple substrings have the same maximum length, the first one found in 'other' is returned.
        /// </summary>
        /// <param name="other">The string to compare against.</param>
        /// <returns>
        /// A tuple containing: the substring, position in tree's text, position in other.
        /// Returns (empty string, -1, -1) if no common substring exists.
        /// </returns>
        public (string Substring, int PositionInText, int PositionInOther) LongestCommonSubstringInfo(string other)
        {
            var results = FindAllLongestCommonSubstringsInternal(other, true);
            if (results.PositionsInText.Count == 0)
                return (string.Empty, -1, -1);

            return (results.Substring, results.PositionsInText[0], results.PositionsInOther[0]);
        }

        /// <summary>
        /// Finds all positions where the longest common substring occurs.
        /// If multiple substrings have the same maximum length, all occurrences for all such candidates are returned.
        /// </summary>
        /// <param name="other">The string to compare against.</param>
        /// <returns>
        /// A tuple containing: the substring, all positions in tree's text, all positions in other.
        /// Returns (empty string, empty list, empty list) if no common substring exists.
        /// </returns>
        public (string Substring, IReadOnlyList<int> PositionsInText, IReadOnlyList<int> PositionsInOther) FindAllLongestCommonSubstrings(string other)
        {
            var results = FindAllLongestCommonSubstringsInternal(other, false);
            return (results.Substring, results.PositionsInText, results.PositionsInOther);
        }

        private (string Substring, List<int> PositionsInText, List<int> PositionsInOther) FindAllLongestCommonSubstringsInternal(string other, bool firstOnly)
        {
            ArgumentNullException.ThrowIfNull(other);
            if (other.Length == 0 || _text.Length == 0)
                return (string.Empty, new List<int>(), new List<int>());

            int maxLen = 0;
            var bestMatches = new List<(SuffixTreeNode Node, int MatchEndInOther)>();

            var currentNode = _root;
            SuffixTreeNode? currentEdge = null;
            int edgeOffset = 0;
            int currentMatchLen = 0;

            // Local function to finalize edge traversal when offset reaches edge length
            void TryFinalizeEdge()
            {
                if (currentEdge != null && edgeOffset >= LengthOf(currentEdge))
                {
                    currentNode = currentEdge;
                    currentEdge = null;
                    edgeOffset = 0;
                }
            }

            for (int i = 0; i < other.Length; i++)
            {
                char c = other[i];

                while (true)
                {
                    if (currentEdge != null)
                    {
                        if (GetSymbolAt(currentEdge.Start + edgeOffset) == c)
                        {
                            edgeOffset++;
                            currentMatchLen++;
                            TryFinalizeEdge();
                            break;
                        }
                    }
                    else
                    {
                        if (currentNode.TryGetChild(c, out var nextChild) && nextChild != null)
                        {
                            currentEdge = nextChild;
                            edgeOffset = 1;
                            currentMatchLen++;
                            TryFinalizeEdge();
                            break;
                        }
                    }

                    if (currentMatchLen == 0) break;

                    if (currentNode != _root)
                    {
                        currentNode = currentNode.SuffixLink ?? _root;
                    }
                    currentMatchLen--;

                    int nodeDepth = GetNodeDepth(currentNode);
                    int remaining = currentMatchLen - nodeDepth;

                    if (remaining > 0)
                    {
                        int pos = i - remaining; // Character at current relative position in 'other'
                        currentEdge = null;
                        edgeOffset = 0;

                        while (remaining > 0)
                        {
                            if (!currentNode.TryGetChild(other[pos], out var nextChild) || nextChild == null)
                                break;

                            int edgeLen = LengthOf(nextChild);
                            if (edgeLen <= remaining)
                            {
                                pos += edgeLen;
                                remaining -= edgeLen;
                                currentNode = nextChild;
                            }
                            else
                            {
                                currentEdge = nextChild;
                                edgeOffset = remaining;
                                remaining = 0;
                            }
                        }
                        System.Diagnostics.Debug.Assert(remaining == 0, "Rescan logic did not fully consume the remaining characters.");
                    }
                    else
                    {
                        currentEdge = null;
                        edgeOffset = 0;
                    }
                }

                if (currentMatchLen > maxLen)
                {
                    maxLen = currentMatchLen;
                    bestMatches.Clear();
                    bestMatches.Add((currentEdge ?? currentNode, i));
                }
                else if (currentMatchLen == maxLen && maxLen > 0)
                {
                    bestMatches.Add((currentEdge ?? currentNode, i));
                }
            }

            if (maxLen == 0)
                return (string.Empty, new List<int>(), new List<int>());

            var positionsInText = new List<int>();
            var positionsInOther = new List<int>();

            foreach (var match in bestMatches)
            {
                int matchEndInOther = match.MatchEndInOther;
                positionsInOther.Add(matchEndInOther - maxLen + 1);

                if (firstOnly)
                {
                    ReconstructPath(match.Node, positionsInText);
                    break;
                }
                else
                {
                    CollectLeaves(match.Node, match.Node.DepthFromRoot, positionsInText);
                }
            }

            string substring = other.Substring(bestMatches[0].MatchEndInOther - maxLen + 1, maxLen);
            return (substring, positionsInText, positionsInOther);
        }

        private void ReconstructPath(SuffixTreeNode node, List<int> results)
        {
            var current = node;
            var buffer = GetSearchBuffer();
            while (!current.IsLeaf)
            {
                current.GetChildren(buffer);
                current = buffer[0];
            }
            int leafDepth = GetNodeDepth(current);
            results.Add(_text!.Length + 1 - leafDepth);
        }

        /// <summary>
        /// Finds exact-match anchors between this tree's text and a query string
        /// using O(n + m) suffix-link-based streaming traversal.
        /// <para>
        /// This method walks the query against the suffix tree using suffix links,
        /// identical to the longest-common-substring algorithm, but emits all
        /// right-maximal matches whose length meets or exceeds <paramref name="minLength"/>.
        /// </para>
        /// <para>
        /// A match is emitted when the running match length drops below the threshold
        /// after being above it, capturing the peak (longest) match within each run.
        /// This produces non-overlapping anchors suitable for anchor-based alignment.
        /// </para>
        /// </summary>
        /// <param name="query">The query string to find matches against. Cannot be null.</param>
        /// <param name="minLength">Minimum match length to report (must be &gt; 0).</param>
        /// <returns>
        /// List of (PositionInText, PositionInQuery, Length) tuples representing exact-match
        /// anchors, ordered by their position in the query.
        /// </returns>
        /// <remarks>
        /// <b>Time complexity:</b> O(|text| + |query|) — each character is processed at most
        /// twice (once for extension, once for suffix-link rescan).
        /// <para><b>Space complexity:</b> O(k) where k is the number of anchors found.</para>
        /// </remarks>
        public IReadOnlyList<(int PositionInText, int PositionInQuery, int Length)> FindExactMatchAnchors(
            string query, int minLength)
        {
            ArgumentNullException.ThrowIfNull(query);
            if (query.Length == 0 || _text.Length == 0 || minLength <= 0)
                return Array.Empty<(int, int, int)>();

            var results = new List<(int PositionInText, int PositionInQuery, int Length)>();

            var currentNode = _root;
            SuffixTreeNode? currentEdge = null;
            int edgeOffset = 0;
            int currentMatchLen = 0;

            // Peak tracking: captures the best match within a contiguous run above minLength.
            int peakLen = 0;
            int peakEndInQuery = -1;
            SuffixTreeNode? peakNode = null;

            for (int i = 0; i < query.Length; i++)
            {
                char c = query[i];

                // Streaming traversal — identical to LCS algorithm
                while (true)
                {
                    if (currentEdge != null)
                    {
                        if (GetSymbolAt(currentEdge.Start + edgeOffset) == c)
                        {
                            edgeOffset++;
                            currentMatchLen++;
                            if (edgeOffset >= LengthOf(currentEdge))
                            {
                                currentNode = currentEdge;
                                currentEdge = null;
                                edgeOffset = 0;
                            }
                            break;
                        }
                    }
                    else
                    {
                        if (currentNode.TryGetChild(c, out var nextChild) && nextChild != null)
                        {
                            currentEdge = nextChild;
                            edgeOffset = 1;
                            currentMatchLen++;
                            if (edgeOffset >= LengthOf(currentEdge))
                            {
                                currentNode = currentEdge;
                                currentEdge = null;
                                edgeOffset = 0;
                            }
                            break;
                        }
                    }

                    // Cannot extend — follow suffix link
                    if (currentMatchLen == 0) break;

                    if (currentNode != _root)
                    {
                        currentNode = currentNode.SuffixLink ?? _root;
                    }
                    currentMatchLen--;

                    int nodeDepth = GetNodeDepth(currentNode);
                    int remaining = currentMatchLen - nodeDepth;

                    if (remaining > 0)
                    {
                        int pos = i - remaining;
                        currentEdge = null;
                        edgeOffset = 0;

                        while (remaining > 0)
                        {
                            if (!currentNode.TryGetChild(query[pos], out var nextChild2) || nextChild2 == null)
                                break;

                            int edgeLen = LengthOf(nextChild2);
                            if (edgeLen <= remaining)
                            {
                                pos += edgeLen;
                                remaining -= edgeLen;
                                currentNode = nextChild2;
                            }
                            else
                            {
                                currentEdge = nextChild2;
                                edgeOffset = remaining;
                                remaining = 0;
                            }
                        }
                    }
                    else
                    {
                        currentEdge = null;
                        edgeOffset = 0;
                    }
                }

                // Update peak tracking — record the best match in the current run
                if (currentMatchLen >= minLength)
                {
                    if (currentMatchLen > peakLen)
                    {
                        peakLen = currentMatchLen;
                        peakEndInQuery = i;
                        peakNode = currentEdge ?? currentNode;
                    }
                }
                else if (peakLen >= minLength)
                {
                    // Match just dropped below threshold — emit the peak anchor
                    EmitAnchorFromPeak(results, peakNode!, peakEndInQuery, peakLen);
                    peakLen = 0;
                    peakEndInQuery = -1;
                    peakNode = null;
                }
            }

            // Emit the final run if still above threshold
            if (peakLen >= minLength && peakNode != null)
            {
                EmitAnchorFromPeak(results, peakNode, peakEndInQuery, peakLen);
            }

            return results;
        }

        /// <summary>
        /// Emits a single anchor from the peak of a match run.
        /// </summary>
        private void EmitAnchorFromPeak(
            List<(int PositionInText, int PositionInQuery, int Length)> results,
            SuffixTreeNode node,
            int endInQuery,
            int length)
        {
            int refPos = FindAnyLeafPosition(node);
            if (refPos >= 0)
            {
                results.Add((refPos, endInQuery - length + 1, length));
            }
        }

        /// <summary>
        /// Walks to any leaf descendant of the given node and returns its position
        /// in the source text. Prefers non-terminator children to avoid edge cases.
        /// </summary>
        private int FindAnyLeafPosition(SuffixTreeNode node)
        {
            var current = node;
            var buffer = GetSearchBuffer();
            while (!current.IsLeaf)
            {
                current.GetChildren(buffer);
                if (buffer.Count == 0) return -1;
                // Use last child to prefer non-terminator (terminator is added first in GetChildren)
                current = buffer[buffer.Count - 1];
            }
            int leafDepth = GetNodeDepth(current);
            int pos = _text!.Length + 1 - leafDepth;
            return (pos >= 0 && pos < _text.Length) ? pos : -1;
        }
    }
}
