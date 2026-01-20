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
    }
}
