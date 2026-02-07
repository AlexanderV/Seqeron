using System;
using System.Collections.Generic;

namespace SuffixTree;

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
        var nav = new SuffixTreeNavigator(this);
        return SuffixTreeAlgorithms.FindAllLcs<SuffixTreeNode, SuffixTreeNavigator>(ref nav, other, firstOnly);
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
        var nav = new SuffixTreeNavigator(this);
        return SuffixTreeAlgorithms.FindExactMatchAnchors<SuffixTreeNode, SuffixTreeNavigator>(ref nav, query, minLength);
    }

    /// <summary>
    /// Walks to any leaf descendant and returns its position in the source text.
    /// Used by <see cref="SuffixTreeNavigator"/>.
    /// </summary>
    private int FindAnyLeafPosition(SuffixTreeNode node)
    {
        var current = node;
        var buffer = GetSearchBuffer();
        while (!current.IsLeaf)
        {
            current.GetChildren(buffer);
            if (buffer.Count == 0) return -1;
            current = buffer[buffer.Count - 1];
        }
        int leafDepth = GetNodeDepth(current);
        int pos = _text!.Length + 1 - leafDepth;
        return (pos >= 0 && pos < _text.Length) ? pos : -1;
    }
}
