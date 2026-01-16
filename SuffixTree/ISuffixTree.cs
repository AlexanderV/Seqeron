using System;
using System.Collections.Generic;

namespace SuffixTree
{
    /// <summary>
    /// Interface for suffix tree operations.
    /// Provides substring search and pattern matching capabilities.
    /// </summary>
    public interface ISuffixTree
    {
        /// <summary>
        /// Gets the original text that this suffix tree was built from.
        /// </summary>
        string Text { get; }

        /// <summary>
        /// Gets the total number of nodes in the tree (including root).
        /// </summary>
        int NodeCount { get; }

        /// <summary>
        /// Gets the number of leaf nodes in the tree.
        /// Equal to the length of the original text (one leaf per suffix).
        /// </summary>
        int LeafCount { get; }

        /// <summary>
        /// Gets the maximum depth of the tree (longest path from root to leaf in characters).
        /// </summary>
        int MaxDepth { get; }

        /// <summary>
        /// Checks if the specified string is a substring of the tree content.
        /// </summary>
        /// <param name="value">The substring to search for.</param>
        /// <returns>True if the substring exists, false otherwise.</returns>
        bool Contains(string value);

        /// <summary>
        /// Checks if the specified character span is a substring of the tree content.
        /// </summary>
        /// <param name="value">The character span to search for.</param>
        /// <returns>True if the substring exists, false otherwise.</returns>
        bool Contains(ReadOnlySpan<char> value);

        /// <summary>
        /// Finds all starting positions where the pattern occurs in the original string.
        /// </summary>
        /// <param name="pattern">The pattern to search for.</param>
        /// <returns>Collection of 0-based starting positions of all occurrences.</returns>
        IReadOnlyList<int> FindAllOccurrences(string pattern);

        /// <summary>
        /// Finds all starting positions where the pattern occurs in the original string.
        /// Zero-allocation overload for performance-critical scenarios.
        /// </summary>
        /// <param name="pattern">The pattern to search for.</param>
        /// <returns>Collection of 0-based starting positions of all occurrences.</returns>
        IReadOnlyList<int> FindAllOccurrences(ReadOnlySpan<char> pattern);

        /// <summary>
        /// Counts the number of occurrences of a pattern in the text.
        /// </summary>
        /// <param name="pattern">The pattern to count.</param>
        /// <returns>Number of times the pattern occurs in the text.</returns>
        int CountOccurrences(string pattern);

        /// <summary>
        /// Counts the number of occurrences of a pattern in the text.
        /// Zero-allocation overload for performance-critical scenarios.
        /// </summary>
        /// <param name="pattern">The pattern to count.</param>
        /// <returns>Number of times the pattern occurs in the text.</returns>
        int CountOccurrences(ReadOnlySpan<char> pattern);

        /// <summary>
        /// Finds the longest substring that appears at least twice in the text.
        /// </summary>
        /// <returns>The longest repeated substring, or empty string if none exists.</returns>
        string LongestRepeatedSubstring();

        /// <summary>
        /// Returns all suffixes of the original string in sorted order.
        /// Useful for debugging and educational purposes.
        /// </summary>
        /// <returns>All suffixes sorted lexicographically.</returns>
        IReadOnlyList<string> GetAllSuffixes();

        /// <summary>
        /// Enumerates all suffixes of the original string in sorted order lazily.
        /// Use this for large strings to avoid O(n²) memory allocation.
        /// </summary>
        /// <returns>Lazy enumerable of suffixes sorted lexicographically.</returns>
        IEnumerable<string> EnumerateSuffixes();

        /// <summary>
        /// Finds the longest common substring between this tree's text and another string.
        /// </summary>
        /// <param name="other">The string to compare against.</param>
        /// <returns>The longest common substring, or empty string if none exists.</returns>
        string LongestCommonSubstring(string other);

        /// <summary>
        /// Finds the longest common substring between this tree's text and another character span.
        /// Zero-allocation overload for performance-critical scenarios.
        /// </summary>
        /// <param name="other">The character span to compare against.</param>
        /// <returns>The longest common substring, or empty string if none exists.</returns>
        string LongestCommonSubstring(ReadOnlySpan<char> other);

        /// <summary>
        /// Finds the longest common substring with position information.
        /// </summary>
        /// <param name="other">The string to compare against.</param>
        /// <returns>
        /// A tuple containing: the substring, position in tree's text, position in other.
        /// Returns (empty string, -1, -1) if no common substring exists.
        /// </returns>
        (string Substring, int PositionInText, int PositionInOther) LongestCommonSubstringInfo(string other);

        /// <summary>
        /// Finds all positions where the longest common substring occurs.
        /// </summary>
        /// <param name="other">The string to compare against.</param>
        /// <returns>
        /// A tuple containing: the substring, all positions in tree's text, all positions in other.
        /// Returns (empty string, empty list, empty list) if no common substring exists.
        /// </returns>
        (string Substring, IReadOnlyList<int> PositionsInText, IReadOnlyList<int> PositionsInOther) FindAllLongestCommonSubstrings(string other);

        /// <summary>
        /// Creates a detailed string representation of the tree structure.
        /// Useful for debugging and visualization.
        /// </summary>
        /// <returns>A multi-line string showing the tree structure.</returns>
        string PrintTree();
    }
}
