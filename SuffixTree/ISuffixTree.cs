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
        /// Finds the longest common substring between this tree's text and another string.
        /// </summary>
        /// <param name="other">The string to compare against.</param>
        /// <returns>The longest common substring, or empty string if none exists.</returns>
        string LongestCommonSubstring(string other);
    }
}
