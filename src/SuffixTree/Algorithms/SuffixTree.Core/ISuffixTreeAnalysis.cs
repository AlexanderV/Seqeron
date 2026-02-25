using System;
using System.Collections.Generic;

namespace SuffixTree;

/// <summary>
/// Analysis and advanced algorithm contract for suffix trees.
/// </summary>
public interface ISuffixTreeAnalysis
{
    /// <summary>
    /// Finds the longest substring that appears at least twice in the text.
    /// </summary>
    string LongestRepeatedSubstring();

    /// <summary>
    /// Returns the longest repeated substring as memory.
    /// Implementations may override for true zero-copy access to the original text.
    /// The default interface implementation wraps <see cref="LongestRepeatedSubstring"/>.
    /// </summary>
    ReadOnlyMemory<char> LongestRepeatedSubstringMemory() => LongestRepeatedSubstring().AsMemory();

    /// <summary>
    /// Returns all suffixes of the original string in sorted order.
    /// </summary>
    IReadOnlyList<string> GetAllSuffixes();

    /// <summary>
    /// Enumerates all suffixes of the original string in sorted order lazily.
    /// </summary>
    IEnumerable<string> EnumerateSuffixes();

    /// <summary>
    /// Finds the longest common substring between this tree's text and another string.
    /// </summary>
    string LongestCommonSubstring(string other);

    /// <summary>
    /// Finds the longest common substring between this tree's text and another character span.
    /// </summary>
    string LongestCommonSubstring(ReadOnlySpan<char> other);

    /// <summary>
    /// Finds the longest common substring with position information.
    /// </summary>
    (string Substring, int PositionInText, int PositionInOther) LongestCommonSubstringInfo(string other);

    /// <summary>
    /// Finds all positions where the longest common substring occurs.
    /// </summary>
    (string Substring, IReadOnlyList<int> PositionsInText, IReadOnlyList<int> PositionsInOther) FindAllLongestCommonSubstrings(string other);

    /// <summary>
    /// Finds exact-match anchors between this tree's text and a query string
    /// using suffix-link-based streaming traversal.
    /// </summary>
    IReadOnlyList<(int PositionInText, int PositionInQuery, int Length)> FindExactMatchAnchors(
        string query, int minLength);
}
