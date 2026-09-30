using System;
using System.Collections.Generic;

namespace SuffixTree;

/// <summary>
/// Analysis and advanced algorithm contract for suffix trees.
/// </summary>
public interface ISuffixTreeAnalysis
{
    /// <summary>
    /// Finds the longest substring that appears at least twice in the text (occurrences may overlap).
    /// On a length tie the returned representative is implementation-specific; its length is always maximal.
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
    /// Finds the longest common substring with position information. On a length tie the substring
    /// first found in <paramref name="other"/> is returned; <c>PositionInOther</c> is its first occurrence
    /// there, <c>PositionInText</c> one (implementation-specific) occurrence in the text; ("", -1, -1) if none.
    /// </summary>
    (string Substring, int PositionInText, int PositionInOther) LongestCommonSubstringInfo(string other);

    /// <summary>
    /// Finds all positions where the longest common substring occurs (the canonical substring of
    /// <see cref="LongestCommonSubstringInfo"/> only); both position lists are 0-based, ascending and duplicate-free.
    /// </summary>
    (string Substring, IReadOnlyList<int> PositionsInText, IReadOnlyList<int> PositionsInOther) FindAllLongestCommonSubstrings(string other);

    /// <summary>
    /// Finds exact-match anchors between this tree's text and a query string
    /// using suffix-link-based streaming traversal: one maximal exact match per maximal run of
    /// matching statistics ≥ <paramref name="minLength"/> (its first peak), ordered by query position.
    /// Anchors of adjacent runs may overlap in the query; minLength ≤ 0 yields an empty list.
    /// </summary>
    IReadOnlyList<(int PositionInText, int PositionInQuery, int Length)> FindExactMatchAnchors(
        string query, int minLength);
}
