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
    /// All tied substrings with all positions: <see cref="ISuffixTree.FindAllLongestRepeatedSubstrings"/>.
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
    /// For every MEM with every text occurrence use <see cref="FindMaximalExactMatches"/>; for MUMs
    /// <see cref="FindMaximalUniqueMatches"/>.
    /// </summary>
    IReadOnlyList<(int PositionInText, int PositionInQuery, int Length)> FindExactMatchAnchors(
        string query, int minLength);

    /// <summary>
    /// Finds all maximal exact matches (MEMs) of length ≥ <paramref name="minLength"/> between this
    /// tree's text (reference) and <paramref name="query"/>: every (PositionInText, PositionInQuery, Length)
    /// that is left- and right-maximal, with every text occurrence — MUMmer 3
    /// <c>mummer -maxmatch -l minLength</c>, forward strand, 0-based. Sorted by query position, then
    /// text position. See <see cref="SuffixTreeAlgorithms.FindMaximalExactMatches{TNode, TNav}"/>.
    /// </summary>
    /// <exception cref="ArgumentNullException"><paramref name="query"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="minLength"/> &lt; 1.</exception>
    /// <remarks>
    /// The default implementation (for implementers outside this library) checks the definition
    /// directly over <see cref="ISuffixTreeSearch.Text"/> in O(|text|·|query|) via
    /// <see cref="SuffixTreeAlgorithms.FindMaximalExactMatchesByDefinition"/>; the library trees
    /// override it with the linear-time suffix-tree algorithm (identical output).
    /// </remarks>
    IReadOnlyList<(int PositionInText, int PositionInQuery, int Length)> FindMaximalExactMatches(
        string query, int minLength)
        => SuffixTreeAlgorithms.FindMaximalExactMatchesByDefinition(DefaultTextOf(this), query, minLength);

    /// <summary>
    /// Finds maximal unique matches (MUMs) of length ≥ <paramref name="minLength"/>: MEMs whose string
    /// occurs exactly once in the text and, for <see cref="MumUniqueness.Both"/>, exactly once in the
    /// query — MUMmer 3 <c>-mum</c> (Both) / <c>-mumreference</c> (Reference), forward strand, 0-based.
    /// Sorted by query position, then text position.
    /// See <see cref="SuffixTreeAlgorithms.FindMaximalUniqueMatches{TNode, TNav}"/>.
    /// </summary>
    /// <exception cref="ArgumentNullException"><paramref name="query"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="minLength"/> &lt; 1 or undefined <paramref name="uniqueness"/>.</exception>
    /// <remarks>
    /// The default implementation checks the definition over <see cref="ISuffixTreeSearch.Text"/>
    /// (<see cref="SuffixTreeAlgorithms.FindMaximalUniqueMatchesByDefinition"/>); the library trees
    /// override it with the linear-time suffix-tree algorithm (identical output).
    /// </remarks>
    IReadOnlyList<(int PositionInText, int PositionInQuery, int Length)> FindMaximalUniqueMatches(
        string query, int minLength, MumUniqueness uniqueness = MumUniqueness.Both)
        => SuffixTreeAlgorithms.FindMaximalUniqueMatchesByDefinition(DefaultTextOf(this), query, minLength, uniqueness);

    private static ITextSource DefaultTextOf(ISuffixTreeAnalysis tree)
        => tree as ISuffixTreeSearch is { } search
            ? search.Text
            : throw new NotSupportedException(
                "The default maximal-match implementation needs ISuffixTreeSearch.Text; override the method.");
}
