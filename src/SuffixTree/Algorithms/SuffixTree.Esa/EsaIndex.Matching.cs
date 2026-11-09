#pragma warning disable S2325, S4144, S1144, S3453, CA1822 // stub of the red run (work package 1b.4): every member throws NotImplementedException; the implementation replaces this file
namespace SuffixTree.Esa;

/// <summary>
/// The matching-statistics operations of the in-memory ESA (work package 1b.4): the longest-common-substring family, the anchors, the maximal exact and unique matches, all
/// by the generic algorithms of <c>SuffixTree.Core</c> instantiated on <see cref="EsaNavigator"/>. Signatures and results are those of <c>ISuffixTreeAnalysis</c> (the facade of 1b.6
/// forwards to them). Stub of the red run.
/// </summary>
internal sealed partial class EsaIndex
{
    internal string LongestCommonSubstring(string other) => throw new NotImplementedException();

    internal string LongestCommonSubstring(ReadOnlySpan<char> other) => throw new NotImplementedException();

    internal (string Substring, int PositionInText, int PositionInOther) LongestCommonSubstringInfo(string other) => throw new NotImplementedException();

    internal (string Substring, int PositionInText, int PositionInOther) LongestCommonSubstringInfo(ReadOnlySpan<char> other) => throw new NotImplementedException();

    internal (string Substring, IReadOnlyList<int> PositionsInText, IReadOnlyList<int> PositionsInOther) FindAllLongestCommonSubstrings(string other) => throw new NotImplementedException();

    internal IReadOnlyList<(string Substring, IReadOnlyList<int> PositionsInText, IReadOnlyList<int> PositionsInOther)> FindAllDistinctLongestCommonSubstrings(string other)
        => throw new NotImplementedException();

    internal IReadOnlyList<(int PositionInText, int PositionInQuery, int Length)> FindExactMatchAnchors(string query, int minLength) => throw new NotImplementedException();

    internal IReadOnlyList<(int PositionInText, int PositionInQuery, int Length)> FindMaximalExactMatches(string query, int minLength) => throw new NotImplementedException();

    internal IReadOnlyList<(int PositionInText, int PositionInQuery, int Length)> FindMaximalUniqueMatches(string query, int minLength, MumUniqueness uniqueness = MumUniqueness.Both)
        => throw new NotImplementedException();
}
