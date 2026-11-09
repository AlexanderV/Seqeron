namespace SuffixTree.Esa;

/// <summary>
/// The matching-statistics operations of the in-memory ESA (work package 1b.4): the longest-common-substring family, the anchors, the maximal exact and unique matches, all by the generic
/// algorithms of <c>SuffixTree.Core</c> (<see cref="SuffixTreeAlgorithms"/>) instantiated on <see cref="EsaNavigator"/>, as <c>PersistentSuffixTree.Analysis.cs</c> does for the v6 tree. Signatures
/// and results are those of <c>ISuffixTreeAnalysis</c> (the facade of 1b.6 forwards to them). They need suffix links, so the first call builds the navigation pack (<see cref="Navigation"/>); the
/// arguments are validated BEFORE that, so an invalid call never pays for it. Nothing is rewritten: MEM, MUM, the anchors and the LCS lists are the code that runs on the Ukkonen tree and on
/// the v6 tree, with the same output and the same order; the occurrence reported by the anchors and by <see cref="LongestCommonSubstringInfo(string)"/> is the one the navigator offers first
/// (<c>SA[Lb]</c> of the match locus), which the contract leaves open.
/// </summary>
internal sealed partial class EsaIndex
{
    /// <summary>The longest common substring of the text and <paramref name="other"/> (on a length tie the one that occurs first in <paramref name="other"/>); empty when they share no character.</summary>
    /// <exception cref="ArgumentNullException"><paramref name="other"/> is null.</exception>
    internal string LongestCommonSubstring(string other) => LongestCommonSubstringInfo(other).Substring;

    /// <summary>The same for a character span.</summary>
    internal string LongestCommonSubstring(ReadOnlySpan<char> other) => LongestCommonSubstringInfo(new string(other)).Substring;

    /// <summary>
    /// The longest common substring with its positions: <c>PositionInOther</c> is its first occurrence in <paramref name="other"/>, <c>PositionInText</c> one occurrence in the text; <c>("", -1, -1)</c>
    /// when there is none.
    /// </summary>
    /// <exception cref="ArgumentNullException"><paramref name="other"/> is null.</exception>
    internal (string Substring, int PositionInText, int PositionInOther) LongestCommonSubstringInfo(string other)
    {
        var results = FindAllLcsInternal(other, firstOnly: true);
        if (results.PositionsInText.Count == 0) return (string.Empty, -1, -1);
        return (results.Substring, results.PositionsInText[0], results.PositionsInOther[0]);
    }

    /// <summary>The same for a character span.</summary>
    internal (string Substring, int PositionInText, int PositionInOther) LongestCommonSubstringInfo(ReadOnlySpan<char> other) => LongestCommonSubstringInfo(new string(other));

    /// <summary>All positions of the canonical longest common substring in the text and in <paramref name="other"/> (both ascending and duplicate-free).</summary>
    /// <exception cref="ArgumentNullException"><paramref name="other"/> is null.</exception>
    internal (string Substring, IReadOnlyList<int> PositionsInText, IReadOnlyList<int> PositionsInOther) FindAllLongestCommonSubstrings(string other)
    {
        var results = FindAllLcsInternal(other, firstOnly: false);
        return (results.Substring, results.PositionsInText, results.PositionsInOther);
    }

    /// <summary>Every distinct longest common substring (all length ties) with all its positions, ordered by first occurrence in <paramref name="other"/>.</summary>
    /// <exception cref="ArgumentNullException"><paramref name="other"/> is null.</exception>
    internal IReadOnlyList<(string Substring, IReadOnlyList<int> PositionsInText, IReadOnlyList<int> PositionsInOther)> FindAllDistinctLongestCommonSubstrings(string other)
    {
        ArgumentNullException.ThrowIfNull(other);
        EsaNavigator nav = CreateNavigator();
        return SuffixTreeAlgorithms.FindAllDistinctLongestCommonSubstrings<EsaNode, EsaNavigator>(ref nav, other);
    }

    /// <summary>One maximal exact match per maximal run of matching statistics of at least <paramref name="minLength"/> (its first peak), ordered by query position; empty for <c>minLength &lt;= 0</c>.</summary>
    /// <exception cref="ArgumentNullException"><paramref name="query"/> is null.</exception>
    internal IReadOnlyList<(int PositionInText, int PositionInQuery, int Length)> FindExactMatchAnchors(string query, int minLength)
    {
        ArgumentNullException.ThrowIfNull(query);
        EsaNavigator nav = CreateNavigator();
        return SuffixTreeAlgorithms.FindExactMatchAnchors<EsaNode, EsaNavigator>(ref nav, query, minLength);
    }

    /// <summary>All maximal exact matches of length at least <paramref name="minLength"/> (MUMmer 3 <c>-maxmatch</c>), ordered by query position, then text position.</summary>
    /// <exception cref="ArgumentNullException"><paramref name="query"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="minLength"/> is less than 1.</exception>
    internal IReadOnlyList<(int PositionInText, int PositionInQuery, int Length)> FindMaximalExactMatches(string query, int minLength)
    {
        ArgumentNullException.ThrowIfNull(query);
        ArgumentOutOfRangeException.ThrowIfLessThan(minLength, 1);
        EsaNavigator nav = CreateNavigator();
        return SuffixTreeAlgorithms.FindMaximalExactMatches<EsaNode, EsaNavigator>(ref nav, query, minLength);
    }

    /// <summary>The maximal unique matches of length at least <paramref name="minLength"/> (MUMmer 3 <c>-mum</c> for <see cref="MumUniqueness.Both"/>, <c>-mumreference</c> for <see cref="MumUniqueness.Reference"/>).</summary>
    /// <exception cref="ArgumentNullException"><paramref name="query"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="minLength"/> is less than 1, or <paramref name="uniqueness"/> is not a defined value.</exception>
    internal IReadOnlyList<(int PositionInText, int PositionInQuery, int Length)> FindMaximalUniqueMatches(string query, int minLength, MumUniqueness uniqueness = MumUniqueness.Both)
    {
        ArgumentNullException.ThrowIfNull(query);
        ArgumentOutOfRangeException.ThrowIfLessThan(minLength, 1);
        if (uniqueness is not (MumUniqueness.Both or MumUniqueness.Reference))
            throw new ArgumentOutOfRangeException(nameof(uniqueness), uniqueness, "Undefined MumUniqueness value.");
        EsaNavigator nav = CreateNavigator();
        return SuffixTreeAlgorithms.FindMaximalUniqueMatches<EsaNode, EsaNavigator>(ref nav, query, minLength, uniqueness);
    }

    private (string Substring, List<int> PositionsInText, List<int> PositionsInOther) FindAllLcsInternal(string other, bool firstOnly)
    {
        ArgumentNullException.ThrowIfNull(other);
        EsaNavigator nav = CreateNavigator();
        return SuffixTreeAlgorithms.FindAllLcs<EsaNode, EsaNavigator>(ref nav, other, firstOnly);
    }
}
