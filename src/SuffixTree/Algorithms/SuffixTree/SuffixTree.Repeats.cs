using System;
using System.Collections.Generic;

namespace SuffixTree;

public partial class SuffixTree
{
    /// <inheritdoc />
    public IReadOnlyList<(string Substring, IReadOnlyList<int> PositionsInText, IReadOnlyList<int> PositionsInOther)>
        FindAllDistinctLongestCommonSubstrings(string other)
    {
        var nav = new SuffixTreeNavigator(this);
        return SuffixTreeAlgorithms.FindAllDistinctLongestCommonSubstrings<SuffixTreeNode, SuffixTreeNavigator>(ref nav, other);
    }

    /// <inheritdoc cref="ISuffixTree.FindMaximalRepeatedPairs"/>
    public IReadOnlyList<(int FirstPosition, int SecondPosition, int Length)> FindMaximalRepeatedPairs(
        int minLength, Func<char, bool>? isUniqueSymbol = null)
        => SuffixTreeAlgorithms.FindMaximalRepeatedPairs(this, minLength, isUniqueSymbol);

    /// <summary>
    /// Every distinct longest substring occurring in at least <paramref name="minSupport"/> of
    /// <paramref name="texts"/> (default: all of them — the longest common substring of k strings,
    /// Rosalind LCSM; Gusfield 1997 §7.6 k-common substring), sorted ordinally, via a generalized suffix
    /// tree. See <see cref="SuffixTreeAlgorithms.FindLongestCommonSubstrings"/>.
    /// </summary>
    public static IReadOnlyList<string> FindLongestCommonSubstrings(IReadOnlyList<string> texts, int? minSupport = null)
    {
        ArgumentNullException.ThrowIfNull(texts);
        return SuffixTreeAlgorithms.FindLongestCommonSubstrings(texts, minSupport ?? texts.Count, static s => Build(s));
    }

    /// <summary>
    /// Gusfield's l(q) for every q = 1..k: element q is the length of the longest substring occurring in at
    /// least q of <paramref name="texts"/> (element 0 = element 1). See
    /// <see cref="SuffixTreeAlgorithms.LongestCommonSubstringLengthsBySupport"/>.
    /// </summary>
    public static int[] LongestCommonSubstringLengthsBySupport(IReadOnlyList<string> texts)
        => SuffixTreeAlgorithms.LongestCommonSubstringLengthsBySupport(texts, static s => Build(s));
}
