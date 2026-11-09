using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.Esa.Tests.Analysis;

/// <summary>The texts that the analysis tests of work package 1b.3 sweep (every text is deterministic; the random ones come from <see cref="SplitMix64"/>).</summary>
internal static class AnalysisCorpus
{
    /// <summary>Hand-picked texts: the literature examples, n = 0..3, ties, no repeat, all equal, code units at the ends of the alphabet.</summary>
    internal static IEnumerable<string> Edge()
    {
        foreach (string s in new[]
        {
            "", "a", "aa", "ab", "aaa", "aab", "aba", "abc", "baa", "aaaa", "abab", "abcd", "aabb", "banana", "mississippi", "abracadabra", "ATATAT", "GATTACAGATTTACA",
            "abcabcxyzxyz", "xyxzyzxz", "abcdbcdacdab", "abcxbcaxcab", "xyzzyxxz", "xabcyabcwabcyz", "acgtacgtacgt", "gattacagattacca", "aaaaaaaabaaaaaaaab",
            "\0", "\0\0", "￿\0￿\0", "\0￿\0￿a", "é中😀é中😀", "𐀀𐀀\uD800", "\uDC00𐀀\uD800", "NNNNACGTNNNN",
        })
            yield return s;
    }

    /// <summary>Every string over <c>sigma</c> letters up to the given length for the three alphabet sizes that the default suite affords.</summary>
    internal static IEnumerable<string> ExhaustiveSmall()
    {
        foreach (var (sigma, maxLen) in new[] { (2, 12), (3, 8), (4, 6) })
            foreach (string s in ExhaustiveStrings.All(sigma, maxLen))
                yield return s;
    }

    /// <summary>The adversarial families at each of the given lengths (see <see cref="StringFamilies.Adversarial"/>); names carry the parameters.</summary>
    internal static IEnumerable<(string Name, string Text)> Adversarial(IEnumerable<int> sizes, ulong seed)
    {
        foreach (int n in sizes)
            foreach (var (name, text) in StringFamilies.Adversarial(n, seed))
                yield return ($"{name}-n{n}", text);
    }

    /// <summary>The sizes of the "every n up to 200" style sweep that is affordable in the default suite: all up to 24, then a sparse ladder to 200.</summary>
    internal static readonly int[] SweepSizes = [.. Enumerable.Range(0, 25), 31, 32, 33, 47, 64, 65, 100, 144, 200];

    /// <summary>Everything the default suite checks against the definition-level references: edge, exhaustive-small, adversarial sweep.</summary>
    internal static IEnumerable<(string Name, string Text)> Definitional()
    {
        foreach (string s in Edge()) yield return ($"edge '{s}'", s);
        foreach (string s in ExhaustiveSmall()) yield return ($"all '{s}'", s);
        foreach (var item in Adversarial(SweepSizes, 0x1B03)) yield return item;
    }
}
