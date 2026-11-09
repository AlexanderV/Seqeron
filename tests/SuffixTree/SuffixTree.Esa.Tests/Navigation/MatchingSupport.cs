using System.Text;
using SuffixTree.Esa.Tests.TestSupport;
using UkkonenTree = global::SuffixTree.SuffixTree;

namespace SuffixTree.Esa.Tests.Navigation;

using Triple = (int PositionInText, int PositionInQuery, int Length);

/// <summary>Definition-level oracles and the shared checks of the matching operations of work package 1b.4 (MEM, MUM, anchors, the LCS family).</summary>
internal static class MatchingChecks
{
    internal static List<Triple> Triples(IReadOnlyList<Triple> list) => [.. list];

    /// <summary>ms[i]: the length of the longest suffix of <c>query[0..i]</c> that occurs in <paramref name="text"/> (brute force; a suffix of an occurring string occurs, so the lengths that work are a prefix of 1, 2, 3, ...).</summary>
    internal static int[] MatchingStatistics(string text, string query)
    {
        var ms = new int[query.Length];
        for (int i = 0; i < query.Length; i++)
        {
            int best = 0;
            for (int length = 1; length <= i + 1; length++)
            {
                if (text.AsSpan().IndexOf(query.AsSpan(i - length + 1, length), StringComparison.Ordinal) < 0) break;
                best = length;
            }
            ms[i] = best;
        }
        return ms;
    }

    /// <summary>The (query position, length) pairs of the anchors by their definition (<c>FindExactMatchAnchors</c> remarks): one per maximal run of ms &gt;= minLength, the first peak.</summary>
    internal static List<(int PositionInQuery, int Length)> AnchorShapes(string text, string query, int minLength)
    {
        var result = new List<(int, int)>();
        if (query.Length == 0 || text.Length == 0 || minLength <= 0) return result;
        int[] ms = MatchingStatistics(text, query);
        int peakLength = 0, peakEnd = -1;
        for (int i = 0; i < ms.Length; i++)
        {
            if (ms[i] >= minLength)
            {
                if (ms[i] > peakLength)
                {
                    peakLength = ms[i];
                    peakEnd = i;
                }
            }
            else if (peakLength >= minLength)
            {
                result.Add((peakEnd - peakLength + 1, peakLength));
                peakLength = 0;
                peakEnd = -1;
            }
        }
        if (peakLength >= minLength) result.Add((peakEnd - peakLength + 1, peakLength));
        return result;
    }

    internal static string Short(string s) => s.Length > 50 ? s[..50] + $"...({s.Length})" : s;

    /// <summary>Anchors: the (query, length) pairs equal the definition and the Ukkonen tree's; the text position is a valid occurrence (the contract leaves WHICH occurrence open).</summary>
    internal static void CheckAnchors(EsaIndex esa, UkkonenTree? ukkonen, string text, string query, int minLength, string context)
    {
        IReadOnlyList<Triple> actual = esa.FindExactMatchAnchors(query, minLength);
        var shapes = AnchorShapes(text, query, minLength);
        if (actual.Count != shapes.Count || !actual.Select(a => (a.PositionInQuery, a.Length)).SequenceEqual(shapes))
            Assert.Fail($"{context} '{Short(text)}' / '{Short(query)}' L {minLength}: anchors {string.Join(" ", actual)} against the definition {string.Join(" ", shapes)}");
        foreach (Triple a in actual)
            if (string.CompareOrdinal(text, a.PositionInText, query, a.PositionInQuery, a.Length) != 0)
                Assert.Fail($"{context} '{Short(text)}' / '{Short(query)}' L {minLength}: anchor {a} is not an occurrence");
        if (ukkonen is not null)
        {
            var reference = ukkonen.FindExactMatchAnchors(query, minLength).Select(a => (a.PositionInQuery, a.Length)).ToList();
            if (!actual.Select(a => (a.PositionInQuery, a.Length)).SequenceEqual(reference))
                Assert.Fail($"{context} '{Short(text)}' / '{Short(query)}' L {minLength}: anchors differ from the Ukkonen tree's");
        }
    }

    internal static void CheckMemMum(EsaIndex esa, UkkonenTree? ukkonen, string text, string query, int minLength, string context)
    {
        var source = new StringTextSource(text);
        var mems = Triples(esa.FindMaximalExactMatches(query, minLength));
        var memsExpected = Triples(SuffixTreeAlgorithms.FindMaximalExactMatchesByDefinition(source, query, minLength));
        if (!mems.SequenceEqual(memsExpected))
            Assert.Fail($"{context} '{Short(text)}' / '{Short(query)}' L {minLength}: MEM {string.Join(" ", mems)} against the definition {string.Join(" ", memsExpected)}");
        if (ukkonen is not null && !mems.SequenceEqual(Triples(ukkonen.FindMaximalExactMatches(query, minLength))))
            Assert.Fail($"{context} '{Short(text)}' / '{Short(query)}' L {minLength}: MEM differs from the Ukkonen tree's");
        foreach (MumUniqueness mode in new[] { MumUniqueness.Both, MumUniqueness.Reference })
        {
            var mums = Triples(esa.FindMaximalUniqueMatches(query, minLength, mode));
            var mumsExpected = Triples(SuffixTreeAlgorithms.FindMaximalUniqueMatchesByDefinition(source, query, minLength, mode));
            if (!mums.SequenceEqual(mumsExpected))
                Assert.Fail($"{context} '{Short(text)}' / '{Short(query)}' L {minLength} {mode}: MUM {string.Join(" ", mums)} against the definition {string.Join(" ", mumsExpected)}");
            if (ukkonen is not null && !mums.SequenceEqual(Triples(ukkonen.FindMaximalUniqueMatches(query, minLength, mode))))
                Assert.Fail($"{context} '{Short(text)}' / '{Short(query)}' L {minLength} {mode}: MUM differs from the Ukkonen tree's");
        }
    }

    /// <summary>The LCS family: the distinct list against the definition, the canonical list, the info and the strings against the first distinct entry; Ukkonen's lists when given.</summary>
    internal static void CheckLcs(EsaIndex esa, UkkonenTree? ukkonen, string text, string query, string context)
    {
        string where = $"{context} '{Short(text)}' / '{Short(query)}'";
        var source = new StringTextSource(text);
        var expected = SuffixTreeAlgorithms.FindAllDistinctLongestCommonSubstringsByDefinition(source, query);
        var distinct = esa.FindAllDistinctLongestCommonSubstrings(query);
        if (distinct.Count != expected.Count) Assert.Fail($"{where}: {distinct.Count} distinct longest common substrings, expected {expected.Count}");
        for (int i = 0; i < expected.Count; i++)
        {
            if (distinct[i].Substring != expected[i].Substring
                || !distinct[i].PositionsInText.SequenceEqual(expected[i].PositionsInText)
                || !distinct[i].PositionsInOther.SequenceEqual(expected[i].PositionsInOther))
                Assert.Fail($"{where}: distinct entry {i}: '{distinct[i].Substring}' [{string.Join(",", distinct[i].PositionsInText)}] [{string.Join(",", distinct[i].PositionsInOther)}] against '{expected[i].Substring}' [{string.Join(",", expected[i].PositionsInText)}] [{string.Join(",", expected[i].PositionsInOther)}]");
        }

        var all = esa.FindAllLongestCommonSubstrings(query);
        if (expected.Count == 0)
        {
            if (all.Substring.Length != 0 || all.PositionsInText.Count != 0 || all.PositionsInOther.Count != 0) Assert.Fail($"{where}: no common character, got '{all.Substring}'");
            var none = esa.LongestCommonSubstringInfo(query);
            if (none != ("", -1, -1)) Assert.Fail($"{where}: info of no match is {none}");
            if (esa.LongestCommonSubstring(query).Length != 0 || esa.LongestCommonSubstring(query.AsSpan()).Length != 0) Assert.Fail($"{where}: LCS of no match is not empty");
            if (esa.LongestCommonSubstringInfo(query.AsSpan()) != ("", -1, -1)) Assert.Fail($"{where}: span info of no match");
            return;
        }

        var first = expected[0];
        if (all.Substring != first.Substring || !all.PositionsInText.SequenceEqual(first.PositionsInText) || !all.PositionsInOther.SequenceEqual(first.PositionsInOther))
            Assert.Fail($"{where}: FindAllLongestCommonSubstrings differs from the first distinct entry");
        foreach (var info in new[] { esa.LongestCommonSubstringInfo(query), esa.LongestCommonSubstringInfo(query.AsSpan()) })
        {
            if (info.Substring != first.Substring || info.PositionInOther != first.PositionsInOther[0] || !first.PositionsInText.Contains(info.PositionInText))
                Assert.Fail($"{where}: info {info}: the substring and its first position in the other string are fixed, the text position must be an occurrence (of {string.Join(",", first.PositionsInText)})");
        }
        if (esa.LongestCommonSubstring(query) != first.Substring || esa.LongestCommonSubstring(query.AsSpan()) != first.Substring) Assert.Fail($"{where}: LongestCommonSubstring");

        if (ukkonen is not null)
        {
            var uAll = ukkonen.FindAllLongestCommonSubstrings(query);
            if (uAll.Substring != all.Substring || !uAll.PositionsInText.SequenceEqual(all.PositionsInText) || !uAll.PositionsInOther.SequenceEqual(all.PositionsInOther))
                Assert.Fail($"{where}: FindAllLongestCommonSubstrings differs from the Ukkonen tree's");
            var uInfo = ukkonen.LongestCommonSubstringInfo(query);
            var info = esa.LongestCommonSubstringInfo(query);
            if (uInfo.Substring != info.Substring || uInfo.PositionInOther != info.PositionInOther) Assert.Fail($"{where}: LongestCommonSubstringInfo differs from the Ukkonen tree's");
        }
    }

    internal static void CheckAll(EsaIndex esa, UkkonenTree? ukkonen, string text, string query, int[] minLengths, string context)
    {
        CheckLcs(esa, ukkonen, text, query, context);
        foreach (int minLength in minLengths)
        {
            CheckMemMum(esa, ukkonen, text, query, minLength, context);
            CheckAnchors(esa, ukkonen, text, query, minLength, context);
        }
    }

    /// <summary>A copy of <paramref name="source"/> with substitutions (from its own letters plus <paramref name="extra"/>), deletions and insertions at the given rate.</summary>
    internal static string Mutate(string source, double rate, SplitMix64 rng, string extra = "")
    {
        if (source.Length == 0) return source;
        string letters = new string(source.Distinct().ToArray()) + extra;
        var sb = new StringBuilder(source.Length + 8);
        foreach (char c in source)
        {
            if (rng.NextInt(1_000_000) >= rate * 1_000_000)
            {
                sb.Append(c);
                continue;
            }
            switch (rng.NextInt(3))
            {
                case 0: sb.Append(letters[rng.NextInt(letters.Length)]); break;
                case 1: break;
                default: sb.Append(c).Append(letters[rng.NextInt(letters.Length)]); break;
            }
        }
        return sb.ToString();
    }
}
