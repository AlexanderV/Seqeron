using SuffixTree.Esa.Tests.Analysis;
using SuffixTree.Esa.Tests.TestSupport;
using SuffixTree.Indexing;
using UkkonenTree = global::SuffixTree.SuffixTree;

namespace SuffixTree.Esa.Tests.Navigation;

using Triple = (int PositionInText, int PositionInQuery, int Length);

/// <summary>
/// The matching-statistics operations of the in-memory ESA (work package 1b.4): the generic algorithms of <c>SuffixTree.Core</c> instantiated on <see cref="EsaNavigator"/> (LCS family, anchors,
/// MEM, MUM). Oracles: the definition-level implementations of <c>SuffixTree.Core</c> (<c>FindMaximalExactMatchesByDefinition</c>, <c>FindMaximalUniqueMatchesByDefinition</c>,
/// <c>FindAllDistinctLongestCommonSubstringsByDefinition</c>), a brute-force matching-statistics oracle for the anchors, and the Ukkonen tree. MEM, MUM and the lists of the LCS family are compared for
/// equality (their order is part of the contract); the occurrence that <c>FindExactMatchAnchors</c> and <c>LongestCommonSubstringInfo</c> report is only checked to be one (the contract leaves it open).
/// </summary>
[TestFixture]
public class EsaMatchingTests
{
    private static readonly int[] SmallMinLengths = [1, 2, 3];

    // ---- MUMmer 3.23 literal outputs (MaximalMatchTests of the Ukkonen suite, 0-based) ----

    [Test]
    public void MummerParity_GattacaRepeat()
    {
        EsaIndex esa = EsaIndex.Build("GATTACAGATTACA");
        Assert.That(esa.FindMaximalExactMatches("ATTACAGGATTACAT", 3), Is.EqualTo(new Triple[] { (1, 0, 7), (8, 0, 6), (0, 7, 7), (7, 7, 7) }));
        Assert.That(esa.FindMaximalUniqueMatches("ATTACAGGATTACAT", 3, MumUniqueness.Both), Is.EqualTo(new Triple[] { (1, 0, 7) }));
        Assert.That(esa.FindMaximalUniqueMatches("ATTACAGGATTACAT", 3, MumUniqueness.Reference), Is.EqualTo(new Triple[] { (1, 0, 7) }));
    }

    [Test]
    public void MummerParity_TandemRepeat()
    {
        EsaIndex esa = EsaIndex.Build("ACGTACGTTACGT");
        Assert.That(esa.FindMaximalExactMatches("TACGTACGTT", 3), Is.EqualTo(new Triple[] { (3, 0, 5), (8, 0, 5), (0, 1, 9), (8, 4, 5), (0, 5, 4) }));
        Assert.That(esa.FindMaximalUniqueMatches("TACGTACGTT", 3, MumUniqueness.Both), Is.EqualTo(new Triple[] { (0, 1, 9) }));
        Assert.That(esa.FindMaximalUniqueMatches("TACGTACGTT", 3, MumUniqueness.Reference), Is.EqualTo(new Triple[] { (0, 1, 9) }));
    }

    [Test]
    public void MummerParity_ThreeSetsDiffer()
    {
        EsaIndex esa = EsaIndex.Build("TTTCCTCATGCAATTCAAAACCATGTCCGTAATGTAGGCG");
        const string Query = "TTTCCTCATGCAAATAGCCTCATGCATCCGTAATGTAG";
        Assert.That(esa.FindMaximalExactMatches(Query, 4), Is.EqualTo(new Triple[] { (0, 0, 13), (21, 6, 4), (15, 10, 4), (3, 17, 9), (21, 20, 4), (25, 26, 12), (22, 32, 4) }));
        Assert.That(esa.FindMaximalUniqueMatches(Query, 4, MumUniqueness.Reference), Is.EqualTo(new Triple[] { (0, 0, 13), (15, 10, 4), (3, 17, 9), (25, 26, 12) }));
        Assert.That(esa.FindMaximalUniqueMatches(Query, 4, MumUniqueness.Both), Is.EqualTo(new Triple[] { (0, 0, 13), (15, 10, 4), (25, 26, 12) }));
    }

    [Test]
    public void MummerParity_LengthOne_MumAtReferenceStartIsKept()
    {
        EsaIndex esa = EsaIndex.Build("ACGT");
        Assert.That(esa.FindMaximalExactMatches("TTAG", 1), Is.EqualTo(new Triple[] { (3, 0, 1), (3, 1, 1), (0, 2, 1), (2, 3, 1) }));
        Assert.That(esa.FindMaximalUniqueMatches("TTAG", 1, MumUniqueness.Reference), Is.EqualTo(new Triple[] { (3, 0, 1), (3, 1, 1), (0, 2, 1), (2, 3, 1) }));
        Assert.That(esa.FindMaximalUniqueMatches("TTAG", 1, MumUniqueness.Both), Is.EqualTo(new Triple[] { (0, 2, 1), (2, 3, 1) }));
    }

    [Test]
    public void Definitions_RepeatsAndHomopolymers()
    {
        Assert.That(EsaIndex.Build("ACGTACGTACGT").FindMaximalExactMatches("ACG", 3), Is.EqualTo(new Triple[] { (0, 0, 3), (4, 0, 3), (8, 0, 3) }));
        Assert.That(EsaIndex.Build("ACGTACGTACGT").FindMaximalUniqueMatches("ACG", 3, MumUniqueness.Reference), Is.Empty);
        Assert.That(EsaIndex.Build("AAAA").FindMaximalExactMatches("AA", 1), Is.EqualTo(new Triple[] { (0, 0, 2), (1, 0, 2), (2, 0, 2), (3, 0, 1), (0, 1, 1) }));
        Assert.That(EsaIndex.Build("CCGATTACACC").FindMaximalUniqueMatches("GATTACATTGATTACA", 5, MumUniqueness.Reference), Is.EqualTo(new Triple[] { (2, 0, 7), (2, 9, 7) }));
        Assert.That(EsaIndex.Build("CCGATTACACC").FindMaximalUniqueMatches("GATTACATTGATTACA", 5, MumUniqueness.Both), Is.Empty);
        Assert.That(EsaIndex.Build("CCGATTACACC").FindMaximalUniqueMatches("GATTACATTGATTACA", 5), Is.Empty, "the default mode is Both");
    }

    [Test]
    public void Lcs_LiteratureAndDistinctTies()
    {
        EsaIndex abc = EsaIndex.Build("abcxyzabc");
        var all = abc.FindAllDistinctLongestCommonSubstrings("xyzqabcpxyz");
        Assert.That(all.Select(r => (r.Substring, r.PositionsInText.ToArray(), r.PositionsInOther.ToArray())).ToArray(), Is.EqualTo(new[]
        {
            ("xyz", new[] { 3 }, new[] { 0, 8 }),
            ("abc", new[] { 0, 6 }, new[] { 4 }),
        }));
        var gattaca = EsaIndex.Build("GATTACA").FindAllDistinctLongestCommonSubstrings("TAGACCA");
        Assert.That(gattaca.Select(r => (r.Substring, r.PositionsInText.ToArray(), r.PositionsInOther.ToArray())).ToArray(), Is.EqualTo(new[]
        {
            ("TA", new[] { 3 }, new[] { 0 }),
            ("GA", new[] { 0 }, new[] { 2 }),
            ("AC", new[] { 4 }, new[] { 3 }),
            ("CA", new[] { 5 }, new[] { 5 }),
        }));
        Assert.That(EsaIndex.Build("banana").LongestCommonSubstring("ananas"), Is.EqualTo("anana"));
        Assert.That(EsaIndex.Build("hello").LongestCommonSubstring("hello"), Is.EqualTo("hello"));
        Assert.That(EsaIndex.Build("abc").LongestCommonSubstring("xay"), Is.EqualTo("a"));
        Assert.That(EsaIndex.Build("xyzdef").LongestCommonSubstring("abcdef"), Is.EqualTo("def"));
        Assert.That(EsaIndex.Build("xxxabcyyy").LongestCommonSubstringInfo("zzzabcwww"), Is.EqualTo(("abc", 3, 3)));
    }

    [Test]
    public void Anchors_Literature()
    {
        EsaIndex esa = EsaIndex.Build("abcdefghij");
        var anchors = esa.FindExactMatchAnchors("defgh", 3);
        Assert.That(anchors, Has.Count.EqualTo(1));
        Assert.That(anchors[0], Is.EqualTo((3, 0, 5)));
        Assert.That(EsaIndex.Build("abcdef").FindExactMatchAnchors("abcdef", 3), Is.EqualTo(new Triple[] { (0, 0, 6) }));
        // text "aba", query "ababa", minLength 3: adjacent runs overlap by fewer than minLength characters (FindExactMatchAnchors remarks)
        Assert.That(EsaIndex.Build("aba").FindExactMatchAnchors("ababa", 3).Select(a => (a.PositionInQuery, a.Length)).ToArray(), Is.EqualTo(new[] { (0, 3), (2, 3) }));
        Assert.That(EsaIndex.Build("ACGTACGTTACGT").FindExactMatchAnchors("TACGTACGTT", 3), Is.SubsetOf(EsaIndex.Build("ACGTACGTTACGT").FindMaximalExactMatches("TACGTACGTT", 3)));
    }

    // ---- exhaustive small ----

    [Test]
    public void ExhaustiveBinaryPairs_TextsUpTo6_QueriesUpTo6()
    {
        string[] queries = ExhaustiveStrings.All(2, 6).ToArray();
        foreach (string text in ExhaustiveStrings.All(2, 6))
        {
            EsaIndex esa = EsaIndex.Build(text);
            UkkonenTree ukkonen = UkkonenTree.Build(text);
            foreach (string query in queries) MatchingChecks.CheckAll(esa, ukkonen, text, query, SmallMinLengths, "bin");
        }
    }

    [Test]
    public void ExhaustiveTernaryPairs_TextsUpTo4_QueriesUpTo4()
    {
        string[] queries = ExhaustiveStrings.All(3, 4).ToArray();
        foreach (string text in ExhaustiveStrings.All(3, 4))
        {
            EsaIndex esa = EsaIndex.Build(text);
            UkkonenTree ukkonen = UkkonenTree.Build(text);
            foreach (string query in queries) MatchingChecks.CheckAll(esa, ukkonen, text, query, SmallMinLengths, "ter");
        }
    }

    [Test]
    public void ForeignSymbolsInTheQuery()
    {
        // the text uses {a, b}, the query {a, b, c}: c never matches, a match ends in front of it
        string[] queries = ExhaustiveStrings.All(3, 5).ToArray();
        foreach (string text in ExhaustiveStrings.All(2, 5))
        {
            EsaIndex esa = EsaIndex.Build(text);
            UkkonenTree ukkonen = UkkonenTree.Build(text);
            foreach (string query in queries) MatchingChecks.CheckAll(esa, ukkonen, text, query, [1, 2], "foreign");
        }
    }

    [Test]
    public void EdgeCases_EmptyAndShortInputs()
    {
        foreach (string text in new[] { "", "a", "ab", "aa" })
        {
            EsaIndex esa = EsaIndex.Build(text);
            foreach (string query in new[] { "", "a", "b", "ab", "aaa", "ba" })
                MatchingChecks.CheckAll(esa, UkkonenTree.Build(text), text, query, [1, 2, 5], $"edge '{text}'");
        }
        foreach (string text in AnalysisCorpus.Edge())
        {
            EsaIndex esa = EsaIndex.Build(text);
            UkkonenTree ukkonen = UkkonenTree.Build(text);
            foreach (string query in new[] { text, text.Length > 2 ? text[1..^1] : "x", "\0", "￿a" })
                MatchingChecks.CheckAll(esa, ukkonen, text, query, [1, 3], "edge corpus");
        }
    }

    [Test]
    public void Unicode_SurrogatesAndLargeCodeUnits()
    {
        foreach (string text in new[] { "é中😀é中😀", "𐀀𐀀\uD800", "\uDC00𐀀\uD800", "\0￿\0￿a", StringFamilies.SurrogatePairs(120, 5), StringFamilies.Bmp(200, 6), StringFamilies.NulAndFfff(150, 7) })
        {
            EsaIndex esa = EsaIndex.Build(text);
            UkkonenTree ukkonen = UkkonenTree.Build(text);
            string query = text.Length > 8 ? string.Concat(text.AsSpan(2, text.Length - 5), "é") : text;
            MatchingChecks.CheckAll(esa, ukkonen, text, query, [1, 2, 4], "unicode");
        }
    }

    // ---- adversarial families, random, two-text ----

    [Test]
    public void AdversarialFamilies_WithMutatedAndForeignQueries()
    {
        var rng = new SplitMix64(0x1B4B30);
        foreach (int n in new[] { 3, 8, 21, 55, 120 })
            foreach (var (name, text) in StringFamilies.Adversarial(n, 41))
            {
                EsaIndex esa = EsaIndex.Build(text);
                UkkonenTree ukkonen = UkkonenTree.Build(text);
                var queries = new List<string>
                {
                    text,
                    MatchingChecks.Mutate(text, 0.05, rng),
                    MatchingChecks.Mutate(text, 0.2, rng, "xyz"),
                    text.Length > 4 ? text.Substring(1, text.Length / 2) : text,
                    new string(text.Reverse().ToArray()),
                    StringFamilies.Adversarial(n, 42).ElementAt(rng.NextInt(StringFamilies.Adversarial(n, 42).Count())).Text,
                };
                foreach (string query in queries) MatchingChecks.CheckAll(esa, ukkonen, text, query, [1, 3, 8], $"{name} n {n}");
            }
    }

    [Test]
    public void RandomPairs_AllAlphabets()
    {
        var rng = new SplitMix64(0x1B4B31);
        foreach (int sigma in new[] { 1, 2, 3, 4, 20, 256 })
            for (int rep = 0; rep < 20; rep++)
            {
                int n = 1 + rng.NextInt(300);
                string text = StringFamilies.Random(n, sigma, rng.Next());
                string query = (rep % 4) switch
                {
                    0 => StringFamilies.Random(1 + rng.NextInt(200), sigma, rng.Next()),
                    1 => MatchingChecks.Mutate(text, 0.03, rng),
                    2 => MatchingChecks.Mutate(Window(text, rng), 0.1, rng),
                    _ => StringFamilies.Random(1 + rng.NextInt(60), sigma + 1, rng.Next()),
                };
                MatchingChecks.CheckAll(EsaIndex.Build(text), UkkonenTree.Build(text), text, query, [1, 2, 4, 9], $"random sigma {sigma} #{rep}");
            }
    }

    private static string Window(string text, SplitMix64 rng)
    {
        int start = rng.NextInt(text.Length);
        return text.Substring(start, 1 + rng.NextInt(Math.Min(text.Length - start, 150)));
    }

    [Test]
    public void LargerTexts_AgainstUkkonen()
    {
        var rng = new SplitMix64(0x1B4B32);
        foreach (var (name, text) in new[]
        {
            ("dna-20000", StringFamilies.RandomDna(20000, 51)),
            ("protein-10000", StringFamilies.RandomProtein(10000, 52)),
            ("fibonacci-8000", StringFamilies.Fibonacci(8000)),
            ("periodic-10000", StringFamilies.Periodic(10000, 61, 53)),
            ("a^3000", StringFamilies.AllSame(3000)),
            ("runs-10000", StringFamilies.Runs(10000, 40, 54)),
        })
        {
            EsaIndex esa = EsaIndex.Build(text);
            UkkonenTree ukkonen = UkkonenTree.Build(text);
            string window = text.Substring(rng.NextInt(text.Length / 2), Math.Min(text.Length / 3, 3000));
            string query = MatchingChecks.Mutate(window, 0.02, rng);
            foreach (int minLength in new[] { 8, 12, 20, 60 })
            {
                var mems = MatchingChecks.Triples(esa.FindMaximalExactMatches(query, minLength));
                Assert.That(mems, Is.EqualTo(MatchingChecks.Triples(ukkonen.FindMaximalExactMatches(query, minLength))), $"{name} MEM L {minLength}");
                foreach (MumUniqueness mode in new[] { MumUniqueness.Both, MumUniqueness.Reference })
                    Assert.That(MatchingChecks.Triples(esa.FindMaximalUniqueMatches(query, minLength, mode)), Is.EqualTo(MatchingChecks.Triples(ukkonen.FindMaximalUniqueMatches(query, minLength, mode))), $"{name} MUM {mode} L {minLength}");
                var anchors = esa.FindExactMatchAnchors(query, minLength);
                Assert.That(anchors.Select(a => (a.PositionInQuery, a.Length)).ToList(), Is.EqualTo(ukkonen.FindExactMatchAnchors(query, minLength).Select(a => (a.PositionInQuery, a.Length)).ToList()), $"{name} anchors L {minLength}");
                foreach (Triple a in anchors) Assert.That(string.CompareOrdinal(text, a.PositionInText, query, a.PositionInQuery, a.Length), Is.EqualTo(0), $"{name}: anchor {a} is an occurrence");
            }
            var lcs = esa.FindAllDistinctLongestCommonSubstrings(query);
            var uLcs = ukkonen.FindAllDistinctLongestCommonSubstrings(query);
            Assert.That(lcs.Count, Is.EqualTo(uLcs.Count), $"{name}: LCS ties");
            for (int i = 0; i < lcs.Count; i++)
            {
                Assert.That(lcs[i].Substring, Is.EqualTo(uLcs[i].Substring));
                Assert.That(lcs[i].PositionsInText, Is.EqualTo(uLcs[i].PositionsInText));
                Assert.That(lcs[i].PositionsInOther, Is.EqualTo(uLcs[i].PositionsInOther));
            }
            var info = esa.LongestCommonSubstringInfo(query);
            Assert.That(info.Substring, Is.EqualTo(ukkonen.LongestCommonSubstringInfo(query).Substring), $"{name}: LCS");
            Assert.That(string.CompareOrdinal(text, info.PositionInText, info.Substring, 0, info.Substring.Length), Is.EqualTo(0), $"{name}: LCS text position is an occurrence");
        }
    }

    // ---- contract ----

    [Test]
    public void Contract_ArgumentChecks()
    {
        EsaIndex esa = EsaIndex.Build("banana");
        Assert.Throws<ArgumentNullException>(() => esa.FindMaximalExactMatches(null!, 3));
        Assert.Throws<ArgumentNullException>(() => esa.FindMaximalUniqueMatches(null!, 3));
        Assert.Throws<ArgumentNullException>(() => esa.FindExactMatchAnchors(null!, 3));
        Assert.Throws<ArgumentNullException>(() => esa.LongestCommonSubstring((string)null!));
        Assert.Throws<ArgumentNullException>(() => esa.LongestCommonSubstringInfo((string)null!));
        Assert.Throws<ArgumentNullException>(() => esa.FindAllLongestCommonSubstrings(null!));
        Assert.Throws<ArgumentNullException>(() => esa.FindAllDistinctLongestCommonSubstrings(null!));
        foreach (int bad in new[] { 0, -1, int.MinValue })
        {
            Assert.Throws<ArgumentOutOfRangeException>(() => esa.FindMaximalExactMatches("ana", bad), $"MEM minLength {bad}");
            Assert.Throws<ArgumentOutOfRangeException>(() => esa.FindMaximalUniqueMatches("ana", bad), $"MUM minLength {bad}");
            Assert.That(esa.FindExactMatchAnchors("ana", bad), Is.Empty, $"anchors minLength {bad}: empty, as the Ukkonen tree");
        }
        Assert.Throws<ArgumentOutOfRangeException>(() => esa.FindMaximalUniqueMatches("ana", 2, (MumUniqueness)99));
        Assert.That(esa.FindMaximalExactMatches("an", 3), Is.Empty, "a query shorter than minLength");
        Assert.That(EsaIndex.Build("ab").FindMaximalExactMatches("banana", 3), Is.Empty, "a text shorter than minLength");
    }

    [Test]
    public void Life_TheMatchingOperationsBuildThePackOnFirstUse()
    {
        EsaIndex index = EsaIndex.Build(StringFamilies.RandomDna(4000, 9));
        Assert.That(index.SuffixLinkSupport, Is.EqualTo(IndexSupport.Degraded));
        string query = index.Text.Substring(10, 100);
        _ = index.FindMaximalExactMatches(query, 12);
        Assert.That(index.SuffixLinkSupport, Is.EqualTo(IndexSupport.Native));
        Assert.That(index.NavigationBuildCount, Is.EqualTo(1));
        _ = index.FindMaximalUniqueMatches(query, 12);
        _ = index.FindExactMatchAnchors(query, 12);
        _ = index.LongestCommonSubstring(query);
        Assert.That(index.NavigationBuildCount, Is.EqualTo(1), "one build serves every operation");
    }

    [Test]
    public void Concurrency_ParallelQueriesOnOneIndexAgree()
    {
        string text = StringFamilies.RandomDna(30_000, 10);
        EsaIndex index = EsaIndex.Build(text);
        var rng = new SplitMix64(77);
        string[] queries = Enumerable.Range(0, 24).Select(_ => MatchingChecks.Mutate(text.Substring(rng.NextInt(20_000), 2000), 0.02, rng)).ToArray();
        UkkonenTree ukkonen = UkkonenTree.Build(text);
        var expected = queries.Select(q => MatchingChecks.Triples(ukkonen.FindMaximalExactMatches(q, 14))).ToArray();
        var actual = new List<Triple>[queries.Length];
        Parallel.For(0, queries.Length, new ParallelOptions { MaxDegreeOfParallelism = 8 }, i => actual[i] = MatchingChecks.Triples(index.FindMaximalExactMatches(queries[i], 14)));
        for (int i = 0; i < queries.Length; i++) Assert.That(actual[i], Is.EqualTo(expected[i]), $"query {i}");
        Assert.That(index.NavigationBuildCount, Is.EqualTo(1));
    }

    [Test]
    public void ChildLayouts_GiveTheSameAnswers()
    {
        string text = StringFamilies.Periodic(3000, 29, 12);
        string query = MatchingChecks.Mutate(text.Substring(100, 1500), 0.03, new SplitMix64(5));
        var expected = MatchingChecks.Triples(EsaIndex.Build(text).FindMaximalExactMatches(query, 10));
        foreach (ChildLayoutKind layout in Enum.GetValues<ChildLayoutKind>())
        {
            EsaIndex index = EsaIndex.Build(text, layout);
            Assert.That(MatchingChecks.Triples(index.FindMaximalExactMatches(query, 10)), Is.EqualTo(expected), $"{layout}");
            string small = text[..300];
            MatchingChecks.CheckAll(EsaIndex.Build(small, layout), null, small, MatchingChecks.Mutate(text.Substring(20, 200), 0.05, new SplitMix64(6)), [2, 6], $"{layout}");
        }
    }
}

/// <summary>Longer enumerations of <see cref="EsaMatchingTests"/>: all binary text/query pairs up to length 8 and 9, ternary up to 5/6, minLength 1-4. Heavy; explicit, category <c>Exhaustive</c>.</summary>
[TestFixture]
[Category("Exhaustive")]
[Explicit("Heavy enumeration; run with --filter Category=Exhaustive (see PROGRESS-S1.md, Log 1b.4)")]
public class EsaMatchingExhaustiveTests
{
    private static void Pairs(int sigma, int textLength, int queryLength, int querySigma)
    {
        string[] queries = ExhaustiveStrings.All(querySigma, queryLength).ToArray();
        foreach (string text in ExhaustiveStrings.All(sigma, textLength))
        {
            EsaIndex esa = EsaIndex.Build(text);
            UkkonenTree ukkonen = UkkonenTree.Build(text);
            foreach (string query in queries) MatchingChecks.CheckAll(esa, ukkonen, text, query, [1, 2, 3, 4], $"sigma {sigma}");
        }
    }

    [Test]
    public void BinaryPairs_TextsUpTo9_QueriesUpTo8() => Pairs(2, 9, 8, 2);

    [Test]
    public void TernaryPairs_TextsUpTo6_QueriesUpTo5() => Pairs(3, 6, 5, 3);

    [Test]
    public void FourLetterPairs_TextsUpTo4_QueriesUpTo4() => Pairs(4, 4, 4, 4);

    [Test]
    public void BinaryTexts_ForeignQueries_UpTo7() => Pairs(2, 7, 6, 3);
}
