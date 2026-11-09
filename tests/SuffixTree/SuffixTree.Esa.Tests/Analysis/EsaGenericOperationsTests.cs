using SuffixTree.Esa.Tests.TestSupport;
using UkkonenTree = global::SuffixTree.SuffixTree;

namespace SuffixTree.Esa.Tests.Analysis;

/// <summary>
/// The operations that the ESA serves through <c>Traverse</c> and the generic algorithms of <c>SuffixTree.Core</c>, without a specialised implementation (work package 1b.3):
/// <c>FindMaximalRepeatedPairs</c> and the k-common-substring operations (<c>FindLongestCommonSubstrings</c>, <c>LongestCommonSubstringLengthsBySupport</c>). The tests
/// wire them through the test view <see cref="EsaNativeTreeShim"/> (a call through <see cref="ISuffixTree"/> reaches the interface default, which runs on <c>Traverse</c> of the ESA)
/// and check them against the definition (L0), against the Ukkonen tree (L2) and against each other.
/// </summary>
[TestFixture]
public class EsaGenericOperationsTests
{
    private static readonly (string Name, Func<char, bool>? Predicate)[] Predicates =
    [
        ("none", null),
        ("N", c => c == 'N'),
        ("not-ACGT", c => "ACGT".IndexOf(c) < 0),
        ("a", c => c == 'a'),
        ("a-or-NUL", c => c is 'a' or '\0'),
    ];

    /// <summary>The ESA behind the interface, as the facade of work package 1b.6 will be (the default members run on <c>Traverse</c>).</summary>
    private static ISuffixTree Tree(string text) => new EsaNativeTreeShim(EsaIndex.Build(text));

    private static List<(int, int, int)> Pairs(IReadOnlyList<(int FirstPosition, int SecondPosition, int Length)> list) => list.Select(p => (p.FirstPosition, p.SecondPosition, p.Length)).ToList();

    private static void CheckPairs(string text, int minLength, (string Name, Func<char, bool>? Predicate) unique, string context)
    {
        ISuffixTree viaInterface = Tree(text);
        var actual = Pairs(viaInterface.FindMaximalRepeatedPairs(minLength, unique.Predicate));
        var expected = AnalysisReference.MaximalRepeatedPairs(text, minLength, unique.Predicate);
        Assert.That(actual, Is.EqualTo(expected), $"{context} L {minLength} unique {unique.Name}: against the definition");
        Assert.That(actual, Is.EqualTo(Pairs(UkkonenTree.Build(text).FindMaximalRepeatedPairs(minLength, unique.Predicate))), $"{context} L {minLength} unique {unique.Name}: against Ukkonen");
        Assert.That(actual, Is.EqualTo(Pairs(SuffixTreeAlgorithms.FindMaximalRepeatedPairs(new EsaTreeShim(EsaIndex.Build(text)), minLength, unique.Predicate))), $"{context}: the Traverse-default shim takes the same path");
    }

    [Test]
    public void MaximalRepeatedPairs_LiteratureValues()
    {
        // repeat-match -f of MUMmer 3.23, converted to 0-based (MaximalRepeatedPairsTests)
        ISuffixTree gusfield = Tree("xabcyabcwabcyz");
        Assert.That(Pairs(gusfield.FindMaximalRepeatedPairs(1)), Is.EqualTo(new[] { (1, 5, 3), (1, 9, 4), (5, 9, 3) }));
        Assert.That(Pairs(gusfield.FindMaximalRepeatedPairs(4)), Is.EqualTo(new[] { (1, 9, 4) }));
        Assert.That(Pairs(Tree("acgtacgtacgt").FindMaximalRepeatedPairs(3)), Is.EqualTo(new[] { (0, 4, 8), (0, 8, 4) }));
        Assert.That(Pairs(Tree("aaaaa").FindMaximalRepeatedPairs(1)), Is.EqualTo(new[] { (0, 1, 4), (0, 2, 3), (0, 3, 2), (0, 4, 1) }));
        Assert.That(Pairs(Tree("gattacagattacca").FindMaximalRepeatedPairs(2)), Is.EqualTo(new[] { (0, 7, 6), (5, 13, 2) }));
    }

    [Test]
    public void MaximalRepeatedPairs_EdgeCases()
    {
        foreach (string text in new[] { "", "a", "ab", "aa", "aaa", "abc", "\0\0", "￿￿" })
            foreach (int minLength in new[] { 1, 2, 3, 10 })
                CheckPairs(text, minLength, Predicates[0], $"'{text}'");
        foreach (int bad in new[] { 0, -1, int.MinValue })
            Assert.Throws<ArgumentOutOfRangeException>(() => Tree("banana").FindMaximalRepeatedPairs(bad), $"minLength {bad}");
    }

    [Test]
    public void MaximalRepeatedPairs_ExhaustiveSmallStrings_AllThresholdsAndUniqueSymbols()
    {
        foreach (var (sigma, maxLen) in new[] { (2, 9), (3, 6), (4, 5) })
            foreach (string text in ExhaustiveStrings.All(sigma, maxLen).Where(s => s.Length > 0))
                foreach (int minLength in new[] { 1, 2, 3 })
                    CheckPairs(text, minLength, Predicates[text.Length % 2 == 0 ? 0 : 3], $"'{text}'");
    }

    [Test]
    public void MaximalRepeatedPairs_AdversarialFamilies()
    {
        foreach (var (name, text) in AnalysisCorpus.Adversarial([8, 21, 55, 120], 0x1B3D))
            foreach (int minLength in new[] { 1, 3, 8 })
                CheckPairs(text, minLength, Predicates[(text.Length + minLength) % Predicates.Length], name);
    }

    [Test]
    public void MaximalRepeatedPairs_LargerTexts_AgainstUkkonen()
    {
        var rng = new SplitMix64(0x1B3E);
        for (int i = 0; i < 6; i++)
        {
            int n = 1000 + rng.NextInt(4001);
            ulong seed = rng.Next();
            foreach (var (name, text) in StringFamilies.Adversarial(n, seed).Where((_, k) => (k + i) % 7 == 0))
            {
                // keep the output bounded: the random-like families have few long pairs, the structured ones are cut at a quarter of the text
                bool randomLike = new[] { "random-", "bmp", "skewed", "all-distinct", "surrogate", "lone-", "nul-" }.Any(prefix => name.StartsWith(prefix, StringComparison.Ordinal));
                int minLength = randomLike ? 24 : n / 4;
                ISuffixTree viaInterface = Tree(text);
                var actual = Pairs(viaInterface.FindMaximalRepeatedPairs(minLength));
                var expected = Pairs(UkkonenTree.Build(text).FindMaximalRepeatedPairs(minLength));
                Assert.That(actual, Is.EqualTo(expected), $"{name} n {n} L {minLength}");
            }
        }
    }

    private static IEnumerable<string[]> TextSets()
    {
        yield return ["GATTACA", "TAGACCA", "ATACA"];                           // Rosalind LCSM
        yield return ["banana"];
        yield return ["banana", "ananas"];
        yield return ["AB", "CDE", "XYZ", "AB"];
        yield return ["", ""];
        yield return ["ACGT", "", "ACGT"];
        yield return ["AAAAAAA", "CAAC", "GGG"];
        yield return ["abab", "baba", "abba"];
        yield return ["ab", "ab", "b"];                      // texts that use the private-use separators
        yield return ["aaaa", "aaaa", "aaaa", "aaaa"];
        var rng = new SplitMix64(0x1B3F);
        for (int i = 0; i < 120; i++)
        {
            int k = 1 + rng.NextInt(5);
            string letters = (i % 3) switch { 0 => "ab", 1 => "acg", _ => "acgtn" };
            yield return Enumerable.Range(0, k).Select(_ => new string(Enumerable.Range(0, rng.NextInt(25)).Select(_ => letters[rng.NextInt(letters.Length)]).ToArray())).ToArray();
        }
    }

    [Test]
    public void KCommonSubstrings_FromTheEsa_AgreeWithTheDefinitionAndUkkonen()
    {
        static ISuffixTree BuildEsa(string s) => Tree(s);
        foreach (string[] texts in TextSets())
        {
            string context = "{" + string.Join("|", texts) + "}";
            var (lengths, bySupport) = AnalysisReference.CommonSubstrings(texts);
            Assert.That(SuffixTreeAlgorithms.LongestCommonSubstringLengthsBySupport(texts, BuildEsa), Is.EqualTo(lengths), $"{context}: lengths against the definition");
            Assert.That(SuffixTreeAlgorithms.LongestCommonSubstringLengthsBySupport(texts, BuildEsa), Is.EqualTo(UkkonenTree.LongestCommonSubstringLengthsBySupport(texts)), $"{context}: lengths against Ukkonen");
            for (int q = 1; q <= texts.Length; q++)
            {
                IReadOnlyList<string> actual = SuffixTreeAlgorithms.FindLongestCommonSubstrings(texts, q, BuildEsa);
                Assert.That(actual, Is.EqualTo(bySupport[q]), $"{context} q {q}: against the definition");
                Assert.That(actual, Is.EqualTo(UkkonenTree.FindLongestCommonSubstrings(texts, q)), $"{context} q {q}: against Ukkonen");
            }
        }
    }

    [Test]
    public void KCommonSubstrings_Validation_IsTheGenericOnes()
    {
        static ISuffixTree BuildEsa(string s) => Tree(s);
        Assert.Throws<ArgumentException>(() => SuffixTreeAlgorithms.FindLongestCommonSubstrings(Array.Empty<string>(), 1, BuildEsa));
        Assert.Throws<ArgumentOutOfRangeException>(() => SuffixTreeAlgorithms.FindLongestCommonSubstrings(["a", "b"], 3, BuildEsa));
        Assert.Throws<ArgumentOutOfRangeException>(() => SuffixTreeAlgorithms.FindLongestCommonSubstrings(["a", "b"], 0, BuildEsa));
    }

    [Test]
    public void NativeAndTraverseDefaultViews_AreTheSameTreeForTheGenericAlgorithms()
    {
        // the two test views differ only in who answers the five default-implemented members; the event stream they hand to the generic algorithms is the same
        foreach (string text in new[] { "banana", "mississippi", "abracadabra", StringFamilies.Fibonacci(80) })
        {
            EsaIndex esa = EsaIndex.Build(text);
            var native = new TraverseRecorder();
            var generic = new TraverseRecorder();
            new EsaNativeTreeShim(esa).Traverse(native);
            new EsaTreeShim(esa).Traverse(generic);
            Assert.That(native.Events, Is.EqualTo(generic.Events), text);
        }
    }
}
