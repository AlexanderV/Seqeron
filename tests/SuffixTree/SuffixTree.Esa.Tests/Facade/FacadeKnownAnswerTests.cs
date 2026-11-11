namespace SuffixTree.Esa.Tests.Facade;

/// <summary>
/// The literal answers that the Mcp.Core tests carry for the library operations (hand-derived or produced by MUMmer 3.23 / Gusfield / Rosalind, see the originals), re-asserted
/// through the backend factory for every backend. The Mcp.Core tests themselves check the argument validation and the DTO mapping of the wrapper, which is not backend agnostic.
/// </summary>
[TestFixtureSource(typeof(SuffixTreeBackends), nameof(SuffixTreeBackends.AllFixtures))]
public class FacadeKnownAnswerTests(SuffixTreeBackend backend) : PortedBackendTestBase(backend)
{
    [TestCase("banana", "ana", true, 2)]
    [TestCase("banana", "xyz", false, 0)]
    [TestCase("banana", "", true, 6)]
    [TestCase("aaaaa", "aa", true, 4)]
    [TestCase("abc", "", true, 3)]
    [TestCase("AbCd", "abcd", false, 0)]
    [TestCase("AbCd", "AbC", true, 1)]
    [TestCase("banana", "b", true, 1)]
    public void ContainsCountFindAll_LiteralAnswers(string text, string pattern, bool contains, int count)
    {
        ISuffixTree tree = Tree(text);
        Assert.Multiple(() =>
        {
            Assert.That(tree.Contains(pattern), Is.EqualTo(contains));
            Assert.That(tree.CountOccurrences(pattern), Is.EqualTo(count));
            Assert.That(tree.FindAllOccurrences(pattern), Has.Count.EqualTo(count));
        });
    }

    [Test]
    public void ContainsCountFindAll_LiteralAnswers_Positions()
    {
        Assert.Multiple(() =>
        {
            Assert.That(Tree("banana").FindAllOccurrences("ana").Order(), Is.EqualTo(new[] { 1, 3 }));
            Assert.That(Tree("aaaaa").FindAllOccurrences("aa").Order(), Is.EqualTo(new[] { 0, 1, 2, 3 }));
            Assert.That(Tree("abc").FindAllOccurrences("").Order(), Is.EqualTo(new[] { 0, 1, 2 }));
            Assert.That(Tree("abc").FindAllOccurrences("xyz"), Is.Empty);
            Assert.That(Tree("GATATATGCATATACTT").FindAllOccurrences("ATAT").Order(), Is.EqualTo(new[] { 1, 3, 9 }));
            Assert.That(Tree("mississippi").FindAllOccurrences("i").Order(), Is.EqualTo(new[] { 1, 4, 7, 10 }));
            Assert.That(Tree("abracadabra").FindAllOccurrences("a").Order(), Is.EqualTo(new[] { 0, 3, 5, 7, 10 }));
        });
    }

    [Test]
    public void LongestCommonSubstring_LiteralAnswers()
    {
        Assert.Multiple(() =>
        {
            Assert.That(Tree("xabxac").LongestCommonSubstring("abcabxabcd"), Is.EqualTo("abxa"));
            Assert.That(Tree("ab").LongestCommonSubstring("bandana"), Is.EqualTo("b"));
            Assert.That(Tree("banana").LongestCommonSubstring("panama"), Is.EqualTo("ana"));
            Assert.That(Tree("abcdef").LongestCommonSubstring("xyzcdeww"), Is.EqualTo("cde"));
            Assert.That(Tree("abc").LongestCommonSubstring("xyz"), Is.EqualTo(""));
        });
        const string text1 = "abcabxabcd";
        const string text2 = "zzabcabczz";
        string lcs = Tree(text1).LongestCommonSubstring(text2);
        Assert.Multiple(() =>
        {
            Assert.That(text1.Contains(lcs, StringComparison.Ordinal), Is.True);
            Assert.That(text2.Contains(lcs, StringComparison.Ordinal), Is.True);
            Assert.That(lcs.Length, Is.EqualTo(LongestCommonSubstringLength(text1, text2)));
        });
    }

    private static int LongestCommonSubstringLength(string a, string b)
    {
        var dp = new int[b.Length + 1];
        int best = 0;
        for (int i = 1; i <= a.Length; i++)
            for (int j = b.Length; j >= 1; j--)
            {
                dp[j] = a[i - 1] == b[j - 1] ? dp[j - 1] + 1 : 0;
                best = Math.Max(best, dp[j]);
            }
        return best;
    }

    [TestCase("mississippi", "issi")]
    [TestCase("abracadabra", "abra")]
    [TestCase("ATATAT", "ATAT")]
    [TestCase("banana", "ana")]
    [TestCase("aaaaa", "aaaa")]
    [TestCase("abcdef", "")]
    public void LongestRepeatedSubstring_LiteralAnswers(string text, string expected)
    {
        ISuffixTree tree = Tree(text);
        Assert.That(tree.LongestRepeatedSubstring(), Is.EqualTo(expected));
        Assert.That(tree.LongestRepeatedSubstringMemory().ToString(), Is.EqualTo(expected));
    }

    [Test]
    public void AllLongestRepeatedSubstrings_LiteralAnswers()
    {
        var ties = Tree("abcxbcaxcab").FindAllLongestRepeatedSubstrings();
        Assert.That(ties.Select(t => t.Substring), Is.EqualTo(new[] { "ab", "bc", "ca" }));
        Assert.That(ties.Select(t => t.Positions.ToArray()), Is.EqualTo(new[] { new[] { 0, 9 }, new[] { 1, 4 }, new[] { 5, 8 } }));
        foreach (var (text, substring, positions) in new[] { ("banana", "ana", new[] { 1, 3 }), ("mississippi", "issi", new[] { 1, 4 }), ("aaaa", "aaa", new[] { 0, 1 }) })
        {
            var single = Tree(text).FindAllLongestRepeatedSubstrings();
            Assert.That(single, Has.Count.EqualTo(1), text);
            Assert.That(single[0].Substring, Is.EqualTo(substring), text);
            Assert.That(single[0].Positions.ToArray(), Is.EqualTo(positions), text);
        }
        Assert.That(Tree("abcd").FindAllLongestRepeatedSubstrings(), Is.Empty);
        const string xyz = "xyzzyxxz";
        var expected = global::SuffixTree.SuffixTree.Build(xyz).FindAllLongestRepeatedSubstrings();
        var actual = Tree(xyz).FindAllLongestRepeatedSubstrings();
        Assert.That(actual.Select(l => l.Substring), Is.EqualTo(expected.Select(l => l.Substring)));
        Assert.That(actual.Select(l => l.Positions.ToArray()), Is.EqualTo(expected.Select(l => l.Positions.ToArray())));
    }

    [Test]
    public void AllDistinctLongestCommonSubstrings_LiteralAnswers()
    {
        var ties = Tree("abcxyzabc").FindAllDistinctLongestCommonSubstrings("xyzqabcpxyz");
        Assert.That(ties.Select(s => s.Substring), Is.EqualTo(new[] { "xyz", "abc" }));
        Assert.That(ties.Select(s => s.PositionsInText.ToArray()), Is.EqualTo(new[] { new[] { 3 }, new[] { 0, 6 } }));
        Assert.That(ties.Select(s => s.PositionsInOther.ToArray()), Is.EqualTo(new[] { new[] { 0, 8 }, new[] { 4 } }));
        Assert.That(Tree("AAAA").FindAllDistinctLongestCommonSubstrings("CCC"), Is.Empty);
        var gattaca = Tree("GATTACA").FindAllDistinctLongestCommonSubstrings("TAGACCA");
        Assert.That(gattaca.Select(s => s.Substring), Is.EqualTo(new[] { "TA", "GA", "AC", "CA" }));
    }

    [Test]
    public void MaximalMatches_LiteralAnswers()
    {
        const string reference = "TTTCCTCATGCAATTCAAAACCATGTCCGTAATGTAGGCG";
        const string query = "TTTCCTCATGCAAATAGCCTCATGCATCCGTAATGTAG";
        ISuffixTree gattaca = Tree("GATTACAGATTACA");
        ISuffixTree three = Tree(reference);
        Assert.Multiple(() =>
        {
            Assert.That(gattaca.FindMaximalExactMatches("ATTACAGGATTACAT", 3), Is.EqualTo(new[] { (1, 0, 7), (8, 0, 6), (0, 7, 7), (7, 7, 7) }));
            Assert.That(gattaca.FindMaximalExactMatches("ATTACAGGATTACAT", 20), Is.Empty);
            Assert.That(three.FindMaximalExactMatches(query, 4), Is.EqualTo(new[] { (0, 0, 13), (21, 6, 4), (15, 10, 4), (3, 17, 9), (21, 20, 4), (25, 26, 12), (22, 32, 4) }));
            Assert.That(three.FindMaximalUniqueMatches(query, 4, MumUniqueness.Both), Is.EqualTo(new[] { (0, 0, 13), (15, 10, 4), (25, 26, 12) }));
            Assert.That(three.FindMaximalUniqueMatches(query, 4), Is.EqualTo(new[] { (0, 0, 13), (15, 10, 4), (25, 26, 12) }), "the default is Both");
            Assert.That(three.FindMaximalUniqueMatches(query, 4, MumUniqueness.Reference), Is.EqualTo(new[] { (0, 0, 13), (15, 10, 4), (3, 17, 9), (25, 26, 12) }));
            Assert.That(Tree("ACGT").FindMaximalUniqueMatches("TTAG", 1, MumUniqueness.Both), Is.EqualTo(new[] { (0, 2, 1), (2, 3, 1) }));
        });
    }

    [Test]
    public void KCommonSubstrings_LiteralAnswers()
    {
        string[] lcsm = ["GATTACA", "TAGACCA", "ATACA"];
        Assert.Multiple(() =>
        {
            Assert.That(FindLongestCommonSubstrings(lcsm), Is.EqualTo(new[] { "AC", "CA", "TA" }));
            Assert.That(FindLongestCommonSubstrings(lcsm, 2), Is.EqualTo(new[] { "TACA" }));
            Assert.That(LongestCommonSubstringLengthsBySupport(lcsm), Is.EqualTo(new[] { 7, 7, 4, 2 }));
            Assert.That(FindLongestCommonSubstrings(["AAA", "CCC"]), Is.Empty);
        });
    }

    [Test]
    public void MaximalRepeatedPairs_LiteralAnswers()
    {
        Assert.Multiple(() =>
        {
            Assert.That(Tree("xabcyabcwabcyz").FindMaximalRepeatedPairs(3).Select(p => (p.FirstPosition, p.SecondPosition, p.Length)), Is.EqualTo(new[] { (1, 5, 3), (1, 9, 4), (5, 9, 3) }));
            ISuffixTree acgt = Tree("ACGTNACGTNACGT");
            Assert.That(acgt.FindMaximalRepeatedPairs(3).Select(p => (p.FirstPosition, p.SecondPosition, p.Length)), Is.EqualTo(new[] { (0, 5, 9), (0, 10, 4) }));
            Assert.That(acgt.FindMaximalRepeatedPairs(3, c => c == 'N').Select(p => (p.FirstPosition, p.SecondPosition, p.Length)), Is.EqualTo(new[] { (0, 5, 4), (0, 10, 4), (5, 10, 4) }));
        });
    }

    [TestCase("banana")]
    [TestCase("abcdef")]
    [TestCase("aaaaa")]
    [TestCase("mississippi")]
    public void Statistics_FundamentalInvariants(string text)
    {
        ISuffixTree tree = Tree(text);
        Assert.Multiple(() =>
        {
            Assert.That(tree.Text.Length, Is.EqualTo(text.Length));
            Assert.That(tree.LeafCount, Is.EqualTo(text.Length));
            Assert.That(tree.MaxDepth, Is.EqualTo(text.Length));
            Assert.That(tree.NodeCount, Is.GreaterThanOrEqualTo(text.Length + 1));
            Assert.That(tree.NodeCount, Is.LessThanOrEqualTo((2 * text.Length) + 1));
            Assert.That(tree.MaxDepth, Is.GreaterThanOrEqualTo(tree.LongestRepeatedSubstring().Length));
            Assert.That(tree.CountOccurrences(""), Is.EqualTo(tree.LeafCount));
        });
    }
}
