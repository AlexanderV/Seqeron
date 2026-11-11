using SuffixTree.Esa.Tests.Fuzz;
using SuffixTree.Esa.Tests.TestSupport;
using SuffixTree.Persistent;
using UkkonenTree = global::SuffixTree.SuffixTree;

namespace SuffixTree.Esa.Tests.Facade;

/// <summary>
/// The facade of every child-table layout against the definitions, the Ukkonen tree and the persistent v6 tree, which are reached through the same backend seam
/// (<see cref="SuffixTreeBackends"/>), on named texts, the adversarial families and seeded random strings. Every member of <see cref="ISuffixTree"/> is compared
/// (<see cref="FacadeChecks.CheckEverything"/>); the results the contract leaves open are checked by property.
/// </summary>
[TestFixtureSource(typeof(SuffixTreeBackends), nameof(SuffixTreeBackends.EsaFixtures))]
public class FacadeDifferentialTests(SuffixTreeBackend backend) : PortedBackendTestBase(backend)
{
    private static readonly string[] Named =
    [
        "", "a", "aa", "ab", "ba", "aba", "abab", "banana", "mississippi", "abracadabra", "aaaaaaaaaa", "abababababab", "abcabxabcd", "xabxac", "GATTACAGATTACA", "ACGTACGTTACGT",
        "\uD800a\uD800", "\0￿\0", "a$b$", "🧬αβγ🧪$", "ab😀ab😀", "the quick brown fox jumps over the lazy dog", "aabaabaaabaabaaab",
    ];

    private void Compare(string text, ulong seed, string context, int maxPatternLength = 6)
    {
        ISuffixTree candidate = Tree(text);
        using var v6 = PersistentSuffixTreeFactory.CreatePersistent(new StringTextSource(text));
        FacadeChecks.CheckEverything(candidate, UkkonenTree.Build(text), v6, text, seed, $"[{Backend.Name}] {context}", maxPatternLength);
    }

    [Test]
    public void NamedTexts_EveryMember()
    {
        ulong seed = 1;
        foreach (string text in Named) Compare(text, seed++, "named");
    }

    [Test]
    public void AdversarialFamilies_EveryMember_N150()
    {
        foreach (var (name, text) in StringFamilies.Adversarial(150, 3))
            Compare(text, 7, name, 4);
    }

    [Test]
    public void AdversarialFamilies_SmallSizes_EveryMember()
    {
        foreach (int n in new[] { 1, 2, 3, 4, 5, 8, 13 })
            foreach (var (name, text) in StringFamilies.Adversarial(n, 11))
                Compare(text, 5, $"{name} n={n}", 4);
    }

    [Test]
    public void RandomStrings_EveryMember()
    {
        var rng = new SplitMix64(2024);
        for (int i = 0; i < 40; i++)
        {
            int n = rng.NextInt(120);
            int sigma = new[] { 1, 2, 3, 4, 4, 20 }[rng.NextInt(6)];
            Compare(StringFamilies.Random(n, sigma, rng.Next()), rng.Next(), $"random {i} sigma {sigma} n {n}", 5);
        }
    }

    [Test]
    public void LargerTexts_SearchAndDiagnostics_AgainstUkkonen()
    {
        foreach (var (name, text) in new (string, string)[]
        {
            ("dna-5000", StringFamilies.RandomDna(5000, 17)),
            ("fibonacci-4000", StringFamilies.Fibonacci(4000)),
            ("thue-morse-4096", StringFamilies.ThueMorse(4096)),
            ("a^3000", new string('a', 3000)),
            ("periodic-5", StringFamilies.Periodic(3000, 5, 3)),
            ("protein-3000", StringFamilies.RandomProtein(3000, 9)),
        })
        {
            ISuffixTree candidate = Tree(text);
            UkkonenTree ukkonen = UkkonenTree.Build(text);
            FacadeChecks.CheckAgainstTree(candidate, ukkonen, "ukkonen", text, $"[{Backend.Name}] {name}");
            var rng = new SplitMix64(41);
            for (int k = 0; k < 40; k++)
            {
                int start = rng.NextInt(text.Length);
                string pattern = text.Substring(start, Math.Min(1 + rng.NextInt(40), text.Length - start));
                Assert.That(candidate.CountOccurrences(pattern), Is.EqualTo(ukkonen.CountOccurrences(pattern)), $"{name} count '{pattern}'");
                Assert.That(candidate.FindAllOccurrences(pattern).Order(), Is.EqualTo(ukkonen.FindAllOccurrences(pattern).Order()), $"{name} positions '{pattern}'");
            }
            Assert.That(candidate.CountDistinctSubstrings(), Is.EqualTo(ukkonen.CountDistinctSubstrings()), name);
            Assert.That(candidate.CountDistinctSubstringsByLength(50), Is.EqualTo(ukkonen.CountDistinctSubstringsByLength(50)), name);
            Assert.That(candidate.LongestRepeatedSubstring().Length, Is.EqualTo(ukkonen.LongestRepeatedSubstring().Length), name);
            var allCandidate = candidate.FindAllLongestRepeatedSubstrings();
            var allUkkonen = ukkonen.FindAllLongestRepeatedSubstrings();
            Assert.That(allCandidate.Select(t => t.Substring), Is.EqualTo(allUkkonen.Select(t => t.Substring)), name);
            Assert.That(allCandidate.Select(t => t.Positions), Is.EqualTo(allUkkonen.Select(t => t.Positions)), name);
            string query = string.Concat(StringFamilies.RandomDna(300, 5), text.AsSpan(text.Length / 3, Math.Min(200, text.Length - (text.Length / 3))));
            Assert.That(candidate.FindExactMatchAnchors(query, 6).Select(a => (a.PositionInQuery, a.Length)), Is.EqualTo(ukkonen.FindExactMatchAnchors(query, 6).Select(a => (a.PositionInQuery, a.Length))), $"{name} anchors");
            Assert.That(candidate.FindMaximalExactMatches(query, 8), Is.EqualTo(ukkonen.FindMaximalExactMatches(query, 8)), $"{name} MEM");
            Assert.That(candidate.FindMaximalUniqueMatches(query, 8), Is.EqualTo(ukkonen.FindMaximalUniqueMatches(query, 8)), $"{name} MUM");
            Assert.That(candidate.FindAllDistinctLongestCommonSubstrings(query).Select(t => (t.Substring, t.PositionsInText.ToArray(), t.PositionsInOther.ToArray())),
                Is.EqualTo(ukkonen.FindAllDistinctLongestCommonSubstrings(query).Select(t => (t.Substring, t.PositionsInText.ToArray(), t.PositionsInOther.ToArray()))), $"{name} LCS");
        }
    }

    [Test]
    public void LogicalHash_EqualsThatOfTheOtherBackends()
    {
        foreach (string text in Named.Concat(StringFamilies.Adversarial(90, 2).Select(f => f.Text)))
        {
            ISuffixTree candidate = Tree(text);
            using var v6 = PersistentSuffixTreeFactory.CreatePersistent(new StringTextSource(text));
            byte[] expected = SuffixTreeSerializer.CalculateLogicalHash(UkkonenTree.Build(text));
            Assert.That(SuffixTreeSerializer.CalculateLogicalHash(candidate), Is.EqualTo(expected), $"[{Backend.Name}] '{FacadeChecks.Short(text)}' against Ukkonen");
            Assert.That(SuffixTreeSerializer.CalculateLogicalHash(v6), Is.EqualTo(expected), $"v6 '{FacadeChecks.Short(text)}'");
        }
    }

    [Test]
    public void TheFacade_IsIndependentOfTheChildTableLayout()
    {
        // the three layouts are three backends; the same text must give byte-identical diagnostics
        foreach (string text in Named.Concat(StringFamilies.Adversarial(70, 8).Select(f => f.Text)))
        {
            string reference = SuffixTreeBackends.Esa.FromString(text).PrintTree();
            Assert.That(Tree(text).PrintTree(), Is.EqualTo(reference), $"[{Backend.Name}] PrintTree of '{FacadeChecks.Short(text)}'");
            var a = new FacadeChecks.EventCollector();
            var b = new FacadeChecks.EventCollector();
            SuffixTreeBackends.Esa.FromString(text).Traverse(a);
            Tree(text).Traverse(b);
            Assert.That(b.Events, Is.EqualTo(a.Events), $"[{Backend.Name}] Traverse of '{FacadeChecks.Short(text)}'");
        }
    }

    [Test]
    public void TheTextSourceOverload_BuildsTheSameTree()
    {
        const string text = "abracadabra abracadabra";
        ISuffixTree fromString = Tree(text);
        ISuffixTree fromStringSource = Tree(new StringTextSource(text));
        ISuffixTree fromCustomSource = Tree(new ArraySource(text));
        foreach (ISuffixTree tree in new[] { fromStringSource, fromCustomSource })
        {
            Assert.That(tree.Text.ToString(), Is.EqualTo(text));
            Assert.That(tree.NodeCount, Is.EqualTo(fromString.NodeCount));
            Assert.That(FacadeChecks.Events(tree), Is.EqualTo(FacadeChecks.Events(fromString)));
            Assert.That(tree.FindAllOccurrences("abra").Order(), Is.EqualTo(fromString.FindAllOccurrences("abra").Order()));
        }
    }

    /// <summary>A text source that is not a <see cref="StringTextSource"/>.</summary>
    private sealed class ArraySource(string text) : ITextSource
    {
        private readonly char[] _chars = text.ToCharArray();

        public int Length => _chars.Length;

        public char this[int index] => _chars[index];

        public string Substring(int start, int length) => new(_chars, start, length);

        public ReadOnlySpan<char> Slice(int start, int length) => _chars.AsSpan(start, length);

        public IEnumerator<char> GetEnumerator() => ((IEnumerable<char>)_chars).GetEnumerator();

        System.Collections.IEnumerator System.Collections.IEnumerable.GetEnumerator() => GetEnumerator();

        public override string ToString() => new(_chars);
    }
}
