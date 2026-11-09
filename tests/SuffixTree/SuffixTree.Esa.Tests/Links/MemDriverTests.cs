using SuffixTree.Esa.Bench.Links;
using SuffixTree.Esa.Tests.TestSupport;
using UkkonenTree = global::SuffixTree.SuffixTree;

namespace SuffixTree.Esa.Tests.Links;

/// <summary>
/// The throwaway matching-statistics driver of prototype P2 (work package 1b.4a): its MEM set must equal the definition
/// (<c>SuffixTreeAlgorithms.FindMaximalExactMatchesByDefinition</c>) and the Ukkonen tree's <c>FindMaximalExactMatches</c> for every link kernel, both endpoint sources and
/// one or two loci, on exhaustive small texts and queries, adversarial families with mutated queries, foreign symbols, empty and over-long queries.
/// </summary>
[TestFixture]
public class MemDriverTests
{
    private sealed class Setup
    {
        internal required EsaIndex Esa { get; init; }

        internal required LinkContext Links { get; init; }
    }

    private static Setup Prepare(string text, int fanout = 4)
    {
        EsaIndex esa = EsaIndex.Build(text, ChildLayoutKind.Tagged);
        int[] isa = LinkTables.BuildIsa(esa.SuffixArray);
        return new Setup { Esa = esa, Links = new LinkContext(esa, isa, LinkTables.BuildPsi(esa.SuffixArray, isa), new LcpRangeHelper(esa.Lcp, fanout)) };
    }

    private static List<(int, int, int)> Sorted(IEnumerable<(int PositionInText, int PositionInQuery, int Length)> mems)
        => mems.Select(m => (m.PositionInText, m.PositionInQuery, m.Length)).OrderBy(m => m.Item2).ThenBy(m => m.Item1).ThenBy(m => m.Item3).ToList();

    private static void Check(Setup s, UkkonenTree? ukkonen, string text, string query, int minLength, LinkKernel kernel, EndpointSource source, int loci, string label)
    {
        var expected = Sorted(SuffixTreeAlgorithms.FindMaximalExactMatchesByDefinition(new StringTextSource(text), query, minLength));
        var driver = new MemDriver(s.Esa, s.Links, kernel, source, loci);
        var actual = Sorted(driver.Run(query, minLength));
        if (!expected.SequenceEqual(actual))
            Assert.Fail($"{label} text '{Clip(text)}' query '{Clip(query)}' minLength {minLength} {kernel}/{source}/loci {loci}: {actual.Count} matches, expected {expected.Count}");
        if (ukkonen is not null && !Sorted(ukkonen.FindMaximalExactMatches(query, minLength)).SequenceEqual(actual))
            Assert.Fail($"{label} text '{Clip(text)}' query '{Clip(query)}' minLength {minLength}: the driver differs from the Ukkonen tree (the two oracles agree with each other?)");
    }

    private static string Clip(string s) => s.Length > 50 ? s[..50] + "..." : s;

    private static readonly (LinkKernel Kernel, EndpointSource Source, int Loci)[] AllCombos =
        (from k in Enum.GetValues<LinkKernel>() from src in Enum.GetValues<EndpointSource>() from loci in new[] { 1, 2 } select (k, src, loci)).ToArray();

    [Test]
    public void ExhaustiveTinyTextsAndQueries()
    {
        foreach (string text in ExhaustiveStrings.All(2, 6))
        {
            UkkonenTree tree = UkkonenTree.Build(text);
            Setup s = Prepare(text, 2);
            foreach (string query in ExhaustiveStrings.All(2, 6).Where(q => q.Length > 0))
                for (int minLength = 1; minLength <= 4; minLength++)
                {
                    var combo = AllCombos[(text.Length * 7 + query.Length * 3 + minLength) % AllCombos.Length];
                    Check(s, tree, text, query, minLength, combo.Kernel, combo.Source, combo.Loci, "exhaustive");
                }
        }
    }

    [Test]
    public void ExhaustiveTernaryTexts_AllCombos()
    {
        foreach (string text in ExhaustiveStrings.All(3, 4).Where(t => t.Length > 0))
        {
            Setup s = Prepare(text, 2);
            foreach (string query in ExhaustiveStrings.All(3, 4).Where(q => q.Length > 0))
                foreach (var (kernel, source, loci) in AllCombos)
                    Check(s, null, text, query, 1 + (text.Length + query.Length) % 3, kernel, source, loci, "ternary");
        }
    }

    private static string Mutate(string text, SplitMix64 rng, int length, double substitutionRate, int indels)
    {
        int start = text.Length > length ? rng.NextInt(text.Length - length + 1) : 0;
        var chars = text.Substring(start, Math.Min(length, text.Length)).ToList();
        char[] alphabet = text.Distinct().ToArray();
        for (int i = 0; i < chars.Count; i++)
            if (rng.NextInt(1000) < (int)(substitutionRate * 1000)) chars[i] = alphabet[rng.NextInt(alphabet.Length)];
        for (int k = 0; k < indels && chars.Count > 2; k++)
        {
            int at = rng.NextInt(chars.Count);
            if (rng.NextInt(2) == 0) chars.RemoveAt(at);
            else chars.Insert(at, alphabet[rng.NextInt(alphabet.Length)]);
        }
        return new string(chars.ToArray());
    }

    [Test]
    public void AdversarialFamilies_MutatedQueries_AllCombos()
    {
        var rng = new SplitMix64(0x1B4A02);
        foreach (int n in new[] { 20, 100, 400 })
            foreach (var (name, text) in StringFamilies.Adversarial(n, 13).Where(f => f.Text.Distinct().Count() <= 256))
            {
                UkkonenTree tree = UkkonenTree.Build(text);
                Setup s = Prepare(text, 4);
                for (int rep = 0; rep < 3; rep++)
                {
                    string query = Mutate(text, rng, 10 + rng.NextInt(n), 0.05 * rep, rep);
                    int minLength = 1 + rng.NextInt(rep == 0 ? 3 : 12);
                    var combo = AllCombos[rng.NextInt(AllCombos.Length)];
                    Check(s, tree, text, query, minLength, combo.Kernel, combo.Source, combo.Loci, $"{name} n {n}");
                }
            }
    }

    [Test]
    public void RandomTexts_AllCombosOnAFewQueries()
    {
        var rng = new SplitMix64(0x1B4A03);
        foreach (int sigma in new[] { 2, 4, 20 })
            for (int rep = 0; rep < 12; rep++)
            {
                string text = StringFamilies.Random(50 + rng.NextInt(500), sigma, rng.Next());
                UkkonenTree tree = UkkonenTree.Build(text);
                Setup s = Prepare(text, 16);
                string query = Mutate(text, rng, 30 + rng.NextInt(300), 0.02 * (rep % 5), rep % 3);
                foreach (var (kernel, source, loci) in AllCombos)
                    Check(s, tree, text, query, 1 + rng.NextInt(10), kernel, source, loci, $"random sigma {sigma} rep {rep}");
            }
    }

    [Test]
    public void ForeignSymbols_EmptyAndOverlongQueries_AndDegenerateTexts()
    {
        string text = "abracadabraabracadabra";
        Setup s = Prepare(text);
        UkkonenTree tree = UkkonenTree.Build(text);
        foreach (string query in new[] { "", "z", "zzzz", "abrazabra", "abra", "cadabrazzabra", text, text + text, "ra", "xaxbxrxa", "aaaa", "\0abra", "abéra" })
            foreach (int minLength in new[] { 1, 2, 4, 30 })
                foreach (var (kernel, source, loci) in AllCombos)
                    Check(s, tree, text, query, minLength, kernel, source, loci, "special");

        // a^n against a^m, a single-symbol text, the empty text
        string a = new('a', 200);
        Setup sa = Prepare(a, 4);
        foreach (var (kernel, source, loci) in AllCombos)
        {
            Check(sa, null, a, new string('a', 77), 5, kernel, source, loci, "a^n");
            Check(sa, null, a, new string('a', 300), 100, kernel, source, loci, "a^n long query");
            Check(sa, null, a, "baaab", 1, kernel, source, loci, "a^n with foreign");
        }
        Setup one = Prepare("a");
        Check(one, null, "a", "aaa", 1, LinkKernel.RmqPsvNsv, EndpointSource.Psi, 2, "n=1");
        Setup empty = Prepare("");
        Assert.That(new MemDriver(empty.Esa, empty.Links, LinkKernel.Scan, EndpointSource.Isa, 1).Run("abc", 1), Is.Empty);
    }

    [Test]
    public void MinLength_AboveTextOrQuery_GivesNothing_AndIsValidated()
    {
        Setup s = Prepare("acgtacgt");
        var driver = new MemDriver(s.Esa, s.Links, LinkKernel.PsvNsv, EndpointSource.Psi, 1);
        Assert.That(driver.Run("acgt", 9), Is.Empty);
        Assert.That(driver.Run("acgtacgtacgt", 9), Is.Empty);
        Assert.Throws<ArgumentOutOfRangeException>(() => driver.Run("acgt", 0));
        Assert.Throws<ArgumentNullException>(() => driver.Run(null!, 1));
        Assert.Throws<ArgumentOutOfRangeException>(() => _ = new MemDriver(s.Esa, s.Links, LinkKernel.Scan, EndpointSource.Isa, 3));
    }

    [Test]
    public void Stats_CountLinksAndSteps_AndTwoLociFollowMoreLinksThanOne()
    {
        string text = StringFamilies.RandomDna(5000, 77);
        var rng = new SplitMix64(5);
        string query = Mutate(text, rng, 1500, 0.03, 0);
        Setup s = Prepare(text, 16);
        var one = new MemDriver(s.Esa, s.Links, LinkKernel.RmqPsvNsv, EndpointSource.Psi, 1);
        var two = new MemDriver(s.Esa, s.Links, LinkKernel.RmqPsvNsv, EndpointSource.Psi, 2);
        var r1 = Sorted(one.Run(query, 12));
        var r2 = Sorted(two.Run(query, 12));
        Assert.That(r2, Is.EqualTo(r1));
        Assert.That(one.Stats.Starts, Is.EqualTo(query.Length));
        Assert.That(one.Stats.Links, Is.GreaterThan(0));
        Assert.That(one.Stats.ExtendSteps, Is.GreaterThan(query.Length / 2));
        Assert.That(one.Stats.Emitted, Is.EqualTo(r1.Count));
        Assert.That(two.Stats.Links, Is.GreaterThan(one.Stats.Links), "the second locus follows its own links");
        Assert.That(two.Stats.Links, Is.LessThanOrEqualTo(2 * one.Stats.Links));
    }

    [Test]
    public void ARepeatedRun_GivesTheSameResult()
    {
        string text = StringFamilies.Fibonacci(3000);
        Setup s = Prepare(text);
        var driver = new MemDriver(s.Esa, s.Links, LinkKernel.RmqPsvNsv, EndpointSource.Isa, 2);
        string query = text.Substring(100, 500);
        var first = Sorted(driver.Run(query, 8));
        var again = Sorted(driver.Run(query, 8));
        Assert.That(again, Is.EqualTo(first));
        Assert.That(first, Is.EqualTo(Sorted(SuffixTreeAlgorithms.FindMaximalExactMatchesByDefinition(new StringTextSource(text), query, 8))));
    }
}
