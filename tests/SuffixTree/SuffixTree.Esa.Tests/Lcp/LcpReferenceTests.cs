using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.Esa.Tests.Lcp;

/// <summary>1a.4: Phi / PLCP LCP construction against the naive and Kasai oracles on every corpus of the strategy.</summary>
[TestFixture]
public class LcpReferenceTests
{
    private static CompactLcp BuildFromSa(string text, int[] sa)
    {
        var p = TextPreparation.Prepare(text);
        return p.Width == SymbolWidth.Byte ? LcpReference.Build(p.ByteReader(), sa) : LcpReference.Build(p.UShortReader(), sa);
    }

    private static void AssertLcp(string text, int[] sa, string context)
    {
        int[] expected = NaiveLcp.Build(text, sa);
        CompactLcp actual = BuildFromSa(text, sa);
        int[] dense = actual.ToArray();
        if (!expected.AsSpan().SequenceEqual(dense))
            Assert.Fail($"{context}: LCP differs from the naive oracle (n = {text.Length})");
        var cursor = actual.CreateCursor();
        for (int r = 0; r < dense.Length; r++)
            if (cursor.Get(r) != expected[r]) Assert.Fail($"{context}: cursor at rank {r}");
    }

    [Test]
    public void KnownExamples()
    {
        Assert.That(BuildFromSa("banana", NaiveSuffixArray.Build("banana")).ToArray(), Is.EqualTo(new[] { 0, 0, 1, 3, 0, 0, 2 }));
        Assert.That(BuildFromSa("", NaiveSuffixArray.Build("")).ToArray(), Is.EqualTo(new[] { 0 }));
        Assert.That(BuildFromSa("x", NaiveSuffixArray.Build("x")).ToArray(), Is.EqualTo(new[] { 0, 0 }));
        Assert.That(BuildFromSa("aaaa", NaiveSuffixArray.Build("aaaa")).ToArray(), Is.EqualTo(new[] { 0, 0, 1, 2, 3 }));
        Assert.That(BuildFromSa("abab", NaiveSuffixArray.Build("abab")).ToArray(), Is.EqualTo(new[] { 0, 0, 2, 0, 1 }));
    }

    [Test]
    public void Exhaustive_AllStrings_Sigma2UpTo14_Sigma3UpTo9_Sigma4UpTo7_EqualNaiveAndKasai()
    {
        long count = 0;
        foreach (var (sigma, maxLen) in new[] { (2, 14), (3, 9), (4, 7) })
            foreach (string s in ExhaustiveStrings.All(sigma, maxLen))
            {
                int[] sa = NaiveSuffixArray.Build(s);
                AssertLcp(s, sa, $"sigma {sigma} '{s}'");
                if (!KasaiLcp.Build(s, sa).AsSpan().SequenceEqual(BuildFromSa(s, sa).ToArray()))
                    Assert.Fail($"sigma {sigma}: PLCP LCP differs from Kasai on '{s}'");
                count++;
            }
        Assert.That(count, Is.EqualTo(32767 + 29524 + 21845));
    }

    [TestCase(1)]
    [TestCase(2)]
    [TestCase(3)]
    [TestCase(4)]
    [TestCase(20)]
    [TestCase(256)]
    [TestCase(1000)]
    public void RandomTexts_UpTo2000_EqualNaive(int sigma)
    {
        var rng = new SplitMix64(2000UL + (ulong)sigma);
        for (int i = 0; i < 300; i++)
        {
            int n = rng.NextInt(2001);
            ulong seed = rng.Next();
            string s = StringFamilies.Random(n, sigma, seed);
            AssertLcp(s, NaiveSuffixArray.Build(s), $"random sigma {sigma} seed {seed} n {n}");
        }
    }

    [Test]
    public void AdversarialFamilies_UpTo2000_EqualNaive_FromBothTheNaiveAndTheSaisSuffixArray()
    {
        foreach (int n in new[] { 0, 1, 2, 3, 4, 5, 7, 8, 16, 17, 33, 64, 65, 128, 129, 255, 256, 257, 500, 1000, 1024, 2000 })
            foreach (var (name, text) in StringFamilies.Adversarial(n, 4))
            {
                AssertLcp(text, NaiveSuffixArray.Build(text), $"{name} n {n} (naive SA)");
                AssertLcp(text, SaisReference.BuildSuffixArray(text), $"{name} n {n} (SA-IS)");
            }
    }

    [Test]
    public void AllSameSymbol_70000_ForcesEscapesAtBoth255And65535()
    {
        const int n = 70_000;
        string s = StringFamilies.AllSame(n);
        int[] sa = SaisReference.BuildSuffixArray(s);
        CompactLcp lcp = BuildFromSa(s, sa);

        // LCP[r] = r - 1 for r >= 1: the values >= 255 are r = 256..n (69,745 of them), the values >= 65,535 are r = 65,536..n (4,465)
        for (int r = 1; r <= n; r++)
            if (lcp[r] != r - 1) Assert.Fail($"LCP[{r}] = {lcp[r]}");
        Assert.That(lcp[0], Is.EqualTo(0));
        Assert.That(lcp.Kind, Is.EqualTo(LcpKind.U16)); // 69,745 escapes rule out U8; 8 * 4,465 <= 70,000 keeps U16
        Assert.That(lcp.OverflowCount, Is.EqualTo(4465));
        Assert.That(lcp.ToArray().AsSpan().SequenceEqual(KasaiLcp.Build(s, sa)), Is.True);

        // the same values forced into every representation round-trip exactly
        int[] dense = lcp.ToArray();
        var u8 = CompactLcp.FromDense(dense, LcpKind.U8);
        Assert.That(u8.OverflowCount, Is.EqualTo(69_745));
        Assert.That(u8.ToArray().AsSpan().SequenceEqual(dense), Is.True);
        Assert.That(CompactLcp.FromDense(dense, LcpKind.U32).ToArray().AsSpan().SequenceEqual(dense), Is.True);
    }

    [TestCase(254, 0)]
    [TestCase(255, 1)]
    [TestCase(256, 2)]
    public void TwoCopiesOfARandomBlock_PlaceTheLargestLcpAroundTheU8EscapeValue(int block, int expectedU8Escapes)
    {
        // text = R R: LCP(suffix i, suffix i + L) = L - i, so the LCP values are L, L-1, ..., 1 (plus the small random ones)
        string r = StringFamilies.RandomDna(block, 77);
        string s = r + r;
        int[] sa = SaisReference.BuildSuffixArray(s);
        CompactLcp lcp = BuildFromSa(s, sa);
        Assert.That(lcp.Kind, Is.EqualTo(LcpKind.U8));
        Assert.That(lcp.OverflowCount, Is.EqualTo(expectedU8Escapes));
        Assert.That(lcp.ToArray().AsSpan().SequenceEqual(KasaiLcp.Build(s, sa)), Is.True);
        Assert.That(lcp.ToArray().Max(), Is.EqualTo(block));
    }

    [TestCase(65534, 0)]
    [TestCase(65535, 1)]
    [TestCase(65536, 2)]
    public void TwoCopiesOfARandomBlock_PlaceTheLargestLcpAroundTheU16EscapeValue(int block, int expectedU16Escapes)
    {
        string r = StringFamilies.RandomDna(block, 78);
        string s = r + r;
        int[] sa = SaisReference.BuildSuffixArray(s);
        CompactLcp lcp = BuildFromSa(s, sa);
        Assert.That(lcp.Kind, Is.EqualTo(LcpKind.U16)); // ~65,000 values >= 255 rule out U8; at most 2 values >= 65,535
        Assert.That(lcp.OverflowCount, Is.EqualTo(expectedU16Escapes));
        Assert.That(lcp.ToArray().AsSpan().SequenceEqual(KasaiLcp.Build(s, sa)), Is.True);
        Assert.That(lcp.ToArray().Max(), Is.EqualTo(block));
    }

    [Test]
    public void RandomDna_UsesU8WithoutEscapes_AndAShortRepeatKeepsU8WithFewEscapes()
    {
        string dna = StringFamilies.RandomDna(100_000, 11);
        var plain = BuildFromSa(dna, SaisReference.BuildSuffixArray(dna));
        Assert.That(plain.Kind, Is.EqualTo(LcpKind.U8));
        Assert.That(plain.OverflowCount, Is.EqualTo(0));

        string withRun = StringFamilies.AllSame(400) + StringFamilies.RandomDna(20_000, 12);
        var lcp = BuildFromSa(withRun, SaisReference.BuildSuffixArray(withRun));
        Assert.That(lcp.Kind, Is.EqualTo(LcpKind.U8));
        Assert.That(lcp.OverflowCount, Is.InRange(100, 400));
        Assert.That(lcp.ToArray().AsSpan().SequenceEqual(KasaiLcp.Build(withRun, SaisReference.BuildSuffixArray(withRun))), Is.True);
    }

    [Test]
    public void Utf16Texts_EqualNaive()
    {
        foreach (ulong seed in new ulong[] { 1, 2, 3 })
            foreach (string s in new[]
            {
                StringFamilies.Bmp(1500, seed),
                StringFamilies.SurrogatePairs(1500, seed),
                StringFamilies.LoneSurrogates(1500, seed),
                StringFamilies.NulAndFfff(1500, seed),
                StringFamilies.AllDistinct(1500, seed),
            })
                AssertLcp(s, NaiveSuffixArray.Build(s), "utf-16");
    }

    [Test]
    public void LargeStructuredTexts_EqualKasai()
    {
        foreach (var (name, text) in new[]
        {
            ("fibonacci", StringFamilies.Fibonacci(200_000)),
            ("thue-morse", StringFamilies.ThueMorse(200_000)),
            ("abab", StringFamilies.AbRepeat(200_000)),
            ("random-protein", StringFamilies.RandomProtein(200_000, 5)),
            ("periodic-aab-p65", StringFamilies.PeriodicAaB(200_000, 65)),
        })
        {
            int[] sa = SaisReference.BuildSuffixArray(text);
            Assert.That(BuildFromSa(text, sa).ToArray().AsSpan().SequenceEqual(KasaiLcp.Build(text, sa)), Is.True, name);
        }
    }

    [Test]
    public void Plcp_IsTheLcpPermutedIntoTextOrder()
    {
        foreach (var (name, text) in StringFamilies.Adversarial(300, 8))
        {
            int[] sa = NaiveSuffixArray.Build(text);
            int[] naive = NaiveLcp.Build(text, sa);
            var p = TextPreparation.Prepare(text);
            int[] plcp = p.Width == SymbolWidth.Byte ? LcpReference.ComputePlcp(p.ByteReader(), sa) : LcpReference.ComputePlcp(p.UShortReader(), sa);
            Assert.That(plcp, Has.Length.EqualTo(text.Length), name);
            for (int r = 1; r < sa.Length; r++)
                if (plcp[sa[r]] != naive[r]) Assert.Fail($"{name}: PLCP[SA[{r}]] = {plcp[sa[r]]}, naive LCP[{r}] = {naive[r]}");
        }
    }

    [Test]
    public void InvalidSuffixArrays_AreRejected()
    {
        byte[] symbols = { 0, 1, 0 };
        Assert.Throws<ArgumentException>(() => LcpReference.Build(new ByteSymbols(symbols), new[] { 3, 2, 0 }));      // wrong length
        Assert.Throws<ArgumentException>(() => LcpReference.Build(new ByteSymbols(symbols), new[] { 0, 1, 2, 3 }));   // sa[0] != n
    }
}
