using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.Esa.Tests.Sais;

/// <summary>
/// 1a.5 (decision 2026-11-03, review MAJ-5): only the safe reference is generic over the index type; the <c>long</c> instantiation is an
/// oracle for the index arithmetic (it must give the identical SA and LCP). The optimized kernels stay <c>int</c>-only.
/// </summary>
[TestFixture]
public class IndexWidthTests
{
    private static void AssertWidthsAgree(string text, string context)
    {
        int[] narrow = SaisReference<int>.BuildSuffixArray(text);
        long[] wide = SaisReference<long>.BuildSuffixArray(text);
        if (narrow.Length != wide.Length || !narrow.Select(x => (long)x).SequenceEqual(wide))
            Assert.Fail($"{context}: SaisReference<long> differs from SaisReference<int> (n = {text.Length})");

        var p = TextPreparation.Prepare(text);
        CompactLcp lcpNarrow = p.Width == SymbolWidth.Byte ? LcpReference<int>.Build(p.ByteReader(), narrow) : LcpReference<int>.Build(p.UShortReader(), narrow);
        CompactLcp lcpWide = p.Width == SymbolWidth.Byte ? LcpReference<long>.Build(p.ByteReader(), wide) : LcpReference<long>.Build(p.UShortReader(), wide);
        if (lcpNarrow.Kind != lcpWide.Kind || lcpNarrow.OverflowCount != lcpWide.OverflowCount || !lcpNarrow.ToArray().AsSpan().SequenceEqual(lcpWide.ToArray()))
            Assert.Fail($"{context}: LcpReference<long> differs from LcpReference<int> (n = {text.Length})");
    }

    [Test]
    public void KnownExample_Banana()
    {
        Assert.That(SaisReference<long>.BuildSuffixArray("banana"), Is.EqualTo(new long[] { 6, 5, 3, 1, 0, 4, 2 }));
        Assert.That(SaisReference<long>.BuildSuffixArray(""), Is.EqualTo(new long[] { 0 }));
        Assert.That(SaisReference<int>.BuildSuffixArray("banana"), Is.EqualTo(SaisReference.BuildSuffixArray("banana")));
    }

    [Test]
    public void Exhaustive_AllStrings_Sigma2UpTo14_Sigma3UpTo9_Sigma4UpTo7()
    {
        long count = 0;
        foreach (var (sigma, maxLen) in new[] { (2, 14), (3, 9), (4, 7) })
            foreach (string s in ExhaustiveStrings.All(sigma, maxLen))
            {
                AssertWidthsAgree(s, $"sigma {sigma} '{s}'");
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
    public void RandomTexts_UpTo2000(int sigma)
    {
        var rng = new SplitMix64(3000UL + (ulong)sigma);
        for (int i = 0; i < 100; i++)
        {
            int n = rng.NextInt(2001);
            ulong seed = rng.Next();
            AssertWidthsAgree(StringFamilies.Random(n, sigma, seed), $"random sigma {sigma} seed {seed} n {n}");
        }
    }

    [Test]
    public void AdversarialFamilies_UpTo2000()
    {
        foreach (int n in new[] { 0, 1, 2, 3, 4, 5, 8, 17, 64, 65, 129, 256, 257, 1000, 2000 })
            foreach (var (name, text) in StringFamilies.Adversarial(n, 5))
                AssertWidthsAgree(text, $"{name} n {n}");
    }

    [Test]
    public void LongTexts_AllSameAndFibonacci_AgreeThroughTheVerifier()
    {
        foreach (var (name, text) in new[] { ("a^n", StringFamilies.AllSame(70_000)), ("fibonacci", StringFamilies.Fibonacci(70_000)) })
        {
            long[] wide = SaisReference<long>.BuildSuffixArray(text);
            Assert.That(SuffixArrayVerifier.Check(text, wide.Select(x => (int)x).ToArray()), Is.Null, name);
            Assert.That(wide.Select(x => (int)x), Is.EqualTo(SaisReference.BuildSuffixArray(text)), name);
        }
    }

    [Test]
    public void RecursionTraces_AreIdentical()
    {
        foreach (var (name, text) in StringFamilies.Adversarial(1500, 6))
        {
            var narrow = new SaisTrace();
            var wide = new SaisTrace();
            SaisReference<int>.BuildSuffixArray(text, narrow);
            SaisReference<long>.BuildSuffixArray(text, wide);
            Assert.That(wide.Levels, Is.EqualTo(narrow.Levels), name);
        }
    }

    [Test]
    public void PlcpAndEscapeCounts_AreIdentical()
    {
        string text = StringFamilies.AllSame(300) + StringFamilies.RandomDna(2000, 3);
        var p = TextPreparation.Prepare(text);
        int[] narrow = SaisReference<int>.BuildSuffixArray(text);
        long[] wide = SaisReference<long>.BuildSuffixArray(text);
        int[] plcpNarrow = LcpReference<int>.ComputePlcp(p.ByteReader(), narrow);
        long[] plcpWide = LcpReference<long>.ComputePlcp(p.ByteReader(), wide);
        Assert.That(plcpWide, Is.EqualTo(plcpNarrow.Select(x => (long)x)));
    }

    [Test]
    public void TheFacadesDelegateToTheIntInstantiation()
    {
        string text = StringFamilies.RandomDna(500, 9);
        Assert.That(SaisReference.BuildSuffixArray(text), Is.EqualTo(SaisReference<int>.BuildSuffixArray(text)));
        var p = TextPreparation.Prepare(text);
        int[] sa = SaisReference.BuildSuffixArray(text);
        Assert.That(LcpReference.Build(p.ByteReader(), sa).ToArray(), Is.EqualTo(LcpReference<int>.Build(p.ByteReader(), sa).ToArray()));
    }

    [Test]
    public void ConversionsAreChecked_ANarrowIndexTypeThrowsInsteadOfWrapping()
    {
        // short is not a supported index type; it only demonstrates that every conversion is checked
        string small = StringFamilies.RandomDna(200, 1);
        Assert.That(SaisReference<short>.BuildSuffixArray(small).Select(x => (int)x), Is.EqualTo(SaisReference.BuildSuffixArray(small)));
        Assert.Throws<OverflowException>(() => SaisReference<short>.BuildSuffixArray(StringFamilies.AllSame(40_000)));
    }

    [Test]
    public void GuardsAndArgumentChecks_AreTheSameForBothWidths()
    {
        Assert.Throws<NotSupportedException>(() => SaisReference<long>.Build(new FakeReader(int.MaxValue), 4));
        Assert.Throws<NotSupportedException>(() => SaisReference<long>.Build(new FakeReader(Array.MaxLength), 4, new long[4]));
        Assert.Throws<ArgumentException>(() => SaisReference<long>.Build(new IntSymbols(new[] { 0, 2 }), 2));
        Assert.Throws<ArgumentException>(() => SaisReference<long>.Build(new IntSymbols(new[] { 0, 1 }), 2, new long[2]));
        Assert.Throws<ArgumentOutOfRangeException>(() => SaisReference<long>.Build(new IntSymbols(new[] { 0, 1 }), -1));
        byte[] symbols = { 0, 1, 0 };
        Assert.Throws<ArgumentException>(() => LcpReference<long>.Build(new ByteSymbols(symbols), new long[] { 3, 2, 0 }));
        Assert.Throws<ArgumentException>(() => LcpReference<long>.Build(new ByteSymbols(symbols), new long[] { 0, 1, 2, 3 }));
    }

    private readonly struct FakeReader(int length) : ISymbolReader
    {
        public int Length { get; } = length;

        public int this[int index] => throw new InvalidOperationException("the size guard must fire before any symbol is read");
    }
}
