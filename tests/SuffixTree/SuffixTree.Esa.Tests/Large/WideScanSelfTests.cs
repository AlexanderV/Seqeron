using SuffixTree.Esa.Tests.External;
using SuffixTree.Esa.Tests.Search;
using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.Esa.Tests.Large;

/// <summary>
/// Work package S2.4, the passes over the LCP array at 64 bits at a size where the 32-bit scan of S1 is the oracle (the default set; the passes over 2^31 + 10 generated values are <see cref="WideEngineLargeTests"/>): the wide
/// <see cref="LcpIntervalScan"/> gives the events and the statistics of the 32-bit scan on real texts and the closed forms on the generated ones, <see cref="LeftmostOccurrence"/> reads the wide views and refuses a text whose table
/// cannot exist, and <see cref="CompactLcp.LongCursor"/> reads values above <see cref="int.MaxValue"/> where a 32-bit cursor would refuse or wrap them.
/// </summary>
[TestFixture]
public unsafe class WideScanSelfTests
{
    private struct Recorder32 : IIntervalSink
    {
        internal List<string> Events;

        public readonly void LIndex(int previous, int k) => Events.Add($"L {previous} {k}");

        public readonly void Close(int ell, int lb, int rb, int first, bool lastChild) => Events.Add($"C {ell} {lb} {rb} {first} {lastChild}");
    }

    private struct Recorder64 : IWideIntervalSink
    {
        internal List<string> Events;

        public readonly void LIndex(long previous, long k) => Events.Add($"L {previous} {k}");

        public readonly void Close(long ell, long lb, long rb, long first, bool lastChild) => Events.Add($"C {ell} {lb} {rb} {first} {lastChild}");
    }

    public static IEnumerable<string> Texts()
    {
        var rng = new SplitMix64(0x5CA4);
        yield return "a";
        yield return "ab";
        yield return "abc";
        yield return "mississippi";
        yield return "abracadabra";
        yield return new string('a', 600);
        yield return string.Concat(Enumerable.Repeat("ab", 400));
        for (int i = 0; i < 6; i++) yield return new string([.. Enumerable.Range(0, 200 + (50 * i)).Select(_ => (char)('a' + rng.NextInt(2 + i)))]);
    }

    [Test]
    public void TheWideScan_GivesTheEventsAndTheStatisticsOfTheThirtyTwoBitScan([ValueSource(nameof(Texts))] string text)
    {
        EsaIndex index = EsaIndex.Build(text);
        CompactLcp lcp = index.Lcp;
        var sink32 = new Recorder32 { Events = [] };
        var stream32 = new CompactLcpStream(lcp);
        int peak32 = LcpIntervalScan.Run(text.Length, ref stream32, ref sink32);
        var sink64 = new Recorder64 { Events = [] };
        var stream64 = new WideCompactLcpStream(lcp);
        long peak64 = LcpIntervalScan.RunWide(text.Length, ref stream64, ref sink64);
        Assert.That(sink64.Events, Is.EqualTo(sink32.Events));
        Assert.That(peak64, Is.EqualTo(peak32));
        Assert.That(LcpIntervalScan.StatisticsOf(lcp), Is.EqualTo(index.Statistics));
    }

    [Test]
    public void TheWideScan_OfTheGeneratedLcpOfAClosedForm_GivesTheStatisticsOfTheClosedForm([Values(1L, 2L, 3L, 17L, 300L, 1500L)] long n, [Values(0, 1, 3, 11)] int shape)
    {
        ClosedFormText f = shape == 0 ? new PowerOfA(n) : new PeriodicAkB(n, shape);
        CompactLcp lcp = CompactLcp.FromReaders(LcpKind.U64, new PosGen<W8>(f, GenSection.LcpBase), new PosGen<W5>(f, GenSection.LcpPairsOfPowerOfA, long.MaxValue), isChecked: true);
        EsaStatistics st = LcpIntervalScan.StatisticsOf(lcp);
        Assert.That((st.Length, st.NodeCount, st.InternalCount, st.LeafCount, st.MaxDepth, st.LrsDepth, st.LrsRank), Is.EqualTo((f.Length, f.NodeCount, f.InternalCount, f.Length, f.Length, f.LrsDepth, f.LrsRank == -1 ? 0 : f.LrsRank)), f.ToString());
        Assert.That((Int128)EsaIndex.CountDistinctSubstrings(lcp, n), Is.EqualTo(((Int128)n * (n + 1) / 2) - f.LcpSum), f.ToString());
        long[] byLength = EsaIndex.CountDistinctSubstringsByLength(lcp, n, 5);
        Assert.That(byLength[0], Is.EqualTo(1));
        Assert.That(byLength, Has.Length.EqualTo(Math.Min(5, n) + 1));
    }

    [Test]
    public void TheStackPeak_OfTheWideScan_IsTheHeightOfTheStack()
    {
        Assert.That(LcpIntervalScan.StatisticsOf(EsaIndex.Build("abc").Lcp).StackPeak, Is.EqualTo(1), "a text with no repeat has the root only");
        Assert.That(LcpIntervalScan.StatisticsOf(EsaIndex.Build("aaaa").Lcp).StackPeak, Is.EqualTo(EsaIndex.Build("aaaa").Statistics.StackPeak));
        Assert.That(LcpIntervalScan.StatisticsOf(EsaIndex.Build("aaaa").Lcp).StackPeak, Is.EqualTo(4), "a^4: the root and the nodes of a, aa, aaa");
    }

    [Test]
    public void TheWideScan_RefusesALcpArrayThatNoTextHas()
    {
        var negative = new FixedStream([0L, 1, -1, 0]);
        Assert.That(() => RunStream(negative), Throws.TypeOf<ArgumentException>().With.Message.Contains("LCP[3]"));
        var startsAtOne = new FixedStream([1L, 0, 0, 0]);
        Assert.That(() => RunStream(startsAtOne), Throws.TypeOf<ArgumentException>().With.Message.Contains("LCP[1]"));
        Assert.That(() => RunStream(new FixedStream([]), -1), Throws.TypeOf<ArgumentOutOfRangeException>());
        var empty = new FixedStream([]);
        var sink = new WideStatisticsSink();
        Assert.That(LcpIntervalScan.RunWide(0, ref empty, ref sink), Is.EqualTo(0));
    }

    private static void RunStream(FixedStream stream, long n = 4)
    {
        var sink = new WideStatisticsSink();
        LcpIntervalScan.RunWide(n, ref stream, ref sink);
    }

    private struct FixedStream(long[] values) : IWideLcpStream
    {
        private int _next;

        public long Next() => values[_next++];
    }

    // ------------------------------------------------------------------------------------------ LeftmostOccurrence over the wide views

    [Test]
    public void LeftmostOccurrence_OverTheWideViews_EqualsTheInMemoryTable([ValueSource(nameof(Texts))] string text)
    {
        EsaIndex memory = EsaIndex.Build(text);
        int[] expected = LeftmostOccurrence.Compute(memory.SuffixArray, memory.Lcp);
        LcpChildCase c = LcpChildCase.Of(text);
        char[] alphabet = [.. text.Distinct().Order()];
        foreach ((int width, LcpKind kind, bool tagged) in WideImage.Layouts)
        {
            using var image = new WideImage(c, width, kind, tagged);
            EsaIndex wide = image.Index(isChecked: true, alphabet);
            Assert.That(LeftmostOccurrence.Compute(wide.Sa, wide.Lcp), Is.EqualTo(expected), $"width {width} {kind} {(tagged ? "tagged" : string.Empty)}");
        }
    }

    [Test]
    public void LeftmostOccurrence_RefusesATextWhoseTableCannotExist_BeforeAnyAllocation()
    {
        foreach (long n in new[] { Array.MaxLength, Array.MaxLength + 1L, (1L << 31) + 10 })
        {
            var f = new PowerOfA(n);
            CompactLcp lcp = CompactLcp.FromReaders(LcpKind.U64, new PosGen<W8>(f, GenSection.LcpBase), new PosGen<W5>(f, GenSection.LcpPairsOfPowerOfA, long.MaxValue), isChecked: true);
            var sa = new SaView((int*)null, n + 1, isChecked: true);
            NotSupportedException e = Assert.Throws<NotSupportedException>(() => LeftmostOccurrence.Compute(sa, lcp), $"n = {n}")!;
            Assert.That(e.Message, Does.Contain("traversal").And.Contain(n.ToString(System.Globalization.CultureInfo.InvariantCulture)));
        }
    }

    [Test]
    public void LeftmostOccurrence_RefusesAnLcpArrayOfAnotherLength()
    {
        EsaIndex a = EsaIndex.Build("banana");
        EsaIndex b = EsaIndex.Build("bananas");
        Assert.That(() => LeftmostOccurrence.Compute(a.Sa, b.Lcp), Throws.TypeOf<ArgumentException>());
        Assert.That(() => LeftmostOccurrence.Compute(b.Sa, a.Lcp), Throws.TypeOf<ArgumentException>());
    }

    // ------------------------------------------------------------------------------------------ the long cursor

    /// <summary>A reader over an array of 64-bit values: a small LCP array that holds values above 2^32.</summary>
    private readonly struct ArrayReader(ulong[] values) : IPosReader
    {
        public static int ElementBytes => 8;

        public long Count => values.Length;

        public ulong Get(long index) => (ulong)index < (ulong)values.Length ? values[index] : throw new InvalidDataException("outside");
    }

    private readonly struct PairReader(ulong[] values) : IPosReader
    {
        public static int ElementBytes => 5;

        public long Count => values.Length;

        public ulong Get(long index) => (ulong)index < (ulong)values.Length ? values[index] : throw new InvalidDataException("outside");
    }

    [Test]
    public void TheLongCursor_ReadsValuesAboveInt32_OfASmallWideArray_AndTheStoreAgrees()
    {
        ulong[] values = [0, 0, 5_000_000_000, 7, 1UL << 33, 3, (1UL << 32) + 1];
        CompactLcp lcp = CompactLcp.FromReaders(LcpKind.U64, new ArrayReader(values), new PairReader([]), isChecked: false);
        CompactLcp.LongCursor cursor = lcp.CreateLongCursor();
        for (long r = 0; r < values.Length; r++)
        {
            Assert.That(cursor.Get(r), Is.EqualTo((long)values[r]), $"rank {r}");
            Assert.That(lcp[r], Is.EqualTo((long)values[r]));
        }
        Assert.That(() => cursor.Get(values.Length), Throws.TypeOf<InvalidDataException>());
        Assert.That(() => cursor.Get(-1), Throws.TypeOf<InvalidDataException>());
        Assert.That(() => lcp.CreateCursor().Get(2), Throws.TypeOf<InvalidDataException>(), "the 32-bit cursor refuses a value that is no int, it does not wrap it");
    }

    [Test]
    public void TheLongCursor_ResolvesTheEscapesOfASmallWideArray_InAnyOrder()
    {
        // an 8-bit array with escapes at ranks 2 and 5: the pairs (2, 700) and (5, 90000)
        ulong[] baseValues = [0, 0, 255, 1, 2, 255];
        CompactLcp lcp = CompactLcp.FromReaders(LcpKind.U8, new Pos1(baseValues), new PairReader([2, 700, 5, 90_000]), isChecked: false);
        CompactLcp.LongCursor cursor = lcp.CreateLongCursor();
        foreach (long r in new long[] { 0, 1, 2, 3, 4, 5, 5, 2, 0 })
            Assert.That(cursor.Get(r), Is.EqualTo(r == 2 ? 700 : r == 5 ? 90_000 : (long)baseValues[r]), $"rank {r}");
        Assert.That(() => cursor.Get(6), Throws.TypeOf<InvalidDataException>());
    }

    private readonly struct Pos1(ulong[] values) : IPosReader
    {
        public static int ElementBytes => 1;

        public long Count => values.Length;

        public ulong Get(long index) => (ulong)index < (ulong)values.Length ? values[index] : throw new InvalidDataException("outside");
    }

    [Test]
    public void TheLongCursor_OfTheThirtyTwoBitForms_EqualsTheIntegerCursor_AndRefusesARankOutsideTheArray([ValueSource(nameof(Texts))] string text)
    {
        EsaIndex index = EsaIndex.Build(text);
        CompactLcp.Cursor narrow = index.Lcp.CreateCursor();
        CompactLcp.LongCursor wide = index.Lcp.CreateLongCursor();
        for (int r = 0; r <= text.Length; r++) Assert.That(wide.Get(r), Is.EqualTo((long)narrow.Get(r)), $"rank {r}");
        Assert.That(() => wide.Get(text.Length + 1), Throws.InstanceOf<ArgumentOutOfRangeException>().Or.InstanceOf<InvalidDataException>());
        Assert.That(() => wide.Get(-1), Throws.InstanceOf<ArgumentOutOfRangeException>().Or.InstanceOf<InvalidDataException>());
        Assert.That(() => wide.Get((1L << 32) + 1), Throws.InstanceOf<ArgumentOutOfRangeException>().Or.InstanceOf<InvalidDataException>(), "a rank that wraps to a valid one in 32 bits");
    }
}
