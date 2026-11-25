using System.IO.MemoryMappedFiles;
using SuffixTree.Esa.External;
using SuffixTree.Esa.Tests.Views;
using SuffixTree.Indexing;

namespace SuffixTree.Esa.Tests.Large;

/// <summary>
/// Work package S2.3, the generator-backed views at n = 2^31 + 10, 2^32 + 10 and 2^33 (ADR-005 section 3.7, layers 1 and 2; category <c>Large</c>, explicit): the views over <see cref="PosGen{TWidth}"/> read the closed forms of S2.1 (<c>a^n</c>,
/// <c>(a^k b)*</c>) with no memory, and the views over the sparse shells (<see cref="SparseContainer"/>; the real ends of <c>a^n</c> in a file of 24 - 78 GB nominal) read the same numbers through a mapping. Every read at the ends, at powers of two
/// and at random ranks is compared with the closed form; a bulk copy across 2^31 and 2^32; the LCP of 64 bits (a depth of <c>a^n</c> above 2^32), the overflow list searched at 2^33 pairs; a descent through the generated child cells; the pyramid
/// over 2^31 + 10 entries; the text of 2^33 symbols. A <c>long</c> to <c>int</c> narrowing anywhere on these paths shows as a wrong number or as the wrong exception.
/// </summary>
[TestFixture]
[Category("Large")]
[Explicit(ExplicitReason)]
public unsafe class GeneratedViewLargeTests : LargeFixtureBase
{
    private const long TwoTo31 = 1L << 31;
    private const long TwoTo32 = 1L << 32;
    private const long TwoTo33 = 1L << 33;

    internal static IEnumerable<ClosedFormText> Huge()
    {
        foreach (long n in new[] { TwoTo31 + 10, TwoTo32 + 10, TwoTo33 })
        {
            yield return new PowerOfA(n);
            yield return new PeriodicAkB(n, 1);
            yield return new PeriodicAkB(n, 3);
        }
    }

    /// <summary>The ends, the powers of two and their neighbours, and 200 random ranks.</summary>
    internal static long[] Ranks(long n)
    {
        var set = new SortedSet<long> { 0, 1, 2, 3, 4, n - 3, n - 2, n - 1, n };
        for (int k = 8; k < 41; k++)
        {
            long p = 1L << k;
            if (p > n) break;
            foreach (long d in new long[] { -1, 0, 1 })
            {
                if (p + d >= 0 && p + d <= n) set.Add(p + d);
            }
        }
        ulong x = (ulong)n * 0x9E3779B97F4A7C15UL;
        for (int i = 0; i < 200; i++)
        {
            x = (x ^ (x >> 29)) * 0xBF58476D1CE4E5B9UL + 1;
            set.Add((long)(x % (ulong)(n + 1)));
        }
        return [.. set];
    }

    // ---------------------------------------------------------------------------------------------- the suffix array

    [Test]
    public void SaView_AtTheEndsAtPowersOfTwoAndAtRandomRanks_EqualsTheClosedForm([ValueSource(nameof(Huge))] ClosedFormText f)
    {
        var view = new SaView<PosGen<W5>>(new PosGen<W5>(f, GenSection.Sa));
        Assert.That(view.Length, Is.EqualTo(f.Length + 1));
        foreach (long r in Ranks(f.Length))
        {
            Assert.That(view[r], Is.EqualTo(f.SuffixAt(r)), $"{f} rank {r}");
            Assert.That(view.At(r), Is.EqualTo(f.SuffixAt(r)), $"{f} At({r})");
        }
        Assert.That(() => view.At(f.Length + 1), Throws.TypeOf<InvalidDataException>());
        Assert.That(() => view.At(long.MaxValue), Throws.TypeOf<InvalidDataException>());
    }

    [Test]
    public void SaView_CopiesBlocksAcross2To31And2To32_AndNeverNarrowsAPosition([ValueSource(nameof(Huge))] ClosedFormText f)
    {
        var view = new SaView<PosGen<W5>>(new PosGen<W5>(f, GenSection.Sa));
        foreach (long boundary in new[] { TwoTo31, TwoTo32, f.Length - 5 })
        {
            if (boundary + 6 > f.Length) continue;
            long start = boundary - 5;
            var longs = new long[12];
            view.CopyTo(start, longs);
            for (int i = 0; i < longs.Length; i++) Assert.That(longs[i], Is.EqualTo(f.SuffixAt(start + i)), $"{f} block at {start}, element {i}");
            Assert.That(view.ToArray(start, 12), Is.EqualTo(longs));
            var ints = new int[12];
            if (longs.All(v => v <= int.MaxValue))
            {
                view.CopyTo(start, ints);
                Assert.That(ints, Is.EqualTo(longs.Select(v => (int)v)));
            }
            else
            {
                Assert.That(() => view.CopyTo(start, ints), Throws.TypeOf<InvalidDataException>(), $"{f}: a position above int.MaxValue in a block of 32-bit positions");
            }
        }
    }

    // ---------------------------------------------------------------------------------------------- the LCP array

    [Test]
    public void Lcp_U64_AtTheEndsAtPowersOfTwoAndAtRandomRanks_EqualsTheClosedForm([ValueSource(nameof(Huge))] ClosedFormText f)
    {
        CompactLcp lcp = CompactLcp.FromReaders(LcpKind.U64, new PosGen<W8>(f, GenSection.LcpBase), new PosGen<W5>(f, GenSection.LcpPairsOfPowerOfA, long.MaxValue), isChecked: true);
        Assert.That(lcp.Length, Is.EqualTo(f.Length + 1));
        Assert.That(lcp.OverflowCount, Is.EqualTo(0));
        int above = 0;
        foreach (long r in Ranks(f.Length))
        {
            long expected = f.Lcp(r);
            Assert.That(lcp[r], Is.EqualTo(expected), $"{f} rank {r}");
            if (expected > int.MaxValue)
            {
                above++;
                // the 32-bit indexer at a rank that it can address says that the value is no int, it does not wrap it
                if (r <= int.MaxValue) Assert.That(() => lcp[(int)r], Throws.TypeOf<InvalidDataException>(), $"{f} rank {r}");
            }
        }
        Assert.That(above, Is.GreaterThan(0), "the test must read values above int.MaxValue");
        Assert.That(() => lcp[f.Length + 1], Throws.TypeOf<InvalidDataException>());
    }

    [Test]
    public void Lcp_AValueAboveTheIntegers_IsNotNarrowedByTheIntegerIndexer()
    {
        var f = new PowerOfA(TwoTo33);
        CompactLcp lcp = CompactLcp.FromReaders(LcpKind.U64, new PosGen<W8>(f, GenSection.LcpBase), new PosGen<W5>(f, GenSection.LcpPairsOfPowerOfA, long.MaxValue), isChecked: true);
        Assert.That(lcp[TwoTo32 + 1], Is.EqualTo(TwoTo32));
        Assert.That(lcp[TwoTo31 + 1], Is.EqualTo(TwoTo31));
        Assert.That(lcp[TwoTo31], Is.EqualTo(TwoTo31 - 1));
        Assert.That(lcp[int.MaxValue], Is.EqualTo(int.MaxValue - 1));
        Assert.That(() => lcp.CreateCursor().Get(int.MaxValue), Throws.Nothing);
        Assert.That(() => lcp.ToArray(), Throws.TypeOf<NotSupportedException>());
    }

    [Test]
    public void Lcp_U32_OfTheTextOf2To31Plus10_EqualsTheClosedForm([Values(1, 3)] int k)
    {
        ClosedFormText f = new PeriodicAkB(TwoTo31 + 10, k);
        CompactLcp lcp = CompactLcp.FromReaders(LcpKind.U32, new PosGen<W4>(f, GenSection.LcpBase), new PosGen<W4>(f, GenSection.LcpPairsOfPowerOfA, long.MaxValue), isChecked: true);
        foreach (long r in Ranks(f.Length)) Assert.That(lcp[r], Is.EqualTo(f.Lcp(r)), $"{f} rank {r}");
    }

    [Test]
    public void Lcp_TheOverflowList_OfPowerOfA_IsSearchedOverAHugeNumberOfPairs([Values(255, 65_535)] int escape, [Values(2147483658L, 4294967306L, 8589934592L)] long n)
    {
        var f = new PowerOfA(n);
        LcpKind kind = escape == 255 ? LcpKind.U8 : LcpKind.U16;
        CompactLcp lcp = kind == LcpKind.U8
            ? CompactLcp.FromReaders(kind, new PosGen<W1>(f, GenSection.LcpBase, escape), new PosGen<W5>(f, GenSection.LcpPairsOfPowerOfA, escape), isChecked: true)
            : CompactLcp.FromReaders(kind, new PosGen<W2>(f, GenSection.LcpBase, escape), new PosGen<W5>(f, GenSection.LcpPairsOfPowerOfA, escape), isChecked: true);
        Assert.That(lcp.OverflowCount, Is.EqualTo(n - escape));
        var ranks = new List<long>(Ranks(n)) { escape - 1, escape, escape + 1, escape + 2, escape + 3 };
        foreach (long r in ranks) Assert.That(lcp[r], Is.EqualTo(f.Lcp(r)), $"a^{n} rank {r} escape {escape}");
        Assert.That(lcp[n], Is.EqualTo(n - 1));
        if (n > TwoTo32) Assert.That(lcp[TwoTo32], Is.EqualTo(TwoTo32 - 1));
    }

    // ---------------------------------------------------------------------------------------------- the child cells and a descent over them

    [Test]
    public void ChildCells_Packed_ReadTheCellsOfTheClosedForm([ValueSource(nameof(Huge))] ClosedFormText f)
    {
        CompactLcp lcp = CompactLcp.FromReaders(LcpKind.U64, new PosGen<W8>(f, GenSection.LcpBase), new PosGen<W5>(f, GenSection.LcpPairsOfPowerOfA, long.MaxValue), isChecked: true);
        var cells = new ChildCells<PosGen<W5>>(new PosGen<W5>(f, GenSection.ChildPacked), lcp, tagged: false);
        Assert.That(cells.Length, Is.EqualTo(f.Length + 1));
        long[] ranks = Ranks(f.Length);
        foreach (long lb in ranks)
        {
            Assert.That(cells.FirstLIndex(lb, 0, true), Is.EqualTo(f.PackedCell(lb)), $"{f} cell {lb} (left bound)");
            Assert.That(cells.FirstLIndex(0, lb, false), Is.EqualTo(f.PackedCell(lb)), $"{f} cell {lb} (right bound)");
        }
        Assert.That(() => cells.FirstLIndex(f.Length + 1, 0, true), Throws.TypeOf<InvalidDataException>());
    }

    /// <summary>The pattern descent of <see cref="MappedSearch"/> over the generated views, with 64-bit ranks: the interval of a factor of the text, or empty.</summary>
    private static RankRange Descend<TPos>(ClosedFormText f, SaView<TPos> sa, ChildCells<TPos> cells, CompactLcp lcp, long start, long length)
        where TPos : struct, IPosReader
    {
        long n = f.Length;
        long lb = 0, rb = n, ell = 0, c = 0;
        bool isLast = true;
        while (c < length)
        {
            int symbol = f.CodeAt(start + c);
            long a = lb;
            long boundary = cells.FirstLIndex(lb, rb, isLast);
            long b;
            while (true)
            {
                b = boundary >= 0 ? boundary - 1 : rb;
                long p = sa[a] + ell;
                int key = p < n ? f.CodeAt(p) : -1;
                if (key == symbol) break;
                if (key > symbol || boundary < 0) return RankRange.Empty;
                a = boundary;
                boundary = cells.NextLIndex(a, rb);
            }
            bool childIsLast = boundary < 0;
            if (a == b)
            {
                long suffix = sa[a] + c;
                if (length - c > n - suffix) return RankRange.Empty;
                for (long i = c; i < length; i++)
                {
                    if (f.CodeAt(suffix + (i - c)) != f.CodeAt(start + i)) return RankRange.Empty;
                }
                return new RankRange(a, a);
            }
            long childEll = lcp[cells.FirstLIndex(a, b, childIsLast)];
            long limit = Math.Min(childEll, length);
            long from = sa[a] + c;
            for (long i = c + 1; i < limit; i++)
            {
                if (f.CodeAt(from + (i - c)) != f.CodeAt(start + i)) return RankRange.Empty;
            }
            c = limit;
            lb = a;
            rb = b;
            ell = childEll;
            isLast = childIsLast;
        }
        return new RankRange(lb, rb);
    }

    [Test]
    public void TheDescent_OverTheGeneratedViews_FindsTheFactorsOfTheText([ValueSource(nameof(Huge))] ClosedFormText f)
    {
        CompactLcp lcp = CompactLcp.FromReaders(LcpKind.U64, new PosGen<W8>(f, GenSection.LcpBase), new PosGen<W5>(f, GenSection.LcpPairsOfPowerOfA, long.MaxValue), isChecked: true);
        var sa = new SaView<PosGen<W5>>(new PosGen<W5>(f, GenSection.Sa));
        var cells = new ChildCells<PosGen<W5>>(new PosGen<W5>(f, GenSection.ChildPacked), lcp, tagged: false);
        foreach (long start in new long[] { 0, 1, 2, TwoTo31 - 7, TwoTo31, TwoTo32 + 1, f.Length / 2, f.Length - 600 })
        {
            foreach (long length in new long[] { 1, 2, 3, 5, 17, 100, 600 })
            {
                if (start + length > f.Length) continue;
                RankRange expected = f.RangeOfFactor(start, length);
                RankRange actual = Descend(f, sa, cells, lcp, start, length);
                Assert.That(actual.SameAs(expected), Is.True, $"{f} factor [{start}, +{length}): the views found {actual}, the closed form says {expected}");
            }
        }
    }

    [Test]
    public void TheRuns_OfPowerOfA_HaveTheirBoundsInTheViews([Values(1L, 2L, 10_000L, 2147483647L, 2147483648L, 2147483657L, 4294967295L, 4294967296L)] long k, [Values(2147483658L, 4294967306L, 8589934592L)] long n)
    {
        if (k > n) return;
        var f = new PowerOfA(n);
        CompactLcp lcp = CompactLcp.FromReaders(LcpKind.U64, new PosGen<W8>(f, GenSection.LcpBase), new PosGen<W5>(f, GenSection.LcpPairsOfPowerOfA, long.MaxValue), isChecked: true);
        var sa = new SaView<PosGen<W5>>(new PosGen<W5>(f, GenSection.Sa));
        RankRange run = f.RangeOfRun(k);
        Assert.That(run.Count, Is.EqualTo(n - k + 1), $"the count of a^{k} in a^{n}");
        // the interval [first..last] of the pattern a^k: the suffix at the first rank has k symbols or more, the one before has fewer; no value of the LCP inside is below k, the one at the left border is
        Assert.That(sa[run.First], Is.EqualTo(n - run.First));
        Assert.That(n - sa[run.First], Is.GreaterThanOrEqualTo(k));
        Assert.That(run.Last, Is.EqualTo(n));
        Assert.That(sa[run.Last], Is.EqualTo(0));
        if (run.First > 0) Assert.That(lcp[run.First], Is.LessThan(k));
        if (run.First + 1 <= n) Assert.That(lcp[run.First + 1], Is.GreaterThanOrEqualTo(k));
        Assert.That(lcp[n], Is.EqualTo(n - 1));
    }

    // ---------------------------------------------------------------------------------------------- the pyramid and the text

    [Test]
    public void ThePyramid_OverAnLcpArrayOf2To31Plus10Entries_AnswersAtTheBoundary()
    {
        // the structure of an LCP array is what the pyramid sees: a hash of the rank below 50 (a text's LCP at this size would hold a value above int.MaxValue, which the 32-bit levels refuse: see the bridge)
        var f = new PowerOfA(TwoTo31 + 10);
        CompactLcp lcp = CompactLcp.FromReaders(LcpKind.U32, new PosGen<W4>(f, GenSection.Synthetic), new PosGen<W4>(f, GenSection.LcpPairsOfPowerOfA, long.MaxValue), isChecked: false);
        using ILcpRanges pyramid = LcpPyramid.Create(lcp, 32);
        Assert.That(pyramid.Length, Is.EqualTo(TwoTo31 + 11));
        foreach (long r in new long[] { 0, 1, 31, 32, TwoTo31 - 1, TwoTo31, TwoTo31 + 1, TwoTo31 + 10, 1L << 30, 3L << 29 })
        {
            long value = (long)(PosGen<W4>.Hash(r) % 50);
            Assert.That(pyramid.ValueAt(r), Is.EqualTo(value), $"ValueAt({r})");
            foreach (long bound in new long[] { 1, 3, 25, 49, 50, 51 })
            {
                long prev = r;
                while (prev >= 0 && (long)(PosGen<W4>.Hash(prev) % 50) >= bound) prev--;
                long next = r;
                while (next <= TwoTo31 + 10 && (long)(PosGen<W4>.Hash(next) % 50) >= bound) next++;
                Assert.That(pyramid.PrevLess(r, bound), Is.EqualTo(prev), $"PrevLess({r}, {bound})");
                Assert.That(pyramid.NextLess(r, bound), Is.EqualTo(next), $"NextLess({r}, {bound})");
            }
        }
        long lo = TwoTo31 - 3000, hi = TwoTo31 + 10;
        long min = long.MaxValue;
        for (long j = lo; j <= hi; j++) min = Math.Min(min, (long)(PosGen<W4>.Hash(j) % 50));
        Assert.That(pyramid.Min(lo, hi), Is.EqualTo(min));
        Assert.That(pyramid.Min(0, TwoTo31 + 10), Is.EqualTo(0));
        Assert.That(() => pyramid.ValueAt(TwoTo31 + 11), Throws.TypeOf<ArgumentOutOfRangeException>());
    }

    [Test]
    public void ThePyramid_HoldsAnLcpValueAboveInt32_InItsLongLevels_AtHugeSize()
    {
        // a^n at 2^31 + 70: the last block (ranks 2^31 + 64 and up) has the minimum 2^31 + 63, above int.MaxValue: the 32-bit levels of S2.3 refused it by name, the 64-bit levels of S2.5 hold it
        var f = new PowerOfA(TwoTo31 + 70);
        CompactLcp lcp = CompactLcp.FromReaders(LcpKind.U64, new PosGen<W8>(f, GenSection.LcpBase), new PosGen<W5>(f, GenSection.LcpPairsOfPowerOfA, long.MaxValue), isChecked: false);
        using ILcpRanges pyramid = LcpPyramid.Create(lcp, 32);
        Assert.That(pyramid.Min(TwoTo31 + 64, TwoTo31 + 70), Is.EqualTo(TwoTo31 + 63));
        Assert.That(pyramid.Min(0, TwoTo31 + 70), Is.EqualTo(0));
        Assert.That(pyramid.PrevLess(TwoTo31 + 70, TwoTo31 + 10), Is.EqualTo(TwoTo31 + 10), "a^n: LCP[r] = r - 1 is below r - 1 + 1 at r");
        Assert.That(pyramid.NextLess(TwoTo31 + 64, TwoTo31 + 64), Is.EqualTo(TwoTo31 + 64), "LCP[r] = r - 1 is below r");
    }

    [TestCase(TwoTo31 + 10)]
    [TestCase(TwoTo32 + 10)]
    [TestCase(TwoTo33)]
    public void TheText_AtN_IsReadAt64BitPositions(long n)
    {
        const int bits = 1;
        long total = PackedTextFormat.TotalBytes(n, bits);
        using var section = new VirtualSection(total);
        long[] marks = [0, 1, TwoTo31 - 1, TwoTo31, TwoTo31 + 3, TwoTo32 - 2, TwoTo32 + 7, n - 1];
        section.Open(0, 4096);
        section.Open(total - 4096, 4096);
        foreach (long m in marks)
        {
            if (m >= n) continue;
            section.Open(m / 8, 8);
            section.Data[m / 8] |= (byte)(0x80 >> (int)(m % 8));
        }
        EsaText text = EsaText.FromPacked(MappedTextReader.FromMemory(section.Data, total, n, bits), ['a', 'b'], n);
        Assert.That(text.Length, Is.EqualTo(n));
        foreach (long m in marks)
        {
            if (m >= n) continue;
            Assert.That(text[m], Is.EqualTo('b'), $"position {m}");
            Assert.That(text.AsIndexText().SymbolAt(m), Is.EqualTo((int)'b'));
        }
        Assert.That(text[2], Is.EqualTo('a'));
        Assert.That(text.AsIndexText().SymbolAt(n), Is.EqualTo(-1));
        var block = new char[16];
        text.CopyTo(n - 16, block);
        Assert.That(block[15], Is.EqualTo('b'));
        Assert.That(text.AsIndexText().CopyTo(n - 5, block), Is.EqualTo(5));
        ITextSource source = text.AsTextSource();
        Assert.That(() => source.Length, Throws.TypeOf<NotSupportedException>());
        Assert.That(() => source[-1], Throws.TypeOf<IndexOutOfRangeException>(), "a caller's negative index is an index error whatever the length of the text");
        Assert.That(() => source[int.MinValue], Throws.TypeOf<IndexOutOfRangeException>());
        Assert.That(source[2], Is.EqualTo('a'));
    }

    // ---------------------------------------------------------------------------------------------- the 32-bit members over sections of more than 2^31 entries (the checked bridges)

    [Test]
    public void TheSuffixArrayView_OfMoreThan2To31FourByteEntries_CopiesAtAStartAbove2To31_AndTheIntegerReadsAreTheBridge()
    {
        long count = TwoTo31 + 100;
        long guard = TwoTo33; // address space in front of the section: where a start narrowed to int would read
        using var section = new VirtualSection(guard + (4 * count));
        byte* sa = section.Data + guard;
        long start = TwoTo31 + 5;
        long narrowed = guard + (4 * (long)(int)start); // the byte that the pointer of a copy reaches if the start is narrowed
        section.Open(guard, 4096);
        section.Open(guard + (4 * start), 64);
        section.Open(narrowed, 64);
        var words = (uint*)sa;
        words[50] = (uint)(TwoTo31 + 9);
        words[51] = 12;
        for (int i = 0; i < 4; i++)
        {
            words[start + i] = (uint)(7 + i);
            ((uint*)(section.Data + narrowed))[i] = 0xEEEE_0000u + (uint)i;
        }
        foreach (bool isChecked in new[] { true, false })
        {
            SaView view = SaView.FromMapped(sa, count, 4, isChecked);
            Assert.That(view.Length, Is.EqualTo(count));
            var ints = new int[4];
            var longs = new long[4];
            view.CopyTo(start, ints);
            view.CopyTo(start, longs);
            Assert.That(ints, Is.EqualTo(new[] { 7, 8, 9, 10 }));
            Assert.That(longs, Is.EqualTo(new long[] { 7, 8, 9, 10 }));
            Assert.That(view[start], Is.EqualTo(7));
            // an entry above int.MaxValue is a position of a text that the 32-bit engine cannot hold: the long members read it, the 32-bit members refuse it
            Assert.That(view[50L], Is.EqualTo(TwoTo31 + 9));
            Assert.That(() => view[50], Throws.TypeOf<InvalidDataException>().With.Message.Contains((TwoTo31 + 9).ToString()));
            Assert.That(() => view.CopyTo(50L, new int[1]), Throws.TypeOf<InvalidDataException>());
            Assert.That(() => view.ToArray(49, 3), Throws.TypeOf<InvalidDataException>());
            Assert.That(view[51], Is.EqualTo(12));
            Assert.That(view.ToArray(51, 1), Is.EqualTo(new[] { 12 }));
            Assert.That(view.ToLongArray(49, 3), Is.EqualTo(new long[] { 0, TwoTo31 + 9, 12 }));
        }
    }

    [Test]
    public void TheChildTable_OfMoreThan2To31Cells_TheIntegerMembersAreTheBridge_OverFourAndFiveByteCells([Values(4, 5)] int width)
    {
        long count = TwoTo31 + 11;
        long target = TwoTo31 + 7;
        using var lcpSection = new VirtualSection(count);
        lcpSection.Open(0, 4096);
        lcpSection.Open(target, 64);
        lcpSection.Data[100] = 3;     // LCP[100] = 3, LCP[101] = 0: the next l-index of 100 is the one that the cell says
        lcpSection.Data[target] = 3;  // LCP[target] = 3: the walk accepts it
        CompactLcp lcp = CompactLcp.FromMappedWide(LcpKind.U8, count, lcpSection.Data, null, 0, width, isChecked: true);
        using var cells = new VirtualSection((width * count) + 3);
        cells.Open(0, 4096);
        // cell 100 is the rank above int.MaxValue, cell 3 as well, cell 4 a small one
        for (int i = 0; i < width; i++)
        {
            cells.Data[(width * 100L) + i] = (byte)(target >> (8 * i));
            cells.Data[(width * 3L) + i] = (byte)(target >> (8 * i));
        }
        cells.Data[width * 4L] = 17;
        MappedChild child = MappedChild.FromMappedWide(cells.Data, count, width, lcp, isChecked: true, tagged: false);
        Assert.That(child.Length, Is.EqualTo(count));
        Assert.That(child.NextLIndex(100L, TwoTo33), Is.EqualTo(target));
        Assert.That(child.FirstLIndex(3L, 3L, true), Is.EqualTo(target));
        Assert.That(() => child.NextLIndex(100, int.MaxValue), Throws.TypeOf<InvalidDataException>().With.Message.Contains(target.ToString()));
        Assert.That(() => child.FirstLIndex(3, 3, true), Throws.TypeOf<InvalidDataException>().With.Message.Contains(target.ToString()));
        Assert.That(child.FirstLIndex(4, 4, true), Is.EqualTo(17));
        Assert.That(child.FirstLIndex(4L, 4L, true), Is.EqualTo(17));
    }

    [Test]
    public void ThePyramid_OverMoreThan2To31Entries_FindsThePreviousSmallerThroughABlockThatStartsAbove2To31()
    {
        // an answer in an earlier block that starts at 2^31 or above: the start of the block is a position, so it is computed in long
        long n = TwoTo31 + 300;
        var f = new PowerOfA(n);
        CompactLcp lcp = CompactLcp.FromReaders(LcpKind.U32, new PosGen<W4>(f, GenSection.Synthetic), new PosGen<W4>(f, GenSection.LcpPairsOfPowerOfA, long.MaxValue), isChecked: false);
        using ILcpRanges pyramid = LcpPyramid.Create(lcp, 32);
        int throughALevel = 0;
        foreach (long bound in new long[] { 1, 2, 3 })
        {
            for (long r = TwoTo31 + 32; r <= n; r++)
            {
                long prev = r;
                while (prev >= 0 && (long)(PosGen<W4>.Hash(prev) % 50) >= bound) prev--;
                Assert.That(pyramid.PrevLess(r, bound), Is.EqualTo(prev), $"PrevLess({r}, {bound})");
                if (prev >= TwoTo31 && (prev >> 5) != (r >> 5)) throughALevel++;
            }
        }
        Assert.That(throughALevel, Is.GreaterThan(50), "the queries that the test is about: the answer lies in an earlier block that starts at 2^31 or above");
        // the 32-bit members: the answer 'no entry' of NextLess is the length of the array, which is no int: a size that the 32-bit engine cannot hold, not a damaged file
        Assert.That(() => pyramid.NextLess(5, 0), Throws.TypeOf<NotSupportedException>().With.Message.Contains("NextLess"));
        Assert.That(() => pyramid.Borders(5, 6, 0, out _, out _, out _, out _), Throws.TypeOf<NotSupportedException>().With.Message.Contains("Borders"));
        Assert.That(pyramid.PrevLess(5, 0), Is.EqualTo(-1));
    }

    // ---------------------------------------------------------------------------------------------- the sparse shells: pointers over a mapping

    private static readonly long[] EndRanks = [0, 1, 2, 17, 255, 256, 1000, 4094, 4095];

    [Test]
    public void TheViewsOfAContainerOfIntMaxValueSymbols_AreTheWideForms_AndTheLongEngineServesThem()
    {
        // 2^31 - 1 symbols: n + 1 = 2^31 entries, one more than a 32-bit count holds: the boundary of the 32-bit forms
        long n = int.MaxValue;
        ClosedFormText f = new PowerOfA(n);
        SparseSpec spec = SparseSpec.For(f);
        string path = WriteContainer("viewsmax", spec, SparseContainer.EndPatches(f, SparseContainer.Plan(spec), 4096));
        using SesaFile file = SesaFile.Open(path, new EsaOpenOptions { Validation = EsaValidationLevel.Header });
        Assert.That(file.Sa.Length, Is.EqualTo(TwoTo31));
        Assert.That(file.Lcp.Length, Is.EqualTo(TwoTo31));
        Assert.That(file.Lcp.IsWide, Is.True);
        Assert.That(file.Child.Mapped.Length, Is.EqualTo(TwoTo31));
        foreach (long r in new long[] { 0, 1, n - 1, n }) Assert.That(file.Sa[r], Is.EqualTo(f.SuffixAt(r)), $"SA[{r}]");
        Assert.That(file.Lcp[n], Is.EqualTo(n - 1));
        Assert.That(EsaIndex.FromFile(file).TextLength, Is.EqualTo(n), "the long engine of S2.4 serves the text that the 32-bit engine refused at S2.3");
    }

    [Test]
    public void TheViewsOfAnOpenedContainerOf2To31Plus10Symbols_ReadTheEndsOfTheText()
    {
        long n = TwoTo31 + 10;
        ClosedFormText f = new PowerOfA(n);
        SparseSpec spec = SparseSpec.For(f);
        string path = WriteContainer("views31", spec, SparseContainer.EndPatches(f, SparseContainer.Plan(spec), 4096));
        using SesaFile file = SesaFile.Open(path, new EsaOpenOptions { Validation = EsaValidationLevel.Header, Access = EsaAccessMode.Checked });
        Assert.That(file.Sa.Length, Is.EqualTo(n + 1));
        Assert.That(file.Lcp.Length, Is.EqualTo(n + 1));
        Assert.That(file.Child.Mapped.Length, Is.EqualTo(n + 1));
        Assert.That(file.TextView.Length, Is.EqualTo(n));
        Assert.That(file.Lcp.Kind, Is.EqualTo(LcpKind.U32));
        Assert.That(file.Lcp.IsMapped, Is.True);
        foreach (long r in EndRanks.Concat(EndRanks.Select(e => n - e)))
        {
            Assert.That(file.Sa[r], Is.EqualTo(f.SuffixAt(r)), $"SA[{r}]");
            Assert.That(file.Sa.At(r), Is.EqualTo(f.SuffixAt(r)), $"SA.At({r})");
            Assert.That(file.Lcp[r], Is.EqualTo(f.Lcp(r)), $"LCP[{r}]");
            Assert.That(file.Child.Mapped.FirstLIndex(r, 0, true), Is.EqualTo(f.PackedCell(r)), $"CHILD[{r}] (left)");
            Assert.That(file.Child.Mapped.FirstLIndex(0, r, false), Is.EqualTo(f.PackedCell(r)), $"CHILD[{r}] (right)");
            if (r < n) Assert.That(file.TextView[r], Is.EqualTo('a'), $"TEXT[{r}]");
        }
        Assert.That(file.Child.FirstLIndex(n - 1, 0, true), Is.EqualTo(f.PackedCell(n - 1)));
        Assert.That(() => file.Sa.At(n + 1), Throws.TypeOf<InvalidDataException>());
        Assert.That(() => file.Lcp[n + 1], Throws.TypeOf<InvalidDataException>());
        Assert.That(() => file.Child.Mapped.FirstLIndex(n + 1, 0, true), Throws.TypeOf<InvalidDataException>());
        Assert.That(() => file.TextView[n], Throws.TypeOf<InvalidDataException>());
        // the int surface says what it cannot hold, it does not read a null
        Assert.That(file.Statistics.Length, Is.EqualTo(n), "the statistics of the header are 64-bit since S2.4");
        SesaSection saSection = file.Parsed.Sa;
        Assert.That(saSection.Length, Is.GreaterThan(int.MaxValue));
        Assert.That(() => file.Payload(saSection), Throws.TypeOf<NotSupportedException>().With.Message.Contains("blocks"), "a section above 2^31 bytes is never one span (risk R8)");
        Assert.That(file.Payload(file.Parsed.Alphabet).Length, Is.EqualTo(file.Parsed.Alphabet.Length));
        Assert.That(EsaIndex.FromFile(file).TextLength, Is.EqualTo(n), "the long engine of S2.4 serves the text (the 32-bit engine of S2.3 refused it by name)");
    }

    /// <summary>A sparse shell mapped read-only without the header reader (which refuses a width of 5 until S2.7): the pointers of its sections, from the layout that the writer made.</summary>
    private sealed unsafe class MappedShell : IDisposable
    {
        private readonly MemoryMappedFile _map;
        private readonly MemoryMappedViewAccessor _view;
        private readonly byte* _base;

        internal MappedShell(string path, SparseLayout layout)
        {
            Layout = layout;
            _map = MemoryMappedFile.CreateFromFile(path, FileMode.Open, null, 0, MemoryMappedFileAccess.Read);
            _view = _map.CreateViewAccessor(0, 0, MemoryMappedFileAccess.Read);
            byte* p = null;
            _view.SafeMemoryMappedViewHandle.AcquirePointer(ref p);
            _base = p + _view.PointerOffset;
        }

        internal SparseLayout Layout { get; }

        internal byte* Section(int kind) => _base + Layout.Section(kind).Offset;

        internal long Count(int kind) => Layout.Section(kind).Count;

        public void Dispose()
        {
            _view.SafeMemoryMappedViewHandle.ReleasePointer();
            _view.Dispose();
            _map.Dispose();
        }
    }

    [Test]
    public void TheViewsOverAFiveByteShellOf2To32Plus10Symbols_ReadTheEndsOfTheText([Values(true, false)] bool isChecked)
    {
        long n = TwoTo32 + 10;
        ClosedFormText f = new PowerOfA(n);
        SparseSpec spec = SparseSpec.For(f);
        SparseLayout layout = SparseContainer.Plan(spec);
        Assert.That((layout.PosWidth, layout.LcpCodec), Is.EqualTo((5, 4)), "width 5 and an LCP of 64 bits");
        string path = WriteContainer("views32", spec, SparseContainer.EndPatches(f, layout, 4096));
        using var shell = new MappedShell(path, layout);
        SaView sa = SaView.FromMapped(shell.Section(SesaFormat.KindSa), shell.Count(SesaFormat.KindSa), 5, isChecked);
        CompactLcp lcp = CompactLcp.FromMappedWide(LcpKind.U64, n + 1, shell.Section(SesaFormat.KindLcp), null, 0, 5, isChecked);
        MappedChild child = MappedChild.FromMappedWide(shell.Section(SesaFormat.KindChild), n + 1, 5, lcp, isChecked, tagged: false);
        Assert.That(sa.Length, Is.EqualTo(n + 1));
        Assert.That(sa.StoredWidth, Is.EqualTo(5));
        foreach (long r in EndRanks.Concat(EndRanks.Select(e => n - e)))
        {
            Assert.That(sa[r], Is.EqualTo(f.SuffixAt(r)), $"SA[{r}]");
            Assert.That(lcp[r], Is.EqualTo(f.Lcp(r)), $"LCP[{r}]");
            Assert.That(child.FirstLIndex(r, 0, true), Is.EqualTo(f.PackedCell(r)), $"CHILD[{r}]");
        }
        Assert.That(lcp[n], Is.EqualTo(n - 1), "a depth above 2^32");
        Assert.That(sa.ToLongArray(n - 10, 11), Is.EqualTo(Enumerable.Range(0, 11).Select(i => f.SuffixAt(n - 10 + i))));
        Assert.That(() => sa.At(n + 1), Throws.TypeOf<InvalidDataException>());
        Assert.That(sa[TwoTo31], Is.EqualTo(0), "a hole of the shell reads as zero");
        if (isChecked)
        {
            Assert.That(() => sa[n + 1], Throws.TypeOf<InvalidDataException>());
            Assert.That(() => child.FirstLIndex(n + 1, 0, true), Throws.TypeOf<InvalidDataException>());
        }
    }

    [Test]
    public void TheSectionsLevel_OfTheLargestFourByteContainer_ReadsEveryByteOfThreeSectionsOf17Gb()
    {
        long n = TwoTo32 - 4;
        ClosedFormText f = new PowerOfA(n);
        SparseSpec spec = SparseSpec.For(f);
        Assert.That(SparseContainer.Plan(spec).PosWidth, Is.EqualTo(4));
        string path = WriteContainer("sections32", spec, SparseContainer.EndPatches(f, SparseContainer.Plan(spec), 4096));
        var clock = System.Diagnostics.Stopwatch.StartNew();
        using SesaFile file = SesaFile.Open(path, new EsaOpenOptions { Validation = EsaValidationLevel.Sections });
        double seconds = clock.Elapsed.TotalSeconds;
        TestContext.Out.WriteLine($"OPEN level=Sections n={n} seconds={seconds:F2} section_bytes={file.Report.SectionBytesChecksummed}");
        Assert.That(file.Report.SectionBytesChecksummed, Is.GreaterThan(3 * 4 * (TwoTo32 - 4)));
        Assert.That(file.Sa.Length, Is.EqualTo(n + 1));
        Assert.That(file.Sa[n], Is.EqualTo(0));
        Assert.That(seconds, Is.LessThan(40), "the checksum pass reads 52 GB of holes at 3.3 - 4.6 GB/s");
    }
}
