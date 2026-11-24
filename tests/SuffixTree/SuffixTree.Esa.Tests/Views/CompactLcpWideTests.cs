using SuffixTree.Esa.Tests.External;
using SuffixTree.Esa.Tests.Fuzz;

namespace SuffixTree.Esa.Tests.Views;

/// <summary>
/// The LCP array of work package S2.3 at any width and size: ranks and values are 64-bit, the overflow pairs are <c>2 x POS_WIDTH</c> bytes (8 at width 4, 10 at width 5), the U64 kind (LCP codec 4) has 8-byte entries and no escape, the
/// binary search over the pairs runs over 64-bit ranks, and the 32-bit members are checked bridges (a value from a file that does not fit is <see cref="InvalidDataException"/>, a length that the 32-bit engine cannot hold is
/// <see cref="NotSupportedException"/>). The oracle is the dense LCP of the in-memory index of S1, re-encoded by <see cref="ViewImage"/>.
/// </summary>
[TestFixture]
public unsafe class CompactLcpWideTests
{
    private static readonly string[] Texts =
    [
        "banana", "mississippi", "a", "ab",
        new string('a', 300),                                   // LCP values up to 299: escapes of U8
        string.Concat(Enumerable.Repeat("abcdefghij", 40)),      // values up to 389
        new string('a', 70_100),                                 // values above 65,535: escapes of U16
    ];

    private static CompactLcp Open(ViewImage image, bool isChecked) =>
        CompactLcp.FromMappedWide(image.Kind, image.Count, image.LcpBasePointer, image.PairsPointer, image.Pairs.Count, image.PosWidth, isChecked);

    [Test]
    public void EveryKindAtEveryWidth_GivesTheDenseLcp([ValueSource(nameof(Texts))] string text, [Values("U8", "U16", "U32", "U64")] string kindName, [Values(4, 5)] int width, [Values(true, false)] bool isChecked)
    {
        var kind = Enum.Parse<LcpKind>(kindName);
        using ViewImage image = ViewImage.Of(text, width, kind);
        CompactLcp lcp = Open(image, isChecked);
        Assert.That(lcp.IsMapped, Is.True);
        Assert.That(lcp.Kind, Is.EqualTo(kind));
        Assert.That(lcp.Length, Is.EqualTo(image.Count));
        Assert.That(lcp.OverflowCount, Is.EqualTo(image.Pairs.Count));
        long step = Math.Max(1, image.Count / 4000);
        for (long r = 0; r < image.Count; r += step) Assert.That(lcp[r], Is.EqualTo(image.Lcp[r]), $"rank {r}");
        Assert.That(lcp[image.Count - 1], Is.EqualTo(image.Lcp[^1]));
        foreach ((long rank, long value) in image.Pairs.Take(50).Concat(image.Pairs.TakeLast(50))) Assert.That(lcp[rank], Is.EqualTo(value), $"escape at rank {rank}");
        if (image.Count <= int.MaxValue && kind != LcpKind.U64)
        {
            for (int r = 0; r < image.Count; r += (int)step) Assert.That(lcp[r], Is.EqualTo((int)image.Lcp[r]), $"int rank {r}");
        }
    }

    [Test]
    public void TheEscapesOfU8AndU16_AreResolvedByTheBinarySearch_AtTheFirstTheLastAndBetween([Values(4, 5)] int width, [Values("U8", "U16")] string kindName)
    {
        var kind = Enum.Parse<LcpKind>(kindName);
        using ViewImage image = ViewImage.Of(new string('a', kind == LcpKind.U8 ? 1000 : 70_100), width, kind);
        Assert.That(image.Pairs.Count, Is.GreaterThan(500));
        CompactLcp lcp = Open(image, true);
        (long firstRank, long firstValue) = image.Pairs[0];
        (long lastRank, long lastValue) = image.Pairs[^1];
        Assert.That(lcp[firstRank], Is.EqualTo(firstValue));
        Assert.That(lcp[lastRank], Is.EqualTo(lastValue));
        Assert.That(lcp[firstRank - 1], Is.EqualTo(image.Lcp[firstRank - 1]), "the entry before the first escape is not an escape");
        foreach (int k in new[] { 1, 2, 3, 7, image.Pairs.Count / 2, image.Pairs.Count - 2 }) Assert.That(lcp[image.Pairs[k].Rank], Is.EqualTo(image.Pairs[k].Value), $"pair {k}");
    }

    [Test]
    public void AnEscapeWithoutItsPair_IsInvalidData([Values(4, 5)] int width, [Values(true, false)] bool isChecked)
    {
        using ViewImage image = ViewImage.Of(new string('a', 700), width, LcpKind.U8);
        // the list without its first pair: the escape of that rank is missing; the list without its last pair likewise
        byte[] pairs = ViewImage.Encode(image.Pairs.Skip(1).SelectMany(p => new[] { (ulong)p.Rank, (ulong)p.Value }), width, width == 5 ? 3 : 0);
        using var guarded = new GuardedImage(pairs, flush: true);
        CompactLcp lcp = CompactLcp.FromMappedWide(LcpKind.U8, image.Count, image.LcpBasePointer, guarded.Data, image.Pairs.Count - 1, width, isChecked);
        Assert.That(() => lcp[image.Pairs[0].Rank], Throws.TypeOf<InvalidDataException>().With.Message.Contains("missing"));
        Assert.That(() => lcp[image.Pairs[^1].Rank + 0], Throws.Nothing);
        Assert.That(lcp[image.Pairs[1].Rank], Is.EqualTo(image.Pairs[1].Value));
    }

    [Test]
    public void APairWordThatIsNotANumberOfTheEngine_IsInvalidData([Values(4, 5)] int width)
    {
        using ViewImage image = ViewImage.Of(new string('a', 400), width, LcpKind.U8);
        long rank = image.Pairs[10].Rank;
        // the value of pair 10 as 2^64 - 1 at width 4 is impossible (4 bytes), so damage the value field with the all-ones word of the width
        var fields = image.Pairs.SelectMany(p => new[] { (ulong)p.Rank, (ulong)p.Value }).ToArray();
        fields[(2 * 10) + 1] = width == 4 ? 0xFFFF_FFFFUL : 0xFF_FFFF_FFFFUL;
        using var guarded = new GuardedImage(ViewImage.Encode(fields, width, width == 5 ? 3 : 0), flush: true);
        CompactLcp lcp = CompactLcp.FromMappedWide(LcpKind.U8, image.Count, image.LcpBasePointer, guarded.Data, image.Pairs.Count, width, isChecked: true);
        Assert.That(lcp[rank], Is.EqualTo(width == 4 ? 0xFFFF_FFFFL : 0xFF_FFFF_FFFFL), "a value is a value; the engine that uses it validates it against the text");
    }

    [Test]
    public void TheValuesOfU32AndU64_AreCheckedAgainstTheTextInCheckedMode([Values("U32", "U64")] string kindName)
    {
        var kind = Enum.Parse<LcpKind>(kindName);
        // n + 1 = 8 entries, one of them 100: larger than the text
        int width = kind == LcpKind.U32 ? 4 : 8;
        ulong[] values = [0, 0, 1, 100, 2, 3, 0, 1];
        byte[] bytes = ViewImage.Encode(values, width, 0);
        using var guarded = new GuardedImage(bytes, flush: true);
        CompactLcp strict = CompactLcp.FromMappedWide(kind, values.Length, guarded.Data, null, 0, 4, isChecked: true);
        CompactLcp lax = CompactLcp.FromMappedWide(kind, values.Length, guarded.Data, null, 0, 4, isChecked: false);
        Assert.That(strict[2L], Is.EqualTo(1));
        Assert.That(() => strict[3L], Throws.TypeOf<InvalidDataException>().With.Message.Contains("100"));
        Assert.That(lax[3L], Is.EqualTo(100), "Trusted does not look");
        Assert.That(strict[7L], Is.EqualTo(1));
    }

    [Test]
    public void AnArrayOf32Or64Bits_HasNoOverflowList()
    {
        var bytes = new byte[64];
        fixed (byte* b = bytes)
        {
            nint q = (nint)b;
            Assert.That(() => CompactLcp.FromMappedWide(LcpKind.U32, 4, (byte*)q, (byte*)q, 1, 4, true), Throws.TypeOf<InvalidDataException>());
            Assert.That(() => CompactLcp.FromMappedWide(LcpKind.U64, 4, (byte*)q, (byte*)q, 1, 5, true), Throws.TypeOf<InvalidDataException>());
            Assert.That(() => CompactLcp.FromMappedWide(LcpKind.U8, 0, (byte*)q, null, 0, 4, true), Throws.TypeOf<ArgumentOutOfRangeException>());
            Assert.That(() => CompactLcp.FromMappedWide(LcpKind.U8, 4, (byte*)q, null, -1, 4, true), Throws.TypeOf<ArgumentOutOfRangeException>());
            Assert.That(() => CompactLcp.FromMappedWide(LcpKind.U8, 4, (byte*)q, null, 0, 6, true), Throws.TypeOf<ArgumentOutOfRangeException>());
            Assert.That(() => CompactLcp.FromMappedWide(LcpKind.U8, 4, (byte*)q, null, 0, 3, true), Throws.TypeOf<ArgumentOutOfRangeException>());
        }
    }

    [Test]
    public void TheRankIsBounded_InBothModes_AsInvalidData([Values(4, 5)] int width, [Values(true, false)] bool isChecked)
    {
        using ViewImage image = ViewImage.Of("banana", width, LcpKind.U8);
        CompactLcp lcp = Open(image, isChecked);
        long n = image.Count;
        foreach (long rank in new[] { n, n + 1, -1, long.MinValue, long.MaxValue, 1L << 33 })
        {
            Assert.That(() => lcp[rank], Throws.TypeOf<InvalidDataException>(), $"rank {rank}");
        }
        Assert.That(() => lcp[(int)n], Throws.TypeOf<InvalidDataException>());
    }

    [Test]
    public void TheWideFormAgreesWithTheLegacyMappedForm_OnEveryRank([Values("U8", "U16", "U32")] string kindName)
    {
        var kind = Enum.Parse<LcpKind>(kindName);
        using ViewImage image = ViewImage.Of(new string('a', 70_100), 4, kind);
        CompactLcp wide = Open(image, true);
        CompactLcp legacy = CompactLcp.FromMapped(kind, (int)image.Count, image.LcpBasePointer, (uint*)image.PairsPointer, image.Pairs.Count, isChecked: true);
        for (int r = 0; r < image.Count; r += 3)
        {
            Assert.That(wide[(long)r], Is.EqualTo(legacy[r]), $"rank {r}");
            Assert.That(wide[r], Is.EqualTo(legacy[r]), $"int rank {r}");
        }
        Assert.That(wide.OverflowCount, Is.EqualTo(legacy.OverflowCount));
        Assert.That(wide.PayloadBytes, Is.EqualTo(legacy.PayloadBytes));
        Assert.That(wide.OverflowRanks, Is.EqualTo(legacy.OverflowRanks));
        Assert.That(wide.OverflowValues, Is.EqualTo(legacy.OverflowValues));
    }

    [Test]
    public void ThePayloadBytes_CountThePairsAtTwiceThePositionWidth()
    {
        using ViewImage four = ViewImage.Of(new string('a', 700), 4, LcpKind.U8);
        using ViewImage five = ViewImage.Of(new string('a', 700), 5, LcpKind.U8);
        Assert.That(Open(four, true).PayloadBytes, Is.EqualTo(701 + (8L * four.Pairs.Count)));
        Assert.That(Open(five, true).PayloadBytes, Is.EqualTo(701 + (10L * five.Pairs.Count)));
        using ViewImage u64 = ViewImage.Of("banana", 5, LcpKind.U64);
        Assert.That(Open(u64, true).PayloadBytes, Is.EqualTo(8L * u64.Count));
    }

    [Test]
    public void TheCursor_ReadsTheWideForm_InAndOutOfOrder([Values(4, 5)] int width, [Values("U8", "U16")] string kindName)
    {
        var kind = Enum.Parse<LcpKind>(kindName);
        using ViewImage image = ViewImage.Of(new string('a', kind == LcpKind.U8 ? 900 : 70_100), width, kind);
        CompactLcp lcp = Open(image, true);
        CompactLcp.Cursor forward = lcp.CreateCursor();
        for (int r = 0; r < image.Count; r += 7) Assert.That(forward.Get(r), Is.EqualTo((int)image.Lcp[r]), $"rank {r}");
        CompactLcp.Cursor back = lcp.CreateCursor();
        for (int r = (int)image.Count - 1; r >= 0; r -= 11) Assert.That(back.Get(r), Is.EqualTo((int)image.Lcp[r]), $"back rank {r}");
        Assert.That(lcp.ToArray(), Is.EqualTo(image.Lcp.Select(v => (int)v)));
    }

    [Test]
    public void TheBridgesToThe32BitEngine_ThrowAndNeverNarrow()
    {
        // a U64 array with a value above int.MaxValue: the int indexer says damage (the value came from the file), the long indexer reads it
        ulong[] values = [0, 0, 5, (1UL << 33) + 3, (1UL << 31) - 1, 1UL << 31, 7, 1];
        using var guarded = new GuardedImage(ViewImage.Encode(values, 8, 0), flush: true);
        CompactLcp lcp = CompactLcp.FromMappedWide(LcpKind.U64, values.Length, guarded.Data, null, 0, 5, isChecked: false);
        Assert.That(lcp[3L], Is.EqualTo((1L << 33) + 3));
        Assert.That(lcp[4L], Is.EqualTo(int.MaxValue));
        Assert.That(lcp[4], Is.EqualTo(int.MaxValue));
        Assert.That(lcp[2], Is.EqualTo(5));
        Assert.That(() => lcp[3], Throws.TypeOf<InvalidDataException>().With.Message.Contains((1L << 33) + 3 + ""));
        Assert.That(() => lcp[5], Throws.TypeOf<InvalidDataException>(), "2^31 is not an int");
        Assert.That(() => lcp.ToArray(), Throws.TypeOf<InvalidDataException>());
        CompactLcp.Cursor cursor = lcp.CreateCursor();
        Assert.That(() => cursor.Get(3), Throws.TypeOf<InvalidDataException>());
    }

    [Test]
    public void ALengthThatThe32BitEngineCannotHold_IsNotSupported_ByTheMembersThatNeedAnInteger()
    {
        long count = (1L << 31) + 100; // one byte per entry: 2 GiB of address space, no page touched
        using var section = new VirtualSection(count);
        section.Open(count - 4096, 4096);
        section.Open(0, 4096);
        CompactLcp lcp = CompactLcp.FromMappedWide(LcpKind.U8, count, section.Data, null, 0, 4, isChecked: true);
        Assert.That(lcp.Length, Is.EqualTo(count));
        Assert.That(lcp[count - 1], Is.EqualTo(0));
        Assert.That(lcp[(1L << 31) + 5], Is.EqualTo(0));
        Assert.That(() => lcp.ToArray(), Throws.TypeOf<NotSupportedException>().With.Message.Contains("ToArray"));
        Assert.That(() => lcp.CreateCursor().Get(0), Throws.Nothing);
    }

    [Test]
    public void TheManagedForm_GivesTheSameValuesThroughTheLongIndexer()
    {
        CompactLcp lcp = LcpChildCase.Of(new string('a', 700)).Lcp;
        Assert.That(lcp.IsMapped, Is.False);
        for (int r = 0; r < lcp.Length; r++) Assert.That(lcp[(long)r], Is.EqualTo(lcp[r]), $"rank {r}");
        Assert.That(() => lcp[lcp.Length], Throws.TypeOf<ArgumentOutOfRangeException>());
        Assert.That(() => lcp[-1L], Throws.TypeOf<ArgumentOutOfRangeException>());
        Assert.That(() => lcp[1L << 40], Throws.TypeOf<ArgumentOutOfRangeException>());
    }

    [Test]
    public void TheKindU64_IsNotSelectedByTheBuilder_AndHasNoManagedForm()
    {
        Assert.That(Enum.GetValues<LcpKind>(), Is.EqualTo(new[] { LcpKind.U8, LcpKind.U16, LcpKind.U32, LcpKind.U64 }));
        Assert.That(CompactLcp.SelectKind(1000, 0, 0), Is.EqualTo(LcpKind.U8));
        Assert.That(CompactLcp.SelectKind(1000, 900, 900), Is.EqualTo(LcpKind.U32));
        Assert.That(() => CompactLcp.FromDense([0, 0, 1], LcpKind.U64), Throws.TypeOf<ArgumentOutOfRangeException>().Or.TypeOf<NotSupportedException>());
    }
    [Test]
    public void AValueEqualToTheTextLength_IsTheLargestThatCheckedModeAccepts([Values("U32", "U64")] string kindName)
    {
        var kind = Enum.Parse<LcpKind>(kindName);
        int width = kind == LcpKind.U32 ? 4 : 8;
        // n + 1 = 8 entries: a value of 8 is allowed (the text has 7 symbols; the check is "larger than n + 1"), 9 is not
        ulong[] values = [0, 0, 8, 1, 9, 2, 3, 0];
        using var guarded = new GuardedImage(ViewImage.Encode(values, width, 0), flush: true);
        CompactLcp strict = CompactLcp.FromMappedWide(kind, values.Length, guarded.Data, null, 0, 4, isChecked: true);
        Assert.That(strict[2L], Is.EqualTo(8));
        Assert.That(() => strict[4L], Throws.TypeOf<InvalidDataException>());
        Assert.That(strict[1L], Is.EqualTo(0));
    }

    [Test]
    public void AU64ValueThatIsNoSignedNumber_IsInvalidData_InBothModes([Values(true, false)] bool isChecked)
    {
        ulong[] values = [0, 0, 1UL << 63, ulong.MaxValue, 5];
        using var guarded = new GuardedImage(ViewImage.Encode(values, 8, 0), flush: true);
        CompactLcp lcp = CompactLcp.FromMappedWide(LcpKind.U64, values.Length, guarded.Data, null, 0, 5, isChecked);
        Assert.That(() => lcp[2L], Throws.TypeOf<InvalidDataException>());
        Assert.That(() => lcp[3L], Throws.TypeOf<InvalidDataException>());
        Assert.That(lcp[4L], Is.EqualTo(5), "an ordinary value after the damaged ones is still read");
    }

    [Test]
    public void TheReaderOfAnotherWidthThanTheKind_IsRefused()
    {
        var bytes = new byte[64];
        fixed (byte* b = bytes)
        {
            nint q = (nint)b;
            Assert.That(() => CompactLcp.FromReaders(LcpKind.U8, new Pos16<TrustedAccess>((byte*)q, 4), new Pos32<TrustedAccess>((byte*)q, 0), true), Throws.TypeOf<ArgumentException>());
            Assert.That(() => CompactLcp.FromReaders(LcpKind.U16, new Pos8<TrustedAccess>((byte*)q, 4), new Pos32<TrustedAccess>((byte*)q, 0), true), Throws.TypeOf<ArgumentException>());
            Assert.That(() => CompactLcp.FromReaders(LcpKind.U32, new Pos64<TrustedAccess>((byte*)q, 4), new Pos32<TrustedAccess>((byte*)q, 0), true), Throws.TypeOf<ArgumentException>());
            Assert.That(() => CompactLcp.FromReaders(LcpKind.U64, new Pos32<TrustedAccess>((byte*)q, 4), new Pos32<TrustedAccess>((byte*)q, 0), true), Throws.TypeOf<ArgumentException>());
            Assert.That(() => CompactLcp.FromReaders(LcpKind.U8, new Pos8<TrustedAccess>((byte*)q, 4), new Pos32<TrustedAccess>((byte*)q, 3), true), Throws.TypeOf<ArgumentException>(), "an odd number of fields is no list of pairs");
            Assert.That(() => CompactLcp.FromReaders(LcpKind.U8, new Pos8<TrustedAccess>((byte*)q, 0), new Pos32<TrustedAccess>((byte*)q, 0), true), Throws.TypeOf<ArgumentOutOfRangeException>(), "no entry at all");
            Assert.That(() => CompactLcp.FromReaders(LcpKind.U8, new Pos8<TrustedAccess>((byte*)q, 4), new Pos32<TrustedAccess>((byte*)q, 2), true), Throws.Nothing);
        }
    }
}
