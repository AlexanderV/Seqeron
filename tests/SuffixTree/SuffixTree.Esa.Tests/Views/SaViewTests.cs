using SuffixTree.Esa.Tests.External;
namespace SuffixTree.Esa.Tests.Views;

/// <summary>
/// The suffix-array views of work package S2.3: the generic <see cref="SaView{TPos}"/> over every reader and the non-generic <see cref="SaView"/> (what a container reader and the 32-bit engine hold) over an array, a 4-byte and a 5-byte
/// section, in both access modes: length, entry, <c>At</c>, bulk copies at a 64-bit start into spans of 32- and 64-bit positions, the bounds, and the values that do not fit a 32-bit position. The oracle is the suffix array of the
/// in-memory index of S1 (<see cref="LcpChildCase"/>) re-encoded by <see cref="ViewImage"/>.
/// </summary>
[TestFixture]
public unsafe class SaViewTests
{
    private static readonly string[] Texts = ["banana", "mississippi", "aaaaaaaaaaaaaaaaaaaaaaaaaaaa", "abracadabraabracadabra", "ACGTACGTTTGACCA", "a", "ab", string.Concat(Enumerable.Range(0, 400).Select(i => (char)('a' + ((i * i) % 5))))];

    // ---------------------------------------------------------------------------------------------- the generic view

    [Test]
    public void Generic_OverEveryReader_GivesTheSuffixArray([ValueSource(nameof(Texts))] string text, [Values(4, 5)] int width, [Values(true, false)] bool isChecked)
    {
        using ViewImage image = ViewImage.Of(text, width, LcpKind.U32);
        long n = image.Count;
        if (width == 4)
        {
            if (isChecked) Verify(new SaView<Pos32<CheckedAccess>>(new Pos32<CheckedAccess>(image.SaPointer, n)), image.Sa);
            else Verify(new SaView<Pos32<TrustedAccess>>(new Pos32<TrustedAccess>(image.SaPointer, n)), image.Sa);
            Verify(new SaView<PosInt>(new PosInt(image.Case.Sa)), image.Sa);
        }
        else if (isChecked) Verify(new SaView<Pos40<CheckedAccess>>(new Pos40<CheckedAccess>(image.SaPointer, n)), image.Sa);
        else Verify(new SaView<Pos40<TrustedAccess>>(new Pos40<TrustedAccess>(image.SaPointer, n)), image.Sa);
    }

    private static void Verify<TPos>(SaView<TPos> view, long[] expected)
        where TPos : struct, IPosReader
    {
        Assert.That(view.Length, Is.EqualTo(expected.Length));
        for (int i = 0; i < expected.Length; i++)
        {
            Assert.That(view[i], Is.EqualTo(expected[i]), $"rank {i}");
            Assert.That(view.At(i), Is.EqualTo(expected[i]), $"At({i})");
        }
        Assert.That(view.ToArray(0, expected.Length), Is.EqualTo(expected));
        for (int start = 0; start < expected.Length; start += 3)
        {
            for (int len = 0; len <= Math.Min(5, expected.Length - start); len++)
            {
                var longs = new long[len];
                var ints = new int[len];
                view.CopyTo(start, longs);
                view.CopyTo(start, ints);
                Assert.That(longs, Is.EqualTo(expected.Skip(start).Take(len)), $"CopyTo(long) {start}+{len}");
                Assert.That(ints, Is.EqualTo(expected.Skip(start).Take(len).Select(v => (int)v)), $"CopyTo(int) {start}+{len}");
            }
        }
        Assert.That(view.ToArray(expected.Length, 0), Is.Empty);
    }

    [Test]
    public void Generic_AtAndTheBulkCopies_AreAlwaysBounded_InvalidData([Values(true, false)] bool isChecked)
    {
        using ViewImage image = ViewImage.Of("mississippi", 4, LcpKind.U32);
        long n = image.Count;
        SaView<Pos32<TrustedAccess>> view = new(new Pos32<TrustedAccess>(image.SaPointer, n));
        _ = isChecked;
        foreach (long rank in new[] { n, n + 1, -1, long.MinValue, long.MaxValue, 1L << 32 })
        {
            Assert.That(() => view.At(rank), Throws.TypeOf<InvalidDataException>(), $"At({rank})");
        }
        Assert.That(() => view.CopyTo(n - 1, new long[2]), Throws.TypeOf<InvalidDataException>());
        Assert.That(() => view.CopyTo(n - 1, new int[2]), Throws.TypeOf<InvalidDataException>());
        Assert.That(() => view.CopyTo(-1, new long[1]), Throws.TypeOf<InvalidDataException>());
        Assert.That(() => view.CopyTo(n + 1, Array.Empty<long>()), Throws.TypeOf<InvalidDataException>());
        Assert.That(() => view.CopyTo(1L << 40, new long[1]), Throws.TypeOf<InvalidDataException>());
        Assert.That(() => view.ToArray(n, 1), Throws.TypeOf<InvalidDataException>());
        Assert.That(() => view.ToArray(0, -1), Throws.TypeOf<ArgumentOutOfRangeException>());
        Assert.That(view.At(n - 1), Is.EqualTo(image.Sa[^1]));
        Assert.That(() => view.CopyTo(n, Array.Empty<long>()), Throws.Nothing, "an empty copy at the end is a copy of nothing");
    }

    [Test]
    public void Generic_CheckedIndexer_IsBounded_TrustedIsNot()
    {
        using ViewImage image = ViewImage.Of("banana", 5, LcpKind.U32);
        long n = image.Count;
        var checkedView = new SaView<Pos40<CheckedAccess>>(new Pos40<CheckedAccess>(image.SaPointer, n));
        Assert.That(() => checkedView[n], Throws.TypeOf<InvalidDataException>());
        Assert.That(() => checkedView[-1], Throws.TypeOf<InvalidDataException>());
        Assert.That(checkedView[n - 1], Is.EqualTo(image.Sa[^1]));
    }

    [Test]
    public void Generic_CopyToIntegers_RefusesAPositionAboveTheIntegers_AsDamage()
    {
        // a 5-byte section whose entry 2 is 2^32 + 5 and entry 3 is int.MaxValue + 1: the 32-bit bridge does not narrow, it says the value
        ulong[] values = [7, 8, (1UL << 32) + 5, 1UL << 31, 3, (1UL << 31) - 1];
        using var guarded = new Fuzz.GuardedImage(ViewImage.Encode(values, 5, 3), flush: true);
        var view = new SaView<Pos40<CheckedAccess>>(new Pos40<CheckedAccess>(guarded.Data, values.Length));
        var longs = new long[values.Length];
        view.CopyTo(0, longs);
        Assert.That(longs, Is.EqualTo(values.Select(v => (long)v)));
        Assert.That(longs[2], Is.EqualTo((1L << 32) + 5));
        var small = new int[2];
        view.CopyTo(0, small);
        Assert.That(small, Is.EqualTo(new[] { 7, 8 }));
        view.CopyTo(4, small);
        Assert.That(small, Is.EqualTo(new[] { 3, int.MaxValue }));
        Assert.That(() => view.CopyTo(2, new int[2]), Throws.TypeOf<InvalidDataException>().With.Message.Contains("4294967301"));
        Assert.That(() => view.CopyTo(3, new int[1]), Throws.TypeOf<InvalidDataException>(), "2^31 is not an int");
    }

    // ---------------------------------------------------------------------------------------------- the non-generic view

    [Test]
    public void NonGeneric_OverEveryForm_GivesTheSuffixArray([ValueSource(nameof(Texts))] string text, [Values(4, 5)] int width, [Values(true, false)] bool isChecked)
    {
        using ViewImage image = ViewImage.Of(text, width, LcpKind.U8);
        SaView view = SaView.FromMapped(image.SaPointer, image.Count, width, isChecked);
        VerifyWide(view, image.Sa, width);
        if (width == 4)
        {
            VerifyWide(new SaView(image.Case.Sa), image.Sa, 4);
            VerifyWide(new SaView((int*)image.SaPointer, image.Count, isChecked), image.Sa, 4);
        }
    }

    private static void VerifyWide(SaView view, long[] expected, int width)
    {
        Assert.That(view.Length, Is.EqualTo(expected.Length));
        Assert.That(view.StoredWidth, Is.EqualTo(width));
        for (int i = 0; i < expected.Length; i++)
        {
            Assert.That(view[(long)i], Is.EqualTo(expected[i]), $"long rank {i}");
            Assert.That(view[i], Is.EqualTo((int)expected[i]), $"int rank {i}");
            Assert.That(view.At(i), Is.EqualTo(expected[i]), $"At({i})");
        }
        Assert.That(view.ToLongArray(0, expected.Length), Is.EqualTo(expected));
        Assert.That(view.ToArray(0, expected.Length), Is.EqualTo(expected.Select(v => (int)v)));
        for (int start = 0; start + 4 <= expected.Length; start += 4)
        {
            var longs = new long[4];
            var ints = new int[4];
            view.CopyTo((long)start, longs);
            view.CopyTo((long)start, ints);
            Assert.That(longs, Is.EqualTo(expected.Skip(start).Take(4)));
            Assert.That(ints, Is.EqualTo(expected.Skip(start).Take(4).Select(v => (int)v)));
        }
    }

    [Test]
    public void NonGeneric_TheBoundsOfTheTwoModes([Values(4, 5)] int width)
    {
        using ViewImage image = ViewImage.Of("mississippi", width, LcpKind.U8);
        long n = image.Count;
        SaView c = SaView.FromMapped(image.SaPointer, n, width, isChecked: true);
        Assert.That(() => c[n], Throws.TypeOf<InvalidDataException>());
        Assert.That(() => c[-1L], Throws.TypeOf<InvalidDataException>());
        Assert.That(() => c[(int)n], Throws.TypeOf<InvalidDataException>());
        Assert.That(() => c[-1], Throws.TypeOf<InvalidDataException>());
        Assert.That(() => c[1L << 40], Throws.TypeOf<InvalidDataException>());
        foreach (bool isChecked in new[] { true, false })
        {
            SaView v = SaView.FromMapped(image.SaPointer, n, width, isChecked);
            Assert.That(() => v.At(n), Throws.TypeOf<InvalidDataException>(), "At is bounded in both modes");
            Assert.That(() => v.At(-1), Throws.TypeOf<InvalidDataException>());
            Assert.That(() => v.CopyTo(n - 1, new long[2]), Throws.TypeOf<InvalidDataException>());
            Assert.That(() => v.CopyTo(n - 1, new int[2]), Throws.TypeOf<InvalidDataException>());
            Assert.That(() => v.ToLongArray(n - 1, 2), Throws.TypeOf<InvalidDataException>());
        }
        Assert.That(SaView.FromMapped(image.SaPointer, n, width, isChecked: false)[n - 1], Is.EqualTo(image.Sa[^1]));
    }

    [Test]
    public void NonGeneric_TheIntegerIndexerRefusesAPositionAboveTheIntegers_AndTheLongOneReadsIt()
    {
        ulong[] values = [7, (1UL << 32) + 5, 1UL << 31, (1UL << 31) - 1];
        using var guarded = new Fuzz.GuardedImage(ViewImage.Encode(values, 5, 3), flush: true);
        SaView view = SaView.FromMapped(guarded.Data, values.Length, 5, isChecked: true);
        Assert.That(view[0], Is.EqualTo(7));
        Assert.That(view[1L], Is.EqualTo((1L << 32) + 5));
        Assert.That(view[2L], Is.EqualTo(1L << 31));
        Assert.That(view[3], Is.EqualTo(int.MaxValue));
        Assert.That(() => view[1], Throws.TypeOf<InvalidDataException>());
        Assert.That(() => view[2], Throws.TypeOf<InvalidDataException>(), "2^31 does not fit an int");
        Assert.That(() => view.ToArray(1, 1), Throws.TypeOf<InvalidDataException>());
        Assert.That(() => view.CopyTo(1L, new int[1]), Throws.TypeOf<InvalidDataException>());
        Assert.That(view.ToLongArray(1, 2), Is.EqualTo(new[] { (1L << 32) + 5, 1L << 31 }));
    }

    [Test]
    public void NonGeneric_APointerOfAFiveByteSection_IsNotAnIntegerPointer()
    {
        using ViewImage image = ViewImage.Of("banana", 5, LcpKind.U8);
        SaView view = SaView.FromMapped(image.SaPointer, image.Count, 5, isChecked: true);
        Assert.That(() => view.Pointer, Throws.TypeOf<NotSupportedException>().With.Message.Contains("5"));
        using ViewImage four = ViewImage.Of("banana", 4, LcpKind.U8);
        Assert.That((long)SaView.FromMapped(four.SaPointer, four.Count, 4, isChecked: true).Pointer, Is.EqualTo((long)four.SaPointer));
        Assert.That(new SaView([1, 2, 3]).Pointer == null, Is.True);
    }

    [Test]
    public void NonGeneric_TheConstructors_RejectWhatIsNoSection()
    {
        using ViewImage image = ViewImage.Of("banana", 4, LcpKind.U8);
        Assert.That(() => SaView.FromMapped(image.SaPointer, image.Count, 3, true), Throws.TypeOf<ArgumentOutOfRangeException>());
        Assert.That(() => SaView.FromMapped(image.SaPointer, image.Count, 8, true), Throws.TypeOf<ArgumentOutOfRangeException>());
        Assert.That(() => SaView.FromMapped(image.SaPointer, -1, 4, true), Throws.TypeOf<ArgumentOutOfRangeException>());
        Assert.That(() => new SaView(null!), Throws.TypeOf<ArgumentNullException>());
        Assert.That(SaView.FromMapped(image.SaPointer, 0, 5, true).Length, Is.EqualTo(0));
    }

    [Test]
    public void NonGeneric_OverAnArray_IsTheArrayOfTheIndex_AndReadsLikeTheMappedForm()
    {
        int[] sa = [5, 3, 1, 0, 4, 2];
        var view = new SaView(sa);
        Assert.That(view.Array, Is.SameAs(sa));
        Assert.That(view[3L], Is.EqualTo(0));
        Assert.That(view.At(5), Is.EqualTo(2));
        Assert.That(() => view.At(6), Throws.TypeOf<InvalidDataException>());
        Assert.That(() => view[6L], Throws.TypeOf<IndexOutOfRangeException>().Or.TypeOf<InvalidDataException>());
        Assert.That(view.ToLongArray(1, 3), Is.EqualTo(new long[] { 3, 1, 0 }));
        var wide = new long[2];
        view.CopyTo(4L, wide);
        Assert.That(wide, Is.EqualTo(new long[] { 4, 2 }));
    }
}
