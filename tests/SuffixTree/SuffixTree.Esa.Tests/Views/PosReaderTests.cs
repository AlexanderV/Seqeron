using System.Runtime.InteropServices;
using SuffixTree.Esa.Tests.Fuzz;

namespace SuffixTree.Esa.Tests.Views;

/// <summary>
/// The position readers of work package S2.3 (ADR-005 sections 3.3 and 6.1): the decode of every stored width at the values that a width must hold (0, 1, 2, 2^31 - 1, 2^31, 2^32 - 1, 2^32, 2^40 - 1), at every alignment of the section and at
/// every 64-byte line and 4,096-byte page straddle of the five-byte entries, on the last element of a section (its three bytes of padding are inside the payload: the section lies against an inaccessible page), in both access modes
/// (<see cref="CheckedAccess"/>: every read is bounded and a read outside the section is an <see cref="InvalidDataException"/>; <see cref="TrustedAccess"/>), and at element indexes above 2^31 and 2^32 over reserved address space
/// (<see cref="VirtualSection"/>: the address arithmetic is 64-bit). The encoder here is written from the byte layout of ADR-005 section 6.1 and shares nothing with the readers.
/// </summary>
[TestFixture]
public unsafe class PosReaderTests
{
    private const ulong Mask40 = 0xFF_FFFF_FFFFUL;
    private static readonly ulong[] Specials40 = [0, 1, 2, (1UL << 31) - 1, 1UL << 31, (1UL << 32) - 1, 1UL << 32, (1UL << 32) + 1, (1UL << 39), Mask40 - 1, Mask40];

    // ------------------------------------------------------------------------------------------------ the reference encoder

    /// <summary>The bytes of <paramref name="values"/> at <paramref name="width"/> bytes per element, little endian, and <paramref name="padding"/> zero bytes after the last element.</summary>
    private static byte[] Encode(ReadOnlySpan<ulong> values, int width, int padding)
    {
        var bytes = new byte[(width * values.Length) + padding];
        for (int i = 0; i < values.Length; i++)
        {
            for (int j = 0; j < width; j++) bytes[(width * i) + j] = (byte)(values[i] >> (8 * j));
        }
        return bytes;
    }

    private static ulong Hash(long i, ulong seed)
    {
        ulong x = (ulong)i + seed + 0x9E3779B97F4A7C15UL;
        x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9UL;
        x = (x ^ (x >> 27)) * 0x94D049BB133111EBUL;
        return x ^ (x >> 31);
    }

    private static ulong[] Random(int count, int bits, ulong seed)
    {
        var values = new ulong[count];
        ulong mask = bits == 64 ? ulong.MaxValue : (1UL << bits) - 1;
        for (int i = 0; i < count; i++) values[i] = Hash(i, seed) & mask;
        return values;
    }

    /// <summary>An aligned block of native memory (page aligned), so that the offset of a section inside it is the address alignment of the section.</summary>
    private sealed class Aligned : IDisposable
    {
        private readonly void* _p;

        internal Aligned(int bytes)
        {
            _p = NativeMemory.AlignedAlloc((nuint)bytes, 4096);
            new Span<byte>(_p, bytes).Clear();
            Data = (byte*)_p;
        }

        internal byte* Data { get; }

        public void Dispose() => NativeMemory.AlignedFree(_p);
    }

    // ------------------------------------------------------------------------------------------------ Pos40: the decode

    [Test]
    public void Pos40_DecodesTheValuesThatTheWidthMustHold_AtEveryAlignmentOfTheSection()
    {
        // each special value at five indexes (the first, the second, the last, and two in the middle) and at every offset of the section inside a page: the load is unaligned
        int count = Specials40.Length * 5;
        var values = new ulong[count];
        for (int i = 0; i < count; i++) values[i] = Specials40[(i + (i / Specials40.Length)) % Specials40.Length];
        byte[] image = Encode(values, 5, 3);
        using var block = new Aligned(image.Length + 128);
        for (int offset = 0; offset < 64; offset++)
        {
            byte* section = block.Data + offset;
            image.CopyTo(new Span<byte>(section, image.Length));
            var checkedReader = new Pos40<CheckedAccess>(section, count);
            var trustedReader = new Pos40<TrustedAccess>(section, count);
            for (int i = 0; i < count; i++)
            {
                Assert.That(checkedReader.Get(i), Is.EqualTo(values[i]), $"checked, offset {offset}, element {i}");
                Assert.That(trustedReader.Get(i), Is.EqualTo(values[i]), $"trusted, offset {offset}, element {i}");
            }
        }
    }

    [Test]
    public void Pos40_DecodesEveryElementOfThreePages_IncludingEveryLineAndPageStraddle()
    {
        const int count = 2500; // 12,503 bytes: three pages and a bit
        ulong[] values = Random(count, 40, 0x40);
        byte[] image = Encode(values, 5, 3);
        using var block = new Aligned(image.Length + 64);
        image.CopyTo(new Span<byte>(block.Data, image.Length));
        var reader = new Pos40<CheckedAccess>(block.Data, count);
        int lineStraddles = 0, pageStraddles = 0, loadStraddles = 0;
        for (int i = 0; i < count; i++)
        {
            long at = 5L * i;
            if ((at & 63) > 59) lineStraddles++; // the five bytes of the entry cross a 64-byte line
            if ((at & 63) > 56) loadStraddles++; // the 8-byte load crosses a 64-byte line
            if ((at & 4095) > 4091) pageStraddles++; // the five bytes of the entry cross a page
            Assert.That(reader.Get(i), Is.EqualTo(values[i]), $"element {i}");
        }
        Assert.That(lineStraddles, Is.GreaterThan(30), "the test must contain line straddles");
        Assert.That(loadStraddles, Is.GreaterThan(lineStraddles), "and loads that cross a line");
        Assert.That(pageStraddles, Is.GreaterThanOrEqualTo(2), "and entries that cross a page");
    }

    [Test]
    public void Pos40_TheEntryThatCrossesAPage_IsDecodedWhole()
    {
        // 4,095 = 5 * 819: the entry 819 has its first byte on the last byte of the first page
        ulong[] values = new ulong[1200];
        values[818] = Mask40;
        values[819] = 0x12_3456_789AUL;
        values[820] = 1UL << 32;
        byte[] image = Encode(values, 5, 3);
        using var block = new Aligned(image.Length + 64);
        image.CopyTo(new Span<byte>(block.Data, image.Length));
        var reader = new Pos40<CheckedAccess>(block.Data, values.Length);
        Assert.That(reader.Get(818), Is.EqualTo(Mask40));
        Assert.That(reader.Get(819), Is.EqualTo(0x12_3456_789AUL));
        Assert.That(reader.Get(820), Is.EqualTo(1UL << 32));
        Assert.That(reader.Get(817), Is.EqualTo(0UL));
    }

    [Test]
    public void Pos40_TheLoadKeepsTheFortyLowBits_AndNeverTheBytesOfTheNeighbours()
    {
        // every element is 0 and every neighbour byte is 0xFF: the three upper bytes of the load belong to the next entry and are discarded
        const int count = 40;
        var bytes = new byte[(5 * count) + 3];
        Array.Fill(bytes, (byte)0xFF);
        for (int i = 0; i < count; i += 2)
        {
            for (int j = 0; j < 5; j++) bytes[(5 * i) + j] = 0;
        }
        bytes[^3] = bytes[^2] = bytes[^1] = 0;
        using var block = new Aligned(bytes.Length + 64);
        bytes.CopyTo(new Span<byte>(block.Data, bytes.Length));
        var reader = new Pos40<TrustedAccess>(block.Data, count);
        for (int i = 0; i < count; i += 2) Assert.That(reader.Get(i), Is.EqualTo(0UL), $"element {i}");
        for (int i = 1; i < count; i += 2) Assert.That(reader.Get(i), Is.EqualTo(Mask40), $"element {i}");
    }

    [Test]
    public void Pos40_TheLastElement_ReadsInsideThePayload_AgainstAnInaccessiblePage([Values(1, 2, 7, 12, 818, 819, 820, 1000)] int count, [Values(true, false)] bool isChecked)
    {
        ulong[] values = Random(count, 40, 0x41);
        values[^1] = Mask40; // the last entry has the high bytes set: a load that read the padding as data, or a mask that kept a bit more, would show
        byte[] payload = Encode(values, 5, PosReader.Padding40);
        Assert.That(payload, Has.Length.EqualTo((5 * count) + 3));
        using var image = new GuardedImage(payload, flush: true);
        Assert.That(image.Guarded || !OperatingSystem.IsLinux(), Is.True);
        ulong last = isChecked ? new Pos40<CheckedAccess>(image.Data, count).Get(count - 1) : new Pos40<TrustedAccess>(image.Data, count).Get(count - 1);
        Assert.That(last, Is.EqualTo(Mask40));
        for (int i = 0; i < count; i++)
        {
            ulong v = isChecked ? new Pos40<CheckedAccess>(image.Data, count).Get(i) : new Pos40<TrustedAccess>(image.Data, count).Get(i);
            Assert.That(v, Is.EqualTo(values[i]), $"element {i}");
        }
    }

    // ------------------------------------------------------------------------------------------------ the other widths

    [Test]
    public void Pos32_DecodesUnsignedWords_WithoutSignExtension()
    {
        ulong[] values = [0, 1, 2, (1UL << 31) - 1, 1UL << 31, 0xFFFF_FFFEUL, 0xFFFF_FFFFUL, 0x8000_0001UL];
        byte[] image = Encode(values, 4, 0);
        fixed (byte* p = image)
        {
            var c = new Pos32<CheckedAccess>(p, values.Length);
            var t = new Pos32<TrustedAccess>(p, values.Length);
            for (int i = 0; i < values.Length; i++)
            {
                Assert.That(c.Get(i), Is.EqualTo(values[i]), $"checked {i}");
                Assert.That(t.Get(i), Is.EqualTo(values[i]), $"trusted {i}");
            }
        }
    }

    [Test]
    public void Pos8_Pos16_Pos64_DecodeTheirWidths([Values(true, false)] bool isChecked)
    {
        ulong[] v8 = [0, 1, 127, 128, 254, 255];
        ulong[] v16 = [0, 1, 255, 256, 32767, 32768, 65534, 65535];
        ulong[] v64 = [0, 1, (1UL << 32) - 1, 1UL << 32, 1UL << 40, (1UL << 63) - 1, 1UL << 63, ulong.MaxValue];
        byte[] b8 = Encode(v8, 1, 0), b16 = Encode(v16, 2, 0), b64 = Encode(v64, 8, 0);
        fixed (byte* p8 = b8, p16 = b16, p64 = b64)
        {
            for (int i = 0; i < v8.Length; i++) Assert.That(isChecked ? new Pos8<CheckedAccess>(p8, v8.Length).Get(i) : new Pos8<TrustedAccess>(p8, v8.Length).Get(i), Is.EqualTo(v8[i]), $"8 bits, {i}");
            for (int i = 0; i < v16.Length; i++) Assert.That(isChecked ? new Pos16<CheckedAccess>(p16, v16.Length).Get(i) : new Pos16<TrustedAccess>(p16, v16.Length).Get(i), Is.EqualTo(v16[i]), $"16 bits, {i}");
            for (int i = 0; i < v64.Length; i++) Assert.That(isChecked ? new Pos64<CheckedAccess>(p64, v64.Length).Get(i) : new Pos64<TrustedAccess>(p64, v64.Length).Get(i), Is.EqualTo(v64[i]), $"64 bits, {i}");
        }
    }

    [Test]
    public void PosInt_ReadsTheArrayAsUnsignedWords_AndAgreesWithPos32OverTheSameBytes()
    {
        int[] array = [0, 1, 2, int.MaxValue, int.MinValue, -1, 12345];
        var reader = new PosInt(array);
        Assert.That(reader.Count, Is.EqualTo(array.Length));
        ulong[] expected = [0, 1, 2, int.MaxValue, 1UL << 31, 0xFFFF_FFFFUL, 12345];
        for (int i = 0; i < array.Length; i++) Assert.That(reader.Get(i), Is.EqualTo(expected[i]), $"element {i}");
        byte[] bytes = new byte[4 * array.Length];
        Buffer.BlockCopy(array, 0, bytes, 0, bytes.Length);
        fixed (byte* p = bytes)
        {
            var wide = new Pos32<CheckedAccess>(p, array.Length);
            for (int i = 0; i < array.Length; i++) Assert.That(wide.Get(i), Is.EqualTo(reader.Get(i)), $"element {i}");
        }
    }

    [Test]
    public void PosInt_AReadOutsideTheArray_IsTheArraysOwnException()
    {
        var reader = new PosInt([1, 2, 3]);
        Assert.That(() => reader.Get(3), Throws.TypeOf<IndexOutOfRangeException>());
        Assert.That(() => reader.Get(-1), Throws.TypeOf<IndexOutOfRangeException>());
        Assert.That(() => reader.Get(1L << 40), Throws.TypeOf<IndexOutOfRangeException>());
        Assert.That(() => new PosInt(null!), Throws.TypeOf<ArgumentNullException>());
    }

    // ------------------------------------------------------------------------------------------------ the bounds of Checked mode

    [Test]
    public void Checked_AReadOutsideTheSection_IsInvalidData_ForEveryWidth()
    {
        var bytes = new byte[8 * 4 + 16];
        fixed (byte* p = bytes)
        {
            nint q = (nint)p;
            long[] outside = [4, 5, 100, 1L << 31, 1L << 32, 1L << 40, long.MaxValue, -1, -2, int.MinValue, long.MinValue];
            foreach (long i in outside)
            {
                Assert.That(() => new Pos8<CheckedAccess>((byte*)q, 4).Get(i), Throws.TypeOf<InvalidDataException>(), $"8 bits, {i}");
                Assert.That(() => new Pos16<CheckedAccess>((byte*)q, 4).Get(i), Throws.TypeOf<InvalidDataException>(), $"16 bits, {i}");
                Assert.That(() => new Pos32<CheckedAccess>((byte*)q, 4).Get(i), Throws.TypeOf<InvalidDataException>(), $"32 bits, {i}");
                Assert.That(() => new Pos40<CheckedAccess>((byte*)q, 4).Get(i), Throws.TypeOf<InvalidDataException>(), $"40 bits, {i}");
                Assert.That(() => new Pos64<CheckedAccess>((byte*)q, 4).Get(i), Throws.TypeOf<InvalidDataException>(), $"64 bits, {i}");
            }
            for (long i = 0; i < 4; i++)
            {
                Assert.That(() => new Pos40<CheckedAccess>((byte*)q, 4).Get(i), Throws.Nothing);
                Assert.That(() => new Pos64<CheckedAccess>((byte*)q, 4).Get(i), Throws.Nothing);
            }
        }
    }

    [Test]
    public void Checked_TheMessageNamesTheElementAndTheCount()
    {
        var bytes = new byte[64];
        fixed (byte* p = bytes)
        {
            nint q = (nint)p;
            InvalidDataException e = Assert.Throws<InvalidDataException>(() => new Pos32<CheckedAccess>((byte*)q, 7).Get(9))!;
            Assert.That(e.Message, Does.Contain("7").And.Contain("9").And.Contain("damaged"));
        }
    }

    [Test]
    public void Checked_TheLastElementIsReadable_AndTheNextOneIsNot_AtTheEdgeOfTheGuardPage()
    {
        // 12 entries of 5 bytes and the 3 bytes of padding against the guard page: element 11 is read, element 12 is refused before any load
        ulong[] values = Random(12, 40, 0x42);
        using var image = new GuardedImage(Encode(values, 5, 3), flush: true);
        var reader = new Pos40<CheckedAccess>(image.Data, 12);
        Assert.That(reader.Get(11), Is.EqualTo(values[11]));
        Assert.That(() => reader.Get(12), Throws.TypeOf<InvalidDataException>());
        var wide = new Pos32<CheckedAccess>(image.Data, 15);
        Assert.That(() => wide.Get(15), Throws.TypeOf<InvalidDataException>());
    }

    [Test]
    public void TheConstructors_RejectANegativeCount_AndTakeZero()
    {
        var bytes = new byte[16];
        fixed (byte* p = bytes)
        {
            nint q = (nint)p;
            Assert.That(() => new Pos8<CheckedAccess>((byte*)q, -1), Throws.TypeOf<ArgumentOutOfRangeException>());
            Assert.That(() => new Pos16<CheckedAccess>((byte*)q, -1), Throws.TypeOf<ArgumentOutOfRangeException>());
            Assert.That(() => new Pos32<TrustedAccess>((byte*)q, -1), Throws.TypeOf<ArgumentOutOfRangeException>());
            Assert.That(() => new Pos40<CheckedAccess>((byte*)q, -1), Throws.TypeOf<ArgumentOutOfRangeException>());
            Assert.That(() => new Pos64<TrustedAccess>((byte*)q, -1), Throws.TypeOf<ArgumentOutOfRangeException>());
            Assert.That(new Pos40<CheckedAccess>(p, 0).Count, Is.EqualTo(0));
            Assert.That(() => new Pos40<CheckedAccess>((byte*)q, 0).Get(0), Throws.TypeOf<InvalidDataException>());
        }
    }

    [Test]
    public void TheReaders_SayTheirCountAndTheirWidth()
    {
        var bytes = new byte[64];
        fixed (byte* p = bytes)
        {
            Assert.That(new Pos8<CheckedAccess>(p, 3).Count, Is.EqualTo(3));
            Assert.That(new Pos16<CheckedAccess>(p, 4).Count, Is.EqualTo(4));
            Assert.That(new Pos32<CheckedAccess>(p, 5).Count, Is.EqualTo(5));
            Assert.That(new Pos40<TrustedAccess>(p, 6).Count, Is.EqualTo(6));
            Assert.That(new Pos64<TrustedAccess>(p, 7).Count, Is.EqualTo(7));
            Assert.That(new Pos40<TrustedAccess>(p, 1L << 40).Count, Is.EqualTo(1L << 40));
        }
        Assert.That(Pos8<CheckedAccess>.ElementBytes, Is.EqualTo(1));
        Assert.That(Pos16<CheckedAccess>.ElementBytes, Is.EqualTo(2));
        Assert.That(Pos32<CheckedAccess>.ElementBytes, Is.EqualTo(4));
        Assert.That(Pos40<TrustedAccess>.ElementBytes, Is.EqualTo(5));
        Assert.That(Pos64<TrustedAccess>.ElementBytes, Is.EqualTo(8));
        Assert.That(PosInt.ElementBytes, Is.EqualTo(4));
        Assert.That(PosReader.Mask40, Is.EqualTo(Mask40));
        Assert.That(PosReader.Padding40, Is.EqualTo(3));
    }

    [Test]
    public void TheReaders_AgreeWithTheReferenceEncoder_OnRandomValuesOfEveryWidth([Values(true, false)] bool isChecked)
    {
        const int count = 3000;
        ulong[] v32 = Random(count, 32, 1), v40 = Random(count, 40, 2), v64 = Random(count, 64, 3), v16 = Random(count, 16, 4), v8 = Random(count, 8, 5);
        byte[] b32 = Encode(v32, 4, 0), b40 = Encode(v40, 5, 3), b64 = Encode(v64, 8, 0), b16 = Encode(v16, 2, 0), b8 = Encode(v8, 1, 0);
        fixed (byte* p32 = b32, p40 = b40, p64 = b64, p16 = b16, p8 = b8)
        {
            for (int i = 0; i < count; i++)
            {
                Assert.That(isChecked ? new Pos32<CheckedAccess>(p32, count).Get(i) : new Pos32<TrustedAccess>(p32, count).Get(i), Is.EqualTo(v32[i]));
                Assert.That(isChecked ? new Pos40<CheckedAccess>(p40, count).Get(i) : new Pos40<TrustedAccess>(p40, count).Get(i), Is.EqualTo(v40[i]));
                Assert.That(isChecked ? new Pos64<CheckedAccess>(p64, count).Get(i) : new Pos64<TrustedAccess>(p64, count).Get(i), Is.EqualTo(v64[i]));
                Assert.That(isChecked ? new Pos16<CheckedAccess>(p16, count).Get(i) : new Pos16<TrustedAccess>(p16, count).Get(i), Is.EqualTo(v16[i]));
                Assert.That(isChecked ? new Pos8<CheckedAccess>(p8, count).Get(i) : new Pos8<TrustedAccess>(p8, count).Get(i), Is.EqualTo(v8[i]));
            }
        }
    }

    // ------------------------------------------------------------------------------------------------ element indexes above 2^31 and 2^32 (reserved address space)

    [Test]
    public void Pos32_ReadsAtIndexesAbove2To31_AndAbove2To32([Values(true, false)] bool isChecked)
    {
        // 2^32 + 16 entries of 4 bytes = 16 GiB of address space; only the pages of the entries that the test writes are open
        long[] at = [(1L << 31) - 1, 1L << 31, (1L << 31) + 1, (1L << 32) - 1, 1L << 32, (1L << 32) + 10];
        long count = (1L << 32) + 16;
        using var section = new VirtualSection(4 * count);
        foreach (long i in at) section.Open(4 * i, 4);
        foreach (long i in at) *(uint*)(section.Data + (4 * i)) = (uint)(0xA000_0000UL + (ulong)(i & 0xFFFF));
        foreach (long i in at)
        {
            ulong expected = 0xA000_0000UL + (ulong)(i & 0xFFFF);
            Assert.That(isChecked ? new Pos32<CheckedAccess>(section.Data, count).Get(i) : new Pos32<TrustedAccess>(section.Data, count).Get(i), Is.EqualTo(expected), $"element {i}");
        }
        Assert.That(() => new Pos32<CheckedAccess>(section.Data, count).Get(count), Throws.TypeOf<InvalidDataException>());
    }

    [Test]
    public void Pos40_ReadsAtIndexesAbove2To31_2To32_And2To33([Values(true, false)] bool isChecked)
    {
        long count = (1L << 33) + 16; // 5 x 2^33 = 43 GB of address space
        long[] at = [(1L << 31) - 1, 1L << 31, (1L << 31) + 1, (1L << 32) - 1, 1L << 32, (1L << 32) + 10, (1L << 33) - 1, 1L << 33, count - 1];
        using var section = new VirtualSection((5 * count) + PosReader.Padding40);
        foreach (long i in at) section.Open(5 * i, 8);
        for (int k = 0; k < at.Length; k++)
        {
            ulong value = Mask40 - (ulong)k;
            long i = at[k];
            for (int j = 0; j < 5; j++) section.Data[(5 * i) + j] = (byte)(value >> (8 * j));
        }
        for (int k = 0; k < at.Length; k++)
        {
            long i = at[k];
            ulong expected = Mask40 - (ulong)k;
            Assert.That(isChecked ? new Pos40<CheckedAccess>(section.Data, count).Get(i) : new Pos40<TrustedAccess>(section.Data, count).Get(i), Is.EqualTo(expected), $"element {i}");
        }
        Assert.That(() => new Pos40<CheckedAccess>(section.Data, count).Get(count), Throws.TypeOf<InvalidDataException>());
        Assert.That(() => new Pos40<CheckedAccess>(section.Data, count).Get(long.MaxValue / 5), Throws.TypeOf<InvalidDataException>());
    }

    [Test]
    public void Pos64_Pos16_Pos8_ReadAtIndexesAbove2To31()
    {
        long count = (1L << 31) + 100;
        long[] at = [(1L << 31) - 1, 1L << 31, (1L << 31) + 5, count - 1];
        using var s64 = new VirtualSection(8 * count);
        using var s16 = new VirtualSection(2 * count);
        using var s8 = new VirtualSection(count);
        foreach (long i in at)
        {
            s64.Open(8 * i, 8);
            s16.Open(2 * i, 2);
            s8.Open(i, 1);
            *(ulong*)(s64.Data + (8 * i)) = (1UL << 40) + (ulong)i;
            *(ushort*)(s16.Data + (2 * i)) = (ushort)(0xFF00 + (i & 0xFF));
            s8.Data[i] = (byte)(0x80 + (i & 0x7F));
        }
        foreach (long i in at)
        {
            Assert.That(new Pos64<CheckedAccess>(s64.Data, count).Get(i), Is.EqualTo((1UL << 40) + (ulong)i));
            Assert.That(new Pos16<TrustedAccess>(s16.Data, count).Get(i), Is.EqualTo((ulong)(0xFF00 + (i & 0xFF))));
            Assert.That(new Pos8<CheckedAccess>(s8.Data, count).Get(i), Is.EqualTo((ulong)(0x80 + (i & 0x7F))));
        }
    }
}
