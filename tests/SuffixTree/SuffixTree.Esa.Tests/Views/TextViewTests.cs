using System.Runtime.InteropServices;
using SuffixTree.Esa.External;
using SuffixTree.Esa.Tests.External;
using SuffixTree.Indexing;

namespace SuffixTree.Esa.Tests.Views;

/// <summary>
/// The text views of work package S2.3 at 64 bits: <see cref="EsaText"/> (length, indexer, substrings, copies in blocks) over a packed section of any size, and the two contracts built on it, <see cref="ITextSource"/> (32-bit:
/// a caller's argument error is the exception of the in-memory string source, a text longer than the integers is <see cref="NotSupportedException"/> naming the long member) and <see cref="IIndexText"/> (64-bit). Review row R06: the
/// caller's errors are <see cref="ArgumentOutOfRangeException"/> / <see cref="IndexOutOfRangeException"/>, never "the file is damaged" - which is what a position that came from the file is.
/// </summary>
[TestFixture]
public unsafe class TextViewTests
{
    private static EsaText Packed(string text, out GCHandle pin)
    {
        char[] letters = [.. text.Distinct().OrderBy(c => c)];
        var symbols = new byte[text.Length];
        for (int i = 0; i < text.Length; i++) symbols[i] = (byte)Array.IndexOf(letters, text[i]);
        int bits = TextOracle.BitsOf(Math.Max(1, letters.Length));
        byte[] image = TextOracle.Pack(symbols, bits);
        pin = GCHandle.Alloc(image, GCHandleType.Pinned);
        return EsaText.FromPacked(MappedTextReader.FromMemory((byte*)pin.AddrOfPinnedObject(), image.Length, text.Length, bits), letters, (long)text.Length);
    }

    private static string Sample(int n) => string.Concat(Enumerable.Range(0, n).Select(i => (char)('a' + (((i * 7) + (i / 5)) % 4))));

    [Test]
    public void TheLongMembers_AgreeWithTheIntegerOnes_OverAPackedSection()
    {
        string text = Sample(10_000);
        EsaText t = Packed(text, out GCHandle pin);
        try
        {
            Assert.That(t.Length, Is.EqualTo(10_000L));
            for (long i = 0; i < 10_000; i += 7) Assert.That(t[i], Is.EqualTo(text[(int)i]), $"position {i}");
            Assert.That(t[9_999L], Is.EqualTo(text[9_999]));
            Assert.That(t.Substring(4_095L, 4_097), Is.EqualTo(text.Substring(4_095, 4_097)));
            Assert.That(t.Substring(0L, 10_000), Is.EqualTo(text));
            var dest = new char[5_000];
            t.CopyTo(4_500L, dest);
            Assert.That(new string(dest), Is.EqualTo(text.Substring(4_500, 5_000)));
            Assert.That(() => t[10_000L], Throws.TypeOf<InvalidDataException>().With.Message.Contains("10000"));
            Assert.That(() => t[-1L], Throws.TypeOf<InvalidDataException>());
            Assert.That(() => t[1L << 33], Throws.TypeOf<InvalidDataException>());
            Assert.That(() => t.Substring(9_999L, 2), Throws.TypeOf<InvalidDataException>());
            Assert.That(() => t.Substring(1L << 33, 0), Throws.TypeOf<InvalidDataException>());
            Assert.That(() => t.CopyTo(9_999L, new char[2]), Throws.TypeOf<InvalidDataException>());
        }
        finally
        {
            pin.Free();
        }
    }

    // ---------------------------------------------------------------------------------------------- R06: the caller's errors

    /// <summary>The exception types of every call of the contract that a caller can get wrong, over the string source and over the packed one: equal, and never the damage of a file.</summary>
    [Test]
    public void ITextSource_TheCallersErrors_AreTheExceptionsOfTheInMemorySource([Values(81, 5, 0)] int n)
    {
        string text = Sample(n);
        ITextSource memory = EsaText.FromString(text).AsTextSource();
        EsaText packed = n == 0 ? EsaText.FromString("") : Packed(text, out _);
        ITextSource source = packed.AsTextSource();
        (int Start, int Length)[] ranges = [(0, 0), (n, 0), (n + 1, 0), (-1, 1), (0, -1), (n, 1), (n - 1, 2), (0, n + 1), (int.MaxValue, 1), (int.MinValue, 1), (1, int.MaxValue)];
        foreach ((int start, int length) in ranges)
        {
            AssertSame(() => memory.Substring(start, length), () => source.Substring(start, length), $"Substring({start}, {length})");
            AssertSame(() => memory.Slice(start, length).ToString(), () => source.Slice(start, length).ToString(), $"Slice({start}, {length})");
        }
        foreach (int index in new[] { n, n + 1, -1, int.MaxValue, int.MinValue })
        {
            AssertSame(() => memory[index], () => source[index], $"[{index}]");
        }
        Assert.That(source.Length, Is.EqualTo(n));
        Assert.That(() => source[n], Throws.Nothing.Or.TypeOf<IndexOutOfRangeException>());
    }

    private static void AssertSame<T>(Func<T> expected, Func<T> actual, string what)
    {
        Exception? e = Capture(expected), a = Capture(actual);
        Assert.That(a?.GetType(), Is.EqualTo(e?.GetType()), $"{what}: {a?.GetType().Name} {a?.Message}");
        Assert.That(a, Is.Not.InstanceOf<InvalidDataException>(), what);
    }

    private static Exception? Capture<T>(Func<T> call)
    {
        try
        {
            call();
            return null;
        }
        catch (Exception e)
        {
            return e;
        }
    }

    [Test]
    public void ITextSource_AValueFromTheFile_IsStillInvalidData()
    {
        // sigma 3 on 2 bits: a symbol 3 in the section is not a letter: the file is damaged, whoever asks
        EsaText t = Packed("abcabcabc", out GCHandle pin);
        try
        {
            byte* data = (byte*)pin.AddrOfPinnedObject();
            data[0] |= 0xC0;
            ITextSource source = t.AsTextSource();
            Assert.That(() => source[0], Throws.TypeOf<InvalidDataException>().With.Message.Contains("symbol 3"));
            Assert.That(() => source.Substring(0, 4), Throws.TypeOf<InvalidDataException>());
            Assert.That(() => t.AsIndexText().SymbolAt(0), Throws.TypeOf<InvalidDataException>());
            Assert.That(source[1], Is.EqualTo('b'));
        }
        finally
        {
            pin.Free();
        }
    }

    [Test]
    public void IIndexText_TheCallersErrors_AreArgumentOutOfRange_AsInMemory()
    {
        string text = Sample(81);
        IIndexText memory = EsaText.FromString(text).AsIndexText();
        IIndexText packed = Packed(text, out _).AsIndexText();
        Assert.That(packed.Length, Is.EqualTo(81L));
        foreach (long position in new long[] { -1, long.MinValue, 80, 81, 82, long.MaxValue, 1L << 40 })
        {
            AssertSame(() => memory.SymbolAt(position), () => packed.SymbolAt(position), $"SymbolAt({position})");
            if (position >= 0) Assert.That(packed.SymbolAt(position), Is.EqualTo(memory.SymbolAt(position)), $"SymbolAt({position})");
        }
        foreach (long start in new long[] { -1, long.MinValue, 0, 40, 81, 82, long.MaxValue, 1L << 40 })
        {
            foreach (int length in new[] { 0, 1, 30, 100 })
            {
                var a = new char[length];
                var b = new char[length];
                Exception? ea = Capture(() => memory.CopyTo(start, a)), eb = Capture(() => packed.CopyTo(start, b));
                Assert.That(eb?.GetType(), Is.EqualTo(ea?.GetType()), $"CopyTo({start}, {length})");
                if (ea is null) Assert.That(b, Is.EqualTo(a), $"CopyTo({start}, {length}) content");
            }
        }
    }

    // ---------------------------------------------------------------------------------------------- above 2^31 symbols

    [Test]
    public void ATextOfMoreThanTheIntegers_IsReadAtThe64BitPositions_AndTheIntegerContractSaysSo()
    {
        // 2^31 + 10 symbols of one bit (a = 0, b = 1): 256 MiB of address space; the bytes around the positions that the test reads are open
        long n = (1L << 31) + 10;
        const int bits = 1;
        long total = PackedTextFormat.TotalBytes(n, bits);
        using var section = new VirtualSection(total);
        section.Open(0, 4096);
        section.Open(total - 4096, 4096);
        section.Open(((1L << 31) / 8) - 4096, 8192);
        // symbols 2^31 - 3 .. 2^31 + 2 are b b a b a b: the bit stream is big endian, symbol i is bit 7 - (i mod 8) of byte i / 8
        long[] ones = [(1L << 31) - 3, (1L << 31) - 2, (1L << 31) , (1L << 31) + 2, n - 1, 0, 7];
        foreach (long i in ones) section.Data[i / 8] |= (byte)(0x80 >> (int)(i % 8));
        var reader = MappedTextReader.FromMemory(section.Data, total, n, bits);
        EsaText text = EsaText.FromPacked(reader, ['a', 'b'], n);
        Assert.That(text.Length, Is.EqualTo(n));
        Assert.That(text.IsPacked, Is.True);
        foreach (long i in ones) Assert.That(text[i], Is.EqualTo('b'), $"position {i}");
        foreach (long i in new long[] { 1, (1L << 31) - 1, (1L << 31) + 1, n - 2 }) Assert.That(text[i], Is.EqualTo('a'), $"position {i}");
        Assert.That(() => text[n], Throws.TypeOf<InvalidDataException>());
        var block = new char[8];
        text.CopyTo((1L << 31) - 3, block);
        Assert.That(new string(block), Is.EqualTo("bbababaa"), "positions 2^31 - 3 .. 2^31 + 4: a block across 2^31");
        Assert.That(text.Substring((1L << 31) - 3, 8), Is.EqualTo(new string(block)));

        IIndexText index = text.AsIndexText();
        Assert.That(index.Length, Is.EqualTo(n));
        Assert.That(index.SymbolAt((1L << 31) + 2), Is.EqualTo((int)'b'));
        Assert.That(index.SymbolAt(n), Is.EqualTo(-1));
        Assert.That(index.SymbolAt(long.MaxValue), Is.EqualTo(-1));
        Assert.That(index.CopyTo((1L << 31) - 3, block), Is.EqualTo(8));
        Assert.That(index.CopyTo(n - 3, block), Is.EqualTo(3));
        Assert.That(index.CopyTo(n, block), Is.EqualTo(0));
        Assert.That(() => index.CopyTo(n + 1, block), Throws.TypeOf<ArgumentOutOfRangeException>());

        ITextSource source = text.AsTextSource();
        Assert.That(() => source.Length, Throws.TypeOf<NotSupportedException>().With.Message.Contains("IIndexText"));
        Assert.That(source[0], Is.EqualTo('b'));
        Assert.That(source[7], Is.EqualTo('b'));
        Assert.That(source.Substring(0, 8), Is.EqualTo(ReadBack(text, 0, 8)));
        Assert.That(() => source[-1], Throws.TypeOf<IndexOutOfRangeException>());
        Assert.That(() => source.Substring(-1, 1), Throws.TypeOf<ArgumentOutOfRangeException>());
    }

    private static string ReadBack(EsaText text, long start, int count)
    {
        var chars = new char[count];
        for (int i = 0; i < count; i++) chars[i] = text[start + i];
        return new string(chars);
    }

    [Test]
    public void TheIntegerSubstring_OfATextAboveTheIntegers_IsNotSupportedWhenItWouldNotFitAString()
    {
        long n = (1L << 31) + 10;
        long total = PackedTextFormat.TotalBytes(n, 1);
        using var section = new VirtualSection(total);
        section.Open(0, 4096);
        section.Open(total - 4096, 4096);
        EsaText text = EsaText.FromPacked(MappedTextReader.FromMemory(section.Data, total, n, 1), ['a', 'b'], n);
        Assert.That(() => text.Substring(5L), Throws.TypeOf<NotSupportedException>().With.Message.Contains("Substring"));
        Assert.That(text.Substring(n - 4), Is.EqualTo("aaaa"));
        Assert.That(text.Substring(n), Is.Empty);
        Assert.That(() => text.Substring(n + 1), Throws.TypeOf<InvalidDataException>());
    }
}
