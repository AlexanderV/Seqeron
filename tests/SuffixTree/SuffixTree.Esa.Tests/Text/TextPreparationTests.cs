using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.Esa.Tests.Text;

/// <summary>1a.2: distinct-code-unit scan, order-preserving dense remap, symbol readers, size guards.</summary>
[TestFixture]
public class TextPreparationTests
{
    [Test]
    public void CountDistinct_CountsCodeUnits()
    {
        Assert.That(TextPreparation.CountDistinct(""), Is.EqualTo(0));
        Assert.That(TextPreparation.CountDistinct("aaaa"), Is.EqualTo(1));
        Assert.That(TextPreparation.CountDistinct("abcabc"), Is.EqualTo(3));
        Assert.That(TextPreparation.CountDistinct("\0￿𐀀"), Is.EqualTo(4));
        Assert.That(TextPreparation.CountDistinct(StringFamilies.AllDistinct(65536, 1)), Is.EqualTo(65536));
    }

    [Test]
    public void Prepare_EmptyAndTinyTexts()
    {
        var empty = TextPreparation.Prepare("");
        Assert.That(empty.Length, Is.EqualTo(0));
        Assert.That(empty.AlphabetSize, Is.EqualTo(0));
        Assert.That(empty.Width, Is.EqualTo(SymbolWidth.Byte));
        Assert.That(empty.Alphabet, Is.Empty);

        var one = TextPreparation.Prepare("x");
        Assert.That(one.Length, Is.EqualTo(1));
        Assert.That(one.AlphabetSize, Is.EqualTo(1));
        Assert.That(one.Alphabet, Is.EqualTo(new[] { 'x' }));
        Assert.That(one.Bytes.ToArray(), Is.EqualTo(new byte[] { 0 }));

        var two = TextPreparation.Prepare("ba");
        Assert.That(two.AlphabetSize, Is.EqualTo(2));
        Assert.That(two.Alphabet, Is.EqualTo(new[] { 'a', 'b' }));
        Assert.That(two.Bytes.ToArray(), Is.EqualTo(new byte[] { 1, 0 }));

        var same = TextPreparation.Prepare("zz");
        Assert.That(same.AlphabetSize, Is.EqualTo(1));
        Assert.That(same.Bytes.ToArray(), Is.EqualTo(new byte[] { 0, 0 }));
    }

    [Test]
    public void Prepare_RemapIsDenseOrderPreservingAndInvertible_OnEveryFamily()
    {
        foreach (var (name, text) in StringFamilies.Adversarial(150, 9))
        {
            var p = TextPreparation.Prepare(text);
            Assert.That(p.Length, Is.EqualTo(text.Length), name);
            Assert.That(p.AlphabetSize, Is.EqualTo(text.Distinct().Count()), name);
            Assert.That(p.Alphabet, Is.EqualTo(text.Distinct().OrderBy(c => c).ToArray()), name);
            for (int i = 0; i < text.Length; i++)
            {
                int ci = p.SymbolAt(i);
                Assert.That(ci, Is.InRange(0, p.AlphabetSize - 1), name);
                Assert.That(p.Alphabet[ci], Is.EqualTo(text[i]), name);
                for (int j = 0; j < text.Length; j += 7)
                    Assert.That(Math.Sign(p.SymbolAt(j).CompareTo(ci)), Is.EqualTo(Math.Sign(text[j].CompareTo(text[i]))), $"{name} {i} {j}");
            }
        }
    }

    [Test]
    public void Prepare_SpecialCodeUnits_KeepCodeUnitOrder()
    {
        // lone surrogates, a valid pair, U+0000 and U+FFFF are ordinary code units; order is by numeric value
        string text = "￿\uDC00\uD800\0a😀􏿿";
        var p = TextPreparation.Prepare(text);
        Assert.That(p.Alphabet, Is.EqualTo(text.Distinct().OrderBy(c => c).ToArray()));
        Assert.That(p.Alphabet[0], Is.EqualTo('\0'));
        Assert.That(p.Alphabet[^1], Is.EqualTo('￿'));
        for (int i = 0; i < text.Length; i++)
            Assert.That(p.Alphabet[p.SymbolAt(i)], Is.EqualTo(text[i]));
    }

    [Test]
    public void Prepare_WidthIsByteUpTo256DistinctUnits_EvenWhenTheUnitsAreAbove255()
    {
        string s256 = new(StringFamilies.Alphabet(256));
        Assert.That(s256.All(c => c > 255), Is.True);
        var p256 = TextPreparation.Prepare(s256);
        Assert.That(p256.AlphabetSize, Is.EqualTo(256));
        Assert.That(p256.Width, Is.EqualTo(SymbolWidth.Byte));
        Assert.That(p256.Bytes.ToArray(), Is.EqualTo(Enumerable.Range(0, 256).Select(i => (byte)i).ToArray()));

        var p257 = TextPreparation.Prepare(new string(StringFamilies.Alphabet(257)));
        Assert.That(p257.AlphabetSize, Is.EqualTo(257));
        Assert.That(p257.Width, Is.EqualTo(SymbolWidth.UShort));

        var pAll = TextPreparation.Prepare(StringFamilies.AllDistinct(65536, 3));
        Assert.That(pAll.AlphabetSize, Is.EqualTo(65536));
        Assert.That(pAll.Width, Is.EqualTo(SymbolWidth.UShort));
        Assert.That(pAll.UShorts.ToArray().OrderBy(x => x), Is.EqualTo(Enumerable.Range(0, 65536).Select(i => (ushort)i)));
    }

    [Test]
    public void Readers_ExposeLengthAndSymbols()
    {
        var b = new ByteSymbols(new byte[] { 3, 1, 2 });
        Assert.That(b.Length, Is.EqualTo(3));
        Assert.That(b[0], Is.EqualTo(3));
        Assert.That(b[2], Is.EqualTo(2));
        Assert.Throws<IndexOutOfRangeException>(() => _ = new ByteSymbols(new byte[] { 3 })[1]);

        var u = new UShortSymbols(new ushort[] { 65535, 0 });
        Assert.That(u.Length, Is.EqualTo(2));
        Assert.That(u[0], Is.EqualTo(65535));

        var i = new IntSymbols(new[] { int.MaxValue, 7 });
        Assert.That(i.Length, Is.EqualTo(2));
        Assert.That(i[0], Is.EqualTo(int.MaxValue));
        Assert.That(i[1], Is.EqualTo(7));

        var p = TextPreparation.Prepare("abca");
        Assert.That(p.ByteReader().Length, Is.EqualTo(4));
    }

    [Test]
    public void TextLengthGuard_AllowsUpToArrayMaxLengthMinusOne_AndThrowsNotSupportedFromThere()
    {
        Assert.DoesNotThrow(() => EsaLimits.EnsureTextLength(0));
        Assert.DoesNotThrow(() => EsaLimits.EnsureTextLength(Array.MaxLength - 1)); // n + 1 == Array.MaxLength
        var ex = Assert.Throws<NotSupportedException>(() => EsaLimits.EnsureTextLength(Array.MaxLength));
        Assert.That(ex!.Message, Does.Contain("64-bit"));
        Assert.Throws<NotSupportedException>(() => EsaLimits.EnsureTextLength(int.MaxValue));
        Assert.Throws<NotSupportedException>(() => EsaLimits.EnsureTextLength(long.MaxValue));
        Assert.Throws<ArgumentOutOfRangeException>(() => EsaLimits.EnsureTextLength(-1));
    }

    [Test]
    public void ReaderGuard_UsesOnlyTheLength_NeverTouchesTheSymbols()
    {
        Assert.Throws<NotSupportedException>(() => EsaLimits.EnsureReaderFits(new FakeReader(int.MaxValue)));
        Assert.Throws<NotSupportedException>(() => EsaLimits.EnsureReaderFits(new FakeReader(Array.MaxLength)));
        Assert.DoesNotThrow(() => EsaLimits.EnsureReaderFits(new FakeReader(Array.MaxLength - 1)));
        Assert.DoesNotThrow(() => EsaLimits.EnsureReaderFits(new FakeReader(0)));
    }

    private readonly struct FakeReader(int length) : ISymbolReader
    {
        public int Length { get; } = length;

        public int this[int index] => throw new InvalidOperationException("the size guard must not read symbols");
    }
}
