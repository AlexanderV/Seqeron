using SuffixTree.Esa.External;

namespace SuffixTree.Esa.Tests.External;

/// <summary>The layout arithmetic of TEXT codec 4 (ADR-003 kind 1, ADR-004 PH0): the width, the byte counts, the packing of one block.</summary>
[TestFixture]
public class PackedTextFormatTests
{
    [Test]
    public void TheConstants_AreTheOnesOfTheLayout()
    {
        Assert.That(PackedTextFormat.CodecId, Is.EqualTo(4));
        Assert.That(PackedTextFormat.PadBytes, Is.EqualTo(64));
        Assert.That(PackedTextFormat.MinBits, Is.EqualTo(1));
        Assert.That(PackedTextFormat.MaxBits, Is.EqualTo(8));
        Assert.That(PackedTextFormat.MaxSymbols, Is.EqualTo(4_294_967_292L), "n <= 2^32 - 4");
    }

    [TestCase(0, 1)]
    [TestCase(1, 1)]
    [TestCase(2, 1)]
    [TestCase(3, 2)]
    [TestCase(4, 2)]
    [TestCase(5, 3)]
    [TestCase(8, 3)]
    [TestCase(9, 4)]
    [TestCase(16, 4)]
    [TestCase(17, 5)]
    [TestCase(32, 5)]
    [TestCase(33, 6)]
    [TestCase(64, 6)]
    [TestCase(65, 7)]
    [TestCase(128, 7)]
    [TestCase(129, 8)]
    [TestCase(255, 8)]
    [TestCase(256, 8)]
    [TestCase(257, 9)]
    [TestCase(65536, 16)]
    [TestCase(65537, 17)]
    public void BitsFor_IsTheSmallestWidthThatHoldsSigma_AndAtLeastOne(long sigma, int bits) => Assert.That(PackedTextFormat.BitsFor(sigma), Is.EqualTo(bits));

    [Test]
    public void BitsFor_AtTheEveryPowerOfTwoAndOneAbove()
    {
        for (int b = 1; b <= 40; b++)
        {
            Assert.That(PackedTextFormat.BitsFor(1L << b), Is.EqualTo(b), $"sigma = 2^{b}");
            Assert.That(PackedTextFormat.BitsFor((1L << b) + 1), Is.EqualTo(b + 1), $"sigma = 2^{b} + 1");
        }
    }

    [Test]
    public void DataBytes_AndTotalBytes_RoundUpToWholeBytes_AndAddThePadding()
    {
        for (int bits = 1; bits <= 8; bits++)
        {
            for (long n = 0; n <= 200; n++)
            {
                long expected = ((n * bits) + 7) / 8;
                Assert.That(PackedTextFormat.DataBytes(n, bits), Is.EqualTo(expected), $"n={n} b={bits}");
                Assert.That(PackedTextFormat.TotalBytes(n, bits), Is.EqualTo(expected + 64), $"n={n} b={bits}");
            }
        }
        Assert.That(PackedTextFormat.DataBytes(4_294_967_292L, 8), Is.EqualTo(4_294_967_292L));
        Assert.That(PackedTextFormat.DataBytes(4_294_967_292L, 2), Is.EqualTo(1_073_741_823L));
        Assert.That(PackedTextFormat.DataBytes(3_000_000_001L, 3), Is.EqualTo(1_125_000_001L), "above 2^31 symbols the bit count exceeds 2^32");
    }

    [Test]
    public void PackInto_EqualsTheBitByBitOracle_ForEveryWidthAndLength()
    {
        for (int bits = 1; bits <= 8; bits++)
        {
            for (int n = 0; n <= 70; n++)
            {
                byte[] s = TextOracle.Random(n, 1 << bits, (ulong)((bits * 1000) + n));
                byte[] expected = TextOracle.Pack(s, bits);
                var dst = new byte[PackedTextFormat.TotalBytes(n, bits)];
                PackedTextFormat.PackInto(s, bits, dst);
                Assert.That(dst, Is.EqualTo(expected), $"n={n} b={bits}");
            }
        }
    }

    [Test]
    public void PackInto_WritesTheWholeLastByte_EvenOnDirtyMemory_AndNothingBeyondTheData()
    {
        foreach (int bits in new[] { 1, 2, 3, 5, 7 })
        {
            byte[] s = TextOracle.Random(13, 1 << bits, 9);
            long data = PackedTextFormat.DataBytes(13, bits);
            var dst = new byte[data + 4];
            Array.Fill(dst, (byte)0xFF);
            PackedTextFormat.PackInto(s, bits, dst);
            Assert.That(dst[..(int)data], Is.EqualTo(TextOracle.Pack(s, bits)[..(int)data]), $"b={bits}: unused bits are zero");
            Assert.That(dst[(int)data..], Is.All.EqualTo((byte)0xFF), $"b={bits}: nothing beyond the data is touched");
        }
    }

    [Test]
    public void PackInto_OfBlocksOfEightSymbols_ConcatenatesToThePackingOfTheWhole()
    {
        foreach (int bits in new[] { 2, 3, 4, 8 })
        {
            byte[] s = TextOracle.Random(8 * 25 + 3, 1 << bits, 5);
            byte[] whole = TextOracle.Pack(s, bits);
            var parts = new byte[whole.Length];
            long at = 0;
            for (int i = 0; i < s.Length; i += 8)
            {
                int take = Math.Min(8, s.Length - i);
                int bytes = (int)PackedTextFormat.DataBytes(take, bits);
                PackedTextFormat.PackInto(s.AsSpan(i, take), bits, parts.AsSpan((int)at, bytes));
                at += bytes;
            }
            Assert.That(parts, Is.EqualTo(whole), $"b={bits}: a block of 8 symbols ends on a byte boundary");
        }
    }

    [Test]
    public void PackInto_RefusesASymbolThatDoesNotFit_AndATooSmallDestination()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => PackedTextFormat.PackInto(new byte[] { 4 }, 2, new byte[8]));
        Assert.Throws<ArgumentOutOfRangeException>(() => PackedTextFormat.PackInto(new byte[] { 0, 8 }, 3, new byte[8]));
        PackedTextFormat.PackInto(new byte[] { 3, 7 }, 3, new byte[8]);
        Assert.Throws<ArgumentException>(() => PackedTextFormat.PackInto(new byte[8], 2, new byte[1]), "2 bytes needed");
        PackedTextFormat.PackInto(new byte[8], 2, new byte[2]);
        Assert.Throws<ArgumentOutOfRangeException>(() => PackedTextFormat.PackInto(new byte[1], 0, new byte[8]));
        Assert.Throws<ArgumentOutOfRangeException>(() => PackedTextFormat.PackInto(new byte[1], 9, new byte[8]));
    }

    [Test]
    public void BytesPerSymbol_Eight_IsTheIdentity_AndTwoBitsPackFourToTheByte()
    {
        byte[] s = TextOracle.Random(100, 256, 3);
        var dst = new byte[100];
        PackedTextFormat.PackInto(s, 8, dst);
        Assert.That(dst, Is.EqualTo(s));
        var four = new byte[] { 0, 1, 2, 3 };
        var one = new byte[1];
        PackedTextFormat.PackInto(four, 2, one);
        Assert.That(one[0], Is.EqualTo(0b00_01_10_11));
    }

    [Test]
    public void PackedTextInfo_DerivesSigmaAndSizes()
    {
        var info = new PackedTextInfo(1000, 3, [1, 2, 3, 4, 5], [200, 200, 200, 200, 200], "00");
        Assert.That(info.Sigma, Is.EqualTo(5));
        Assert.That(info.DataBytes, Is.EqualTo(375));
        Assert.That(info.TotalBytes, Is.EqualTo(375 + 64));
    }
}
