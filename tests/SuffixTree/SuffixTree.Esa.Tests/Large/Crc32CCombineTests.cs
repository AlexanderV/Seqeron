using SuffixTree.Esa.External;
using SuffixTree.Esa.Tests.External;
using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.Esa.Tests.Large;

/// <summary>
/// Work package S2.1, deliverable B (the CRC of a hole): the CRC-32C of N zero bytes by the algebra of the polynomial (square and multiply, <see cref="Crc32CCombine"/>), not by a pass over N bytes. Checked against the
/// definition (bit by bit), against a direct pass of the product (<see cref="Crc32C"/>) for N up to 64 MiB, against the product chained over random zero lengths, and, for the sizes above 2^31 where a pass is too slow for the
/// default set, against values that the <c>Large</c> category computed once by the product's own pass over 2^31, 2^32 and 2^33 zero bytes (see <c>LargeSparseContainerTests</c>, which recomputes them).
/// </summary>
[TestFixture]
public class Crc32CCombineTests
{
    /// <summary>CRC-32C of 2^31, 2^32, 2^33, 2^32 + 10 and 2^31 + 7 zero bytes, by the product's pass (the <c>Large</c> test <c>PinnedZeroCrcs_AreWhatAPassOfTheProductGives</c> recomputes them).</summary>
    private static readonly (long Length, uint Crc)[] Pinned =
    [
        (1L << 31, 0x527D5351u),
        (1L << 32, 0xF16177D2u),
        (1L << 33, 0x48674BC7u),
        ((1L << 32) + 10, 0x2B60B55Du),
        ((1L << 31) + 7, 0x8C28B28Au),
    ];

    public static IEnumerable<long> SmallLengths()
    {
        for (long n = 0; n <= 70; n++) yield return n;
        foreach (long n in new long[] { 127, 128, 129, 255, 256, 257, 1023, 1024, 1025, 4095, 4096, 4097, 65_535, 65_536, 65_537, 1_000_003, 1L << 20, (1L << 24) + 5 }) yield return n;
    }

    public static IEnumerable<long> ShortLengths() => SmallLengths().Where(n => n <= 1_100_000);

    [TestCaseSource(nameof(ShortLengths))]
    public void OfZeros_IsTheDefinitionOfTheCrc(long length)
    {
        Assert.That(Crc32CCombine.OfZeros(length), Is.EqualTo(SesaOracle.Crc32cBitwise(new byte[length])));
    }

    [TestCaseSource(nameof(SmallLengths))]
    public void OfZeros_IsADirectPassOfTheProduct(long length) => Assert.That(Crc32CCombine.OfZeros(length), Is.EqualTo(Crc32C.Compute(new byte[length])));

    [Test]
    public void OfZeros_IsADirectPass_UpTo64MiB()
    {
        byte[] zeros = new byte[64 << 20];
        foreach (long length in new long[] { (1L << 26) - 1, 1L << 26, 12_345_678, 33_554_433, 50_000_000 })
        {
            Assert.That(Crc32CCombine.OfZeros(length), Is.EqualTo(Crc32C.Compute(zeros.AsSpan(0, (int)length))), $"{length} zero bytes");
        }
    }

    [Test]
    public void OfZeros_EqualsTheProductChainedOverRandomZeroLengths()
    {
        byte[] zeros = new byte[1 << 22];
        var rng = new SplitMix64(0xC2C32C);
        for (int c = 0; c < 24; c++)
        {
            long length = (long)(rng.Next() % (ulong)(100L << 20)); // up to 100 MiB, in pieces of 4 MiB
            uint crc = 0;
            for (long done = 0; done < length; done += zeros.Length) crc = Crc32C.Append(crc, zeros.AsSpan(0, (int)Math.Min(zeros.Length, length - done)));
            Assert.That(Crc32CCombine.OfZeros(length), Is.EqualTo(crc), $"{length} zero bytes");
        }
    }

    [Test]
    public void OfZeros_AboveTwoToThe31_IsWhatAPassOfTheProductGave()
    {
        foreach ((long length, uint crc) in Pinned)
        {
            Assert.That(Crc32CCombine.OfZeros(length), Is.EqualTo(crc), $"{length} zero bytes");
        }
    }

    [Test]
    public void OfZeros_IsAdditive_AtEverySizeThatMatters()
    {
        long[] lengths = [0, 1, 7, 4096, 1L << 20, (1L << 31) - 1, 1L << 31, (1L << 31) + 10, (1L << 32) - 1, 1L << 32, (1L << 32) + 10, 1L << 33, (1L << 40) + 3, (1L << 55) + 1];
        foreach (long a in lengths)
            foreach (long b in lengths)
            {
                long sum = a + b;
                Assert.That(Crc32CCombine.Combine(Crc32CCombine.OfZeros(a), Crc32CCombine.OfZeros(b), b), Is.EqualTo(Crc32CCombine.OfZeros(sum)), $"{a} + {b} zero bytes");
            }
    }

    [Test]
    public void Combine_JoinsTheCrcsOfTwoPiecesOfData()
    {
        var rng = new SplitMix64(0xC0B1E);
        byte[] data = new byte[1 << 16];
        for (int i = 0; i < data.Length; i++) data[i] = (byte)rng.Next();
        for (int c = 0; c < 200; c++)
        {
            int cut = rng.NextInt(data.Length + 1);
            uint whole = Crc32C.Compute(data);
            uint a = Crc32C.Compute(data.AsSpan(0, cut)), b = Crc32C.Compute(data.AsSpan(cut));
            Assert.That(Crc32CCombine.Combine(a, b, data.Length - cut), Is.EqualTo(whole), $"cut at {cut}");
        }
        Assert.That(Crc32CCombine.Combine(0x1234ABCD, 0, 0), Is.EqualTo(0x1234ABCDu), "an empty second piece");
        Assert.That(Crc32CCombine.Combine(0, Crc32C.Compute(data), data.Length), Is.EqualTo(Crc32C.Compute(data)), "an empty first piece");
    }

    [Test]
    public void Multiply_IsTheProductOfPolynomials()
    {
        // x^a * x^b = x^(a + b) for a few exponents, one is x^0 = 1, zero annihilates
        foreach ((ulong a, ulong b) in new (ulong, ulong)[] { (0, 0), (1, 1), (3, 5), (31, 1), (32, 0), (255, 1000), (1UL << 20, 1UL << 21), (1UL << 40, 12345) })
            Assert.That(Crc32CCombine.Multiply(Crc32CCombine.XPower(a), Crc32CCombine.XPower(b)), Is.EqualTo(Crc32CCombine.XPower(a + b)), $"x^{a} * x^{b}");
        Assert.That(Crc32CCombine.XPower(0), Is.EqualTo(0x80000000u), "x^0 is the unit");
        Assert.That(Crc32CCombine.XPower(1), Is.EqualTo(0x40000000u), "x^1");
        Assert.That(Crc32CCombine.Multiply(0, 0xDEADBEEF), Is.EqualTo(0u));
        Assert.That(Crc32CCombine.Multiply(0xDEADBEEF, 0), Is.EqualTo(0u));
        Assert.That(Crc32CCombine.Multiply(0x80000000, 0xDEADBEEF), Is.EqualTo(0xDEADBEEFu));
        // x^32 = the reflected polynomial, by the definition of the modulus
        Assert.That(Crc32CCombine.XPower(32), Is.EqualTo(0x82F63B78u));
    }

    [Test]
    public void OfSparse_IsTheCrcOfTheSectionWithItsPatches()
    {
        var rng = new SplitMix64(0x5A55E);
        for (int trial = 0; trial < 300; trial++)
        {
            int length = rng.NextInt(70_000);
            var section = new byte[length];
            var patches = new List<(long Offset, byte[] Data)>();
            int at = 0;
            while (at < length && rng.NextInt(4) != 0)
            {
                int gap = rng.NextInt(Math.Max(1, (length - at) / 2) + 1);
                int size = rng.NextInt(Math.Min(300, length - at - gap) + 1);
                if (size == 0 || at + gap + size > length) break;
                var data = new byte[size];
                for (int i = 0; i < size; i++) data[i] = (byte)rng.Next();
                data.CopyTo(section, at + gap);
                patches.Add((at + gap, data));
                at += gap + size;
            }
            Assert.That(Crc32CCombine.OfSparse(length, patches), Is.EqualTo(Crc32C.Compute(section)), $"trial {trial}: {length} bytes, {patches.Count} patches");
        }
        // a patch at the very start and at the very end, back to back, and none
        byte[] a = [1, 2, 3], b = [9, 9];
        byte[] expected = new byte[1000];
        a.CopyTo(expected, 0);
        b.CopyTo(expected, 998);
        Assert.That(Crc32CCombine.OfSparse(1000, [(0, a), (998, b)]), Is.EqualTo(Crc32C.Compute(expected)));
        Assert.That(Crc32CCombine.OfSparse(1000, []), Is.EqualTo(Crc32C.Compute(new byte[1000])));
        Assert.That(Crc32CCombine.OfSparse(0, []), Is.EqualTo(0u));
    }

    [Test]
    public void OfSparse_AtTheSizeOfTheLargeSections_IsTheCrcOfTheZerosAroundThePatches()
    {
        // a patch of 5 bytes at the offset 2^32 + 1 of a section of 2^33 + 3 bytes: the same as the combine of its three pieces
        byte[] patch = [1, 2, 3, 4, 5];
        long length = (1L << 33) + 3, offset = (1L << 32) + 1;
        uint crc = Crc32CCombine.OfSparse(length, [(offset, patch)]);
        uint expected = Crc32CCombine.Combine(Crc32CCombine.Combine(Crc32CCombine.OfZeros(offset), Crc32C.Compute(patch), patch.Length), Crc32CCombine.OfZeros(length - offset - patch.Length), length - offset - patch.Length);
        Assert.That(crc, Is.EqualTo(expected));
        Assert.That(crc, Is.Not.EqualTo(Crc32CCombine.OfZeros(length)), "the patch changes the CRC");
    }

    [Test]
    public void OfSparse_RefusesPatchesThatOverlapAreUnsortedOrOutside()
    {
        byte[] a = [1, 2, 3];
        Assert.Throws<ArgumentException>(() => Crc32CCombine.OfSparse(100, [(10, a), (11, a)]), "overlap");
        Assert.Throws<ArgumentException>(() => Crc32CCombine.OfSparse(100, [(20, a), (10, a)]), "unsorted");
        Assert.Throws<ArgumentException>(() => Crc32CCombine.OfSparse(100, [(98, a)]), "past the end");
        Assert.Throws<ArgumentException>(() => Crc32CCombine.OfSparse(100, [(-1, a)]), "before the start");
    }
}
