using System.Runtime.InteropServices;
using System.Text;
using SuffixTree.Esa.External;

namespace SuffixTree.Esa.Tests.External;

/// <summary>CRC-32C of the section and header checksums (ADR-003 section 6): known values, the bit-by-bit definition, chaining, and the hardware and table paths against each other.</summary>
[TestFixture]
public sealed unsafe class Crc32CTests
{
    private static readonly byte[] CheckString = Encoding.ASCII.GetBytes("123456789");

    [Test]
    public void TheCheckValueOfTheCatalogue_IsE3069283()
    {
        Assert.That(Crc32C.Compute(CheckString), Is.EqualTo(0xE3069283u));
        Assert.That(Crc32C.AppendTable(0, CheckString), Is.EqualTo(0xE3069283u));
        if (Crc32C.HardwareAccelerated) Assert.That(Crc32C.AppendHardware(0, CheckString), Is.EqualTo(0xE3069283u));
    }

    [Test]
    public void TheValuesOfRfc3720_AreReproduced()
    {
        // iSCSI (RFC 3720 appendix B.4): 32 bytes of zeros, of ones, ascending and descending
        var zeros = new byte[32];
        byte[] ones = [.. Enumerable.Repeat((byte)0xFF, 32)];
        byte[] up = [.. Enumerable.Range(0, 32).Select(i => (byte)i)];
        byte[] down = [.. Enumerable.Range(0, 32).Select(i => (byte)(31 - i))];
        Assert.That(Crc32C.Compute(zeros), Is.EqualTo(0x8A9136AAu));
        Assert.That(Crc32C.Compute(ones), Is.EqualTo(0x62A8AB43u));
        Assert.That(Crc32C.Compute(up), Is.EqualTo(0x46DD794Eu));
        Assert.That(Crc32C.Compute(down), Is.EqualTo(0x113FDB5Cu));
    }

    [Test]
    public void TheEmptyInput_HasCrcZero_AndAppendingNothingKeepsTheValue()
    {
        Assert.That(Crc32C.Compute([]), Is.Zero);
        Assert.That(Crc32C.Append(0xDEADBEEF, []), Is.EqualTo(0xDEADBEEFu));
        Assert.That(Crc32C.AppendTable(0xDEADBEEF, []), Is.EqualTo(0xDEADBEEFu));
    }

    [Test]
    public void EveryLength_AtEveryAlignment_EqualsTheBitwiseDefinition_OnBothPaths()
    {
        var rng = new Random(7);
        var buffer = new byte[400];
        rng.NextBytes(buffer);
        for (int start = 0; start < 9; start++)
        {
            for (int length = 0; length <= 300; length++)
            {
                ReadOnlySpan<byte> data = buffer.AsSpan(start, length);
                uint expected = SesaOracle.Crc32cBitwise(data);
                Assert.That(Crc32C.Append(0, data), Is.EqualTo(expected), $"dispatch, start {start}, length {length}");
                Assert.That(Crc32C.AppendTable(0, data), Is.EqualTo(expected), $"table, start {start}, length {length}");
                if (Crc32C.HardwareAccelerated) Assert.That(Crc32C.AppendHardware(0, data), Is.EqualTo(expected), $"hardware, start {start}, length {length}");
            }
        }
    }

    [Test]
    public void ChunkedAppends_EqualTheOneShotValue_ForEverySplit()
    {
        var data = new byte[257];
        new Random(11).NextBytes(data);
        uint whole = SesaOracle.Crc32cBitwise(data);
        for (int cut = 0; cut <= data.Length; cut++)
        {
            uint crc = Crc32C.Append(0, data.AsSpan(0, cut));
            Assert.That(Crc32C.Append(crc, data.AsSpan(cut)), Is.EqualTo(whole), $"dispatch, cut {cut}");
            uint t = Crc32C.AppendTable(0, data.AsSpan(0, cut));
            Assert.That(Crc32C.AppendTable(t, data.AsSpan(cut)), Is.EqualTo(whole), $"table, cut {cut}");
        }
    }

    [Test]
    public void ABigRandomBuffer_GivesTheSameValueOnBothPaths_AndThroughAPointer()
    {
        var data = new byte[(3 << 20) + 5];
        new Random(13).NextBytes(data);
        uint table = Crc32C.AppendTable(0, data);
        Assert.That(Crc32C.Append(0, data), Is.EqualTo(table));
        if (Crc32C.HardwareAccelerated) Assert.That(Crc32C.AppendHardware(0, data), Is.EqualTo(table));
        fixed (byte* p = data)
        {
            Assert.That(Crc32C.Append(0, p, data.Length), Is.EqualTo(table));
            Assert.That(Crc32C.Append(Crc32C.Append(0, p, 1000), p + 1000, data.Length - 1000), Is.EqualTo(table), "a pointer, in two pieces");
        }
        Assert.That(Crc32C.Append(0, null, 0), Is.Zero, "a null pointer of length 0 is the empty input");
    }

    [Test]
    public void ASingleBitFlip_ChangesTheValue_ForEveryBitOfAShortBuffer()
    {
        var data = new byte[40];
        new Random(17).NextBytes(data);
        uint reference = Crc32C.Compute(data);
        for (int bit = 0; bit < data.Length * 8; bit++)
        {
            data[bit / 8] ^= (byte)(1 << (bit % 8));
            Assert.That(Crc32C.Compute(data), Is.Not.EqualTo(reference), $"bit {bit}");
            data[bit / 8] ^= (byte)(1 << (bit % 8));
        }
    }

    [Test]
    public void TheHardwareFlag_AgreesWithTheCpu_AndTheHardwarePathRefusesWithoutIt()
    {
        bool expected = System.Runtime.Intrinsics.X86.Sse42.X64.IsSupported || System.Runtime.Intrinsics.Arm.Crc32.Arm64.IsSupported;
        Assert.That(Crc32C.HardwareAccelerated, Is.EqualTo(expected));
        if (!expected) Assert.Throws<PlatformNotSupportedException>(() => Crc32C.AppendHardware(0, CheckString));
    }

    [Test]
    public void ALongSequenceOfWordsAsSpans_MatchesMemoryMarshalViews()
    {
        uint[] words = [.. Enumerable.Range(0, 1000).Select(i => (uint)(i * 2654435761u))];
        ReadOnlySpan<byte> bytes = MemoryMarshal.AsBytes(words.AsSpan());
        Assert.That(Crc32C.Compute(bytes), Is.EqualTo(SesaOracle.Crc32cBitwise(bytes)));
    }
}
