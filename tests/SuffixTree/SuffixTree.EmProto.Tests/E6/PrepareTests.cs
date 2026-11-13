using SplitMix64 = SuffixTree.Esa.Tests.TestSupport.SplitMix64;
using System.Security.Cryptography;
using SuffixTree.EmProto.Dcs;
using SuffixTree.EmProto.Harness;
using SuffixTree.EmProto.Pipeline;
using SuffixTree.EmProto.Text;
using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.EmProto.Tests.WholePipeline;

/// <summary>PH0 of ADR-004 section 4.2: alphabet scan, symbol counts, SHA-256, the packed text streamed into its file in blocks.</summary>
[TestFixture]
public class PrepareTests : TempDirFixture
{
    private string Write(byte[] data)
    {
        string p = Path.Combine(Dir, "in-" + Guid.NewGuid().ToString("N")[..6] + ".txt");
        File.WriteAllBytes(p, data);
        return p;
    }

    private static byte[] Random(int n, int sigma, ulong seed, int offset = 'a')
    {
        var rng = new SplitMix64(seed);
        var b = new byte[n];
        for (int i = 0; i < n; i++) b[i] = (byte)(offset + rng.NextInt(sigma));
        return b;
    }

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
    [TestCase(129, 8)]
    [TestCase(256, 8)]
    public void TheWidth_IsTheSmallestBThatHoldsSigma_AndAtLeastOne(int sigma, int bits)
    {
        var data = new byte[3 * sigma];
        for (int i = 0; i < data.Length; i++) data[i] = (byte)((i * 7 % sigma) + (sigma == 256 ? 0 : 3));
        for (int i = 0; i < sigma; i++) data[i] = (byte)(i + (sigma == 256 ? 0 : 3)); // every letter occurs
        PrepareResult r = Prepare.Run(Write(data), Path.Combine(Dir, "p.pk"), 64, 128);
        Assert.That(r.Sigma, Is.EqualTo(sigma));
        Assert.That(r.Bits, Is.EqualTo(bits));
    }

    [Test]
    public void PackedFile_EqualsTheBitByBitReference_AndReadsBackAsTheDenseSymbols()
    {
        byte[] data = Random(1003, 5, 2, 'c');
        string pk = Path.Combine(Dir, "p.pk");
        PrepareResult r = Prepare.Run(Write(data), pk, 24, 100);
        char[] letters = [.. data.Select(b => (char)b).Distinct().OrderBy(c => c)];
        byte[] dense = [.. data.Select(b => (byte)Array.BinarySearch(letters, (char)b))];
        Assert.That(File.ReadAllBytes(pk), Is.EqualTo(PipelineChecks.ReferencePacked(dense, letters.Length)));
        Assert.That(r.Letters, Is.EqualTo(letters.Select(c => (byte)c).ToArray()), "increasing byte order");
        using var mapped = new MappedFile(pk, writable: false);
        unsafe
        {
            using var text = new PackedText(mapped.Pointer, data.Length, r.Bits, mapped, ownsOwner: false);
            for (int i = 0; i < data.Length; i++) Assert.That(text.Symbol(i), Is.EqualTo(dense[i]), $"symbol {i}");
        }
    }

    [Test]
    public void CountsAndHash_AreThoseOfTheInput()
    {
        byte[] data = Random(5000, 4, 7, 'A');
        PrepareResult r = Prepare.Run(Write(data), Path.Combine(Dir, "p.pk"), 1024, 4096);
        Assert.That(r.N, Is.EqualTo(5000));
        Assert.That(r.Counts, Is.EqualTo(r.Letters.Select(l => (long)data.Count(b => b == l)).ToArray()));
        Assert.That(r.Counts.Sum(), Is.EqualTo(5000));
        Assert.That(r.Sha256, Is.EqualTo(Convert.ToHexStringLower(SHA256.HashData(data))));
        string meta = r.MetaText("x.txt");
        Assert.That(meta, Does.StartWith("n=5000\nbits=2\nsigma=4\n"));
        Assert.That(meta, Does.Contain("letters_hex=").And.Contain("counts=").And.Contain($"sha256={r.Sha256}").And.Contain("name=x.txt"));
    }

    [Test]
    public void TheFile_DoesNotDependOnTheBlockSize()
    {
        byte[] data = Random(10_007, 3, 11);
        byte[]? first = null;
        foreach (int block in new[] { 8, 16, 24, 40, 1000, 1 << 20 })
        {
            string pk = Path.Combine(Dir, $"p{block}.pk");
            Prepare.Run(Write(data), pk, block, 256);
            byte[] bytes = File.ReadAllBytes(pk);
            first ??= bytes;
            Assert.That(bytes, Is.EqualTo(first), $"block of {block} symbols");
        }
    }

    [TestCase(0)]
    [TestCase(7)]
    [TestCase(9)]
    [TestCase(100)]
    public void ABlockThatIsNotAMultipleOfEight_IsRefused(int block)
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => Prepare.Run(Write("abc"u8.ToArray()), Path.Combine(Dir, "p.pk"), block, 64));
    }

    [Test]
    public void TheInput_IsReadTwice_TheOutputIsWrittenThroughMappedWindows_NothingIsWrittenExplicitly()
    {
        byte[] data = Random(100_000, 4, 5);
        IoCounters.Reset();
        MappedIo.Reset();
        Prepare.Run(Write(data), Path.Combine(Dir, "p.pk"), 4096, 1000);
        Assert.That(IoCounters.ReadBytes, Is.EqualTo(200_000), "two passes");
        Assert.That(IoCounters.WriteBytes, Is.Zero);
        Assert.That(MappedIo.WriteBytes, Is.EqualTo(25_000), "2 bits per symbol");
        Assert.That(MappedIo.MaxWindowBytes, Is.LessThanOrEqualTo(1000));
    }

    [Test]
    public void SingleSymbol_AndOneSymbolText()
    {
        string pk = Path.Combine(Dir, "p.pk");
        PrepareResult r = Prepare.Run(Write("aaaaaaaaaaaa"u8.ToArray()), pk, 8, 64);
        Assert.That((r.Sigma, r.Bits, r.N), Is.EqualTo((1, 1, 12L)));
        Assert.That(new FileInfo(pk).Length, Is.EqualTo(2 + 64));
        Assert.That(File.ReadAllBytes(pk), Is.All.EqualTo((byte)0), "the only symbol is 0");
        PrepareResult one = Prepare.Run(Write("z"u8.ToArray()), pk, 8, 64);
        Assert.That((one.Sigma, one.Bits, one.N), Is.EqualTo((1, 1, 1L)));
        Assert.That(new FileInfo(pk).Length, Is.EqualTo(1 + 64));
    }

    [Test]
    public void AllByteValues_IncludingNulAndFf_AreSymbolsInIncreasingByteOrder()
    {
        var data = new byte[600];
        for (int i = 0; i < data.Length; i++) data[i] = (byte)(i % 256);
        PrepareResult r = Prepare.Run(Write(data), Path.Combine(Dir, "p.pk"), 32, 64);
        Assert.That(r.Letters, Is.EqualTo(Enumerable.Range(0, 256).Select(i => (byte)i).ToArray()));
        Assert.That(r.Bits, Is.EqualTo(8));
        Assert.That(File.ReadAllBytes(Path.Combine(Dir, "p.pk")).Take(600), Is.EqualTo(data), "8 bits per symbol: the file is the text");
    }

    [Test]
    public void AnInputThatChangesBetweenThePasses_IsDetected()
    {
        byte[] data = Random(50_000, 4, 3);
        string input = Write(data);
        // the file is replaced by another one of the same length (and the same alphabet) after the alphabet scan and before the pack
        void Flip() => File.WriteAllBytes(input, Random(50_000, 4, 99));
        InvalidOperationException? e = Assert.Throws<InvalidOperationException>(() => Prepare.Run(input, Path.Combine(Dir, "p.pk"), 1024, 4096, Flip));
        Assert.That(e!.Message, Does.Contain("changed between the two passes"));
    }

    [Test]
    public void AnInputThatIsTruncatedBetweenThePasses_IsDetected()
    {
        string input = Write(Random(5000, 4, 3));
        void Cut() => File.WriteAllBytes(input, Random(4000, 4, 3));
        Assert.Throws<EndOfStreamException>(() => Prepare.Run(input, Path.Combine(Dir, "p.pk"), 1024, 4096, Cut));
    }
}
