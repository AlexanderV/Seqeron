using System.Diagnostics;
using SuffixTree.Esa.External;

namespace SuffixTree.Esa.Tests.Large;

/// <summary>
/// The CRC of a hole above 2^31 bytes against a real pass of the product's CRC over that many zero bytes (the pass takes 0.3 s per 2^30 bytes with the instruction, 4 s for 2^33 bytes). Category <c>Large</c>, explicit:
/// <c>dotnet test tests/SuffixTree/SuffixTree.Esa.Tests/SuffixTree.Esa.Tests.csproj -c Release --filter Category=Large</c>. The default set pins the values this test recomputes (<c>Crc32CCombineTests.Pinned</c>).
/// </summary>
[TestFixture]
[Category("Large")]
[Explicit("Heavy: a pass of the product's CRC-32C over up to 2^33 zero bytes; run with --filter Category=Large (see tests/SuffixTree/SuffixTree.Esa.Tests/Large/README.md)")]
public class LargeCrcTests
{
    [TestCase(1L << 31)]
    [TestCase(1L << 32)]
    [TestCase((1L << 32) + 10)]
    [TestCase(1L << 33)]
    [TestCase((1L << 31) + 7)]
    public void PinnedZeroCrcs_AreWhatAPassOfTheProductGives(long length)
    {
        byte[] zeros = new byte[1 << 26];
        var clock = Stopwatch.StartNew();
        uint pass = 0;
        for (long done = 0; done < length; done += zeros.Length) pass = Crc32C.Append(pass, zeros.AsSpan(0, (int)Math.Min(zeros.Length, length - done)));
        double passSeconds = clock.Elapsed.TotalSeconds;
        clock.Restart();
        uint combine = Crc32CCombine.OfZeros(length);
        TestContext.Out.WriteLine($"CRC-OF-ZEROS length={length} pass=0x{pass:X8} ({passSeconds:F2} s) combine=0x{combine:X8} ({clock.Elapsed.TotalMilliseconds:F3} ms)");
        Assert.That(combine, Is.EqualTo(pass));
    }
}
