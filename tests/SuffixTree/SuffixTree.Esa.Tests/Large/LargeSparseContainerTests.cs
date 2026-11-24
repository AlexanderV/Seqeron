using System.Diagnostics;
using SuffixTree.Esa.External;
using SuffixTree.Esa.Tests.External;

namespace SuffixTree.Esa.Tests.Large;

/// <summary>
/// Work package S2.1, deliverable B at the sizes that matter (category <c>Large</c>, explicit): the sparse container writer makes a container of 2^31 + 10 and 2^32 + 10 symbols that takes a few MB of disk (a hole is not
/// written: <c>SetLength</c>, the header page and the patches), whose header page and section table are valid for the width, and whose sections carry the CRC-32C of their bytes by the combine, which the S1 reader (it can
/// open a container below 2^31 - 1 symbols) confirms by a real pass over the holes. Run: <c>dotnet test tests/SuffixTree/SuffixTree.Esa.Tests/SuffixTree.Esa.Tests.csproj -c Release --filter FullyQualifiedName~LargeSparseContainerTests</c>.
/// </summary>
[TestFixture]
[Category("Large")]
[Explicit(ExplicitReason)]
public class LargeSparseContainerTests : LargeFixtureBase
{
    private const long TwoTo31 = 1L << 31;
    private const long TwoTo32 = 1L << 32;

    [TestCase(TwoTo31 + 10)]
    [TestCase(TwoTo32 + 10)]
    [TestCase(TwoTo32 - 4)]
    [TestCase(TwoTo32 - 3)]
    public void TheWriter_MakesAHeaderThatDescribesTheFile_AndTakesNoDisk(long n)
    {
        ClosedFormText f = new PowerOfA(n);
        SparseSpec spec = SparseSpec.For(f);
        SparseLayout layout = SparseContainer.Plan(spec);
        var clock = Stopwatch.StartNew();
        string path = WriteContainer("shell", spec, SparseContainer.EndPatches(f, layout, 4096));
        double seconds = clock.Elapsed.TotalSeconds;
        long allocated = LargeDisk.AllocatedBytes(path);
        TestContext.Out.WriteLine($"SPARSE-FILE n={n} width={layout.PosWidth} lcp_codec={layout.LcpCodec} nominal_bytes={layout.FileSize} disk_bytes={allocated} create_seconds={seconds:F3}");

        Assert.That(new FileInfo(path).Length, Is.EqualTo(layout.FileSize), "the file has its nominal size");
        Assert.That(allocated, Is.InRange(0, 4L << 20), "a few blocks: the header page and the patches, never the sections");
        byte[] page = new byte[SesaFormat.PageBytes];
        using (FileStream stream = File.OpenRead(path)) stream.ReadExactly(page);
        SesaOracle o = SesaOracle.FromBytes(page);
        Assert.That((o.TextLen, o.NodeCount, o.FileSize, o.PosWidth, o.Sigma), Is.EqualTo((n, f.NodeCount, layout.FileSize, layout.PosWidth, 1)));
        Assert.That(o.HeaderCrc, Is.EqualTo(o.ExpectedHeaderCrc()), "the header CRC-32C");
        Assert.That((o.VersionMajor, o.VersionMinor), Is.EqualTo(layout.PosWidth == 5 ? (1, 1) : (1, 0)));
        foreach (OracleSection s in o.Sections)
        {
            Assert.That(s.Offset + s.Length, Is.LessThanOrEqualTo(layout.FileSize), $"kind {s.Kind}: inside the file");
            Assert.That(s.Offset % (s.Length >= (1 << 20) ? 4096 : 64), Is.Zero, $"kind {s.Kind}: aligned");
            // the CRC of a hole is the CRC of its bytes: the unpatched middle is zero, the patches are folded in
            Assert.That(s.HasCrc, Is.True, $"kind {s.Kind}");
        }
        OracleSection sa = o.Get(SesaFormat.KindSa);
        Assert.That((sa.Count, sa.ElemWidth, sa.Length), Is.EqualTo((n + 1, layout.PosWidth, ((n + 1) * layout.PosWidth) + (layout.PosWidth == 5 ? 3 : 0))));
    }

    /// <summary>The S1 reader opens a container below 2^31 - 1 symbols: its pass over the holes confirms the CRC values that the combine made, at 0.27 and 8.6 GB per section.</summary>
    [TestCase(1L << 24)]
    [TestCase(1L << 28)]
    [TestCase(TwoTo31 - 3)]
    public void TheShellBelowTwoToThe31_OpensAtTheSectionsLevelOfTheS1Reader_ByARealPassOverTheHoles(long n)
    {
        ClosedFormText f = new PowerOfA(n);
        SparseSpec spec = SparseSpec.For(f);
        string path = WriteContainer("below31", spec, SparseContainer.EndPatches(f, SparseContainer.Plan(spec), 4096));
        var clock = Stopwatch.StartNew();
        using SesaFile file = SesaFile.Open(path, new EsaOpenOptions { Validation = EsaValidationLevel.Sections, Access = EsaAccessMode.Checked });
        double seconds = clock.Elapsed.TotalSeconds;
        long bytes = file.Report.SectionBytesChecksummed;
        TestContext.Out.WriteLine($"SECTIONS-PASS n={n} section_bytes={bytes} open_seconds={seconds:F3} header_seconds={file.Report.HeaderSeconds:F4} sections_seconds={file.Report.SectionsSeconds:F3} gb_per_second={bytes / 1e9 / file.Report.SectionsSeconds:F2}");
        Assert.That(bytes, Is.GreaterThan(3 * 4 * n), "SA, CHILD and LCP were read");
    }

    /// <summary>The CRC-32C of a section of the width-5 container that is a hole is the CRC of that many zero bytes, checked here by the product's own pass over a real file of zeros (a prefix of 8 GiB).</summary>
    [Test]
    public void AHoleOfEightGiB_HasTheCrcThatARealPassOverTheFileGives()
    {
        string path = NewPath("hole");
        const long length = 8L << 30;
        using (var stream = new FileStream(path, FileMode.CreateNew, FileAccess.Write, FileShare.None)) stream.SetLength(length);
        uint pass = 0;
        byte[] buffer = new byte[1 << 26];
        using (FileStream read = File.OpenRead(path))
        {
            int got;
            while ((got = read.Read(buffer, 0, buffer.Length)) > 0) pass = Crc32C.Append(pass, buffer.AsSpan(0, got));
        }
        Assert.That(Crc32CCombine.OfZeros(length), Is.EqualTo(pass));
    }
}
