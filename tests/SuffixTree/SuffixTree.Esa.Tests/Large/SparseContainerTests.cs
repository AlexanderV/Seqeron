using System.Buffers.Binary;
using System.Security.Cryptography;
using System.Text;
using SuffixTree.Esa.External;
using SuffixTree.Esa.Tests.External;

namespace SuffixTree.Esa.Tests.Large;

/// <summary>
/// Work package S2.1, deliverable B at small n (the default set): the sparse container writer is a writer of SESA containers, so with every entry given it must produce the bytes of the product's writer; its shells must open
/// at the <c>Header</c> and <c>Sections</c> levels of the S1 reader and be refused by <c>Structure</c>; its layout at width 5 is the one of ADR-005 section 6.1; its safety rules refuse what they promise to refuse.
/// The files that really are sparse and above 2^31 symbols are the tests of <see cref="LargeSparseContainerTests"/> (category <c>Large</c>).
/// </summary>
[TestFixture]
public class SparseContainerTests
{
    private static SparseSpec ZeroSpec(long n, int width = 0, int sigma = 1, bool psi = false) => new()
    {
        N = n,
        PosWidth = width,
        Sigma = sigma,
        NodeCount = n < 2 ? n + 1 : n + 2 + 1,
        MaxLcp = 0,
        Psi = psi,
    };

    public static IEnumerable<TestCaseData> LayoutShapes()
    {
        foreach (long n in new long[] { 0, 1, 2, 7, 100, 1000, 65_535, 300_000, 1_000_003 })
            foreach (int sigma in new[] { 1, 2, 4, 5, 20 })
                foreach (bool psi in new[] { false, true })
                    if (n > 0 || sigma == 1) yield return new TestCaseData(n, sigma, psi).SetArgDisplayNames($"n={n}", $"sigma={sigma}", $"psi={psi}");
    }

    [TestCaseSource(nameof(LayoutShapes))]
    public void Plan_AtWidth4_IsTheLayoutOfTheProductWriter(long n, int sigma, bool psi)
    {
        if (n == 0) sigma = 0;
        int bits = PackedTextFormat.BitsFor(Math.Max(1, sigma));
        var info = new PackedTextInfo(n, bits, new byte[sigma], new long[sigma], string.Empty);
        SesaPlan plan = SesaLayout.Plan(info, psi);
        SesaTail tail = SesaLayout.PlaceTail(plan, LcpKind.U8, 0);

        SparseLayout layout = SparseContainer.Plan(ZeroSpec(n, 4, sigma, psi));
        Assert.That(layout.Section(SesaFormat.KindText).Offset, Is.EqualTo(plan.TextOffset), "TEXT offset");
        Assert.That(layout.Section(SesaFormat.KindText).Length, Is.EqualTo(plan.TextLength), "TEXT length");
        Assert.That(layout.Section(SesaFormat.KindAlphabet).Offset, Is.EqualTo(plan.AlphabetOffset), "ALPHABET offset");
        Assert.That(layout.Section(SesaFormat.KindCounts).Offset, Is.EqualTo(plan.CountsOffset), "COUNTS offset");
        Assert.That(layout.Section(SesaFormat.KindSa).Offset, Is.EqualTo(plan.SaOffset), "SA offset");
        Assert.That(layout.Section(SesaFormat.KindSa).Length, Is.EqualTo(plan.SaLength), "SA length");
        Assert.That(layout.Section(SesaFormat.KindChild).Offset, Is.EqualTo(plan.ChildOffset), "CHILD offset");
        Assert.That(layout.Section(SesaFormat.KindChild).Length, Is.EqualTo(plan.ChildLength), "CHILD length");
        if (psi)
        {
            Assert.That(layout.Section(SesaFormat.KindPsi).Offset, Is.EqualTo(plan.PsiOffset), "PSI offset");
            Assert.That(layout.Section(SesaFormat.KindPsi).Length, Is.EqualTo(plan.PsiLength), "PSI length");
        }
        else
        {
            Assert.That(layout.Has(SesaFormat.KindPsi), Is.False, "no PSI");
        }
        Assert.That(layout.Section(SesaFormat.KindBuildInfo).Offset, Is.EqualTo(plan.BuildInfoOffset), "BUILD_INFO offset");
        Assert.That(layout.Section(SesaFormat.KindLcp).Offset, Is.EqualTo(plan.LcpOffset), "LCP offset");
        Assert.That(layout.FileSize, Is.EqualTo(tail.FinalLength), "file size");
    }

    public static IEnumerable<TestCaseData> RealContainers()
    {
        foreach (bool psi in new[] { false, true })
        {
            foreach (int n in new[] { 1, 2, 3, 10, 64, 100, 199 })
            {
                int nn = n;
                yield return new TestCaseData((Func<ClosedFormText>)(() => new PowerOfA(nn)), psi).SetArgDisplayNames($"a^{n}", $"psi={psi}");
            }
            foreach (int k in new[] { 1, 2, 3, 6 })
                foreach (int n in new[] { 1, 5, 7, 13, 40, 150, 199 })
                {
                    int kk = k, nn = n;
                    yield return new TestCaseData((Func<ClosedFormText>)(() => new PeriodicAkB(nn, kk)), psi).SetArgDisplayNames($"(a^{k} b)* n={n}", $"psi={psi}");
                }
        }
    }

    [TestCaseSource(nameof(RealContainers))]
    public void EveryEntryGiven_IsTheContainerOfTheProductWriter_ByteForByte(Func<ClosedFormText> make, bool psi)
    {
        ClosedFormText f = make();
        string text = f.Materialize();
        LcpChildCase c = LcpChildCase.Of(text);
        byte[] expected = SesaTestFile.Image(c, psi: psi, counts: true, buildInfo: true, codec: ChildCodec.Packed);
        PackedTextInfo info = TextOracle.Info(c.Symbols, c.Bits, c.Sigma);
        string infoText = SesaBuilder.BuildInfoText(info, c.Lcp.Kind, c.Lcp.OverflowCount, psi, ChildCodec.Packed);
        SparseSpec spec = SparseSpec.For(f, 0, psi) with { BuildInfo = infoText, TextSha256 = SHA256.HashData(Encoding.Latin1.GetBytes(text)) };
        SparseLayout layout = SparseContainer.Plan(spec);

        using var dir = new TempDir();
        string path = dir.File("sparse.sesa");
        SparseContainer.Write(path, spec, SparseContainer.AllPatches(f, layout));
        byte[] actual = File.ReadAllBytes(path);
        if (!actual.AsSpan().SequenceEqual(expected))
        {
            int at = actual.AsSpan().CommonPrefixLength(expected);
            Assert.Fail($"{f}: the container differs from the product's at byte {at} (lengths {actual.Length} and {expected.Length}): {(at < actual.Length ? actual[at] : -1)} against {(at < expected.Length ? expected[at] : -1)}");
        }
        // and it is a container: every level of the S1 reader accepts it
        using SesaFile file = SesaFile.Open(path, new EsaOpenOptions { Validation = EsaValidationLevel.Structure, Access = EsaAccessMode.Checked });
        Assert.That(file.Header.TextLength, Is.EqualTo(f.Length));
    }

    [TestCase(1000L)]
    [TestCase(70_000L)]
    public void AShellOfZeros_OpensAtTheHeaderAndSectionsLevels_AndTheStructureLevelRefusesIt(long n)
    {
        using var dir = new TempDir();
        string path = dir.File("shell.sesa");
        ClosedFormText f = new PowerOfA(n);
        SparseSpec spec = SparseSpec.For(f);
        SparseContainer.Write(path, spec, SparseContainer.EndPatches(f, SparseContainer.Plan(spec), 64));
        foreach (EsaValidationLevel level in new[] { EsaValidationLevel.Header, EsaValidationLevel.Sections })
        {
            using SesaFile file = SesaFile.Open(path, new EsaOpenOptions { Validation = level, Access = EsaAccessMode.Checked });
            Assert.That((file.Header.TextLength, file.Header.NodeCount, file.Header.LrsDepth, file.Header.LrsRank), Is.EqualTo((n, f.NodeCount, f.LrsDepth, f.LrsRank)), $"{level}");
        }
        Assert.Throws<InvalidDataException>(() => SesaFile.Open(path, new EsaOpenOptions { Validation = EsaValidationLevel.Structure, Access = EsaAccessMode.Checked }).Dispose(), "the hole of the middle is not a suffix array");
    }

    [Test]
    public void AShellWithPatches_HasTheCrcOfItsBytes_AtTheSectionsLevelOfTheS1Reader()
    {
        // a patch in every section kind, in the middle too: the CRC (zeros combined with the patch) must equal the product's pass over the real bytes
        using var dir = new TempDir();
        ClosedFormText f = new PeriodicAkB(70_000, 3);
        SparseSpec spec = SparseSpec.For(f, 0, psi: true);
        SparseLayout layout = SparseContainer.Plan(spec);
        var patches = new List<SparsePatch>(SparseContainer.EndPatches(f, layout, 100));
        foreach (int kind in new[] { SesaFormat.KindSa, SesaFormat.KindChild, SesaFormat.KindPsi, SesaFormat.KindLcp })
            patches.Add(new SparsePatch(kind, 4 * 30_000, [1, 2, 3, 4, 5, 6, 7, 8]));
        patches.Add(new SparsePatch(SesaFormat.KindText, 1000, [0xAA, 0x55]));
        string path = dir.File("shell.sesa");
        SparseContainer.Write(path, spec, patches);
        using SesaFile file = SesaFile.Open(path, new EsaOpenOptions { Validation = EsaValidationLevel.Sections, Access = EsaAccessMode.Checked });
        Assert.That(file.Report.SectionBytesChecksummed, Is.GreaterThan(0), "the checksum pass ran");
        // a corrupted byte under a patch is a CRC error (the pass reads the real bytes)
        byte[] bytes = File.ReadAllBytes(path);
        SesaOracle o = SesaOracle.FromBytes(bytes);
        bytes[o.Get(SesaFormat.KindSa).Offset + (4 * 30_000) + 2] ^= 1;
        string damaged = dir.File("damaged.sesa");
        File.WriteAllBytes(damaged, bytes);
        InvalidDataException e = Assert.Throws<InvalidDataException>(() => SesaFile.Open(damaged, new EsaOpenOptions { Validation = EsaValidationLevel.Sections, Access = EsaAccessMode.Checked }).Dispose())!;
        Assert.That(e.Message, Does.Contain("SA").And.Contain("CRC-32C"));
    }

    [TestCase(1000L, 5)]
    [TestCase(1_000_000L, 5)]
    [TestCase(1000L, 4)]
    public void TheLayoutOfAWidth5Container_IsTheOneOfAdr005Section61(long n, int width)
    {
        using var dir = new TempDir();
        ClosedFormText f = new PowerOfA(n);
        SparseSpec spec = SparseSpec.For(f, width, psi: true);
        string path = dir.File("w.sesa");
        SparseContainer.Write(path, spec, SparseContainer.EndPatches(f, SparseContainer.Plan(spec), 100));
        byte[] image = File.ReadAllBytes(path);
        SesaOracle o = SesaOracle.FromBytes(image);
        Assert.That(o.PosWidth, Is.EqualTo(width));
        Assert.That((o.VersionMajor, o.VersionMinor), Is.EqualTo(width == 5 ? (1, 1) : (1, 0)), "a container with 5-byte positions is version 1.1");
        Assert.That(o.FileSize, Is.EqualTo(image.Length));
        Assert.That(o.HeaderCrc, Is.EqualTo(o.ExpectedHeaderCrc()), "the header CRC covers the version word");
        int padding = width == 5 ? 3 : 0;
        foreach (int kind in new[] { SesaFormat.KindSa, SesaFormat.KindChild, SesaFormat.KindPsi })
        {
            OracleSection s = o.Get(kind);
            Assert.That((s.Count, s.ElemWidth, s.Length), Is.EqualTo((n + 1, width, ((n + 1) * width) + padding)), $"kind {kind}: COUNT, ELEM_WIDTH and LENGTH with the padding that belongs to the payload");
        }
        foreach (OracleSection s in o.Sections)
        {
            Assert.That(s.Crc, Is.EqualTo(SesaOracle.Crc32cBitwise(o.Payload(s))), $"kind {s.Kind}: CRC-32C of the payload with its padding");
            Assert.That(s.Offset % (s.Length >= (1 << 20) ? 4096 : 64), Is.Zero, $"kind {s.Kind}: alignment");
        }
    }

    [Test]
    public void TheNominalSize_AtTheSizesThatMatter_IsTheArithmeticOfTheLayout()
    {
        // n = 2^32 + 10, width 5, a^n: LCP values reach 2^32 + 9, so the LCP is U64 (codec 4); no PSI
        long n = (1L << 32) + 10, rows = n + 1;
        SparseLayout w5 = SparseContainer.Plan(SparseSpec.For(new PowerOfA(n)));
        Assert.That(w5.PosWidth, Is.EqualTo(5));
        Assert.That(w5.LcpCodec, Is.EqualTo(4));
        Assert.That(w5.Section(SesaFormat.KindSa).Length, Is.EqualTo((5 * rows) + 3));
        Assert.That(w5.Section(SesaFormat.KindChild).Length, Is.EqualTo((5 * rows) + 3));
        Assert.That(w5.Section(SesaFormat.KindLcp).Length, Is.EqualTo(8 * rows));
        Assert.That(w5.Section(SesaFormat.KindText).Length, Is.EqualTo(((n + 7) / 8) + 64), "one bit per symbol for a one-letter alphabet");
        long lcpEnd = w5.Section(SesaFormat.KindLcp).Offset + w5.Section(SesaFormat.KindLcp).Length;
        Assert.That(w5.FileSize, Is.EqualTo(lcpEnd));
        Assert.That(w5.FileSize, Is.EqualTo(SparseContainer.NominalBytes(SparseSpec.For(new PowerOfA(n)))));
        Assert.That(w5.FileSize, Is.InRange((5 * rows * 2) + (8 * rows), (5 * rows * 2) + (8 * rows) + (1L << 30)), "about 85 GB nominal: the two position sections and the 8-byte LCP");
        // n = 2^31 + 10 at width 4: LCP values reach 2^31 + 9 which U32 holds
        n = (1L << 31) + 10;
        rows = n + 1;
        SparseLayout w4 = SparseContainer.Plan(SparseSpec.For(new PowerOfA(n)));
        Assert.That((w4.PosWidth, w4.LcpCodec), Is.EqualTo((4, 3)));
        Assert.That(w4.Section(SesaFormat.KindSa).Length, Is.EqualTo(4 * rows));
        Assert.That(w4.Section(SesaFormat.KindLcp).Length, Is.EqualTo(4 * rows));
        // the width follows n: 4 up to 2^32 - 4, 5 above
        Assert.That(SparseSpec.WidthFor(SparseSpec.MaxContainer4), Is.EqualTo(4));
        Assert.That(SparseSpec.WidthFor(SparseSpec.MaxContainer4 + 1), Is.EqualTo(5));
        Assert.That(SparseSpec.WidthFor(0), Is.EqualTo(4));
        // the LCP codec follows the largest value
        Assert.That(new[] { 0L, 254, 255, 65_534, 65_535, uint.MaxValue - 1L, uint.MaxValue, uint.MaxValue + 1L }.Select(SparseSpec.LcpCodecFor), Is.EqualTo(new[] { 1, 1, 2, 2, 3, 3, 4, 4 }));
    }

    [Test]
    public void TheWriter_RefusesWhatItPromisesToRefuse()
    {
        using var dir = new TempDir();
        ClosedFormText f = new PowerOfA(1000);
        SparseSpec spec = SparseSpec.For(f);
        string Path(string name) => dir.File(name);
        // more than 64 MiB of real data
        Assert.Throws<ArgumentException>(() => SparseContainer.Write(Path("a"), spec, Big()), "dense region above 64 MiB");
        Assert.That(File.Exists(Path("a")), Is.False);
        // a patch outside its section, into the padding, overlapping, in a section that does not exist
        Assert.Throws<ArgumentException>(() => SparseContainer.Write(Path("b"), spec, [new SparsePatch(SesaFormat.KindSa, 4 * 1001 - 2, new byte[4])]), "past the end of SA");
        Assert.Throws<ArgumentException>(() => SparseContainer.Write(Path("c"), spec, [new SparsePatch(SesaFormat.KindSa, -1, new byte[4])]), "before the start");
        Assert.Throws<ArgumentException>(() => SparseContainer.Write(Path("d"), spec, [new SparsePatch(SesaFormat.KindSa, 0, new byte[8]), new SparsePatch(SesaFormat.KindSa, 4, new byte[8])]), "overlap");
        Assert.Throws<ArgumentException>(() => SparseContainer.Write(Path("e"), spec, [new SparsePatch(SesaFormat.KindPsi, 0, new byte[4])]), "no PSI section");
        SparseSpec w5 = SparseSpec.For(f, 5);
        Assert.Throws<ArgumentException>(() => SparseContainer.Write(Path("f"), w5, [new SparsePatch(SesaFormat.KindSa, (5 * 1001) - 4, new byte[5])]), "the 3 padding bytes of a width-5 section stay zero");
        Assert.Throws<ArgumentException>(() => SparseContainer.Write(Path("g"), spec, [new SparsePatch(SesaFormat.KindText, ((1000 + 7) / 8) - 1, new byte[65])]), "the 64 padding bytes of TEXT stay zero");
        foreach (string name in new[] { "b", "c", "d", "e", "f", "g" }) Assert.That(File.Exists(Path(name)), Is.False, $"nothing is left behind by a refused write ({name})");
        // width 4 holds at most 2^32 - 4 symbols; no width but 4 and 5; at most 2^40 - 1 symbols
        Assert.Throws<ArgumentOutOfRangeException>(() => SparseContainer.Plan(new SparseSpec { N = SparseSpec.MaxContainer4 + 1, PosWidth = 4, NodeCount = 1 }));
        Assert.Throws<ArgumentOutOfRangeException>(() => SparseContainer.Plan(new SparseSpec { N = 10, PosWidth = 8, NodeCount = 1 }));
        Assert.Throws<ArgumentOutOfRangeException>(() => SparseContainer.Plan(new SparseSpec { N = 1L << 40, NodeCount = 1 }));
        // a nominal size above the cap
        Assert.Throws<ArgumentException>(() => SparseContainer.Write(Path("h"), SparseSpec.For(new PowerOfA((1L << 36) + 1)), []), "the cap on the nominal size");
        Assert.That(File.Exists(Path("h")), Is.False);
        // an existing file is neither overwritten nor deleted
        File.WriteAllText(Path("keep"), "keep");
        Assert.Throws<IOException>(() => SparseContainer.Write(Path("keep"), spec, []));
        Assert.That(File.ReadAllText(Path("keep")), Is.EqualTo("keep"));

        static IReadOnlyList<SparsePatch> Big() => [new SparsePatch(SesaFormat.KindSa, 0, new byte[(int)LargeDisk.MaxDenseBytes])];
    }

    [Test]
    public void TheDiskGuard_DecidesBySparseness()
    {
        Assert.That(LargeDisk.Evaluate(0, 1L << 30).Sparse, Is.True);
        Assert.That(LargeDisk.Evaluate(LargeDisk.MaxAllocatedByProbe - 1, 1L << 30).Sparse, Is.True);
        Assert.That(LargeDisk.Evaluate(LargeDisk.MaxAllocatedByProbe, 1L << 30).Sparse, Is.False);
        Assert.That(LargeDisk.Evaluate(1L << 30, 1L << 30).Sparse, Is.False);
        Assert.That(LargeDisk.Evaluate(-1, 1L << 30).Sparse, Is.False, "stat could not be read: not sparse");
        Assert.That(LargeDisk.Evaluate(1L << 30, 1L << 30).Message, Does.Contain("NOT"));
    }

    [Test]
    public void TheProbe_OnTheDirectoryOfTheTests_FindsSparseFiles_AndLeavesNothing()
    {
        using var dir = new TempDir();
        SparseProbe probe = LargeDisk.Probe(dir.Path, 64L << 20);
        if (!probe.Sparse) Assert.Inconclusive(probe.Message);
        Assert.That(probe.AllocatedBytes, Is.LessThan(LargeDisk.MaxAllocatedByProbe));
        Assert.That(Directory.GetFiles(dir.Path), Is.Empty, "the probe file is deleted");
        Assert.That(LargeDisk.Probe(dir.Path, 64L << 20), Is.SameAs(probe), "the answer is remembered");
    }

    [Test]
    public void TheDirectory_IsTheEnvironmentVariableOrASubdirectoryOfTheTemporaryDirectory()
    {
        string? old = Environment.GetEnvironmentVariable("SEQERON_LARGE_DIR");
        try
        {
            Environment.SetEnvironmentVariable("SEQERON_LARGE_DIR", null);
            Assert.That(LargeDisk.Directory, Is.EqualTo(Path.Combine(Path.GetTempPath(), "seqeron-large")));
            Environment.SetEnvironmentVariable("SEQERON_LARGE_DIR", "/some/dir");
            Assert.That(LargeDisk.Directory, Is.EqualTo("/some/dir"));
            Environment.SetEnvironmentVariable("SEQERON_LARGE_DIR", "  ");
            Assert.That(LargeDisk.Directory, Is.EqualTo(Path.Combine(Path.GetTempPath(), "seqeron-large")));
        }
        finally
        {
            Environment.SetEnvironmentVariable("SEQERON_LARGE_DIR", old);
        }
    }

    [Test]
    public void ThePatchesOfAClosedForm_AreItsEntries()
    {
        // the end patches carry the entries of the two ends: SA, LCP, CHILD (packed) in the width of the container, TEXT from the symbols
        ClosedFormText f = new PeriodicAkB(5000, 4);
        SparseLayout layout = SparseContainer.Plan(SparseSpec.For(f, 5));
        IReadOnlyList<SparsePatch> patches = SparseContainer.EndPatches(f, layout, 50);
        SparsePatch headSa = patches.First(p => p.Kind == SesaFormat.KindSa && p.Offset == 0);
        Assert.That(headSa.Data.Length, Is.EqualTo(50 * 5));
        for (int r = 0; r < 50; r++)
        {
            long value = 0;
            for (int b = 0; b < 5; b++) value |= (long)headSa.Data[(5 * r) + b] << (8 * b);
            Assert.That(value, Is.EqualTo(f.SuffixAt(r)), $"SA[{r}]");
        }
        SparsePatch tailSa = patches.First(p => p.Kind == SesaFormat.KindSa && p.Offset > 0);
        Assert.That(tailSa.Offset, Is.EqualTo(5L * (5001 - 50)));
        SparsePatch text = patches.First(p => p.Kind == SesaFormat.KindText && p.Offset == 0);
        for (int s = 0; s < 56; s++)
            Assert.That((text.Data[s / 8] >> (7 - (s % 8))) & 1, Is.EqualTo(f.CodeAt(s)), $"text bit {s}: one bit per symbol, the first bit of a byte is its most significant");
        Assert.That(SparseContainer.PsiOf(f, 0), Is.Zero);
        Assert.That(SparseContainer.PsiOf(f, f.RankOf(4999)), Is.Zero, "the successor of the last suffix is the empty suffix");
        Assert.That(SparseContainer.PsiOf(f, f.RankOf(10)), Is.EqualTo(f.RankOf(11)));
    }
}
