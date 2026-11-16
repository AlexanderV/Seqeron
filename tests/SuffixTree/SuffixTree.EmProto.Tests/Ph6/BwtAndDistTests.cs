using SuffixTree.EmProto.Dcs;
using SuffixTree.EmProto.Harness;
using SuffixTree.EmProto.Ph6;
using SuffixTree.EmProto.Text;
using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.EmProto.Tests.Ph6;

/// <summary>The BWT pass (route (ii) of E9: BWT[j] = T[SA[j] - 1] from the suffix array file and the packed text) and its side file.</summary>
[TestFixture]
public unsafe class BwtPassTests
{
    private static readonly int[] Sigmas = [1, 2, 4, 5, 20, 64, 65, 255, 256];

    private static BwtMeta Run(Ph6Dir dir, Ph6Case c, BwtOptions? o, out string bwtPath, int? sigma = null)
    {
        string sa = c.WriteSa(dir);
        string pk = c.WritePacked(dir);
        bwtPath = dir.File("x.bwt");
        using var mapped = new MappedFile(pk, writable: false);
        using var text = new PackedText(mapped.Pointer, c.N, Math.Max(1, PackedText.BitsFor(c.Sigma)), mapped, ownsOwner: false);
        return BwtPass.Build(sa, text, sigma ?? c.Sigma, bwtPath, o);
    }

    [Test]
    public void BwtEqualsTheDefinition_AllAlphabets_TinyBuffers()
    {
        int v = 0;
        foreach (int sigma in Sigmas)
            for (int rep = 0; rep < 8; rep++, v++)
            {
                int n = sigma + 1 + ((rep * 53) % 200);
                Ph6Case c = Ph6Case.Of(Ph6Case.RandomText(n, sigma, (ulong)(500 + v)));
                using var dir = new Ph6Dir();
                BwtMeta meta = Run(dir, c, new BwtOptions { ReadBytes = new[] { 4, 8, 28, 4096 }[v % 4], WindowBytes = new[] { 1, 5, 64, 1 << 16 }[(v / 4) % 4] }, out string bwt);
                Assert.That(File.ReadAllBytes(bwt), Is.EqualTo(c.ExpectedBwt()), $"sigma {sigma} n {n}");
                Assert.That(meta.N, Is.EqualTo(n));
                Assert.That(meta.PrimaryRow, Is.EqualTo(c.Primary));
                Assert.That(meta.Counts, Is.EqualTo(c.Counts()));
                Assert.That(meta.Sigma, Is.EqualTo(sigma));
            }
    }

    [Test]
    public void MetaFile_RoundTrips_AndIsWrittenNextToTheBwt()
    {
        Ph6Case c = Ph6Case.OfNaive("mississippi");
        using var dir = new Ph6Dir();
        BwtMeta meta = Run(dir, c, null, out string bwt);
        Assert.That(File.Exists(BwtMeta.MetaPath(bwt)), Is.True);
        BwtMeta back = BwtMeta.Read(bwt);
        Assert.That(back.N, Is.EqualTo(meta.N));
        Assert.That(back.PrimaryRow, Is.EqualTo(meta.PrimaryRow));
        Assert.That(back.Counts, Is.EqualTo(meta.Counts));
        BwtMeta parsed = BwtMeta.Parse("n=5\nprimary=3\nsigma=2\ncounts=2,3\n");
        Assert.That((parsed.N, parsed.PrimaryRow), Is.EqualTo((5L, 3L)));
        Assert.That(parsed.Counts, Is.EqualTo(new long[] { 2, 3 }));
        Assert.That(BwtMeta.Parse(new BwtMeta(0, 0, []).Text()).Counts, Is.Empty, "no symbols: an empty counts line");
        Assert.Throws<InvalidDataException>(() => BwtMeta.Parse("n=5\nprimary=3\ncounts=2,3\n"), "sigma line missing");
        Assert.Throws<InvalidDataException>(() => BwtMeta.Parse("n=5\nprimary=3\nsigma=3\ncounts=2,3\n"), "sigma does not match the counts");
        Assert.Throws<InvalidDataException>(() => BwtMeta.Parse("n=5\nsigma=1\ncounts=5\n"), "primary missing");
    }

    [Test]
    public void SentinelRow_IsTheDummyByte_AndTheRowOfSuffixZero()
    {
        foreach (string text in new[] { "abcdefg", "gfedcba", "banana", "aaaa" })
        {
            Ph6Case c = Ph6Case.OfNaive(text);
            using var dir = new Ph6Dir();
            BwtMeta meta = Run(dir, c, null, out string bwt);
            byte[] bytes = File.ReadAllBytes(bwt);
            Assert.That(bytes, Has.Length.EqualTo(c.N + 1));
            Assert.That(meta.PrimaryRow, Is.EqualTo(Array.IndexOf(c.Sa, 0)), text);
            Assert.That(bytes[meta.PrimaryRow], Is.EqualTo(0), "the dummy byte");
            Assert.That(bytes[0], Is.EqualTo(c.Symbols[^1]), "row 0 (the empty suffix) holds the last symbol of the text");
        }
    }

    [Test]
    public void EmptyText_HasOneRow()
    {
        using var dir = new Ph6Dir();
        string sa = dir.File("e.sa");
        File.WriteAllBytes(sa, new byte[4]); // SA = [0] = [n]
        string pk = dir.File("e.pk");
        File.WriteAllBytes(pk, new byte[PackedText.TotalBytes(0, 1)]);
        using var mapped = new MappedFile(pk, writable: false);
        using var text = new PackedText(mapped.Pointer, 0, 1, mapped, ownsOwner: false);
        BwtMeta meta = BwtPass.Build(sa, text, 0, dir.File("e.bwt"));
        Assert.That(meta.N, Is.EqualTo(0));
        Assert.That(meta.PrimaryRow, Is.EqualTo(0));
        Assert.That(File.ReadAllBytes(dir.File("e.bwt")), Is.EqualTo(new byte[1]));
    }

    [Test]
    public void BadSuffixArrays_AreRefused()
    {
        Ph6Case c = Ph6Case.OfNaive("banana");
        // no row holds 0
        using (var dir = new Ph6Dir())
        {
            var bad = (int[])c.Sa.Clone();
            int at = Array.IndexOf(bad, 0);
            bad[at] = 3;
            Assert.Throws<InvalidDataException>(() => RunWithSa(dir, c, bad), "no suffix 0");
        }
        // an entry beyond n
        using (var dir = new Ph6Dir())
        {
            var bad = (int[])c.Sa.Clone();
            bad[2] = c.N + 1;
            Assert.Throws<InvalidDataException>(() => RunWithSa(dir, c, bad), "entry above n");
        }
        // two rows hold 0
        using (var dir = new Ph6Dir())
        {
            var bad = (int[])c.Sa.Clone();
            bad[Array.IndexOf(bad, 3)] = 0;
            Assert.Throws<InvalidDataException>(() => RunWithSa(dir, c, bad), "two rows with suffix 0");
        }
        // SA file of the wrong length
        using (var dir = new Ph6Dir())
        {
            string sa = dir.File("x.sa");
            File.WriteAllBytes(sa, new byte[4 * c.N]);
            string pk = c.WritePacked(dir);
            using var mapped = new MappedFile(pk, writable: false);
            using var text = new PackedText(mapped.Pointer, c.N, 2, mapped, ownsOwner: false);
            Assert.Throws<InvalidDataException>(() => BwtPass.Build(sa, text, c.Sigma, dir.File("x.bwt")));
        }
    }

    [Test]
    public void SymbolOutsideTheDeclaredAlphabet_IsRefused()
    {
        Ph6Case c = Ph6Case.OfNaive("banana"); // sigma 3
        using var dir = new Ph6Dir();
        Assert.Throws<InvalidDataException>(() => Run(dir, c, null, out _, sigma: 2));
    }

    [Test]
    public void ArenaBudget_IsReleased()
    {
        Ph6Case c = Ph6Case.OfNaive("banana");
        var budget = new ArenaBudget(160L << 20);
        using var dir = new Ph6Dir();
        Run(dir, c, new BwtOptions { Budget = budget }, out _);
        Assert.That(budget.LiveBytes, Is.EqualTo(0));
        Assert.That(budget.PeakBytes, Is.GreaterThan(0));
    }

    private static void RunWithSa(Ph6Dir dir, Ph6Case c, int[] sa)
    {
        string saPath = dir.File("x.sa");
        File.WriteAllBytes(saPath, System.Runtime.InteropServices.MemoryMarshal.AsBytes(sa.AsSpan()).ToArray());
        string pk = c.WritePacked(dir);
        using var mapped = new MappedFile(pk, writable: false);
        using var text = new PackedText(mapped.Pointer, c.N, Math.Max(1, PackedText.BitsFor(c.Sigma)), mapped, ownsOwner: false);
        BwtPass.Build(saPath, text, c.Sigma, dir.File("x.bwt"));
    }
}

/// <summary>ISA by one-level radix distribution (the DIST primitive) and Psi from it.</summary>
[TestFixture]
public class IsaByDistTests
{
    private static DistStats Run(Ph6Dir dir, Ph6Case c, long r, int bufferBytes, int readBytes, out string isa, string? scratch = null)
    {
        string sa = c.WriteSa(dir);
        isa = dir.File("x.isa");
        return IsaByDist.Build(sa, c.N, isa, new DistOptions { ScratchDir = scratch ?? dir.File("scratch"), ChunkPositions = r, IoBufferBytes = bufferBytes, ReadBytes = readBytes, WindowBytes = 24 });
    }

    [Test]
    public void IsaEqualsTheInRamIsa_ChunkAndBufferSizes()
    {
        int v = 0;
        foreach (int sigma in new[] { 1, 2, 4, 20, 256 })
            for (int rep = 0; rep < 10; rep++, v++)
            {
                int n = sigma + 1 + ((rep * 37) % 300);
                Ph6Case c = Ph6Case.Of(Ph6Case.RandomText(n, sigma, (ulong)(800 + v)));
                using var dir = new Ph6Dir();
                long r = new long[] { 1, 2, 3, 7, 64, 1000, 1L << 30 }[v % 7];
                DistStats st = Run(dir, c, r, new[] { 8, 16, 24, 65536 }[v % 4], new[] { 4, 12, 4096 }[v % 3], out string isa);
                Assert.That(File.ReadAllBytes(isa), Is.EqualTo(c.ExpectedIsaBytes()), $"sigma {sigma} n {n} R {r}");
                Assert.That(st.Chunks, Is.EqualTo((n + 1 + r - 1) / r), "chunks");
                Assert.That(Directory.GetFileSystemEntries(dir.File("scratch")), Is.Empty, "chunk files are deleted");
            }
    }

    [Test]
    public void Families_AllSameFibonacciThueMorse()
    {
        foreach (string text in new[] { StringFamilies.AllSame(5000), StringFamilies.Fibonacci(5000), StringFamilies.ThueMorse(5000) })
        {
            Ph6Case c = Ph6Case.Of(text);
            using var dir = new Ph6Dir();
            Run(dir, c, 311, 64, 256, out string isa);
            Assert.That(File.ReadAllBytes(isa), Is.EqualTo(c.ExpectedIsaBytes()));
        }
    }

    [Test]
    public void Volume_IsSixWordsPerEntry_WhenEveryTupleGoesThroughAFile()
    {
        Ph6Case c = Ph6Case.Of(Ph6Case.RandomText(1000, 4, 9));
        long rows = c.N + 1;
        using var dir = new Ph6Dir();
        long r0 = IoCounters.ReadBytes, w0 = IoCounters.WriteBytes, m0 = MappedIo.WriteBytes;
        // a buffer of one tuple: every tuple is written to its chunk file and read back
        DistStats st = Run(dir, c, 100, 8, 400, out _);
        Assert.That(st.TupleBytesWritten, Is.EqualTo(8 * rows));
        Assert.That(st.TupleBytesRead, Is.EqualTo(8 * rows));
        Assert.That(IoCounters.ReadBytes - r0, Is.EqualTo(4 * rows + 8 * rows), "SA once, tuples once");
        Assert.That(IoCounters.WriteBytes - w0, Is.EqualTo(8 * rows), "tuples once");
        Assert.That(MappedIo.WriteBytes - m0, Is.EqualTo(4 * rows), "ISA once");
        // a buffer that holds a whole chunk: nothing goes to a file
        using var dir2 = new Ph6Dir();
        r0 = IoCounters.ReadBytes;
        w0 = IoCounters.WriteBytes;
        st = Run(dir2, c, 100, 1 << 20, 400, out _);
        Assert.That(st.TupleBytesWritten, Is.EqualTo(0));
        Assert.That(IoCounters.ReadBytes - r0, Is.EqualTo(4 * rows));
        Assert.That(IoCounters.WriteBytes - w0, Is.EqualTo(0));
    }

    [Test]
    public void NotAPermutation_IsRefused()
    {
        Ph6Case c = Ph6Case.OfNaive("banana");
        foreach (long r in new long[] { 2, 100 })
        {
            using var dir = new Ph6Dir();
            var dup = (int[])c.Sa.Clone();
            dup[3] = dup[4];
            string sa = dir.File("x.sa");
            File.WriteAllBytes(sa, System.Runtime.InteropServices.MemoryMarshal.AsBytes(dup.AsSpan()).ToArray());
            Assert.Throws<InvalidDataException>(() => IsaByDist.Build(sa, c.N, dir.File("x.isa"), new DistOptions { ScratchDir = dir.File("s"), ChunkPositions = r }), $"duplicate entry, R {r}");
            var big = (int[])c.Sa.Clone();
            big[2] = c.N + 1;
            File.WriteAllBytes(sa, System.Runtime.InteropServices.MemoryMarshal.AsBytes(big.AsSpan()).ToArray());
            Assert.Throws<InvalidDataException>(() => IsaByDist.Build(sa, c.N, dir.File("x.isa"), new DistOptions { ScratchDir = dir.File("s2"), ChunkPositions = r }), $"entry above n, R {r}");
        }
    }

    [Test]
    public void SaFileOfTheWrongLength_IsRefused()
    {
        Ph6Case c = Ph6Case.OfNaive("banana");
        using var dir = new Ph6Dir();
        string sa = c.WriteSa(dir);
        Assert.Throws<InvalidDataException>(() => IsaByDist.Build(sa, c.N + 1, dir.File("x.isa"), new DistOptions { ScratchDir = dir.File("s") }));
        Assert.Throws<ArgumentOutOfRangeException>(() => IsaByDist.Build(sa, c.N, dir.File("x.isa"), new DistOptions { ScratchDir = dir.File("s"), ChunkPositions = 0 }));
        Assert.Throws<ArgumentOutOfRangeException>(() => IsaByDist.Build(sa, c.N, dir.File("x.isa"), new DistOptions { ScratchDir = dir.File("s"), IoBufferBytes = 7 }));
    }

    [Test]
    public void ArenaBudget_IsReleased_AndEnforced()
    {
        Ph6Case c = Ph6Case.Of(Ph6Case.RandomText(500, 4, 1));
        using var dir = new Ph6Dir();
        var budget = new ArenaBudget(160L << 20);
        IsaByDist.Build(c.WriteSa(dir), c.N, dir.File("x.isa"), new DistOptions { ScratchDir = dir.File("s"), ChunkPositions = 50, Budget = budget });
        Assert.That(budget.LiveBytes, Is.EqualTo(0));
        Assert.That(budget.PeakBytes, Is.GreaterThan(0));
        var tiny = new ArenaBudget(1000);
        Assert.Throws<ArenaBudgetExceededException>(() => IsaByDist.Build(dir.File("x.sa"), c.N, dir.File("y.isa"), new DistOptions { ScratchDir = dir.File("s"), ChunkPositions = 50, Budget = tiny }));
        Assert.That(tiny.LiveBytes, Is.EqualTo(0));
    }
}

[TestFixture]
public class PsiFromIsaTests
{
    [Test]
    public void PsiEqualsTheInRamPsi_TinyBuffers()
    {
        int v = 0;
        foreach (int sigma in new[] { 1, 2, 5, 256 })
            for (int rep = 0; rep < 8; rep++, v++)
            {
                int n = sigma + 1 + ((rep * 61) % 250);
                Ph6Case c = Ph6Case.Of(Ph6Case.RandomText(n, sigma, (ulong)(300 + v)));
                using var dir = new Ph6Dir();
                string sa = c.WriteSa(dir);
                File.WriteAllBytes(dir.File("x.isa"), c.ExpectedIsaBytes());
                PsiFromIsa.Build(sa, dir.File("x.isa"), c.N, dir.File("x.psi"), readBytes: new[] { 4, 12, 4096 }[v % 3], windowBytes: new[] { 4, 20, 1 << 16 }[v % 3]);
                Assert.That(File.ReadAllBytes(dir.File("x.psi")), Is.EqualTo(c.ExpectedPsiBytes(false)), $"sigma {sigma} n {n}");
            }
    }

    [Test]
    public void IsaThatIsTooShort_IsRefused()
    {
        Ph6Case c = Ph6Case.OfNaive("banana");
        using var dir = new Ph6Dir();
        string sa = c.WriteSa(dir);
        File.WriteAllBytes(dir.File("x.isa"), new byte[4 * c.N]);
        Assert.Throws<InvalidDataException>(() => PsiFromIsa.Build(sa, dir.File("x.isa"), c.N, dir.File("x.psi")));
    }
}

/// <summary>The files and their relation to the tables of P1.</summary>
[TestFixture]
public class PsiFilesTests
{
    [Test]
    public void LoadForP1_MapsEntryZeroToMinusOne_ForBothConventions()
    {
        Ph6Case c = Ph6Case.OfNaive("mississippi");
        using var dir = new Ph6Dir();
        foreach (bool cyclic in new[] { false, true })
        {
            File.WriteAllBytes(dir.File("x.psi"), c.ExpectedPsiBytes(cyclic));
            int[] psi = PsiFiles.LoadForP1(dir.File("x.psi"));
            Assert.That(psi, Is.EqualTo(c.PsiRam), $"cyclic {cyclic}");
            Assert.That(psi[0], Is.EqualTo(-1));
        }
    }

    [Test]
    public void ReadInts_ReadsLittleEndian_AndRefusesARaggedFile()
    {
        using var dir = new Ph6Dir();
        File.WriteAllBytes(dir.File("a.bin"), [1, 0, 0, 0, 0, 1, 0, 0, 255, 255, 255, 127]);
        Assert.That(PsiFiles.ReadInts(dir.File("a.bin")), Is.EqualTo(new[] { 1, 256, int.MaxValue }));
        File.WriteAllBytes(dir.File("b.bin"), new byte[6]);
        Assert.Throws<InvalidDataException>(() => PsiFiles.ReadInts(dir.File("b.bin")));
    }

    [Test]
    public void WriteOracle_WritesIsaAndPsiOfP1()
    {
        Ph6Case c = Ph6Case.Of(Ph6Case.RandomText(400, 5, 4));
        using var dir = new Ph6Dir();
        PsiFiles.WriteOracle(c.WriteSa(dir), dir.File("o.isa"), dir.File("o.psi"));
        Assert.That(File.ReadAllBytes(dir.File("o.isa")), Is.EqualTo(c.ExpectedIsaBytes()));
        Assert.That(File.ReadAllBytes(dir.File("o.psi")), Is.EqualTo(c.ExpectedPsiBytes(false)));
    }
}
