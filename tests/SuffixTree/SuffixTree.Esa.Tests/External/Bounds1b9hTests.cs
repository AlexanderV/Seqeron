using System.Reflection;
using System.Runtime.InteropServices;
using System.Runtime.Intrinsics.X86;
using SuffixTree.Esa.External;
using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.Esa.Tests.External;

/// <summary>
/// The bounds checks of the checked build (<c>-p:EsaBounds=true</c>, <c>ESA_EXPECT_BOUNDS=true</c>) on the two things that PH3 reads through raw pointers, the rank file of the difference-cover sample and the mapped packed text (work package
/// 1b.9h, debt D: 14 of the 28 mutants of 1b.9c killed the test host by an <see cref="AccessViolationException"/>, a read outside a mapping, instead of failing a test). The production build keeps its raw reads, so what a violation does
/// there is not tested (it cannot be: the process would die); what is tested there is that the reads that are right are right, which the existing tests of both readers do for every position.
/// </summary>
[TestFixture]
public class RankTableBounds1b9hTests
{
    private TempDir _dir = null!;

    [SetUp]
    public void SetUp() => _dir = new TempDir();

    [TearDown]
    public void TearDown() => _dir.Dispose();

    private static void RequireChecked()
    {
        if (!EsaBuild.BoundsChecked) Assert.Ignore("only the checked build (-p:EsaBounds=true) bounds-checks the reads of the rank table; the production build reads through the raw pointer");
    }

    private (RankTable Table, DifferenceCover Cover, SaCase Case, uint[] Expected) Open(int period, string text)
    {
        DifferenceCover cover = DifferenceCover.Create(period);
        SaCase c = SaCase.Of(text);
        uint[] expected = DcsOracle.ExpectedRanks(c.Symbols, cover, out _);
        string path = _dir.File("ranks-" + Guid.NewGuid().ToString("N")[..6]);
        File.WriteAllBytes(path, MemoryMarshal.AsBytes(expected.AsSpan()).ToArray());
        return (RankTable.Open(path, cover, c.N), cover, c, expected);
    }

    /// <summary>The first position above <paramref name="from"/> whose residue is in D.</summary>
    private static long NextSamplePosition(DifferenceCover cover, long from)
    {
        long p = from + 1;
        while (cover.IndexOf[(int)(p % cover.Period)] < 0) p++;
        return p;
    }

    [Test]
    public void TheLastPositionsOfTheText_AreRead_InBothBuilds_AndAreTheRanksOfTheOracle()
    {
        foreach (int period in new[] { 7, 13, 31 })
        {
            foreach (string text in new[] { "banana", StringFamilies.RandomDna(100, 3), "a", string.Empty, StringFamilies.Fibonacci(60) })
            {
                (RankTable table, DifferenceCover cover, SaCase c, uint[] expected) = Open(period, text);
                using (table)
                {
                    int found = 0;
                    for (long p = Math.Max(0, c.N - cover.Period); p <= c.N; p++)
                    {
                        int index = cover.IndexOf[(int)(p % cover.Period)];
                        if (index < 0) continue;
                        Assert.That(table.Rank(p), Is.EqualTo(expected[((p / cover.Period) * cover.Residues.Length) + index]), $"X = {period}, n = {c.N}, p = {p}");
                        found++;
                    }
                    if (c.N >= cover.Period) Assert.That(found, Is.GreaterThanOrEqualTo(cover.Residues.Length), "the last X + 1 positions hold every residue of the sample");
                    table.Verify();
                }
            }
        }
    }

    [Test]
    public void APositionBeyondTheText_ThatIsInTheSample_IsAMemoryViolation_ButTheLastPositionIsNot()
    {
        RequireChecked();
        foreach (int period in new[] { 7, 13, 31 })
        {
            (RankTable table, DifferenceCover cover, SaCase c, uint[] expected) = Open(period, StringFamilies.RandomDna(40, 5));
            using (table)
            {
                long last = c.N;
                while (cover.IndexOf[(int)(last % cover.Period)] < 0) last--;
                Assert.That(table.Rank(last), Is.EqualTo(expected[((last / cover.Period) * cover.Residues.Length) + cover.IndexOf[(int)(last % cover.Period)]]), "the last sample position is read");
                long beyond = NextSamplePosition(cover, c.N);
                var e = Assert.Throws<EsaMemoryViolationException>(() => table.Rank(beyond), $"X = {period}, position {beyond} of a text of {c.N} symbols");
                Assert.That(e!.Kind, Is.EqualTo(MemoryViolationKind.OutOfBounds));
                Assert.That(e.Message, Does.Contain("rank table").And.Contain($"position {beyond}").And.Contain($"{c.N} symbols"));
                Assert.Throws<EsaMemoryViolationException>(() => table.Rank(NextSamplePosition(cover, beyond + (10L * cover.Period))), "a block far beyond the file");
                Assert.Throws<EsaMemoryViolationException>(() => table.Rank(long.MaxValue), "the largest position");
                Assert.Throws<EsaMemoryViolationException>(() => table.Rank((long)uint.MaxValue + 1 + cover.Residues[0]), "a position that does not fit 32 bits is not folded into the table");
            }
        }
    }

    [Test]
    public void ANegativePosition_IsAMemoryViolation_NotAReadBeforeTheMapping()
    {
        RequireChecked();
        (RankTable table, DifferenceCover cover, _, _) = Open(7, "banana");
        using (table)
        {
            foreach (long p in new[] { -1L, -(long)cover.Period, -7L * cover.Period, long.MinValue, -1L - uint.MaxValue })
            {
                Assert.Throws<EsaMemoryViolationException>(() => table.Rank(p), $"position {p}");
            }
        }
    }

    [Test]
    public void APositionOfTheTextWhoseResidueIsNotInTheSample_IsAMemoryViolation_NotAnArbitrarySlot()
    {
        RequireChecked();
        foreach (int period in new[] { 7, 13, 31 })
        {
            (RankTable table, DifferenceCover cover, SaCase c, _) = Open(period, StringFamilies.RandomDna(60, 7));
            using (table)
            {
                int outside = Enumerable.Range(0, period).First(r => cover.IndexOf[r] < 0);
                Assert.That(outside, Is.LessThan(c.N), "a text longer than the period");
                var e = Assert.Throws<EsaMemoryViolationException>(() => table.Rank(outside), $"X = {period}, residue {outside}");
                Assert.That(e!.Message, Does.Contain($"residue {outside}"));
                Assert.Throws<EsaMemoryViolationException>(() => table.Rank(outside + period), "the same residue in the next block");
            }
        }
    }

    [Test]
    public void ARankAfterDispose_IsAMemoryViolation_NotAReadOfTheReleasedMapping()
    {
        RequireChecked();
        (RankTable table, DifferenceCover cover, _, _) = Open(7, "banana");
        long p = cover.Residues[0];
        Assert.DoesNotThrow(() => table.Rank(p), "the table is open");
        table.Dispose();
        var e = Assert.Throws<EsaMemoryViolationException>(() => table.Rank(p));
        Assert.That(e!.Message, Does.Contain("disposed"));
    }

    [Test]
    public void TheVerificationOfADisposedTable_IsRefused_AsADisposedObject_NeverReadFromTheReleasedMapping()
    {
        (RankTable table, _, _, _) = Open(7, "banana");
        table.Dispose();
        Assert.Throws<ObjectDisposedException>(table.Verify);
    }
}

/// <summary>The reads of the mapped packed text through its raw pointer in the checked build: a read that would leave the region (the data and the 64 bytes of padding) is an <see cref="EsaMemoryViolationException"/> that names the bytes.</summary>
[TestFixture]
public class MappedTextBounds1b9hTests
{
    private static void RequireChecked()
    {
        if (!EsaBuild.BoundsChecked) Assert.Ignore("only the checked build (-p:EsaBounds=true) bounds-checks the reads of the mapped text; the production build reads through the raw pointer");
    }

    /// <summary>A reader whose idea of the text is longer than the region: the limit and the length that its own bounds test looks at are set to <paramref name="symbols"/> (the situation of a caller that computes a position wrongly).</summary>
    private static void LieAboutTheLength(MappedTextReader reader, long symbols)
    {
        const BindingFlags all = BindingFlags.Instance | BindingFlags.NonPublic;
        typeof(MappedTextReader).GetField("_limit", all)!.SetValue(reader, symbols);
        typeof(MappedTextReader).GetField("_n", all)!.SetValue(reader, symbols);
    }

    private static long TotalBytes(long n, int bits) => PackedTextFormat.TotalBytes(n, bits);

    [Test]
    public void TheWindowOfAPositionInTheText_IsReadInBothBuilds_UpToTheLastSymbol()
    {
        foreach (int bits in new[] { 1, 2, 3, 4, 5, 7, 8 })
        {
            byte[] s = TextOracle.Random(77, 1 << bits, (ulong)(bits * 31));
            using var t = new PinnedText(s, bits);
            for (long i = 0; i < s.Length; i++) Assert.That(t.Reader.Window(i), Is.EqualTo(TextOracle.Window(s, i, bits)), $"b = {bits}, i = {i}");
            Assert.That(t.Reader.Window(s.Length), Is.Zero);
        }
    }

    [Test]
    public void AWindow_ThatIsNotWhollyInsideTheRegion_IsAMemoryViolation_AtTheExactByte()
    {
        RequireChecked();
        foreach (int bits in new[] { 1, 2, 3, 4, 5, 8 })
        {
            byte[] s = TextOracle.Random(50, 1 << bits, (ulong)(bits * 17));
            using var t = new PinnedText(s, bits);
            long total = TotalBytes(s.Length, bits);
            LieAboutTheLength(t.Reader, 1L << 40);
            long lastGood = 0; // the symbol whose window starts at byte total - 9 reads nine bytes, the last of the region
            while (((lastGood + 1) * bits >> 3) + 9 <= total) lastGood++;
            Assert.DoesNotThrow(() => t.Reader.Window(lastGood), $"b = {bits}: the window of symbol {lastGood} ends at the last byte of the region");
            var e = Assert.Throws<EsaMemoryViolationException>(() => t.Reader.Window(lastGood + 1), $"b = {bits}: the window of symbol {lastGood + 1} reads one byte too many");
            Assert.That(e!.Kind, Is.EqualTo(MemoryViolationKind.OutOfBounds));
            Assert.That(e.Message, Does.Contain("mapped text").And.Contain($"of a region of {total} bytes"));
            Assert.Throws<EsaMemoryViolationException>(() => t.Reader.Symbol(1L << 30), "a position far beyond the region");
            Assert.Throws<EsaMemoryViolationException>(() => t.Reader.Key(total * 8), "a key far beyond the region");
        }
    }

    [Test]
    public void TheVectorCommonPrefix_ThatWouldLoadBeyondTheRegion_IsAMemoryViolation_AtTheExactByte()
    {
        RequireChecked();
        if (!Avx2.IsSupported) Assert.Ignore("the vector path needs AVX2");
        foreach (int bits in new[] { 2, 3, 4, 8 })
        {
            byte[] s = TextOracle.Random(300, 1 << bits, (ulong)(bits * 13));
            using var t = new PinnedText(s, bits);
            t.Reader.UseVectors = true;
            long total = TotalBytes(s.Length, bits);
            int vs = t.Reader.VectorSymbols;
            LieAboutTheLength(t.Reader, 1L << 40);
            long steps = 0;
            long lastGood = 0; // a step loads 40 bytes from the byte of its first symbol: the last good start is at byte total - 40
            while (((lastGood + 1) * bits >> 3) + 40 <= total) lastGood++;
            Assert.DoesNotThrow(() => t.Reader.Lce(lastGood, lastGood, vs, ref steps), $"b = {bits}: the step from symbol {lastGood} ends at the last byte of the region");
            var e = Assert.Throws<EsaMemoryViolationException>(() => t.Reader.Lce(lastGood + 1, lastGood + 1, vs, ref steps), $"b = {bits}: the step from symbol {lastGood + 1} loads one byte too many");
            Assert.That(e!.Message, Does.Contain("mapped text"));
        }
    }

    [Test]
    public void TheScalarCommonPrefix_ThatWouldReadBeyondTheRegion_IsAMemoryViolation()
    {
        RequireChecked();
        byte[] s = TextOracle.Random(200, 4, 99);
        using var t = new PinnedText(s, 2);
        t.Reader.UseVectors = false;
        LieAboutTheLength(t.Reader, 1L << 40);
        long steps = 0;
        Assert.Throws<EsaMemoryViolationException>(() => t.Reader.Lce(100, 1L << 30, 64, ref steps), "the second suffix starts far beyond the region");
    }

    [Test]
    public void ABulkRead_ThatWouldLeaveTheRegion_IsAMemoryViolation_AtTheExactByte()
    {
        RequireChecked();
        foreach (int bits in new[] { 3, 8 })
        {
            byte[] s = TextOracle.Random(120, 1 << bits, (ulong)(bits * 7));
            using var t = new PinnedText(s, bits);
            long total = TotalBytes(s.Length, bits);
            LieAboutTheLength(t.Reader, 1L << 40);
            var dest = new byte[64];
            long lastGood = 0; // the read touches the bytes from the first symbol's byte up to two bytes after the last symbol's
            while ((((lastGood + 1) * bits) >> 3) + ((((long)dest.Length * bits) + 7) >> 3) + 2 <= total) lastGood++;
            Assert.DoesNotThrow(() => t.Reader.ReadSymbols(lastGood, dest), $"b = {bits}: the read from symbol {lastGood} ends inside the region");
            var e = Assert.Throws<EsaMemoryViolationException>(() => t.Reader.ReadSymbols(lastGood + 1, dest), $"b = {bits}: the read from symbol {lastGood + 1} ends outside it");
            Assert.That(e!.Message, Does.Contain("mapped text"));
        }
    }

    [Test]
    public void TheBulkReadOfAWholeText_IsTheSymbols_InBothBuilds()
    {
        foreach (int bits in new[] { 1, 3, 5, 8 })
        {
            byte[] s = TextOracle.Random(1000, 1 << bits, (ulong)(bits * 5));
            using var t = new PinnedText(s, bits);
            var got = new byte[s.Length];
            t.Reader.ReadSymbols(0, got);
            Assert.That(got, Is.EqualTo(s), $"b = {bits}");
            var tail = new byte[9];
            t.Reader.ReadSymbols(s.Length - tail.Length, tail);
            Assert.That(tail, Is.EqualTo(s[^tail.Length..]), $"b = {bits}: the last nine symbols");
        }
    }
}

/// <summary>
/// The list of the files that the process has mapped for reading, in the checked build (work package 1b.9h, the mutant <c>L10</c> of 1b.9c): creating a file again truncates it under a live mapping, and the next read of the mapping is a bus error
/// that ends the process; two levels of the sample recursion that were given one file name did exactly that. In the checked build the creation is refused, before the file is touched.
/// </summary>
[TestFixture]
public class MappedFileRegistryTests
{
    private TempDir _dir = null!;

    [SetUp]
    public void SetUp() => _dir = new TempDir();

    [TearDown]
    public void TearDown() => _dir.Dispose();

    private static void RequireChecked()
    {
        if (!EsaBuild.BoundsChecked) Assert.Ignore("only the checked build (-p:EsaBounds=true) keeps the list of the mapped files; the production build creates the file and the next read of the mapping is a bus error");
    }

    private string RankFile(out DifferenceCover cover, out long n)
    {
        cover = DifferenceCover.Create(7);
        SaCase c = SaCase.Of("banana");
        n = c.N;
        string path = _dir.File("ranks");
        File.WriteAllBytes(path, MemoryMarshal.AsBytes(DcsOracle.ExpectedRanks(c.Symbols, cover, out _).AsSpan()).ToArray());
        return path;
    }

    [Test]
    public void TheCreationOfAFileThatIsMapped_IsAMemoryViolation_AndTheFileIsNotTouched()
    {
        RequireChecked();
        string path = RankFile(out DifferenceCover cover, out long n);
        long length = new FileInfo(path).Length;
        using (RankTable.Open(path, cover, n))
        {
            var e = Assert.Throws<EsaMemoryViolationException>(() => new CountedFile(path, FileMode.Create, FileAccess.ReadWrite));
            Assert.That(e!.Kind, Is.EqualTo(MemoryViolationKind.OutOfBounds));
            Assert.That(e.Message, Does.Contain(path).And.Contain("mapped"));
            Assert.Throws<EsaMemoryViolationException>(() => new CountedFile(path, FileMode.Truncate, FileAccess.ReadWrite));
            Assert.Throws<EsaMemoryViolationException>(() => MappedOutputFile.Create(path, 4096));
            Assert.That(new FileInfo(path).Length, Is.EqualTo(length), "the file is as it was: nothing was truncated");
            using var reader = new CountedFile(path, FileMode.Open, FileAccess.Read);
            Assert.That(reader.Length, Is.EqualTo(length), "opening it for reading is not a creation");
            using var append = new CountedFile(path, FileMode.OpenOrCreate, FileAccess.Read);
            Assert.That(append.Length, Is.EqualTo(length), "OpenOrCreate of a file that exists does not truncate it");
        }
        using (var created = new CountedFile(path, FileMode.Create, FileAccess.ReadWrite)) Assert.That(created.Length, Is.Zero, "once the table is disposed the file may be created again");
    }

    [Test]
    public void TheTextReader_ListsItsFile_AndTakesItOffAtDispose()
    {
        RequireChecked();
        byte[] symbols = TextOracle.Random(100, 4, 5);
        string path = _dir.File("text");
        File.WriteAllBytes(path, TextOracle.Pack(symbols, 2));
        int before = MappedFileRegistry.Count;
        MappedTextReader reader = MappedTextReader.Open(path, symbols.Length, 2);
        Assert.That(MappedFileRegistry.Count, Is.EqualTo(before + 1));
        Assert.Throws<EsaMemoryViolationException>(() => new CountedFile(path, FileMode.Create, FileAccess.ReadWrite));
        reader.Dispose();
        reader.Dispose();
        Assert.That(MappedFileRegistry.Count, Is.EqualTo(before), "a second Dispose takes nothing off twice");
        using var created = new CountedFile(path, FileMode.Create, FileAccess.ReadWrite);
        Assert.That(created.Length, Is.Zero);
    }

    [Test]
    public void TwoMappingsOfOneFile_AreCounted_TheFileIsFreeWhenTheLastIsDisposed()
    {
        RequireChecked();
        string path = RankFile(out DifferenceCover cover, out long n);
        RankTable first = RankTable.Open(path, cover, n);
        RankTable second = RankTable.Open(path, cover, n);
        first.Dispose();
        Assert.Throws<EsaMemoryViolationException>(() => new CountedFile(path, FileMode.Create, FileAccess.ReadWrite), "one mapping is still alive");
        second.Dispose();
        using var created = new CountedFile(path, FileMode.Create, FileAccess.ReadWrite);
        Assert.That(created.Length, Is.Zero);
    }

    [Test]
    public void AReaderOverMemory_AndAFailedOpen_ListNothing()
    {
        RequireChecked();
        int before = MappedFileRegistry.Count;
        using (var t = new PinnedText(TextOracle.Random(40, 4, 3), 2)) Assert.That(MappedFileRegistry.Count, Is.EqualTo(before));
        string path = RankFile(out DifferenceCover cover, out long n);
        Assert.Throws<InvalidDataException>(() => RankTable.Open(path, cover, n + 100), "a file of another length");
        Assert.Throws<InvalidDataException>(() => MappedTextReader.Open(path, 1000, 2), "a file that is no packed text of that length");
        Assert.That(MappedFileRegistry.Count, Is.EqualTo(before));
    }

    [Test]
    public void TheCreationOfAFileThatIsNotMapped_IsNotRefused_InBothBuilds()
    {
        string path = _dir.File("free");
        using (var created = new CountedFile(path, FileMode.Create, FileAccess.ReadWrite)) created.WriteAt(new byte[8], 0);
        using (var again = new CountedFile(path, FileMode.Create, FileAccess.ReadWrite)) Assert.That(again.Length, Is.Zero);
        using MappedOutputFile file = MappedOutputFile.Create(path, 4096);
        Assert.That(file.Length, Is.EqualTo(4096));
    }
}
