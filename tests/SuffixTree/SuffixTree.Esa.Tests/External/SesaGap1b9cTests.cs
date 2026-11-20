using System.IO.MemoryMappedFiles;
using SuffixTree.Esa.External;
using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.Esa.Tests.External;

/// <summary>
/// The gaps that the complete mutation audit of work package 1b.9c showed in the code of 1b.7e (PH6, the container writer and the end-to-end builder, <c>External/Sesa</c>, <c>External/Ph6</c>, <c>External/Lcp</c>): the mutants
/// of the lists that the audit of 1b.7e did not run and the "second line of defence" mutants that it only argued. Each test names the mutants it kills (ids of <c>analysis/mutation/1b7e_*_listed.json</c>); every one of them
/// was run against the test before it counted.
/// </summary>
[TestFixture]
public class SesaGap1b9cTests
{
    // ---------------------------------------------------------------- V11: a failed output writer of B3 leaves no writer open on the file of the container

    [Test]
    public void AFailureInTheConstructionOfTheOutputWriterOfB3_LeavesNoWriterOpenOnTheContainer()
    {
        // V11 (the cleanup of the constructor of the output writer): the buffers of B3 are larger than those of a B2 task here, so a cap between the two fails the second allocation of the output writer, after both of its ranges of the
        // container were opened; the file is the caller's, who must be able to complete it (which refuses while a writer is open)
        using var dir = new TempDir();
        LcpChildCase c = LcpChildCase.Of(StringFamilies.RandomDna(300, 3));
        (string sa, string textPath) = LcpChildRun.WriteInputs(dir, c);
        using MappedTextReader text = MappedTextReader.Open(textPath, c.N, c.Bits);
        using MappedOutputFile container = MappedOutputFile.Create(dir.File("container.bin"), 1 << 16, null);
        var options = new ExtLcpOptions
        {
            ScratchDir = dir.Sub("scratch"),
            ChunkPositions = 40,
            Threads = 1,
            Workers = 1,
            ReadBufferBytes = 4096,
            TupleBufferBytes = 48,
            WindowBytes = 128,
            ArenaGuards = ArenaOptions.Guarded(0xCD),
            Budget = new ArenaBudget(8_000),
            Place = _ => new LcpPlacement(container, 0, 1 << 15),
        };
        ArenaBudgetExceededException? e = Assert.Throws<ArenaBudgetExceededException>(() => ExternalLcp.Build(sa, text, string.Empty, string.Empty, options));
        Assert.That(e!.Message, Does.Contain("lcp.pairs"), "the failure is the allocation of the output writer, after B2 and after the first buffer of B3");
        Assert.That(() => container.Complete(false), Throws.Nothing, "no writer of the failed build is open on the file");
        Assert.That(options.Budget!.LiveBytes, Is.Zero);
    }

    // ---------------------------------------------------------------- B104: the file options reach the BWT file

    /// <summary>A reservation source that records the path of every file created through it (the platform side of <see cref="MappedFileOptions"/>).</summary>
    private sealed class RecordingSource : IReservationSource
    {
        private readonly List<string> _created = [];

        public bool CanPreallocate => true;

        internal IReadOnlyList<string> Created
        {
            get
            {
                lock (_created) return [.. _created];
            }
        }

        public FileStream CreatePreallocated(string path, long length)
        {
            lock (_created) _created.Add(path);
            return new FileStream(path, FileMode.Create, FileAccess.ReadWrite, FileShare.ReadWrite);
        }

        public long AvailableFreeBytes(string directory) => -1;
    }

    [Test]
    public void TheFileOptions_ReachTheBwtFile_AndEveryOtherFileOfTheBuildThatIsReserved()
    {
        // B104 (the BWT file gets no file options): every file that the build reserves goes through the source of the options; the BWT scratch is one of them, and only a build with PSI has it
        using var dir = new TempDir();
        var psi = new RecordingSource();
        SesaRun.Build(dir, "mississippimississippi", SesaRun.Tiny(3) with { Files = new MappedFileOptions { Source = psi } }, "-p");
        Assert.That(psi.Created.Count(p => p.EndsWith("bwt.bin", StringComparison.Ordinal)), Is.EqualTo(1), "the BWT scratch file of PH6");
        Assert.That(psi.Created.Count(p => p.EndsWith("index.sesa.tmp", StringComparison.Ordinal)), Is.EqualTo(1), "the container");
        var none = new RecordingSource();
        SesaRun.Build(dir, "mississippimississippi", SesaRun.Tiny(3, psi: false) with { Files = new MappedFileOptions { Source = none } }, "-n");
        Assert.That(none.Created.Count(p => p.EndsWith("bwt.bin", StringComparison.Ordinal)), Is.Zero, "no PSI, no BWT file");
        Assert.That(psi.Created.Count - none.Created.Count, Is.EqualTo(1), "the BWT file is the only difference");
    }

    // ---------------------------------------------------------------- P63, P64: the native blocks of the Psi cursors

    [Test]
    public void ThePsiStream_GivesItsNativeBlocksBack_AfterEveryPass_AndAfterAFailureInThePass()
    {
        // P63 (the cursor array is not freed) and P64 (the end array is not freed): a raw native leak that nothing in the runtime shows; the stream counts its live blocks
        using var dir = new TempDir();
        long before = PsiStream.LiveNativeBlocks;
        byte[] s = TextOracle.Random(300, 4, 9);
        Ph6Run.Outcome o = Ph6Run.Run(dir, s, 4, psi: new PsiOptions { SymbolsPerPass = 1, WindowBytes = 16 }, tag: "p");
        Assert.That(PsiStream.LiveNativeBlocks, Is.EqualTo(before), "four passes, every one gave its two blocks back");
        long[] counts = [.. o.Bwt.Counts];
        counts[0]++;
        counts[1]--;
        using MappedOutputFile file = Ph6Run.CreateContainer(dir.File("container-f.bin"), o.Sa, 64, o.Container.Length);
        Assert.Throws<InvalidDataException>(() => PsiStream.Build(o.BwtPath, s.Length, o.Bwt.PrimaryRow, counts, file, 192, new PsiOptions { SymbolsPerPass = 2, WindowBytes = 16 }), "a count that the BWT does not have");
        Assert.That(PsiStream.LiveNativeBlocks, Is.EqualTo(before), "a failed pass gave its blocks back as well");
    }

    // ---------------------------------------------------------------- B182, B183: the two fsyncs of the publication

    [Test]
    public void ThePublication_SyncsOnceAfterTheHeaderAndOnceAfterTheMagic_AndNotAtAllWithoutTheOption()
    {
        // B182 (the first sync is skipped) and B183 (the second one is skipped): the calls of CountedFile.Flush between the hook points of the publication, one each with the option, none without it
        foreach (bool fsync in new[] { true, false })
        {
            using var dir = new TempDir();
            var counters = new IoCounters();
            var at = new Dictionary<string, long>();
            SesaBuildOptions options = SesaRun.Tiny(3) with { Fsync = fsync, Counters = counters, Hook = point => at[point] = counters.Flushes };
            SesaRun.Build(dir, "mississippimississippi", options);
            Assert.That(at[SesaHooks.HeaderSynced] - at[SesaHooks.HeaderWritten], Is.EqualTo(fsync ? 1 : 0), $"fsync {fsync}: between the header and its sync");
            Assert.That(at[SesaHooks.MagicWritten] - at[SesaHooks.HeaderSynced], Is.EqualTo(fsync ? 1 : 0), $"fsync {fsync}: between the sync of the header and the sync of the magic");
        }
    }

    // ---------------------------------------------------------------- K33: the windows of the child table are flushed when released

    [Test]
    public void TheChildTableWindows_AreFlushedWhenTheyAreReleased()
    {
        // K33 (PH5 does not flush windows): with windows of a few entries and one slot every window that is released is msynced, so the number of flushes is at least that of the evictions
        using var dir = new TempDir();
        SesaRun.Outcome o = SesaRun.Build(dir, StringFamiliesText(), SesaRun.Tiny(2));
        WindowStats w = o.Result.Lcp.Scan.ChildWindows;
        Assert.That(w.Evictions, Is.GreaterThan(0), "the geometry makes the windows move");
        Assert.That(w.Flushes, Is.GreaterThanOrEqualTo(w.Evictions), "every released window is flushed");
        Assert.That(w.FlushedBytes, Is.GreaterThan(0));
    }

    private static string StringFamiliesText() => SaRun.Families(1_500).First(f => f.Name == "repeat-rich").Text;

    // ---------------------------------------------------------------- S10: the seconds of a section of the checksum pass

    [Test]
    public void TheChecksumPass_MeasuresTheSecondsOfEverySection_AndKeepsItsBytesAndKind()
    {
        // S10 (the seconds are not measured): a section of 1 MiB or more takes the CRC pass some tens of microseconds at least, which the clock of the platform resolves
        using var dir = new TempDir();
        string path = dir.File("f.bin");
        var data = new byte[4 << 20];
        new Random(5).NextBytes(data);
        File.WriteAllBytes(path, data);
        (uint Crc, SesaCrcRow Row)[] rows = SesaChecksums.Compute(path, [(7, 0, 1 << 20), (9, 1 << 20, 3 << 20)]);
        Assert.That(rows.Select(r => (r.Row.Kind, r.Row.Bytes)), Is.EqualTo(new[] { (7, 1L << 20), (9, 3L << 20) }));
        foreach ((uint crc, SesaCrcRow row) in rows)
        {
            Assert.That(row.Seconds, Is.GreaterThan(0).And.LessThan(60), $"kind {row.Kind}");
            Assert.That(row.Hardware, Is.EqualTo(Crc32C.HardwareAccelerated), $"kind {row.Kind}");
        }
        Assert.That(rows[0].Crc, Is.EqualTo(SesaOracle.Crc32cBitwise(data.AsSpan(0, 1 << 20))));
    }

    // ---------------------------------------------------------------- K23, K37: the guard option reaches PH4 and PH5 one by one

    [Test]
    public void TheGuardOption_ReachesTheB1B2AndB3PassesOfPh4_AndThePassOfPh5_EachOnItsOwn()
    {
        // K23 (PH4 gets no guards) and K37 (PH5 gets no guards): the passes are told apart by the trace of PH4 (b1, b2, b3) and the hook after PH5; every arena of a pass is guarded when the option is passed
        using var dir = new TempDir();
        var at = new Dictionary<string, long>();
        SesaBuildOptions options = SesaRun.Tiny(3) with
        {
            LcpTrace = name => at.TryAdd(name, NativeArena.GuardedCreated),
            Hook = point => at.TryAdd("hook:" + point, NativeArena.GuardedCreated),
        };
        SesaRun.Build(dir, "mississippimississippi" + "abracadabra" + "banana", options, "-g");
        long Between(string from, string to) => at[to] - at[from];
        Assert.That(Between("b1", "b2"), Is.EqualTo(1), "B1: one guarded arena");
        Assert.That(Between("b2", "b3"), Is.EqualTo(1), "B2: one guarded arena per worker, one worker here");
        Assert.That(Between("b3", "hook:" + SesaHooks.LcpChildBuilt), Is.EqualTo(2), "the emission of B3 and the scan of PH5: one guarded arena each");
    }

    // ---------------------------------------------------------------- B124: a container that is longer than planned is refused

    [Test]
    public void TheLengthOfTheClosedContainer_IsCheckedAgainstThePlan_InBothDirections()
    {
        // B124 and the mirror of it (a shorter container is not refused): the check as a function, which no build can reach with a shorter file
        SesaBuilder.CheckContainerLength(0, 0);
        SesaBuilder.CheckContainerLength(123_456, 123_456);
        foreach ((long actual, long planned) in new[] { (4_095L, 4_096L), (4_097L, 4_096L), (0L, 1L), (1L, 0L), (long.MaxValue, 1L), (1L, long.MaxValue) })
        {
            Assert.That(() => SesaBuilder.CheckContainerLength(actual, planned), Throws.TypeOf<InvalidOperationException>().With.Message.Contains($"the container has {actual} bytes, expected {planned}"), $"{actual} against {planned}");
        }
    }

    [Test]
    public void AContainerThatIsLongerThanThePlan_IsRefusedBeforeTheChecksums()
    {
        // B124 (the length check is "smaller than"): the container is checked against the length that the plan gave it; the hook after BUILD_INFO makes the temporary file one byte longer
        using var dir = new TempDir();
        string input = SesaRun.WriteInput(dir, "mississippimississippi");
        string target = dir.File("index.sesa");
        var points = new List<string>();
        SesaBuildOptions options = SesaRun.Tiny(3) with
        {
            Hook = point =>
            {
                points.Add(point);
                if (point != SesaHooks.InfoWritten) return;
                using var f = new FileStream(target + ".tmp", FileMode.Open, FileAccess.ReadWrite, FileShare.ReadWrite);
                f.SetLength(f.Length + 1);
            },
        };
        InvalidOperationException? e = Assert.Throws<InvalidOperationException>(() => SesaBuilder.Build(input, target, options));
        Assert.That(e!.Message, Does.Contain("the container has").And.Contain("expected"));
        Assert.That(points, Does.Not.Contain(SesaHooks.Checksummed), "refused before the checksums");
        Assert.That(File.Exists(target) || File.Exists(target + ".tmp") || Directory.Exists(target + ".tmp.d"), Is.False, "nothing is left");
    }

    // ---------------------------------------------------------------- B94: the counts of the BWT are compared with those of the text

    [Test]
    public void ASuffixArrayThatDoesNotBelongToTheText_IsRefusedByPh6_BeforePsiIsBuilt()
    {
        // B94: after PH4 / PH5 one entry of the suffix array in the temporary container is replaced by the next one: both are in range, neither is the empty suffix or the suffix 0, so the BWT pass takes the array, but the
        // symbols that precede the two positions differ and the counts of the BWT are not those of the text
        using var dir = new TempDir();
        const string text = "mississippimississippi";
        string input = SesaRun.WriteInput(dir, text);
        string target = dir.File("index.sesa");
        SaCase c = SaCase.Of(text);
        uint[] sa = c.Sa;
        int k = Enumerable.Range(1, text.Length - 1).First(i => sa[i] != 0 && sa[i + 1] != 0 && text[(int)sa[i] - 1] != text[(int)sa[i + 1] - 1]);
        SesaPlan plan = SesaLayout.Plan(PackedTextWriter.Scan(input), true);
        SesaBuildOptions options = SesaRun.Tiny(3) with
        {
            Hook = point =>
            {
                if (point != SesaHooks.LcpChildBuilt) return;
                using var f = new FileStream(target + ".tmp", FileMode.Open, FileAccess.ReadWrite, FileShare.ReadWrite);
                f.Seek(plan.SaOffset + (4L * k), SeekOrigin.Begin);
                f.Write(BitConverter.GetBytes(sa[k + 1]));
            },
        };
        InvalidDataException? e = Assert.Throws<InvalidDataException>(() => SesaBuilder.Build(input, target, options));
        Assert.That(e!.Message, Does.Contain("symbol counts of the BWT"));
        Assert.That(File.Exists(target) || File.Exists(target + ".tmp") || Directory.Exists(target + ".tmp.d"), Is.False, "nothing is left");
        // the control: the same build without the damage succeeds
        Assert.That(() => SesaBuilder.Build(input, target, SesaRun.Tiny(3)), Throws.Nothing);
    }

    // ---------------------------------------------------------------- B90: the option IncludePsi decides, at every thread count

    [Test]
    public void WithoutIncludePsi_NoPsiSectionIsWritten_AtAnyThreadCount_AndWithItTheSectionIsThere()
    {
        // B90 ("IncludePsi || Threads > 1000"): the plan, not another option, decides whether PH6 runs
        const string text = "mississippimississippi";
        foreach (int threads in new[] { 1, 2, 3, 1_001 })
        {
            using var dir = new TempDir();
            var points = new List<string>();
            SesaRun.Outcome without = SesaRun.Build(dir, text, SesaRun.Tiny(3, threads, psi: false) with { Hook = points.Add }, "-n");
            Assert.That(without.Oracle.Find(SesaFormat.KindPsi), Is.Null, $"{threads} threads, no PSI asked: no PSI section");
            Assert.That(without.Result.HasPsi, Is.False, $"{threads} threads: the result says no PSI");
            Assert.That(without.Result.Bwt, Is.Null, $"{threads} threads: no BWT pass");
            Assert.That(without.Result.Psi, Is.Null, $"{threads} threads: no PSI pass");
            Assert.That(points, Does.Not.Contain(SesaHooks.BwtBuilt).And.Not.Contain(SesaHooks.PsiBuilt), $"{threads} threads: no BWT, no PSI");
            SesaRun.Compare(without, $"{threads} threads, no PSI", psi: false);
            SesaRun.Outcome with = SesaRun.Build(dir, text, SesaRun.Tiny(3, threads, psi: true), "-p");
            Assert.That(with.Oracle.Find(SesaFormat.KindPsi), Is.Not.Null, $"{threads} threads, PSI asked: a PSI section");
        }
    }

    // ---------------------------------------------------------------- V10: a failed overflow writer of B3 leaves no base writer open on the file of the container

    [Test]
    public void AnOverflowRangeThatLiesOutsideTheContainer_LeavesNoWriterOpenOnTheFile()
    {
        // V10 (a failed overflow writer leaves the base writer open): the placement puts the base array inside the container and the overflow list behind its end, so the second writer of B3 is refused after the first was opened; the file
        // is the caller's, who must be able to complete it (which refuses while a writer is open)
        using var dir = new TempDir();
        LcpChildCase c = LcpChildCase.Of(StringFamilies.RandomDna(300, 3));
        (string sa, string textPath) = LcpChildRun.WriteInputs(dir, c);
        using MappedTextReader text = MappedTextReader.Open(textPath, c.N, c.Bits);
        using MappedOutputFile container = MappedOutputFile.Create(dir.File("container.bin"), 1 << 16, null);
        var options = new ExtLcpOptions
        {
            ScratchDir = dir.Sub("scratch"),
            ChunkPositions = 40,
            Threads = 1,
            Workers = 1,
            ReadBufferBytes = 4096,
            TupleBufferBytes = 48,
            WindowBytes = 128,
            ArenaGuards = ArenaOptions.Guarded(0xCD),
            Place = _ => new LcpPlacement(container, 0, (1 << 16) + 1),
        };
        Assert.Throws<ArgumentOutOfRangeException>(() => ExternalLcp.Build(sa, text, string.Empty, string.Empty, options));
        Assert.That(() => container.Complete(false), Throws.Nothing, "no writer of the failed build is open on the file");
    }

    // ---------------------------------------------------------------- G16: a growth whose mapping fails is undone

    [Test]
    public void AGrowthWhoseMappingFails_LeavesTheFileAtItsOldLength_AndTheOldMappingUsable()
    {
        // G16 (the failed mapping does not undo the growth): the file is extended before it is mapped again, and a failure of the mapping must give the length back (the doc comment of Grow: "a failure leaves the file at its old length")
        using var dir = new TempDir();
        string path = dir.File("grow.bin");
        using MappedOutputFile file = MappedOutputFile.Create(path, 4096, null);
        var failing = new MappedFileOptions { Mapper = _ => throw new IOException("the mapping failed") };
        var e = Assert.Throws<IOException>(() => file.Grow(40_000, failing));
        Assert.That(e!.Message, Does.Contain("the mapping failed"));
        Assert.That(new FileInfo(path).Length, Is.EqualTo(4096), "the file has its old length");
        Assert.That(file.Length, Is.EqualTo(4096));
        using (MappedWindowWriter w = file.CreateWriter(0, 4096, 256))
        {
            w.Append(new byte[4096]);
        }
        file.Complete(false);
        // the control: without the failure the same growth works
        file.Grow(40_000);
        Assert.That((file.Length, new FileInfo(path).Length), Is.EqualTo((40_000L, 40_000L)));
    }

    // ---------------------------------------------------------------- W06: the window of the BWT pass is checked at the entry

    [Test]
    public void TheBwtPass_RefusesAWindowOfZeroBytes_BeforeItOpensTheSuffixArray()
    {
        // W06 (a window of 0 is accepted): the suffix array that does not exist would be the first error of a pass that went on past its entry; the entry check says what is wrong with the options
        using var dir = new TempDir();
        string textPath = SesaRun.WriteInput(dir, "banana");
        LcpChildCase c = LcpChildCase.Of("banana");
        string packed = TextOracle.WriteFile(dir, "t.pk", c.Symbols, c.Bits);
        using MappedTextReader text = MappedTextReader.Open(packed, c.N, c.Bits);
        string missing = dir.File("no-such.sa");
        string bwt = dir.File("o.bwt");
        Assert.That(File.Exists(textPath), Is.True);
        var e = Assert.Throws<ArgumentOutOfRangeException>(() => BwtPass.Build(missing, 0, text, c.Sigma, bwt, new BwtOptions { WindowBytes = 0 }));
        Assert.That(e!.ParamName, Is.EqualTo("options.WindowBytes"));
        // the control: a window of one byte passes the entry and fails where the suffix array is opened
        Assert.Throws<FileNotFoundException>(() => BwtPass.Build(missing, 0, text, c.Sigma, bwt, new BwtOptions { WindowBytes = 1 }));
        Assert.That(File.Exists(bwt), Is.False);
    }

    // ---------------------------------------------------------------- B258: the default read of the BWT pass is a mebibyte

    [Test]
    public void TheBwtPass_ReadsTheSuffixArrayInOneMebibyteStepsByDefault_AndInTheStepOfTheTuningWhenThereIsOne()
    {
        // B258 (the default read of the BWT pass is smaller): 300,000 symbols, a suffix array of 1,200,004 bytes: two reads of it with the default of 2^20 bytes (262,144 entries), 147 reads of 2,048 entries with a tuning of 8 KiB; the reads
        // of the pass are those between the hooks "lcp-child-built" and "bwt-built" (the pass opens the file of the suffix array once and reads nothing else through the counters)
        string text = StringFamilies.RandomDna(300_000, 3);
        foreach ((int? readBytes, long reads) in new[] { ((int?)null, 2L), (8192, 147L) })
        {
            using var dir = new TempDir();
            var counters = new IoCounters();
            long before = -1, after = -1;
            SesaBuildOptions options = new()
            {
                Threads = 1,
                Fsync = false,
                Counters = counters,
                Ph6Tuning = new SesaPh6Tuning { ReadBytes = readBytes },
                Hook = point =>
                {
                    if (point == SesaHooks.LcpChildBuilt) before = counters.ReadCalls;
                    else if (point == SesaHooks.BwtBuilt) after = counters.ReadCalls;
                },
            };
            SesaRun.Outcome o = SesaRun.Build(dir, text, options, "-r");
            SesaRun.Compare(o, $"read {readBytes}", psi: true);
            Assert.That(after - before, Is.EqualTo(reads), $"reads of the suffix array with ReadBytes = {readBytes?.ToString() ?? "default"}");
        }
    }
}
