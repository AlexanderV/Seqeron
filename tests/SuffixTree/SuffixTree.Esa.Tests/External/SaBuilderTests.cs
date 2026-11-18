using System.Globalization;
using SuffixTree.Esa.External;
using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.Esa.Tests.External;

/// <summary>
/// The stage of work package 1b.7d as a whole (PH1, PH2 and PH3 on a packed text: the suffix array file): the contract against the in-RAM oracles on small texts and on the families, the thread counts, the bucket
/// capacities, the failure modes, the counters, the verdict of the RssAnon gate and the platform without <c>/proc</c> or <c>fallocate</c>.
/// </summary>
[TestFixture]
public class SaBuilderTests
{
    private TempDir _dir = null!;

    [SetUp]
    public void SetUp() => _dir = new TempDir();

    [TearDown]
    public void TearDown() => _dir.Dispose();

    private string Scratch => _dir.Sub("scratch");

    private SaBuildOptions Tiny(int v, int threads = 1) => SaTiny.Build(Scratch, v, threads);

    [Test]
    public void SmallTexts_AreTheSuffixArray_WithRotatingTinyParameters()
    {
        int v = 0;
        foreach (string text in new[] { "a", "ab", "ba", "aa", "banana", "mississippi", "abracadabra", new string('a', 40), "abababababababababab", "abcabcabcabcabcabcabc", "aabaabaabaabaabaab", StringFamilies.Fibonacci(90), StringFamilies.ThueMorse(75), StringFamilies.RandomDna(120, 4), StringFamilies.RandomProtein(100, 5) })
        {
            SaCase c = SaCase.Of(text);
            foreach (int variant in new[] { v, v + 4, v + 9 })
            {
                using var dir = new TempDir();
                SaRun.Outcome o = SaRun.Run(dir, c, SaTiny.Build(dir.Sub("scratch"), variant, 1 + (variant % 4)));
                SaRun.Compare(c, o, $"'{(text.Length > 24 ? text[..24] : text)}' variant {variant}");
            }
            v++;
        }
    }

    [Test]
    public void TheEmptyText_AndTextsOfOneAndTwoSymbols_Work()
    {
        foreach (string text in new[] { string.Empty, "x", "xx", "xy", "yx" })
        {
            SaCase c = SaCase.Of(text);
            SaRun.Outcome o = SaRun.Run(_dir, c, Tiny(1), tag: text.Length.ToString(CultureInfo.InvariantCulture) + text);
            SaRun.Compare(c, o, $"'{text}'");
            Assert.That(new FileInfo(o.SaPath).Length, Is.EqualTo(4 * (text.Length + 1)));
            Assert.That(o.Result.N, Is.EqualTo(text.Length));
            Directory.Delete(o.OutDir, true);
        }
    }

    [Test]
    public void Families_AtEveryTinyVariant_AreTheSuffixArray()
    {
        int v = 0;
        foreach (int n in new[] { 600, 3000 })
        {
            foreach ((string name, string text) in SaRun.Families(n))
            {
                SaCase c = SaCase.Of(text);
                using var dir = new TempDir();
                SaRun.Outcome o = SaRun.Run(dir, c, SaTiny.Build(dir.Sub("scratch"), v, 1 + (v % 4)));
                SaRun.Compare(c, o, $"{name} {n} variant {v}");
                v++;
            }
        }
    }

    [Test]
    public void Families_WithTheDefaultParameters_AndWithTheScaledGeometry_AreTheSuffixArray()
    {
        foreach ((string name, string text) in SaRun.Families(40_000))
        {
            SaCase c = SaCase.Of(text);
            foreach (long scale in new long[] { 0, 3_000_000_000 })
            {
                using var dir = new TempDir();
                SaRun.Outcome o = SaRun.Run(dir, c, new SaBuildOptions { ScratchDir = dir.Sub("scratch"), ScaleN = scale });
                SaRun.Compare(c, o, $"{name} 40000, ScaleN {scale}");
                SaBuildParameters p = o.Result.Parameters;
                if (scale == 0) Assert.That((p.BucketCapacity, p.DcsBaseCase, p.DcsBucketCapacity), Is.EqualTo((1L << 20, 16L << 20, 1L << 20)));
                else Assert.That(p.BucketCapacity, Is.LessThan(100), "the scaled geometry splits a text of 40000 symbols into many buckets");
            }
        }
    }

    [Test]
    public void TheOutput_IsIdenticalForEveryThreadCount_BucketCapacity_AndSeed()
    {
        foreach ((string name, string text) in SaRun.Families(3_000).Where(f => f.Name is "a^n" or "fibonacci" or "tandem" or "random-dna" or "genome-like" or "n-runs"))
        {
            SaCase c = SaCase.Of(text);
            string? reference = null;
            int v = 0;
            foreach (int threads in new[] { 1, 2, 4 })
            {
                foreach (long k in new long[] { 9, 400, 5000 })
                {
                    using var dir = new TempDir();
                    SaBuildOptions o = SaTiny.Build(dir.Sub("scratch"), v++, threads) with { Seed = (ulong)(v * 13) };
                    o = o with { Tuning = o.Tuning with { BucketCapacity = k } };
                    SaRun.Outcome r = SaRun.Run(dir, c, o);
                    string sha = SaRun.Sha(r.SaPath);
                    if (reference is null) SaRun.Compare(c, r, $"{name}: first");
                    reference ??= sha;
                    Assert.That(sha, Is.EqualTo(reference), $"{name}: {threads} threads, K {k}");
                }
            }
        }
    }

    [Test]
    public void TwoRuns_GiveTheSameFile_AndTheSameStatistics()
    {
        SaCase c = SaCase.Of(SaRun.GenomeLike(8000, 5));
        SaRun.Outcome a = SaRun.Run(_dir, c, Tiny(5, 3), "a");
        using var dir2 = new TempDir();
        SaRun.Outcome b = SaRun.Run(dir2, c, SaTiny.Build(dir2.Sub("scratch"), 5, 3), "a");
        Assert.That(SaRun.Sha(a.SaPath), Is.EqualTo(SaRun.Sha(b.SaPath)));
        Assert.That(a.Result.Dcs.Levels.Select(l => (l.TextLength, l.DistinctNames, l.Leaves)), Is.EqualTo(b.Result.Dcs.Levels.Select(l => (l.TextLength, l.DistinctNames, l.Leaves))));
        Assert.That((a.Result.Sort.Leaves, a.Result.Sort.Partitions, a.Result.Sort.MaxDepth), Is.EqualTo((b.Result.Sort.Leaves, b.Result.Sort.Partitions, b.Result.Sort.MaxDepth)));
    }

    [Test]
    public void TheResult_CarriesTheParameters_TheStats_AndThePeaks()
    {
        SaCase c = SaCase.Of(SaRun.GenomeLike(6000, 3));
        SaBuildOptions o = Tiny(4, 2) with { ArenaGuards = null };
        SaRun.Outcome r = SaRun.Run(_dir, c, o);
        SaBuildResult res = r.Result;
        Assert.That(res.N, Is.EqualTo(c.N));
        Assert.That(res.Parameters.N, Is.EqualTo(c.N));
        Assert.That(res.Parameters.BucketCapacity, Is.EqualTo(SaTiny.Tuning(4).BucketCapacity));
        Assert.That(res.Period, Is.EqualTo(SaTiny.Periods[4 % SaTiny.Periods.Length]));
        Assert.That(res.CoverResidues, Is.EqualTo(DifferenceCover.Create(res.Period).Residues.Length));
        Assert.That(res.CoverKind, Is.EqualTo(DifferenceCover.Create(res.Period).Kind));
        Assert.That(res.Ph1Seconds, Is.GreaterThan(0));
        Assert.That(res.Ph23Seconds, Is.GreaterThan(0));
        Assert.That(res.Seconds, Is.EqualTo(res.Ph1Seconds + res.Ph23Seconds));
        Assert.That(res.Ph1ArenaPeakBytes, Is.GreaterThan(0));
        Assert.That(res.Ph23ArenaPeakBytes, Is.EqualTo(24L * SaTiny.Tuning(4).BucketCapacity!.Value).Or.EqualTo(48L * SaTiny.Tuning(4).BucketCapacity!.Value), "keys, keys2, positions, positions2 of K entries, one leaf sorter per worker that sorted a leaf");
        Assert.That(res.RankFileBytes, Is.EqualTo(4L * new SampleGeometry(c.N, DifferenceCover.Create(res.Period)).Slots));
        Assert.That(res.Dcs.Levels, Is.Not.Empty);
        Assert.That(res.RecursionLevels, Is.EqualTo(res.Dcs.MaxLevel));
        Assert.That(res.Sort.Leaves, Is.GreaterThan(1));
        Assert.That(res.Ph1ScratchPeakBytes, Is.GreaterThanOrEqualTo(0));
        Assert.That(res.Ph23ScratchPeakBytes, Is.GreaterThanOrEqualTo(res.Ph1ScratchPeakBytes));
        Assert.That(res.Reservation, Is.EqualTo(ReservationOutcome.Preallocated));
    }

    [Test]
    public void TheCounters_AreTheOutputsTheSpillAndTheRankFile_AndThePhasesAddUp()
    {
        SaCase c = SaCase.Of(SaRun.GenomeLike(5000, 9));
        var counters = new IoCounters();
        SaBuildOptions o = Tiny(6, 2) with { Counters = counters, ArenaGuards = null };
        SaRun.Outcome r = SaRun.Run(_dir, c, o);
        SaBuildResult res = r.Result;
        long slots = res.RankFileBytes;
        Assert.That(res.Ph23Io.MappedWriteBytes, Is.EqualTo(4L * (c.N + 1)), "every entry of the suffix array was written once, through the windows");
        Assert.That(res.Ph1Io.MappedWriteBytes, Is.GreaterThanOrEqualTo(slots), "the rank file went through the windows");
        Assert.That(res.Io.TotalBytes, Is.EqualTo(res.Ph1Io.TotalBytes + res.Ph23Io.TotalBytes));
        Assert.That(res.Io, Is.EqualTo(counters.Snapshot()));
        Assert.That(counters.ScratchLiveBytes, Is.Zero, "the rank file and every scratch file are given back");
        Assert.That(counters.ScratchPeakBytes, Is.GreaterThanOrEqualTo(slots), "the rank file is scratch while PH3 reads it");
        Assert.That(res.Ph23Io.ReadBytes, Is.GreaterThan(0), "the position files are read");
        Assert.That(res.Ph23ScratchPeakBytes, Is.EqualTo(counters.ScratchPeakBytes));
    }

    [Test]
    public void ThePrivateCounters_AreUsed_WhenTheCallerGivesNone()
    {
        SaCase c = SaCase.Of(StringFamilies.RandomDna(2000, 3));
        SaRun.Outcome r = SaRun.Run(_dir, c, Tiny(2) with { Counters = null });
        Assert.That(r.Result.Io.MappedWriteBytes, Is.GreaterThan(0));
    }

    [Test]
    public void OneSharedCounter_ServesSeveralBuilds_AndEachResultIsItsOwnDifference()
    {
        SaCase c = SaCase.Of(StringFamilies.RandomDna(2000, 3));
        var counters = new IoCounters();
        SaRun.Outcome a = SaRun.Run(_dir, c, Tiny(2) with { Counters = counters }, "a");
        using var dir2 = new TempDir();
        SaRun.Outcome b = SaRun.Run(dir2, c, SaTiny.Build(dir2.Sub("scratch"), 2) with { Counters = counters }, "b");
        Assert.That(b.Result.Ph23Io.MappedWriteBytes, Is.EqualTo(a.Result.Ph23Io.MappedWriteBytes), "a phase's figure is a difference, not the running total");
        Assert.That(counters.MappedWriteBytes, Is.EqualTo(a.Result.Io.MappedWriteBytes + b.Result.Ph1Io.MappedWriteBytes + b.Result.Ph23Io.MappedWriteBytes));
    }

    [Test]
    public void TheMappedWindows_AreAtLeastOneMebibyte_WithTheDefaultParameters()
    {
        SaCase c = SaCase.Of(StringFamilies.RandomDna(2_500_000, 11));
        SaRun.Outcome r = SaRun.Run(_dir, c, new SaBuildOptions { ScratchDir = Scratch });
        SaRun.Compare(c, r, "2.5 M");
        Assert.That(r.Result.Io.MappedMaxWindowBytes, Is.EqualTo(1 << 20), "the windows of the suffix array and of the rank file are 1 MiB (the E9 rule)");
        Assert.That(r.Result.Parameters.WindowBytes, Is.EqualTo(1 << 20));
    }

    [Test]
    public void TheRankFile_IsDeletedWhenPh3IsDone_AndLivesInTheScratchDirectory_OrAtTheGivenPath()
    {
        SaCase c = SaCase.Of(StringFamilies.Fibonacci(3000));
        long seenLength = -1;
        string? seenPath = null;
        string scratch = Scratch;
        SaBuildOptions o = Tiny(3) with
        {
            Trace = step =>
            {
                if (step == "ph23" && File.Exists(Path.Combine(scratch, "ranks.bin")))
                {
                    seenPath = Path.Combine(scratch, "ranks.bin");
                    seenLength = new FileInfo(seenPath).Length;
                }
            },
        };
        SaRun.Outcome r = SaRun.Run(_dir, c, o);
        SaRun.Compare(c, r, "rank file");
        Assert.That(seenPath, Is.Not.Null, "the rank file existed when PH3 began");
        Assert.That(seenLength, Is.EqualTo(r.Result.RankFileBytes));
        using var dir2 = new TempDir();
        string custom = dir2.File("custom.rank");
        string customSeen = "";
        SaBuildOptions o2 = SaTiny.Build(dir2.Sub("scratch"), 3) with { RankPath = custom, Trace = step => { if (step == "ph23" && File.Exists(custom)) customSeen = custom; } };
        SaRun.Outcome r2 = SaRun.Run(dir2, c, o2);
        SaRun.Compare(c, r2, "custom rank path");
        Assert.That(customSeen, Is.EqualTo(custom));
        Assert.That(File.Exists(custom), Is.False);
    }

    [Test]
    public void TheRankFileCheck_IsAnOption_AndRunsBetweenTheTwoPhases()
    {
        SaCase c = SaCase.Of(SaRun.GenomeLike(4000, 2));
        var events = new List<string>();
        SaRun.Outcome r = SaRun.Run(_dir, c, Tiny(2) with { VerifyRanks = true, Trace = e => { lock (events) events.Add(e); } });
        SaRun.Compare(c, r, "verified");
        Assert.That(events.Where(e => e is "ph1" or "verify" or "ph23"), Is.EqualTo(new[] { "ph1", "verify", "ph23" }));
        using var dir2 = new TempDir();
        events.Clear();
        SaRun.Outcome q = SaRun.Run(dir2, c, SaTiny.Build(dir2.Sub("scratch"), 2) with { Trace = e => { lock (events) events.Add(e); } });
        SaRun.Compare(c, q, "not verified");
        Assert.That(events, Does.Not.Contain("verify"));
        Assert.That(SaRun.Sha(r.SaPath), Is.EqualTo(SaRun.Sha(q.SaPath)));
    }

    [Test]
    public void TheSuffixArray_IsSyncedWhenAsked_AndIsTheSameFile()
    {
        SaCase c = SaCase.Of(StringFamilies.RandomDna(3000, 5));
        SaRun.Outcome r = SaRun.Run(_dir, c, Tiny(2) with { SyncOutput = true });
        SaRun.Compare(c, r, "synced");
    }

    [Test]
    public void AStaleOutput_IsReplaced_AndAForeignScratchFile_AndTheScratchDirectory_Stay()
    {
        SaCase c = SaCase.Of(StringFamilies.RandomDna(900, 3));
        string outDir = _dir.Sub("out");
        Directory.CreateDirectory(outDir);
        string sa = Path.Combine(outDir, "o.sa");
        File.WriteAllText(sa, "stale and longer than nothing");
        Directory.CreateDirectory(Scratch);
        File.WriteAllText(Path.Combine(Scratch, "foreign.txt"), "mine");
        string textPath = SaRun.WriteText(_dir, c);
        using (MappedTextReader reader = MappedTextReader.Open(textPath, c.N, c.Bits))
        {
            SaBuilder.Build(reader, sa, Tiny(2));
        }
        Assert.That(FileBits.ReadUInts(sa), Is.EqualTo(c.Sa));
        Assert.That(Directory.GetFileSystemEntries(Scratch).Select(Path.GetFileName), Is.EqualTo(new[] { "foreign.txt" }));
    }

    [Test]
    public void TheScratchDirectory_IsCreatedByTheBuild_AndRemovedWithItsSubdirectories()
    {
        SaCase c = SaCase.Of(SaRun.GenomeLike(3000, 4));
        string scratch = _dir.Sub("deep/er/scratch");
        SaRun.Outcome r = SaRun.Run(_dir, c, SaTiny.Build(scratch, 3) with { Tuning = SaTiny.Tuning(3) with { SpillBytes = 16, DcsRunSpillBytes = 0 } });
        SaRun.Compare(c, r, "deep scratch");
        Assert.That(Directory.Exists(scratch), Is.False);
    }

    // ---- failure modes: no output, no rank file, no scratch, no descriptor, no mapping

    [Test]
    public void AnArenaCapThatIsTooSmall_EndsWithTheBudgetException_AtBothPhases()
    {
        SaCase c = SaCase.Of(SaRun.GenomeLike(2000, 3));
        foreach (long cap in new long[] { 0, 64, 400 })
        {
            using var dir = new TempDir();
            SaRun.AssertFailsCleanly(dir, c, SaTiny.Build(dir.Sub("scratch"), 1, 2) with { ArenaCapBytes = cap }, typeof(ArenaBudgetExceededException), $"cap {cap}");
        }
    }

    [Test]
    public void AFullDisk_AtTheCreationOfAnyFile_IsAnIoException_AndNothingIsLeft()
    {
        SaCase c = SaCase.Of(new string('a', 500));
        var counting = new FailingAfter(int.MaxValue);
        using (var dir = new TempDir())
        {
            SaRun.Outcome ok = SaRun.Run(dir, c, SaTiny.Build(dir.Sub("scratch"), 1) with { Files = new MappedFileOptions { Source = counting }, Tuning = SaTiny.Tuning(1) with { DcsBaseCase = 20 } });
            SaRun.Compare(c, ok, "counting");
        }
        int creations = counting.Calls;
        Assert.That(creations, Is.GreaterThanOrEqualTo(3), "the suffix array, the rank file and the files of the levels below");
        for (int succeed = 0; succeed < creations; succeed++)
        {
            using var dir = new TempDir();
            SaBuildOptions o = SaTiny.Build(dir.Sub("scratch"), 1) with { Files = new MappedFileOptions { Source = new FailingAfter(succeed) }, Tuning = SaTiny.Tuning(1) with { DcsBaseCase = 20 } };
            SaRun.AssertFailsCleanly(dir, c, o, typeof(IOException), $"creation {succeed + 1} of {creations} fails");
        }
    }

    [Test]
    public void ASpillThatFails_IsAnIoException_TheBuildEnds_AndNothingIsLeft()
    {
        SaCase c = SaCase.Of(SaRun.GenomeLike(4000, 6));
        foreach (string phase in new[] { "naming", "partition" })
        {
            using var dir = new TempDir();
            string scratch = dir.Sub("scratch");
            bool done = false;
            SaBuildOptions o = SaTiny.Build(scratch, 2, 2) with
            {
                Tuning = SaTiny.Tuning(2) with { SpillBytes = 16, DcsRunSpillBytes = 0 },
                Trace = step =>
                {
                    if (step == phase && !done)
                    {
                        done = true;
                        foreach (string sub in Directory.GetDirectories(scratch)) Directory.Delete(sub, true); // the subdirectory of the phase is gone: no file can be created in it
                    }
                },
            };
            SaRun.AssertFailsCleanly(dir, c, o, typeof(IOException), $"spill failure at {phase}");
        }
    }

    [Test]
    public void AHookThatThrows_ReachesTheCallerWithItsOwnType_AtEveryHookPoint()
    {
        SaCase c = SaCase.Of(SaRun.GenomeLike(3000, 7));
        foreach (string point in new[] { "ph1", "naming", "ph23", "partition", "leaf" })
        {
            foreach (int threads in new[] { 1, 3 })
            {
                using var dir = new TempDir();
                SaBuildOptions o = SaTiny.Build(dir.Sub("scratch"), 3, threads) with { Trace = step => { if (step == point) throw new FormatException(point); } };
                string textPath = SaRun.WriteText(dir, c);
                string outDir = dir.Sub("out");
                Directory.CreateDirectory(outDir);
                using MappedTextReader reader = MappedTextReader.Open(textPath, c.N, c.Bits);
                FormatException e = Assert.Throws<FormatException>(() => SaBuilder.Build(reader, Path.Combine(outDir, "o.sa"), o), $"{point}, {threads} threads")!;
                Assert.That(e.Message, Is.EqualTo(point));
                Assert.That(Directory.GetFileSystemEntries(outDir), Is.Empty, point);
                Assert.That(Directory.Exists(o.ScratchDir), Is.False, point);
                ProcFds.AssertClosed(point, [Path.Combine(outDir, "o.sa")], [o.ScratchDir]);
            }
        }
    }

    [Test]
    public void ADisposedReader_IsAnObjectDisposedException_AndNothingIsLeft()
    {
        SaCase c = SaCase.Of(StringFamilies.RandomDna(500, 3));
        string textPath = SaRun.WriteText(_dir, c);
        string outDir = _dir.Sub("out");
        Directory.CreateDirectory(outDir);
        MappedTextReader reader = MappedTextReader.Open(textPath, c.N, c.Bits);
        reader.Dispose();
        Assert.Throws<ObjectDisposedException>(() => SaBuilder.Build(reader, Path.Combine(outDir, "o.sa"), Tiny(1)));
        Assert.That(Directory.GetFileSystemEntries(outDir), Is.Empty);
        Assert.That(Directory.Exists(Scratch), Is.False);
    }

    [Test]
    public void AShortOrDamagedText_IsInvalidData_AtOpen_BeforeAnythingIsCreated()
    {
        SaCase c = SaCase.Of(StringFamilies.RandomDna(500, 3));
        string textPath = SaRun.WriteText(_dir, c);
        byte[] bytes = File.ReadAllBytes(textPath);
        string shortPath = _dir.File("short.pk");
        File.WriteAllBytes(shortPath, bytes[..^1]);
        Assert.Throws<InvalidDataException>(() => MappedTextReader.Open(shortPath, c.N, c.Bits));
        File.WriteAllBytes(shortPath, []);
        Assert.Throws<InvalidDataException>(() => MappedTextReader.Open(shortPath, c.N, c.Bits));
        Assert.That(Directory.Exists(Scratch), Is.False);
    }

    [Test]
    public void ATextAboveTheDesignLimit_IsRefused_BeforeAnythingIsCreated()
    {
        long n = PackedTextFormat.MaxSymbols + 1;
        string path = _dir.File("huge.pk");
        using (var f = new FileStream(path, FileMode.Create)) f.SetLength(PackedTextFormat.TotalBytes(n, 1));
        string outDir = _dir.Sub("out");
        Directory.CreateDirectory(outDir);
        using MappedTextReader reader = MappedTextReader.Open(path, n, 1);
        Assert.Throws<NotSupportedException>(() => SaBuilder.Build(reader, Path.Combine(outDir, "o.sa"), Tiny(1)));
        Assert.That(Directory.GetFileSystemEntries(outDir), Is.Empty);
        Assert.That(Directory.Exists(Scratch), Is.False);
    }

    [Test]
    public void TheOptionsAndTheArguments_AreChecked()
    {
        SaCase c = SaCase.Of("banana");
        string textPath = SaRun.WriteText(_dir, c);
        string sa = _dir.File("o.sa");
        using MappedTextReader reader = MappedTextReader.Open(textPath, c.N, c.Bits);
        SaBuildOptions ok = Tiny(1);
        Assert.Throws<ArgumentNullException>(() => SaBuilder.Build(null!, sa, ok));
        Assert.Throws<ArgumentNullException>(() => SaBuilder.Build(reader, null!, ok));
        Assert.Throws<ArgumentNullException>(() => SaBuilder.Build(reader, sa, null!));
        foreach ((string what, SaBuildOptions bad) in new (string, SaBuildOptions)[]
        {
            ("threads 0", ok with { Threads = 0 }),
            ("period 1", ok with { Period = 1 }),
            ("arena cap -1", ok with { ArenaCapBytes = -1 }),
            ("scale -1", ok with { ScaleN = -1 }),
            ("rule -1", ok with { RuleN = -1 }),
            ("bucket capacity 1", ok with { Tuning = ok.Tuning with { BucketCapacity = 1 } }),
            ("window 0", ok with { Tuning = ok.Tuning with { WindowBytes = 0 } }),
            ("cap of the verdict -1", ok with { RssAnonCapMb = -1 }),
        })
        {
            Assert.Throws<ArgumentOutOfRangeException>(() => SaBuilder.Build(reader, sa, bad), what);
            Assert.That(File.Exists(sa), Is.False, what);
            Assert.That(Directory.Exists(Scratch), Is.False, what);
        }

        // a reader of 32-bit symbols (a reduced string of the recursion) is not a packed text of the file codec
        string widePath = _dir.File("wide.bin");
        File.WriteAllBytes(widePath, new byte[(4 * 3) + 64]);
        using MappedTextReader wide = MappedTextReader.OpenReduced(widePath, 3);
        Assert.Throws<ArgumentException>(() => SaBuilder.Build(wide, sa, ok), "a reader of 32-bit symbols");
        Assert.That(File.Exists(sa), Is.False);
        Assert.That(Directory.Exists(Scratch), Is.False);
    }

    // ---- the verdict of the RssAnon gate and the platforms

    [Test]
    public void WithoutProc_TheBuildStillWorks_AndSaysTheRssFigureIsUnavailable()
    {
        SaCase c = SaCase.Of(SaRun.GenomeLike(3000, 3));
        SaRun.Outcome r = SaRun.Run(_dir, c, Tiny(2) with { ProcSource = new FakeProcSource() });
        SaRun.Compare(c, r, "no /proc");
        Assert.That(r.Result.RssAnonPeakMb, Is.Null);
        Assert.That(r.Result.Ph1RssAnonPeakMb, Is.Null);
        Assert.That(r.Result.Ph23RssAnonPeakMb, Is.Null);
        Assert.That(r.Result.RssVerdict, Is.EqualTo(RssVerdict.Unavailable));
    }

    [Test]
    public void AHostileProcSource_IsUnavailable_NotAnException()
    {
        SaCase c = SaCase.Of(StringFamilies.RandomDna(1500, 3));
        foreach (Exception e in new Exception[] { new IOException("x"), new UnauthorizedAccessException("y"), new InvalidOperationException("z") })
        {
            using var dir = new TempDir();
            SaRun.Outcome r = SaRun.Run(dir, c, SaTiny.Build(dir.Sub("scratch"), 2) with { ProcSource = new FakeProcSource { Throws = e } });
            Assert.That(r.Result.RssVerdict, Is.EqualTo(RssVerdict.Unavailable), e.GetType().Name);
        }
    }

    [Test]
    public void TheRssVerdict_IsWithin_OrExceeded_AgainstTheCap_WhenTheFigureIsThere()
    {
        if (!OperatingSystem.IsLinux() || !File.Exists("/proc/self/status")) Assert.Ignore("needs /proc");
        SaCase c = SaCase.Of(StringFamilies.RandomDna(3000, 3));
        SaRun.Outcome within = SaRun.Run(_dir, c, Tiny(1) with { RssAnonCapMb = 1e6 });
        Assert.That(within.Result.RssVerdict, Is.EqualTo(RssVerdict.Within));
        Assert.That(within.Result.RssAnonPeakMb, Is.GreaterThan(1.0));
        Assert.That(within.Result.Ph1RssAnonPeakMb, Is.GreaterThan(1.0));
        Assert.That(within.Result.Ph23RssAnonPeakMb, Is.GreaterThan(1.0));
        Assert.That(within.Result.RssAnonPeakMb, Is.GreaterThanOrEqualTo(Math.Max(within.Result.Ph1RssAnonPeakMb!.Value, within.Result.Ph23RssAnonPeakMb!.Value) - 10), "the whole run is at least as high as a phase, up to the sampler's grain");
        using var dir2 = new TempDir();
        SaRun.Outcome exceeded = SaRun.Run(dir2, c, SaTiny.Build(dir2.Sub("scratch"), 1) with { RssAnonCapMb = 0.001 });
        Assert.That(exceeded.Result.RssVerdict, Is.EqualTo(RssVerdict.Exceeded), "a verdict, not an exception: the build finished");
        SaRun.Compare(c, exceeded, "exceeded");
    }

    [Test]
    public void WithRssSamplingOff_TheVerdictIsUnavailable()
    {
        SaCase c = SaCase.Of(StringFamilies.RandomDna(1500, 3));
        SaRun.Outcome r = SaRun.Run(_dir, c, Tiny(2) with { MeasureRss = false });
        Assert.That(r.Result.RssVerdict, Is.EqualTo(RssVerdict.Unavailable));
        Assert.That(r.Result.RssAnonPeakMb, Is.Null);
    }

    [Test]
    public void OnTheRealLinux_TheReservationOfTheSuffixArrayWasMade()
    {
        if (!OperatingSystem.IsLinux()) Assert.Ignore("needs Linux");
        SaCase c = SaCase.Of(StringFamilies.RandomDna(600, 3));
        SaRun.Outcome r = SaRun.Run(_dir, c, Tiny(1));
        Assert.That(r.Result.Reservation, Is.EqualTo(ReservationOutcome.Preallocated));
    }

    [Test]
    public void OnASimulatedPlatformWithoutFallocateOrProc_TheStageStillRuns_AndSaysWhatItDidNotCheck()
    {
        SaCase c = SaCase.Of(StringFamilies.RandomDna(600, 3));
        foreach ((long free, ReservationOutcome expected) in new[] { (-1L, ReservationOutcome.Unchecked), (1L << 40, ReservationOutcome.FreeSpaceChecked) })
        {
            using var dir = new TempDir();
            var reservation = new FakeReservationSource { CanPreallocate = false, Free = free };
            SaBuildOptions o = SaTiny.Build(dir.Sub("scratch"), 2, 2) with { Files = new MappedFileOptions { Source = reservation }, ProcSource = new FakeProcSource() };
            SaRun.Outcome r = SaRun.Run(dir, c, o);
            SaRun.Compare(c, r, $"free {free}");
            Assert.That(r.Result.Reservation, Is.EqualTo(expected));
            Assert.That(reservation.PreallocateCalls, Is.Zero);
            Assert.That(reservation.FreeCalls, Is.GreaterThanOrEqualTo(2), "the free space was read for the suffix array and the rank file");
            Assert.That(r.Result.RssVerdict, Is.EqualTo(RssVerdict.Unavailable));
        }
    }

    [Test]
    public void OnASimulatedPlatform_AFullDiskIsStillAnIoException_FromTheFreeSpaceCheck()
    {
        SaCase c = SaCase.Of(StringFamilies.RandomDna(600, 3));
        var reservation = new FakeReservationSource { CanPreallocate = false, Free = 100 };
        SaRun.AssertFailsCleanly(_dir, c, Tiny(2) with { Files = new MappedFileOptions { Source = reservation } }, typeof(IOException), "simulated full disk");
    }

    [Test]
    public void TheScalarCommonPrefix_GivesTheSameFile_AsTheVectorOne()
    {
        string text = SaRun.GenomeLike(7000, 5) + StringFamilies.RandomDna(900, 12) + StringFamilies.RandomDna(77, 6) + StringFamilies.RandomDna(900, 12);
        SaCase c = SaCase.Of(text);
        string? reference = null;
        foreach (bool vectors in new[] { true, false })
        {
            using var dir = new TempDir();
            SaRun.Outcome r = SaRun.Run(dir, c, SaTiny.Build(dir.Sub("scratch"), 3, 2), configure: reader => reader.UseVectors = vectors);
            SaRun.Compare(c, r, $"vectors {vectors}");
            string h = SaRun.Sha(r.SaPath);
            reference ??= h;
            Assert.That(h, Is.EqualTo(reference), $"vectors {vectors}");
        }
    }

    [Test]
    public void EveryParameter_ReachesItsStage()
    {
        SaCase c = SaCase.Of(new string('a', 3000));
        var t = new SaBuildTuning
        {
            DcsBucketCapacity = 40,
            DcsBaseCase = 50,
            DcsChunkSlots = 100,
            DcsSortCapacity = 30,
            DcsRunSpillBytes = 0,
            DcsFanOut = 3,
            DcsReservoir = 12,
            DcsIoBytes = 64,
            DcsMergeBudget = 100,
            BucketCapacity = 70,
            FanOut = 5,
            Oversample = 3,
            SpillBytes = 32,
            WindowBytes = 96,
        };
        var counters = new IoCounters();
        SaBuildOptions o = new SaBuildOptions { ScratchDir = Scratch, Threads = 1, Period = 7, ExactPeriod = true, Tuning = t, Counters = counters };
        SaRun.Outcome r = SaRun.Run(_dir, c, o);
        SaRun.Compare(c, r, "plumbing");
        SaBuildResult res = r.Result;
        Assert.That(res.Parameters, Is.EqualTo(SaBuildRules.Resolve(c.N, o)));
        Assert.That(res.Sort.MaxLeaf, Is.LessThanOrEqualTo(70), "K reaches PH3");
        Assert.That(res.Ph23ArenaPeakBytes, Is.EqualTo(24L * 70), "K sizes the leaf buffers");
        Assert.That(res.Dcs.Levels[0].FileRuns, Is.GreaterThan(0), "the run spill size reaches PH1");
        Assert.That(res.Dcs.MaxLevel, Is.GreaterThanOrEqualTo(2), "N_base reaches PH1");
        Assert.That(res.Io.MappedMaxWindowBytes, Is.EqualTo(96), "the window reaches both writers");
        Assert.That(res.Dcs.Levels[0].SplitScans, Is.GreaterThan(0), "K_s reaches PH1");
    }

    [Test]
    public void TheOutputFile_IsCreatedBeforePh1_SoAFullDiskIsFoundFirst()
    {
        SaCase c = SaCase.Of(StringFamilies.RandomDna(800, 3));
        var order = new List<string>();
        string outDir = _dir.Sub("out");
        Directory.CreateDirectory(outDir);
        string sa = Path.Combine(outDir, "o.sa");
        SaBuildOptions o = Tiny(2) with { Trace = e => { lock (order) order.Add(e + ":" + (File.Exists(sa) ? "sa" : "-")); } };
        string textPath = SaRun.WriteText(_dir, c);
        using MappedTextReader reader = MappedTextReader.Open(textPath, c.N, c.Bits);
        SaBuilder.Build(reader, sa, o);
        Assert.That(order.First(e => e.StartsWith("ph1:", StringComparison.Ordinal)), Is.EqualTo("ph1:sa"));
    }
}

/// <summary>The scale-down rule of ADR-004 (E3, E3b, E6) for PH1 and PH2 / PH3.</summary>
[TestFixture]
public class SaBuildRulesTests
{
    private static SaBuildOptions Options(long scale = 0, long rule = 0, SaBuildTuning? tuning = null, int threads = 4) => new() { ScratchDir = "s", ScaleN = scale, RuleN = rule, Tuning = tuning ?? new SaBuildTuning(), Threads = threads };

    [Test]
    public void Unscaled_EveryStageHasItsReferenceValue()
    {
        SaBuildParameters p = SaBuildRules.Resolve(12_345_678, Options());
        Assert.That(p.N, Is.EqualTo(12_345_678));
        Assert.That(p.Scaled, Is.False);
        Assert.That(p.Factor, Is.EqualTo(1.0));
        Assert.That((p.DcsBucketCapacity, p.DcsBaseCase, p.DcsChunkSlots, p.DcsSortCapacity), Is.EqualTo((1048576L, 16_777_216L, 4194304L, 1048576L)));
        Assert.That((p.DcsRunSpillBytes, p.DcsFanOut, p.DcsReservoir, p.DcsIoBytes, p.DcsMergeBudget), Is.EqualTo((65536, 1024, 16384, 65536, 16L << 20)));
        Assert.That((p.BucketCapacity, p.FanOut, p.Oversample, p.SpillBytes, p.WindowBytes), Is.EqualTo((1048576L, 64, 16, 65536, 1 << 20)));
    }

    [Test]
    public void Scaled_At100Mbp_IsTheGeometryOfThreeTimesTenToTheNinth()
    {
        SaBuildParameters p = SaBuildRules.Resolve(100_000_000, Options(scale: 3_000_000_000));
        Assert.That(p.Scaled, Is.True);
        Assert.That(p.Factor, Is.EqualTo(1e8 / 3e9).Within(1e-12));
        Assert.That((p.DcsBucketCapacity, p.DcsBaseCase), Is.EqualTo((34953L, 559241L)), "the E6 parameters row");
        Assert.That((p.DcsChunkSlots, p.DcsSortCapacity, p.DcsRunSpillBytes), Is.EqualTo((139810L, 34953L, 2185)));
        Assert.That(p.BucketCapacity, Is.EqualTo(33_333), "K = n / 3000");
        Assert.That((p.FanOut, p.Oversample, p.SpillBytes, p.WindowBytes), Is.EqualTo((64, 16, 65536, 1 << 20)), "the buffers keep their real sizes");
    }

    [Test]
    public void TheRuleN_FixesTheGeometry_WhateverTheTextLength()
    {
        SaBuildParameters at1e8 = SaBuildRules.Resolve(100_000_000, Options(scale: 3_000_000_000));
        SaBuildParameters at5e8 = SaBuildRules.Resolve(500_000_000, Options(scale: 3_000_000_000, rule: 100_000_000));
        Assert.That(at5e8 with { N = at1e8.N }, Is.EqualTo(at1e8));
        Assert.That(at5e8.N, Is.EqualTo(500_000_000));
        Assert.That(SaBuildRules.Resolve(500_000_000, Options(scale: 0, rule: 100_000_000)).BucketCapacity, Is.EqualTo(1L << 20), "the rule is ignored when unscaled");
    }

    [Test]
    public void TheFloors_KeepATinyTextUsable()
    {
        SaBuildParameters p = SaBuildRules.Resolve(10, Options(scale: 3_000_000_000));
        Assert.That(p.DcsBucketCapacity, Is.EqualTo(3));
        Assert.That(p.DcsBaseCase, Is.EqualTo(4));
        Assert.That(p.DcsChunkSlots, Is.EqualTo(16));
        Assert.That(p.DcsSortCapacity, Is.EqualTo(3));
        Assert.That(p.DcsRunSpillBytes, Is.Zero);
        Assert.That(p.BucketCapacity, Is.EqualTo(2));
    }

    [Test]
    public void AnExplicitValue_WinsOverTheRule()
    {
        var t = new SaBuildTuning { DcsBucketCapacity = 11, DcsBaseCase = 12, DcsChunkSlots = 13, DcsSortCapacity = 14, DcsRunSpillBytes = 15, DcsFanOut = 16, DcsReservoir = 17, DcsIoBytes = 18, DcsMergeBudget = 19, BucketCapacity = 20, FanOut = 21, Oversample = 22, SpillBytes = 23, WindowBytes = 24 };
        foreach (long scale in new long[] { 0, 3_000_000_000 })
        {
            SaBuildParameters p = SaBuildRules.Resolve(100_000, Options(scale, tuning: t));
            Assert.That(p, Is.EqualTo(new SaBuildParameters { N = 100_000, Factor = p.Factor, Scaled = scale > 0, DcsBucketCapacity = 11, DcsBaseCase = 12, DcsChunkSlots = 13, DcsSortCapacity = 14, DcsRunSpillBytes = 15, DcsFanOut = 16, DcsReservoir = 17, DcsIoBytes = 18, DcsMergeBudget = 19, BucketCapacity = 20, FanOut = 21, Oversample = 22, SpillBytes = 23, WindowBytes = 24 }));
        }
    }

    [Test]
    public void BadValues_AreRefused_WithTheirName()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => SaBuildRules.Resolve(-1, Options()));
        Assert.Throws<ArgumentOutOfRangeException>(() => SaBuildRules.Resolve(10, Options(scale: -1)));
        Assert.Throws<ArgumentOutOfRangeException>(() => SaBuildRules.Resolve(10, Options(rule: -1)));
        Assert.Throws<ArgumentOutOfRangeException>(() => SaBuildRules.Resolve(10, Options(threads: 0)));
        Assert.Throws<ArgumentNullException>(() => SaBuildRules.Resolve(10, null!));
        foreach ((string what, SaBuildTuning t) in new (string, SaBuildTuning)[]
        {
            ("DcsBucketCapacity 1", new() { DcsBucketCapacity = 1 }),
            ("DcsBucketCapacity too large", new() { DcsBucketCapacity = (int.MaxValue / 2) + 1L }),
            ("DcsBaseCase 0", new() { DcsBaseCase = 0 }),
            ("DcsChunkSlots 0", new() { DcsChunkSlots = 0 }),
            ("DcsSortCapacity 1", new() { DcsSortCapacity = 1 }),
            ("DcsRunSpillBytes -1", new() { DcsRunSpillBytes = -1 }),
            ("DcsFanOut 1", new() { DcsFanOut = 1 }),
            ("DcsReservoir 2", new() { DcsReservoir = 2 }),
            ("DcsIoBytes 15", new() { DcsIoBytes = 15 }),
            ("DcsMergeBudget 7", new() { DcsMergeBudget = 7 }),
            ("BucketCapacity 1", new() { BucketCapacity = 1 }),
            ("BucketCapacity too large", new() { BucketCapacity = (int.MaxValue / 2) + 1L }),
            ("FanOut 1", new() { FanOut = 1 }),
            ("Oversample 0", new() { Oversample = 0 }),
            ("SpillBytes 0", new() { SpillBytes = 0 }),
            ("WindowBytes 0", new() { WindowBytes = 0 }),
        })
        {
            Assert.Throws<ArgumentOutOfRangeException>(() => SaBuildRules.Resolve(10, Options(tuning: t)), what);
        }
    }

    [Test]
    public void TheReferenceConstants_AreThoseOfTheAdr()
    {
        Assert.That(SaBuildRules.ReferenceN, Is.EqualTo(3_000_000_000));
        Assert.That(SaBuildRules.BucketReference, Is.EqualTo(1_000_000));
        Assert.That(SaBuildRules.BucketUnscaled, Is.EqualTo(1 << 20));
        Assert.That(SaBuildRules.DcsBucketReference, Is.EqualTo(1 << 20));
        Assert.That(SaBuildRules.DcsBaseReference, Is.EqualTo(16 << 20), "N_base = 16 Mi symbols, the owner's decision of 2026-11-18");
        Assert.That(SaBuildRules.DcsChunkReference, Is.EqualTo(4 << 20));
        Assert.That(SaBuildRules.DcsSortReference, Is.EqualTo(1 << 20));
        Assert.That(SaBuildRules.RunSpillReference, Is.EqualTo(64 * 1024));
    }
}

/// <summary>The target of the DIST primitive can be created through the reservation of <see cref="MappedOutputFile"/> (a change of work package 1b.7d to 1b.7a: PH1's rank file is one of its targets).</summary>
[TestFixture]
public class DistTargetCreationTests
{
    private TempDir _dir = null!;

    [SetUp]
    public void SetUp() => _dir = new TempDir();

    [TearDown]
    public void TearDown() => _dir.Dispose();

    private (string Sa, string Isa) Inputs(int n)
    {
        SaCase c = SaCase.Of(StringFamilies.RandomDna(n, 3));
        string sa = _dir.File("sa.bin");
        FileBits.WriteUInts(sa, c.Sa);
        return (sa, _dir.File("isa.bin"));
    }

    [Test]
    public void WithoutFiles_TheTargetIsCreatedSparse_AsBefore()
    {
        (string sa, string isa) = Inputs(500);
        Dist.BuildIsa(sa, 500, isa, new DistOptions { ScratchDir = _dir.Sub("s"), ChunkEntries = 100 });
        SaCase c = SaCase.Of(StringFamilies.RandomDna(500, 3));
        uint[] got = FileBits.ReadUInts(isa);
        for (int i = 0; i <= 500; i++) Assert.That(got[c.Sa[i]], Is.EqualTo((uint)i));
    }

    [Test]
    public void WithFiles_TheTargetIsCreatedThroughTheReservation_AndHasTheSameBytes()
    {
        (string sa, string isa) = Inputs(500);
        var reservation = new FakeReservationSource { CanPreallocate = true };
        var counters = new IoCounters();
        Dist.BuildIsa(sa, 500, isa, new DistOptions { ScratchDir = _dir.Sub("s"), ChunkEntries = 100, Files = new MappedFileOptions { Source = reservation }, Counters = counters, WindowBytes = 256 });
        Assert.That(counters.Snapshot().MappedMaxWindowBytes, Is.LessThanOrEqualTo(256), "the target written through the reservation uses the window of the options");
        Assert.That(reservation.PreallocateCalls, Is.EqualTo(1));
        string plain = _dir.File("plain.bin");
        Dist.BuildIsa(sa, 500, plain, new DistOptions { ScratchDir = _dir.Sub("s2"), ChunkEntries = 100 });
        Assert.That(File.ReadAllBytes(isa), Is.EqualTo(File.ReadAllBytes(plain)));
        Assert.That(counters.MappedWriteBytes, Is.EqualTo(4 * 501));
    }

    [Test]
    public void AFullDisk_AtTheCreationOfTheTarget_IsAnIoException_AndNoTargetIsLeft()
    {
        (string sa, string isa) = Inputs(500);
        Assert.Throws<IOException>(() => Dist.BuildIsa(sa, 500, isa, new DistOptions { ScratchDir = _dir.Sub("s"), ChunkEntries = 100, Files = new MappedFileOptions { Source = new FailingAfter(0) } }));
        Assert.That(File.Exists(isa), Is.False);
        Assert.That(Directory.Exists(_dir.Sub("s")) ? Directory.GetFileSystemEntries(_dir.Sub("s")) : [], Is.Empty);
        ProcFds.AssertClosed("full disk", [isa, sa], [_dir.Sub("s")]);
    }

    [Test]
    public void TheMappingOfTheTarget_IsReleasedWhenTheScatterIsDone()
    {
        (string sa, string isa) = Inputs(500);
        Dist.BuildIsa(sa, 500, isa, new DistOptions { ScratchDir = _dir.Sub("s"), ChunkEntries = 100, Files = new MappedFileOptions() });
        Assert.That(ProcMaps.Contains(isa) ?? false, Is.False);
        ProcFds.AssertClosed("scatter", [isa, sa], [_dir.Sub("s")]);
    }
}
