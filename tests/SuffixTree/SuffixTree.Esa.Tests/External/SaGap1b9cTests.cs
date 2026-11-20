using SuffixTree.Esa.External;
using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.Esa.Tests.External;

/// <summary>
/// The gaps that the complete mutation audit of work package 1b.9c showed in the code of 1b.7d (PH1-PH3, <c>External/Dcs</c>, <c>External/Sort</c>, <c>External/Sa</c>): the mutants of the lists that the audit of 1b.7d did not run
/// and the "second line of defence" mutants that it only argued. Each test names the mutants it kills (ids of <c>analysis/mutation/1b7d_*_listed.json</c>); every one of them was run against the test before it counted.
/// </summary>
[TestFixture]
public class SaGap1b9cTests
{
    private TempDir _dir = null!;

    [SetUp]
    public void SetUp() => _dir = new TempDir();

    [TearDown]
    public void TearDown() => _dir.Dispose();

    /// <summary>A sink that fails the test if it is ever asked for a writer: the entry of the stage must refuse first.</summary>
    private sealed class NoSink : ISaSliceSink
    {
        public ISaSliceWriter CreateWriter() => throw new InvalidOperationException("the stage went on after the entry checks");
    }

    // ---------------------------------------------------------------- B04: a write of a run that fails after its file was opened

    [Test]
    public void AWriteOfARun_ThatFailsWithADiskFull_GivesTheBytesBack_AndLeavesNoFile()
    {
        // B04 (RunStore: the failed write leaves its file): the path of the run is a link to /dev/full, which opens and then refuses every write with "no space left on device", so the file of the run exists when the write fails
        // and the failure path must remove it (File.Delete removes the link, never the device)
        if (!OperatingSystem.IsLinux() || !File.Exists("/dev/full")) Assert.Ignore("there is no /dev/full");
        string dir = _dir.Sub("runs");
        Directory.CreateDirectory(dir);
        var counters = new IoCounters();
        var ledger = new ScratchLedger(counters);
        var store = new RunStore(dir, "t", 0, ledger, counters);
        string path = Path.Combine(dir, "t-run-00000001.bin");
        File.CreateSymbolicLink(path, "/dev/full");
        Assert.That(Directory.GetFileSystemEntries(dir), Is.EqualTo(new[] { path }), "the file of the run exists before the write");
        Assert.Throws(Is.InstanceOf<IOException>(), () => store.Create([(1UL << 32) | 1, (2UL << 32) | 2]));
        Assert.That(Directory.GetFileSystemEntries(dir), Is.Empty, "the file of the failed run is removed");
        Assert.That(ledger.Live, Is.Zero, "the bytes are given back");
        Assert.That(counters.ScratchLiveBytes, Is.Zero);
        Assert.That(File.Exists("/dev/full"), Is.True, "the device is still there");
    }

    // ---------------------------------------------------------------- H64, H66: a failure in the middle of a partition and in the hand-over of its children

    private void SortWithFailureAt(string point, int occurrence, ScratchLedger ledger, string scratch)
    {
        SaCase c = SaCase.Of(StringFamilies.RandomDna(40_000, 3));
        DifferenceCover cover = DifferenceCover.Create(7);
        string textPath = SaRun.WriteText(_dir, c);
        using MappedTextReader text = MappedTextReader.Open(textPath, c.N, c.Bits);
        using RankTable ranks = RankTable.Open(SaInputs.WriteRankFile(_dir, c, cover), cover, c.N);
        int seen = 0;
        SortOptions options = SaTiny.Sort(scratch, 1, 1) with
        {
            Ledger = ledger,
            Trace = step =>
            {
                if (step == point && ++seen == occurrence) throw new InjectedFailureException($"{point} #{occurrence}");
            },
        };
        Assert.Throws<InjectedFailureException>(() => BucketSorter.Run(text, options, new MemorySaSink(c.N + 1), ranks), $"{point} #{occurrence}");
        Assert.That(seen, Is.EqualTo(occurrence), "the hook point was reached that often");
    }

    [Test]
    public void AFailureInTheMiddleOfTheScanOfAPartition_ClosesEveryWriterThatHoldsAFile()
    {
        // H64 (the writers of a partition are not disposed): the scan of the root bucket has three chunks; the failure is at the second, when every writer holds a spill file; the descriptors must be closed when the exception arrives
        string scratch = _dir.Sub("scratch");
        var ledger = new ScratchLedger(new IoCounters());
        SortWithFailureAt("scan", 2, ledger, scratch);
        ProcFds.AssertClosed("after the failure in the scan", [], [scratch]);
        Assert.That(ledger.Live, Is.Zero, "the bytes of the position files are given back");
        Assert.That(Directory.Exists(scratch), Is.False, "nothing is left");
    }

    [Test]
    public void AFailureInTheHandOverOfTheChildrenOfAPartition_DisposesTheChildrenThatWereCompleted()
    {
        // H66 (the children of a failed hand-over are not disposed): the failure is before the third child is completed; the first two are position files that the ledger holds, which their disposal gives back
        string scratch = _dir.Sub("scratch");
        var ledger = new ScratchLedger(new IoCounters());
        SortWithFailureAt("handover", 3, ledger, scratch);
        Assert.That(ledger.Live, Is.Zero, "the first two children were disposed, the third and the writers were disposed uncompleted");
        ProcFds.AssertClosed("after the failure in the hand-over", [], [scratch]);
        Assert.That(Directory.Exists(scratch), Is.False, "nothing is left");
    }

    // ---------------------------------------------------------------- L29, A08: a failure after the rank file was written removes it

    [Test]
    public void AFailureAfterTheRankFileOfLevelZeroWasWritten_RemovesTheFile_InEveryWayThatLevelZeroEnds()
    {
        // L29 (the inversion of a level above the base does not mark the rank file) and A08 (the base case does not): a rank file at a path outside the scratch directory belongs to the caller from the moment it was asked for; the hook
        // after the last level to write one (level 0) fails, and the file must not be left
        (string Label, string Text, int BaseCase)[] cases =
        [
            ("unique names, written by the merge", StringFamilies.RandomDna(900, 4), 40),
            ("repeated names, base case", new string('a', 300), 1_000),
            ("repeated names, the inversion of the level below", StringFamilies.Periodic(1_200, 37, 5), 5),
        ];
        foreach ((string label, string text, int baseCase) in cases)
        {
            SaCase c = SaCase.Of(text);
            DifferenceCover cover = DifferenceCover.Create(7);
            string textPath = SaRun.WriteText(_dir, c, "-" + label[..5]);
            using MappedTextReader reader = MappedTextReader.Open(textPath, c.N, c.Bits);
            string rank = _dir.File("rank-" + label[..5] + ".bin");
            string scratch = _dir.Sub("scratch-" + label[..5]);
            ExtDcsOptions plain = SaTiny.Dcs(scratch, 1, 1) with { BaseCase = baseCase };
            int events = 0;
            ExtDcsStats stats = ExternalDcs.BuildRankFile(reader, cover, rank, plain with { Trace = step => events += step == "rank-written" ? 1 : 0 });
            ExtDcsLevelStats top = stats.Levels[0];
            Assert.That((top.NamesUnique, top.BaseCase, stats.Levels.Count > 1), Is.EqualTo(label[0] == 'u' ? (true, false, false) : label.Contains("base", StringComparison.Ordinal) ? (false, true, false) : (false, false, true)), $"{label}: the case takes the way it is named for");
            Assert.That(events, Is.EqualTo(stats.Levels.Count), $"{label}: every level writes its rank file");
            File.Delete(rank);
            int seen = 0;
            Assert.Throws<InjectedFailureException>(() => ExternalDcs.BuildRankFile(reader, cover, rank, plain with { Trace = step => { if (step == "rank-written" && ++seen == events) throw new InjectedFailureException(label); } }), label);
            Assert.That(File.Exists(rank), Is.False, $"{label}: the rank file of the failed build is left");
            Assert.That(Directory.Exists(scratch), Is.False, $"{label}: the scratch directory is left");
        }
    }

    // ---------------------------------------------------------------- L15: the position files of the sort of a level carry the tag of the build

    [Test]
    public void ThePositionFilesOfTheSortOfALevel_CarryTheTagOfTheBuild_AtEveryPointOfTheBuild()
    {
        // L15 (the sort files carry no tag): the nested sort of a level below the first is run with a file prefix made of the tag of the build, so that a failed build finds and removes exactly its own files
        SaCase c = SaCase.Of(StringFamilies.Periodic(12_000, 37, 5));
        DifferenceCover cover = DifferenceCover.Create(7);
        string scratch = _dir.Sub("scratch");
        var seen = new SortedSet<string>(StringComparer.Ordinal);
        ExtDcsOptions options = SaTiny.Dcs(scratch, 1, 1) with
        {
            Trace = _ =>
            {
                foreach (string f in Directory.GetFiles(scratch)) seen.Add(Path.GetFileName(f));
            },
        };
        string textPath = SaRun.WriteText(_dir, c);
        using MappedTextReader text = MappedTextReader.Open(textPath, c.N, c.Bits);
        ExtDcsStats stats = ExternalDcs.BuildRankFile(text, cover, _dir.File("r.bin"), options);
        Assert.That(stats.Levels.Count, Is.GreaterThan(1), "the text goes through a level below the first");
        string[] positionFiles = [.. seen.Where(f => f.Contains("-pos-", StringComparison.Ordinal))];
        Assert.That(positionFiles, Is.Not.Empty, "the sort of the level below spilled positions to files");
        Assert.That(positionFiles, Has.All.Matches<string>(f => System.Text.RegularExpressions.Regex.IsMatch(f, @"^e3b-[0-9a-f]{8}-L[1-9][0-9]*-pos-[0-9]{8}\.bin$")), string.Join(", ", positionFiles));
        Assert.That(seen.Select(f => f[..12]).Distinct().Count(), Is.EqualTo(1), "every file of the scratch directory carries the one tag of the build: " + string.Join(", ", seen.Take(5)));
        Assert.That(Directory.Exists(scratch), Is.False, "nothing is left");
    }

    // ---------------------------------------------------------------- L14: the floor of the spill of the sort of a level

    [Test]
    public void TheSortOfALevel_SpillsAtLeast64BytesWhateverTheIoBufferIs_SoTinyBuffersAndTheFloorGiveTheSameFiles()
    {
        // L14 (the sort spill has no floor): the sort of a level takes the buffer of the options but not below 64 bytes, so an I/O buffer of 16 bytes and one of 64 spill the same buckets to the same number of position files
        SaCase c = SaCase.Of(StringFamilies.Periodic(12_000, 37, 5));
        DifferenceCover cover = DifferenceCover.Create(7);
        string textPath = SaRun.WriteText(_dir, c);
        using MappedTextReader text = MappedTextReader.Open(textPath, c.N, c.Bits);
        int FilesWithIoBuffer(int bytes)
        {
            string scratch = _dir.Sub("scratch" + bytes);
            var seen = new HashSet<string>(StringComparer.Ordinal);
            ExtDcsOptions options = SaTiny.Dcs(scratch, 1, 1) with { IoBufferBytes = bytes, Trace = _ => seen.UnionWith(Directory.GetFiles(scratch).Select(f => Path.GetFileName(f)).Where(f => f.Contains("-pos-", StringComparison.Ordinal)).Select(f => f[13..])) };
            ExternalDcs.BuildRankFile(text, cover, _dir.File($"r{bytes}.bin"), options);
            return seen.Count;
        }
        int tiny = FilesWithIoBuffer(16), floor = FilesWithIoBuffer(64);
        Assert.That(floor, Is.GreaterThan(0), "the sort of the level below spills");
        Assert.That(tiny, Is.EqualTo(floor), "16 bytes are raised to the floor of 64");
        Assert.That(FilesWithIoBuffer(32), Is.EqualTo(floor));
    }

    // ---------------------------------------------------------------- Q64: the thread count reaches PH1 and PH3

    /// <summary>
    /// Holds the first thread that calls <see cref="Hit"/> until a second one calls it as well (a rendezvous of two): two tasks of a phase are running at the same time only if the phase has two threads. A phase with one
    /// thread times out, which is what the test looks for.
    /// </summary>
    private sealed class Rendezvous(TimeSpan wait, int skip)
    {
        private readonly ManualResetEventSlim _both = new(false);
        private int _arrived;
        private int _seen;

        internal bool TimedOut { get; private set; }

        /// <summary>Two tasks were in <see cref="Hit"/> at the same time.</summary>
        internal bool Met => _both.IsSet;

        internal void Hit()
        {
            if (_both.IsSet || Interlocked.Increment(ref _seen) <= skip) return; // the first task of a phase is the only one (the root): it cannot meet another
            if (Interlocked.Increment(ref _arrived) >= 2)
            {
                _both.Set();
                return;
            }
            if (!_both.Wait(wait)) TimedOut = true;
        }
    }

    private bool RunsConcurrently(int threads, TimeSpan wait, out bool ph1, out bool ph3)
    {
        SaCase c = SaCase.Of(StringFamilies.RandomDna(3000, 5));
        var inPh1 = new Rendezvous(wait, 1);
        var inPh3 = new Rendezvous(wait, 0);
        bool phase3 = false;
        SaBuildOptions options = SaTiny.Build(_dir.Sub("scratch" + threads), 1, threads) with
        {
            Trace = step =>
            {
                if (step == "ph23") phase3 = true;
                else if (step == "naming" && !phase3) inPh1.Hit();
                else if (step == "leaf" && phase3) inPh3.Hit();
            },
        };
        SaRun.Run(_dir, c, options, "-" + threads);
        ph1 = inPh1.Met && !inPh1.TimedOut;
        ph3 = inPh3.Met && !inPh3.TimedOut;
        return ph1 && ph3;
    }

    [Test]
    public void TheThreadCountOfTheBuild_ReachesPh1AndPh3_TwoTasksOfEachRunAtTheSameTime()
    {
        // Q64 (PH1 gets one thread): the build output does not depend on the thread count, the concurrency does; a task of PH1 (naming) and a task of PH3 (leaf) wait for a second task of their phase
        Assert.That(RunsConcurrently(4, TimeSpan.FromSeconds(20), out bool ph1, out bool ph3), Is.True, $"four threads: PH1 {ph1}, PH3 {ph3}");
        // the control: one thread cannot meet itself, the instrument sees that
        Assert.That(RunsConcurrently(1, TimeSpan.FromMilliseconds(300), out bool c1, out bool c3), Is.False);
        Assert.That((c1, c3), Is.EqualTo((false, false)), "one thread: neither phase has two tasks at the same time");
    }

    // ---------------------------------------------------------------- d_Q76: the guard option reaches PH1 and PH3 each on its own

    [Test]
    public void TheGuardOption_ReachesPh1_AndPh3_EachOnItsOwn()
    {
        // d_Q76 (1b.7d, PH1 gets no arena guards): the trace of the stage tells the two phases apart; every arena of a phase is guarded when the option is passed
        SaCase c = SaCase.Of(StringFamilies.RandomDna(600, 7));
        var at = new Dictionary<string, long>();
        SaBuildOptions options = SaTiny.Build(_dir.Sub("scratch"), 1, 1) with { Trace = step => at.TryAdd(step, NativeArena.GuardedCreated) };
        SaRun.Run(_dir, c, options);
        long end = NativeArena.GuardedCreated;
        Assert.That(at["ph23"] - at["ph1"], Is.GreaterThan(0), "PH1 creates guarded arenas");
        Assert.That(end - at["ph23"], Is.GreaterThan(0), "PH2 / PH3 creates guarded arenas");
    }

    // ---------------------------------------------------------------- M18: the merge window check, as a function

    [Test]
    public void TheMergeWindowCheck_RefusesFewerAndMoreRecordsThanSlots_AndAcceptsExactlyTheNumberOfSlots()
    {
        // M18 (a short count of records is accepted): the check as a function; a window that received fewer records than it has slots lost the name of a sample position, one that received more has two records for a slot
        SampleNamer.CheckWindow(0, 0, 0);
        SampleNamer.CheckWindow(100, 164, 64);
        foreach ((long start, long end, long got) in new[] { (100L, 164L, 63L), (100L, 164L, 65L), (0L, 1L, 0L), (0L, 1L, 2L), (7L, 7L, 1L) })
        {
            Assert.That(() => SampleNamer.CheckWindow(start, end, got), Throws.TypeOf<InvalidDataException>().With.Message.Contains($"merge window [{start}, {end}) received {got} records"), $"[{start}, {end}) with {got}");
        }
    }

    // ---------------------------------------------------------------- M52: a fan-out cap above the number of buckets needed changes nothing

    [Test]
    public void TheFanOutCapOfTheNaming_DoesNotChangeTheShapeWhenFewerBucketsAreNeeded()
    {
        // M52 (the fan-out is the cap only): a root of about 1,300 samples in buckets of 600 needs 5 children; a cap of 64 and a cap of 4,096 give the same tree, so the same leaves and scans
        SaCase c = SaCase.Of(StringFamilies.RandomDna(3_000, 31));
        DifferenceCover cover = DifferenceCover.Create(7);
        string textPath = SaRun.WriteText(_dir, c);
        using MappedTextReader reader = MappedTextReader.Open(textPath, c.N, c.Bits);
        var shapes = new List<(long Samples, long Leaves, long Scans, long SplitScans)>();
        foreach (int cap in new[] { 64, 4_096 })
        {
            ExtDcsStats stats = ExternalDcs.BuildRankFile(reader, cover, _dir.File($"r{cap}.bin"), SaTiny.Dcs(_dir.Sub("scratch" + cap), 1, 1) with { BucketCapacity = 600, FanOut = cap });
            ExtDcsLevelStats l0 = stats.Levels[0];
            shapes.Add((l0.SampleSize, l0.Leaves, l0.Scans, l0.SplitScans));
        }
        Assert.That(shapes[0].Samples, Is.GreaterThan(600).And.LessThan(3 * 600 / 2 * 2), "the root is split into a few buckets");
        Assert.That(shapes[1], Is.EqualTo(shapes[0]));
        Assert.That(shapes[0].Leaves, Is.InRange(3, 12), "a handful of leaves, not hundreds");
    }

    // ---------------------------------------------------------------- L24: every file that PH1 reserves goes through the file options

    [Test]
    public void TheFileOptions_ReachTheRankFileAndTheSuffixArrayFileOfTheLevelBelow()
    {
        // L24 (the inversion has no reservation options): the rank file of a level that is inverted from the sort of the level below, and the suffix array file of that sort, are created through the source of the options
        SaCase c = SaCase.Of(StringFamilies.Periodic(1_200, 37, 5));
        string textPath = SaRun.WriteText(_dir, c);
        using MappedTextReader reader = MappedTextReader.Open(textPath, c.N, c.Bits);
        var source = new RecordingSource();
        ExtDcsOptions options = SaTiny.Dcs(_dir.Sub("scratch"), 1, 1) with { BaseCase = 5, Files = new MappedFileOptions { Source = source } };
        ExtDcsStats stats = ExternalDcs.BuildRankFile(reader, DifferenceCover.Create(7), _dir.File("rank.bin"), options);
        Assert.That(stats.Levels.Count, Is.GreaterThan(1), "the text goes through the inversion");
        Assert.That(source.Created.Count(p => p.EndsWith("rank.bin", StringComparison.Ordinal)), Is.EqualTo(1), "the rank file of the caller");
        Assert.That(source.Created.Count(p => p.EndsWith(".sa", StringComparison.Ordinal)), Is.EqualTo(stats.Levels.Count - 1), "one suffix array file per level below the first");
        Assert.That(source.Created.Count(p => p.EndsWith(".rank", StringComparison.Ordinal)), Is.EqualTo(stats.Levels.Count - 1), "one rank file per level below the first (the file that its level above is inverted from)");
    }

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

    // ---------------------------------------------------------------- H37: the seconds of the leaves and of the partitions

    [Test]
    public void TheSortCountsTheSecondsOfItsLeavesAndOfItsPartitions()
    {
        // H37 (the leaf time is not summed): both are sums of timer ticks of tasks that sort or split thousands of positions
        SaCase c = SaCase.Of(StringFamilies.RandomDna(20_000, 8));
        DifferenceCover cover = DifferenceCover.Create(7);
        string textPath = SaRun.WriteText(_dir, c);
        using MappedTextReader reader = MappedTextReader.Open(textPath, c.N, c.Bits);
        using RankTable ranks = RankTable.Open(SaInputs.WriteRankFile(_dir, c, cover), cover, c.N);
        SortStats stats = BucketSorter.Run(reader, SaTiny.Sort(_dir.Sub("scratch"), 1, 1), new MemorySaSink(c.N + 1), ranks);
        Assert.That(stats.Leaves, Is.GreaterThan(100));
        Assert.That(stats.LeafSeconds, Is.GreaterThan(0).And.LessThan(120));
        Assert.That(stats.PartitionSeconds, Is.GreaterThan(0).And.LessThan(120));
    }

    // ---------------------------------------------------------------- H02: the upper bound of the leaf capacity, with the allocation turned into an exception

    [Test]
    public void ALeafCapacityAboveTheBound_IsRefusedAtTheEntry_AndTheBoundItselfIsNot()
    {
        // H02 (no upper bound on K): without the entry check the next step allocates a leaf sorter of 8 GB or more, which killed the test host of the audit; a small budget turns that step into an exception, so that the mutant
        // fails this test instead of the machine. The bound itself passes the entry and is stopped by the budget.
        SaCase c = SaCase.Of("banana");
        DifferenceCover cover = DifferenceCover.Create(7);
        string textPath = SaRun.WriteText(_dir, c);
        using MappedTextReader reader = MappedTextReader.Open(textPath, c.N, c.Bits);
        using RankTable ranks = RankTable.Open(SaInputs.WriteRankFile(_dir, c, cover), cover, c.N);
        SortOptions Options(long k) => SaTiny.Sort(_dir.Sub("scratch"), 1, 1) with { LeafCapacity = k, Budget = new ArenaBudget(1 << 20) };
        Assert.Throws<ArgumentOutOfRangeException>(() => BucketSorter.Run(reader, Options((int.MaxValue / 2) + 1L), new NoSink(), ranks));
        Assert.Throws<ArgumentOutOfRangeException>(() => BucketSorter.Run(reader, Options(int.MaxValue), new NoSink(), ranks));
        Assert.That(() => BucketSorter.Run(reader, Options(int.MaxValue / 2), new MemorySaSink(c.N + 1), ranks), Throws.TypeOf<ArenaBudgetExceededException>(), "the bound is allowed; the budget refuses its buffers");
    }

    // ---------------------------------------------------------------- M06: a sample that exactly fits the bucket is a leaf

    [Test]
    public void ASampleThatExactlyFitsTheBucketCapacity_IsOneLeaf_AndOneMoreRecordIsSplit()
    {
        // M06 (the root is a leaf one record early): the root of the naming is split when the sample is larger than the capacity of a bucket, not when it is as large
        DifferenceCover cover = DifferenceCover.Create(7);
        SaCase c = SaCase.Of(StringFamilies.RandomDna(100, 12));
        long samples = cover.SampleCount(c.N);
        string textPath = SaRun.WriteText(_dir, c);
        using MappedTextReader reader = MappedTextReader.Open(textPath, c.N, c.Bits);
        ExtDcsLevelStats Level0(long capacity)
        {
            ExtDcsStats stats = ExternalDcs.BuildRankFile(reader, cover, _dir.File("r-" + capacity + ".bin"), SaTiny.Dcs(_dir.Sub("scratch-" + capacity), 1, 1) with { BucketCapacity = capacity });
            return stats.Levels[0];
        }
        ExtDcsLevelStats fits = Level0(samples);
        Assert.That((fits.SampleSize, fits.Leaves, fits.SplitScans), Is.EqualTo((samples, 1L, 0L)), "capacity = sample size: one leaf, no split scan");
        ExtDcsLevelStats over = Level0(samples - 1);
        Assert.That(over.SplitScans, Is.GreaterThanOrEqualTo(1), "capacity = sample size - 1: the root is split");
        Assert.That(over.Leaves, Is.GreaterThan(1));
    }

    // ---------------------------------------------------------------- N04: the sink of the rank file releases its window

    [Test]
    public void ThePlainRankFileSink_ReleasesItsWindowAndItsFile_WhenTheBuildIsDone()
    {
        // N04 (RankFileSink: the writer is not disposed): the sample names are unique here, so the merge writes the rank file through the sink; its last window must not stay mapped
        SaCase c = SaCase.Of(StringFamilies.RandomDna(900, 4));
        string textPath = SaRun.WriteText(_dir, c);
        using MappedTextReader reader = MappedTextReader.Open(textPath, c.N, c.Bits);
        string rank = _dir.File("rank.bin");
        ExtDcsStats stats = ExternalDcs.BuildRankFile(reader, DifferenceCover.Create(7), rank, SaTiny.Dcs(_dir.Sub("scratch"), 1, 1));
        Assert.That(stats.Levels[0].NamesUnique, Is.True, "the case takes the way of the sink");
        Assert.That(ProcMaps.Contains(rank) ?? false, Is.False, "the rank file is still mapped");
        ProcFds.AssertClosed("after the rank file sink", [rank], [_dir.Sub("scratch")]);
    }

    // ---------------------------------------------------------------- Q82: the optional check of the rank file reads every slot

    [Test]
    public void WithVerifyRanks_ASlotOfTheRankFileThatWasDamagedAfterPh1_IsRefusedBeforePh3()
    {
        // Q82 (the check of the rank file is replaced by reading its length): the hook "verify" fires right before the check; it damages slot 0 (a rank of 0 is outside [1, S + |D|]) and the check must say so
        SaCase c = SaCase.Of(StringFamilies.RandomDna(500, 9));
        string scratch = _dir.Sub("scratch");
        string rankPath = Path.Combine(scratch, "ranks.bin");
        var steps = new List<string>();
        SaBuildOptions options = SaTiny.Build(scratch, 1, 1) with
        {
            VerifyRanks = true,
            Trace = step =>
            {
                steps.Add(step);
                if (step != "verify") return;
                using Microsoft.Win32.SafeHandles.SafeFileHandle handle = File.OpenHandle(rankPath, FileMode.Open, FileAccess.Write, FileShare.ReadWrite);
                RandomAccess.Write(handle, new byte[4], 0);
            },
        };
        string textPath = SaRun.WriteText(_dir, c);
        using MappedTextReader reader = MappedTextReader.Open(textPath, c.N, c.Bits);
        InvalidDataException? e = Assert.Throws<InvalidDataException>(() => SaBuilder.Build(reader, _dir.File("o.sa"), options));
        Assert.That(e!.Message, Does.Contain("slot 0 holds the rank 0"));
        Assert.That(steps, Does.Contain("verify").And.Not.Contain("ph23"), "the build stopped before PH2 / PH3");
        Assert.That(File.Exists(_dir.File("o.sa")) || Directory.Exists(scratch), Is.False, "nothing is left");
        // without the option there is no check, and no hook
        steps.Clear();
        SaBuildOptions plain = SaTiny.Build(scratch, 1, 1) with { Trace = steps.Add };
        SaBuilder.Build(reader, _dir.File("o2.sa"), plain);
        Assert.That(steps, Does.Not.Contain("verify"));
    }

    // ---------------------------------------------------------------- X04: the shape of a text in memory is checked before its layout

    private static unsafe void Reduced(byte[] data, long n)
    {
        fixed (byte* p = data)
        {
            using MappedTextReader r = MappedTextReader.FromMemoryReduced(p, data.Length, n);
        }
    }

    private static unsafe void Plain(byte[] data, long n, int bits)
    {
        fixed (byte* p = data)
        {
            using MappedTextReader r = MappedTextReader.FromMemory(p, data.Length, n, bits);
        }
    }

    [Test]
    public void ATextInMemory_WithNoLengthOfATextOrNoWidthOfTheCodec_IsRefusedByTheShapeCheck_NotByTheLayoutCheck()
    {
        // X04 (FromMemoryReduced without CheckShape): the layout check of the same call refuses a wrong n as well, but with the words of a wrong length of the data; the shape check says that the number is no length of a text
        byte[] data = new byte[64 + 16];
        foreach (long n in new long[] { -1, -7, long.MaxValue / 8 })
        {
            InvalidDataException? reduced = Assert.Throws<InvalidDataException>(() => Reduced(data, n), $"reduced, n = {n}");
            Assert.That(reduced!.Message, Does.Contain("not a length of a text"), $"reduced, n = {n}");
            InvalidDataException? plain = Assert.Throws<InvalidDataException>(() => Plain(data, n, 3), $"plain, n = {n}");
            Assert.That(plain!.Message, Does.Contain("not a length of a text"), $"plain, n = {n}");
        }
        foreach (int bits in new[] { 0, -1, 9, 33 })
        {
            InvalidDataException? e = Assert.Throws<InvalidDataException>(() => Plain(data, 1, bits), $"{bits} bits");
            Assert.That(e!.Message, Does.Contain("bits per symbol"), $"{bits} bits");
        }
    }

    // ---------------------------------------------------------------- E13, H07: the limit of the positions at the entry of PH1 and of PH3

    [Test]
    public void ATextAboveTheLimitOfThePositions_IsRefusedAtTheEntryOfPh1AndPh3_BeforeTheScratchIsTouched()
    {
        // E13 (ExternalDcs.BuildRankFile) and H07 (BucketSorter.Run): the limit is checked at the entry of each stage, not only by SaBuilder. A directory cannot be created below a regular file, so a stage that goes on past its entry
        // check fails with an IOException when it makes its scratch directory: only the refusal at the entry gives the NotSupportedException.
        string blocker = _dir.File("blocker");
        File.WriteAllText(blocker, "x");
        string scratch = Path.Combine(blocker, "scratch");
        DifferenceCover cover = DifferenceCover.Create(7);
        SaCase small = SaCase.Of("banana");
        using RankTable ranks = RankTable.Open(SaInputs.WriteRankFile(_dir, small, cover), cover, small.N);
        foreach (long over in new long[] { 1, 2, 1_000_000 })
        {
            long n = PackedTextFormat.MaxSymbols + over;
            string path = _dir.File($"huge{over}.pk");
            using (var f = new FileStream(path, FileMode.Create)) f.SetLength(PackedTextFormat.TotalBytes(n, 1));
            using MappedTextReader text = MappedTextReader.Open(path, n, 1);
            Assert.Throws<NotSupportedException>(() => ExternalDcs.BuildRankFile(text, cover, _dir.File("r.bin"), SaTiny.Dcs(scratch, 1)), $"PH1, {n} symbols");
            Assert.Throws<NotSupportedException>(() => BucketSorter.Run(text, SaTiny.Sort(scratch, 1), new NoSink(), ranks), $"PH3, {n} symbols");
            Assert.Throws<NotSupportedException>(() => SaBuilder.Build(text, _dir.File("o.sa"), SaTiny.Build(scratch, 1)), $"the build, {n} symbols");
            File.Delete(path);
        }
        Assert.That(File.Exists(_dir.File("r.bin")), Is.False);
    }

    [Test]
    public void ATextAtTheLimitOfThePositions_PassesTheEntryOfPh1AndPh3()
    {
        // the control of the test above: the limit itself is not refused at the entry (the stage goes on and fails at the scratch directory, which cannot be made below a regular file)
        string blocker = _dir.File("blocker");
        File.WriteAllText(blocker, "x");
        string scratch = Path.Combine(blocker, "scratch");
        DifferenceCover cover = DifferenceCover.Create(7);
        SaCase small = SaCase.Of("banana");
        using RankTable ranks = RankTable.Open(SaInputs.WriteRankFile(_dir, small, cover), cover, small.N);
        long n = PackedTextFormat.MaxSymbols;
        string path = _dir.File("edge.pk");
        using (var f = new FileStream(path, FileMode.Create)) f.SetLength(PackedTextFormat.TotalBytes(n, 1));
        using MappedTextReader text = MappedTextReader.Open(path, n, 1);
        Assert.That(() => ExternalDcs.BuildRankFile(text, cover, _dir.File("r.bin"), SaTiny.Dcs(scratch, 1)), Throws.Exception.Not.TypeOf<NotSupportedException>(), "PH1");
        Assert.That(() => BucketSorter.Run(text, SaTiny.Sort(scratch, 1), new NoSink(), ranks), Throws.Exception.Not.TypeOf<NotSupportedException>(), "PH3");
        Assert.That(() => SaBuilder.Build(text, _dir.File("o.sa"), SaTiny.Build(scratch, 1)), Throws.Exception.Not.TypeOf<NotSupportedException>(), "the build (F46: the entry of the build is at the limit as well)");
    }
}
