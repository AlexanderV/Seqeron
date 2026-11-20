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
        ph1 = !inPh1.TimedOut;
        ph3 = !inPh3.TimedOut;
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
    }
}
