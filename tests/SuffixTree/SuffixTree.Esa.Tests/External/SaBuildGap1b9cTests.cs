using System.Diagnostics;
using SuffixTree.Esa.External;
using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.Esa.Tests.External;

/// <summary>
/// The gaps that the audit of work package 1b.9c showed in the report and in the wiring of <see cref="SaBuilder"/> (the mutants Q88-Q126 of <c>analysis/mutation/1b7d_sa_sabuilder_listed.json</c>): which directory PH3 spills to,
/// whether it sees the hook of the build, whether the output is synced when it should be, the length check of the output, the figures and the verdict of the RssAnon gate, the seconds of the phases, the period that was used,
/// and the scratch bytes of a build that failed. Each test names the mutants it kills.
/// </summary>
[TestFixture]
public class SaBuildGap1b9cTests
{
    private TempDir _dir = null!;

    [SetUp]
    public void SetUp() => _dir = new TempDir();

    [TearDown]
    public void TearDown() => _dir.Dispose();

    private string Scratch => _dir.Sub("scratch");

    // ---------------------------------------------------------------- Q88, Q94: PH3 spills to its own directory and sees the hook of the build

    [Test]
    public void ThePositionFilesOfPh3_AreInTheDirectoryOfPh3_AndThePhaseSeesTheHookOfTheBuild()
    {
        // Q88 (PH3 gets the scratch directory of PH1) and Q94 (PH3 gets no hook): the hook "handover" of the partition fires when the writers of the root bucket hold their spill files; after "ph23" nothing of PH1 runs any more
        SaCase c = SaCase.Of(StringFamilies.RandomDna(3000, 3));
        string scratch = Scratch;
        string ph1 = Path.Combine(scratch, "ph1"), ph3 = Path.Combine(scratch, "ph3");
        bool inPh3 = false;
        var events = new List<string>();
        string[]? inPh3Dir = null, inPh1Dir = null;
        SaBuildOptions options = SaTiny.Build(scratch, 1) with
        {
            Trace = step =>
            {
                if (step == "ph23") inPh3 = true;
                if (inPh3) events.Add(step);
                if (inPh3 && step == "handover" && inPh3Dir is null)
                {
                    inPh3Dir = Directory.GetFiles(ph3);
                    inPh1Dir = Directory.GetFiles(ph1);
                }
            },
        };
        SaRun.Outcome o = SaRun.Run(_dir, c, options);
        SaRun.Compare(c, o, "traced");
        Assert.That(events, Does.Contain("partition").And.Contain("leaf"), "the phase calls the hook of the build");
        Assert.That(inPh3Dir, Is.Not.Null.And.Not.Empty, "the position files of the partition are in the directory of PH3");
        Assert.That(inPh1Dir, Is.Empty, "and nothing of them in the directory of PH1");
    }

    // ---------------------------------------------------------------- Q98, Q99: the output is synced when the option says so, and only then

    [Test]
    public void TheOutputIsSyncedOnce_WithTheOption_AndNeverWithout()
    {
        // Q98 (the output is always synced) and Q99 (the output is never synced): the fsync of the file of the suffix array is counted as a flush (MappedOutputFile.Complete(fsync, counters))
        SaCase c = SaCase.Of(StringFamilies.RandomDna(800, 3));
        foreach ((bool sync, long flushes) in new[] { (false, 0L), (true, 1L) })
        {
            using var dir = new TempDir();
            var counters = new IoCounters();
            SaRun.Outcome o = SaRun.Run(dir, c, SaTiny.Build(dir.Sub("scratch"), 1) with { SyncOutput = sync, Counters = counters });
            SaRun.Compare(c, o, $"sync {sync}");
            Assert.That(counters.Flushes, Is.EqualTo(flushes), $"SyncOutput = {sync}");
        }
    }

    // ---------------------------------------------------------------- Q100: the output has the length of the plan when the build ends

    [Test]
    public void AnOutputThatChangedItsLengthDuringTheBuild_IsRefused_AndRemoved()
    {
        // Q100 (the final length is not checked): the hook of PH3 lengthens the file of the suffix array (the mapping of the build stays valid); the check at the end must say so, and the failed build removes the file
        SaCase c = SaCase.Of(StringFamilies.RandomDna(600, 3));
        string textPath = SaRun.WriteText(_dir, c);
        string outDir = _dir.Sub("out");
        Directory.CreateDirectory(outDir);
        string saPath = Path.Combine(outDir, "o.sa");
        using MappedTextReader reader = MappedTextReader.Open(textPath, c.N, c.Bits);
        SaBuildOptions options = SaTiny.Build(Scratch, 1) with
        {
            Trace = step =>
            {
                if (step != "ph23") return;
                using var stream = new FileStream(saPath, FileMode.Open, FileAccess.Write, FileShare.ReadWrite);
                stream.SetLength(stream.Length + 100);
            },
        };
        var e = Assert.Throws<InvalidOperationException>(() => SaBuilder.Build(reader, saPath, options));
        Assert.That(e!.Message, Does.Contain($"the suffix array has {(4 * (c.N + 1)) + 100} bytes, expected {4 * (c.N + 1)}"));
        Assert.That(File.Exists(saPath), Is.False, "the output of a failed build is removed");
    }

    // ---------------------------------------------------------------- Q102, Q104, Q106: the figures and the verdict of the RssAnon gate

    private static string Status(long anonKb) => $"Name:\tx\nVmHWM:\t{anonKb + 1000} kB\nRssAnon:\t{anonKb} kB\nRssFile:\t5 kB\nThreads:\t4\n";

    [Test]
    public void TheRssAnonFigures_AreInMegabytesOfTenToTheSixth_AndTheVerdictIsWithinAtTheCapItself()
    {
        // Q102 (the verdict is Within only below the cap), Q104 (the figure is in KiB): a source that always reads the same RssAnon; 312,500 KiB is 320 MB exactly
        SaCase c = SaCase.Of(StringFamilies.RandomDna(900, 3));
        foreach ((long anonKb, double cap, RssVerdict verdict, double mb) in new[]
        {
            (100_000L, 320.0, RssVerdict.Within, 102.4),
            (312_500L, 320.0, RssVerdict.Within, 320.0),
            (312_501L, 320.0, RssVerdict.Exceeded, 320.00102),
            (400_000L, 320.0, RssVerdict.Exceeded, 409.6),
            (400_000L, 500.0, RssVerdict.Within, 409.6),
        })
        {
            using var dir = new TempDir();
            var source = new FakeProcSource().With("/proc/self/status", Status(anonKb));
            SaRun.Outcome o = SaRun.Run(dir, c, SaTiny.Build(dir.Sub("scratch"), 1) with { ProcSource = source, RssAnonCapMb = cap });
            Assert.That(o.Result.RssVerdict, Is.EqualTo(verdict), $"{anonKb} KiB against {cap} MB");
            Assert.That(o.Result.RssAnonPeakMb, Is.EqualTo(mb).Within(1e-4), $"{anonKb} KiB: the whole run");
            Assert.That(o.Result.Ph1RssAnonPeakMb, Is.EqualTo(mb).Within(1e-4), $"{anonKb} KiB: PH1");
            Assert.That(o.Result.Ph23RssAnonPeakMb, Is.EqualTo(mb).Within(1e-4), $"{anonKb} KiB: PH2 / PH3");
        }
    }

    [Test]
    public void TheRssAnonFigureOfPh3_IsThatOfPh3_NotThatOfPh1()
    {
        // Q106 (the PH3 figure is taken from the PH1 sampler): the source rises with every reading, and PH3 begins after PH1 has ended, so every reading of PH3 is above every reading of PH1
        SaCase c = SaCase.Of(StringFamilies.RandomDna(900, 3));
        var source = new FakeProcSource().WithSeries("/proc/self/status", call => Status(10_000 + (1000L * call)));
        SaRun.Outcome o = SaRun.Run(_dir, c, SaTiny.Build(Scratch, 1) with { ProcSource = source });
        Assert.That(o.Result.Ph1RssAnonPeakMb, Is.Not.Null);
        Assert.That(o.Result.Ph23RssAnonPeakMb, Is.Not.Null);
        Assert.That(o.Result.Ph23RssAnonPeakMb, Is.GreaterThan(o.Result.Ph1RssAnonPeakMb), "the phase that began later saw a higher figure");
        Assert.That(o.Result.RssAnonPeakMb, Is.GreaterThanOrEqualTo(o.Result.Ph23RssAnonPeakMb), "the whole run saw at least what PH3 saw");
    }

    // ---------------------------------------------------------------- Q112: the seconds of the phases

    [Test]
    public void TheSecondsOfThePhases_AreMeasuredWithinTheBuild_AndAddUpToNoMoreThanItsWallTime()
    {
        // Q112 (the PH1 time includes PH3: here a second more): the phases are consecutive intervals of the build, so their sum cannot exceed the time that the call took
        SaCase c = SaCase.Of(StringFamilies.RandomDna(1500, 3));
        string textPath = SaRun.WriteText(_dir, c);
        string saPath = Path.Combine(_dir.Sub("out"), "o.sa");
        Directory.CreateDirectory(Path.GetDirectoryName(saPath)!);
        using MappedTextReader reader = MappedTextReader.Open(textPath, c.N, c.Bits);
        long t0 = Stopwatch.GetTimestamp();
        SaBuildResult r = SaBuilder.Build(reader, saPath, SaTiny.Build(Scratch, 1));
        double wall = (Stopwatch.GetTimestamp() - t0) / (double)Stopwatch.Frequency;
        Assert.That(r.Ph1Seconds, Is.GreaterThan(0));
        Assert.That(r.Ph23Seconds, Is.GreaterThan(0));
        Assert.That(r.Seconds, Is.LessThanOrEqualTo(wall), "PH1 and PH2 / PH3 together took no longer than the call");
    }

    // ---------------------------------------------------------------- Q126: the period that was used

    [Test]
    public void TheResultNamesThePeriodOfTheCoverThatWasUsed_NotTheOneThatWasAskedFor()
    {
        // Q126 (the period is the nominal one): the nominal 1024 becomes the Singer period 993 with 32 residues unless the exact period is asked for
        SaCase c = SaCase.Of(StringFamilies.RandomDna(700, 3));
        SaRun.Outcome nominal = SaRun.Run(_dir, c, SaTiny.Build(Scratch, 1) with { ExactPeriod = false, Period = 1024 });
        Assert.That((nominal.Result.Period, nominal.Result.CoverResidues, nominal.Result.CoverKind), Is.EqualTo((993, 32, CoverKind.Singer)));
        using var dir2 = new TempDir();
        SaRun.Outcome exact = SaRun.Run(dir2, c, SaTiny.Build(dir2.Sub("scratch"), 1) with { ExactPeriod = true, Period = 1024 });
        Assert.That(exact.Result.Period, Is.EqualTo(1024));
        Assert.That(exact.Result.CoverKind, Is.EqualTo(CoverKind.Pruned));
    }

    // ---------------------------------------------------------------- Q116: the scratch bytes of a failed build are given back

    [Test]
    public void ABuildThatFailsAfterTheRankFile_GivesItsScratchBytesBack_AtThePointWhereThePhaseBegins_AndInTheMiddleOfPh3()
    {
        // Q116 (the scratch bytes of a failed build are not given back): from the end of PH1 the ledger of the build holds the bytes of the rank file, and in PH3 those of the position files; the account of the counters must be empty again
        SaCase c = SaCase.Of(StringFamilies.RandomDna(3000, 3));
        foreach ((string point, int occurrence) in new[] { ("ph23", 1), ("leaf", 2), ("handover", 2) })
        {
            using var dir = new TempDir();
            int seen = 0;
            bool inPh3 = false;
            SaBuildOptions options = SaTiny.Build(dir.Sub("scratch"), 1) with
            {
                Counters = new IoCounters(),
                Trace = step =>
                {
                    if (step == "ph23") inPh3 = true;
                    if (inPh3 && step == point && ++seen == occurrence) throw new InjectedFailureException($"{point} #{occurrence}");
                },
            };
            SaRun.AssertFailsCleanly(dir, c, options, typeof(InjectedFailureException), $"{point} #{occurrence}");
            Assert.That(seen, Is.EqualTo(occurrence), $"{point} #{occurrence}: the hook point was reached");
        }
    }

    // ---------------------------------------------------------------- F45: each phase has its own seconds

    [Test]
    public void AWaitInOnePhase_ShowsInTheSecondsOfThatPhaseAlone()
    {
        // F45 (the PH1 time is the time of PH3): a hook that waits 400 ms once, inside PH1 (the first "rank-written") or inside PH3 (the first "leaf"), must show in the figure of that phase and not in the other
        SaCase c = SaCase.Of(StringFamilies.RandomDna(800, 3));
        foreach (string phase in new[] { "ph1", "ph3" })
        {
            using var dir = new TempDir();
            bool inPh3 = false, waited = false;
            SaBuildOptions options = SaTiny.Build(dir.Sub("scratch"), 1) with
            {
                Trace = step =>
                {
                    if (step == "ph23") inPh3 = true;
                    bool here = phase == "ph1" ? step == "rank-written" : inPh3 && step == "leaf";
                    if (here && !waited)
                    {
                        waited = true;
                        using var gate = new ManualResetEventSlim(false);
                        gate.Wait(400); // nothing sets it: a wait of 400 ms
                    }
                },
            };
            SaRun.Outcome o = SaRun.Run(dir, c, options);
            Assert.That(waited, Is.True, phase);
            double waiting = phase == "ph1" ? o.Result.Ph1Seconds : o.Result.Ph23Seconds;
            double other = phase == "ph1" ? o.Result.Ph23Seconds : o.Result.Ph1Seconds;
            Assert.That(waiting, Is.GreaterThanOrEqualTo(0.4), $"the phase that waited ({phase})");
            Assert.That(other, Is.LessThan(0.4), $"the other phase ({phase} waited)");
        }
    }

    // ---------------------------------------------------------------- F49, F62: the scratch bytes of a build that succeeded

    [Test]
    public void ABuildThatSucceeded_LeavesNoScratchByteInTheAccount_AndPh3ReportsItsOwnPeakWhichHoldsTheRankFile()
    {
        // F49 (the rank file is not booked in the ledger: the account ends at minus the size of the rank file) and F62 (the PH3 peak is that of PH1, which does not count the rank file yet)
        SaCase c = SaCase.Of(StringFamilies.RandomDna(2500, 3));
        var counters = new IoCounters();
        SaBuildOptions options = SaTiny.Build(Scratch, 1) with { Counters = counters };
        SaRun.Outcome o = SaRun.Run(_dir, c, options);
        SaRun.Compare(c, o, "counted");
        Assert.That(counters.ScratchLiveBytes, Is.Zero, "every scratch byte that the build took was given back");
        Assert.That(o.Result.Ph23ScratchPeakBytes, Is.EqualTo(counters.ScratchPeakBytes), "the peak of the whole build is reached in PH3 (the rank file and the position files)");
        Assert.That(o.Result.Ph23ScratchPeakBytes, Is.GreaterThan(o.Result.Ph1ScratchPeakBytes), "the rank file is booked from the end of PH1 on, so the peak of PH3 is above the peak of PH1");
        Assert.That(o.Result.Ph23ScratchPeakBytes, Is.GreaterThanOrEqualTo(o.Result.RankFileBytes));
    }
}
