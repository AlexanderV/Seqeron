using SuffixTree.Esa.External;
using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.Esa.Tests.External;

/// <summary>
/// What PH1-PH3 owe to the resumable container build (work package 1b.9e, <see cref="ISaCheckpoint"/>), tested on the stage itself with a checkpoint of the test's own: the order of the calls (the reuse is asked first, the rank
/// file is recorded when it is complete, the suffix array when it is complete and the rank file is still there), the rank file that a failure must leave when the state vouches for it and must not touch when this run did
/// not build it, the cover that the build made once and the seconds of the recording that are not the seconds of a phase.
/// </summary>
[TestFixture]
public sealed class SaCheckpoint1b9eTests
{
    private const string Text = "mississippimississippimississippiabracadabraabracadabra";

    private sealed class Fake : ISaCheckpoint
    {
        internal List<string> Calls { get; } = [];

        public DifferenceCover? Cover { get; init; }

        internal bool Reuse { get; init; }

        public bool KeepRanks { get; init; }

        public double CheckpointSeconds { get; private set; }

        internal double RanksJump { get; init; }

        internal double SaJump { get; init; }

        internal Action? OnRanks { get; init; }

        internal Action? OnSa { get; init; }

        public bool TryReuseRanks(DifferenceCover cover)
        {
            Calls.Add("reuse");
            return Reuse;
        }

        public void RanksBuilt(DifferenceCover cover)
        {
            Calls.Add("ranks");
            OnRanks?.Invoke();
            CheckpointSeconds += RanksJump;
        }

        public void SaWritten()
        {
            Calls.Add("sa");
            OnSa?.Invoke();
            CheckpointSeconds += SaJump;
        }
    }

    private static SaBuildOptions Options(TempDir dir, ISaCheckpoint cp, Action<string>? trace = null) => SaTiny.Build(dir.Sub("scratch"), 3, 2) with { Checkpoint = cp, Trace = trace };

    [Test]
    public void TheCallsComeInTheOrderOfThePhases_AndEachSeesWhatItShould()
    {
        SaCase c = SaCase.Of(Text);
        using var dir = new TempDir();
        string rankPath = Path.Combine(dir.Sub("scratch"), "ranks.bin");
        string? saPath = null;
        long ranksBytes = -1;
        bool saComplete = false, rankFileThere = false;
        var fake = new Fake
        {
            OnRanks = () => ranksBytes = File.Exists(rankPath) ? new FileInfo(rankPath).Length : -1,
            OnSa = () =>
            {
                rankFileThere = File.Exists(rankPath);
                uint[] sa = FileBits.ReadUInts(Directory.GetFiles(dir.Sub("out"), "o.sa").Single());
                saComplete = sa.SequenceEqual(c.Sa);
                saPath = dir.Sub("out");
            },
        };
        SaRun.Outcome o = SaRun.Run(dir, c, Options(dir, fake));
        SaRun.Compare(c, o, "with a checkpoint");
        Assert.That(fake.Calls, Is.EqualTo(new[] { "reuse", "ranks", "sa" }));
        Assert.That(ranksBytes, Is.EqualTo(o.Result.RankFileBytes), "the rank file is complete when it is recorded");
        Assert.That(ranksBytes, Is.GreaterThan(0));
        Assert.That(saComplete, Is.True, "the suffix array is complete when it is recorded");
        Assert.That(rankFileThere, Is.True, "and the rank file is still there: a kill between the recording and the removal leaves a file that the next run removes, a kill before leaves nothing to trust");
        Assert.That(saPath, Is.Not.Null);
        Assert.That(o.Result.Dcs.Levels, Is.Not.Empty, "PH1 ran");
    }

    [Test]
    public void ARankFileThatIsVouchedFor_SurvivesAFailureInPh3_AndOneThatIsNot_DoesNot()
    {
        foreach (bool keep in new[] { true, false })
        {
            SaCase c = SaCase.Of(Text);
            using var dir = new TempDir();
            string rankPath = Path.Combine(dir.Sub("scratch"), "ranks.bin");
            var fake = new Fake { KeepRanks = keep };
            int seen = 0;
            bool ph3 = false;
            Assert.Throws<InjectedFailureException>(() => SaRun.Run(dir, c, Options(dir, fake, p =>
            {
                if (p == "ph23") ph3 = true;
                if (ph3 && p == "leaf" && ++seen == 1) throw new InjectedFailureException(p);
            })));
            Assert.That(File.Exists(rankPath), Is.EqualTo(keep), $"vouched for: {keep}");
            Assert.That(fake.Calls, Is.EqualTo(new[] { "reuse", "ranks" }), "the failure came after the rank file and before the suffix array");
        }
    }

    [Test]
    public void ARankFileOfAnEarlierRun_IsReused_NotBuiltAgain_NotRecordedAgain_AndNotDeletedByAFailureOfThisRun()
    {
        SaCase c = SaCase.Of(Text);
        using var dir = new TempDir();
        string rankPath = Path.Combine(dir.Sub("scratch"), "ranks.bin");
        int seen = 0;
        bool ph3 = false;
        Assert.Throws<InjectedFailureException>(() => SaRun.Run(dir, c, Options(dir, new Fake { KeepRanks = true }, p =>
        {
            if (p == "ph23") ph3 = true;
            if (ph3 && p == "leaf" && ++seen == 1) throw new InjectedFailureException(p);
        })));
        long length = new FileInfo(rankPath).Length;
        byte[] bytes = File.ReadAllBytes(rankPath);

        // the second run reuses it and fails in PH3 itself, with a checkpoint that does not vouch for the file: a file that this run did not build is not deleted by its failure
        var reuseAndFail = new Fake { Reuse = true };
        seen = 0;
        ph3 = false;
        Assert.Throws<InjectedFailureException>(() => SaRun.Run(dir, c, Options(dir, reuseAndFail, p =>
        {
            if (p == "ph23") ph3 = true;
            if (ph3 && p == "leaf" && ++seen == 1) throw new InjectedFailureException(p);
        })));
        Assert.That(File.Exists(rankPath), Is.True);
        Assert.That(File.ReadAllBytes(rankPath), Is.EqualTo(bytes), "untouched");
        Assert.That(reuseAndFail.Calls, Is.EqualTo(new[] { "reuse" }));

        // the third run reuses it and finishes
        var finish = new Fake { Reuse = true };
        var traces = new List<string>();
        SaRun.Outcome o = SaRun.Run(dir, c, Options(dir, finish, traces.Add));
        Assert.That(FileBits.ReadUInts(o.SaPath), Is.EqualTo(c.Sa), "the suffix array of a build that reused the rank file");
        Assert.That(Directory.Exists(o.ScratchDir) && Directory.GetFileSystemEntries(o.ScratchDir).Length == 0, Is.True, "the scratch root was not made by this run, so it stays, empty: the owner of it removes it");
        Assert.That(File.Exists(o.RankPath), Is.False, "the rank file has served its consumer");
        Assert.That(finish.Calls, Is.EqualTo(new[] { "reuse", "sa" }), "the rank file is not recorded again");
        Assert.That(o.Result.Dcs.Levels, Is.Empty, "PH1 did not run");
        Assert.That(traces, Does.Not.Contain("naming").And.Not.Contain("rank-written"), "no sample bucket was named");
        Assert.That(o.Result.RankFileBytes, Is.EqualTo(length), "the rank file is on the ledger of this run");
        Assert.That(o.Result.Ph23ScratchPeakBytes, Is.GreaterThanOrEqualTo(length));
    }

    [Test]
    public void TheCoverThatTheCheckpointHolds_IsTheCoverOfTheBuild()
    {
        SaCase c = SaCase.Of(Text);
        using var dir = new TempDir();
        DifferenceCover other = DifferenceCover.Create(9);
        SaBuildOptions options = Options(dir, new Fake { Cover = other });
        Assert.That(options.Period, Is.Not.EqualTo(9), "the options ask for another period");
        SaRun.Outcome o = SaRun.Run(dir, c, options);
        SaRun.Compare(c, o, "another cover");
        Assert.That((o.Result.Period, o.Result.CoverResidues, o.Result.CoverKind), Is.EqualTo((other.Period, other.Residues.Length, other.Kind)), "the build used the cover that it was handed instead of making its own");
        using var dir2 = new TempDir();
        SaRun.Outcome plain = SaRun.Run(dir2, c, Options(dir2, new Fake()));
        Assert.That(plain.Result.Period, Is.Not.EqualTo(9), "without one it makes the cover of its options");
    }

    [Test]
    public void TheSecondsOfTheRecording_AreNotTheSecondsOfThePhases()
    {
        SaCase c = SaCase.Of(Text);
        using var dir = new TempDir();
        SaRun.Outcome o = SaRun.Run(dir, c, Options(dir, new Fake { RanksJump = 1000, SaJump = 2000 }));
        Assert.That(o.Result.Ph1Seconds, Is.InRange(-1000.0, -990.0), "PH1 without the 1000 s that its recording took");
        Assert.That(o.Result.Ph23Seconds, Is.InRange(-2000.0, -1990.0), "PH2 and PH3 without the 2000 s of theirs");
    }

    [Test]
    public void WithoutACheckpoint_NothingChanges_TheSecondsAreTheSecondsOfThePhases()
    {
        SaCase c = SaCase.Of(Text);
        using var dir = new TempDir();
        SaRun.Outcome o = SaRun.Run(dir, c, SaTiny.Build(dir.Sub("scratch"), 3, 2));
        SaRun.Compare(c, o, "no checkpoint");
        Assert.That(o.Result.Ph1Seconds, Is.GreaterThanOrEqualTo(0));
        Assert.That(o.Result.Ph23Seconds, Is.GreaterThanOrEqualTo(0));
    }
}
