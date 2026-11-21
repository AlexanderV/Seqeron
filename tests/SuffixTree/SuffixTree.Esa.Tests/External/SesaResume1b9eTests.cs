using System.Collections.Concurrent;
using System.Security.Cryptography;
using SuffixTree.Esa.External;
using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.Esa.Tests.External;

/// <summary>A reservation source that reports a fixed free space and otherwise behaves like the system's.</summary>
internal sealed class FixedFreeSource(long free) : IReservationSource
{
    public bool CanPreallocate => true;

    public FileStream CreatePreallocated(string path, long length) => SystemReservationSource.Instance.CreatePreallocated(path, length);

    public long AvailableFreeBytes(string directory) => free;
}

/// <summary>The directory of one build of a resumable container: the input, the output directory, the paths the build derives from the target.</summary>
internal sealed class ResumeCase : IDisposable
{
    private readonly TempDir _dir = new();

    internal ResumeCase(string text)
    {
        Text = text;
        Input = SesaRun.WriteInput(_dir, text);
        OutDir = _dir.Sub("out");
        Directory.CreateDirectory(OutDir);
        Target = Path.Combine(OutDir, "x.sesa");
    }

    internal string Text { get; }

    internal string Input { get; }

    internal string OutDir { get; }

    internal string Target { get; }

    internal string Tmp => Target + ".tmp";

    internal string Scratch => Tmp + ".d";

    internal string State => Path.Combine(Scratch, SesaResumeFormat.FileName);

    internal string Lock => Target + ".lock";

    internal string[] Entries => [.. Directory.GetFileSystemEntries(OutDir).Select(Path.GetFileName).Order(StringComparer.Ordinal)!];

    public void Dispose() => _dir.Dispose();
}

/// <summary>Resume of the production build (work package 1b.9e): the failure at every boundary and inside every phase, the files that a stage vouches for and what is done when they are not what they were, the state file, the lock.</summary>
[TestFixture]
public sealed class SesaResume1b9eTests
{
    private static readonly string s_text = StringFamilies.RandomDna(1500, 11) + new string('A', 300) + StringFamilies.RandomDna(900, 12);

    private static readonly ConcurrentDictionary<string, byte[]> s_reference = new();

    /// <summary>An LCP of two bytes per entry (the container grows behind its reserved length).</summary>
    private static string TextU16 => StringFamilies.RandomDna(30_000, 5) + new string('A', 2_600);

    private static string TextU32 => new('a', 100_000);

    private static string TextOverflow => StringFamilies.RandomDna(50_000, 6) + new string('A', 700);

    internal static readonly string[] HookNames = [.. SesaHooks.All];

    private static byte[] Read(string path)
    {
        using var fs = new FileStream(path, FileMode.Open, FileAccess.Read, FileShare.ReadWrite | FileShare.Delete);
        var bytes = new byte[fs.Length];
        fs.ReadExactly(bytes);
        return bytes;
    }

    /// <summary>The container of an uninterrupted build without Resume (the oracle of every test here: the same bytes whatever the options).</summary>
    private static byte[] Reference(string text, bool psi = true, ChildCodec? codec = null)
    {
        string key = $"{text.GetHashCode(StringComparison.Ordinal)}-{text.Length}-{psi}-{codec}";
        return s_reference.GetOrAdd(key, _ =>
        {
            using var c = new ResumeCase(text);
            SesaBuilder.Build(c.Input, c.Target, new SesaBuildOptions { Threads = 2, IncludePsi = psi, ChildCodec = codec, Fsync = false, MeasureRss = false });
            return Read(c.Target);
        });
    }

    private static SesaBuildOptions Opts(int v = 3, int threads = 2, bool psi = true, ChildCodec? codec = null) => SesaRun.Tiny(v, threads, psi) with { Resume = true, ChildCodec = codec, MeasureRss = false };

    private static Action<string> FailAt(string point, int nth = 1)
    {
        int seen = 0;
        return p =>
        {
            if (p == point && ++seen == nth) throw new InjectedFailureException(p);
        };
    }

    private static Action<string> FailAtTrace(string point)
    {
        int seen = 0;
        return p =>
        {
            if (p == point && Interlocked.Increment(ref seen) == 1) throw new InjectedFailureException(p);
        };
    }

    private static void Fails(ResumeCase c, SesaBuildOptions options, string point)
    {
        InjectedFailureException? e = Assert.Throws<InjectedFailureException>(() => SesaBuilder.Build(c.Input, c.Target, options with { Hook = FailAt(point) }));
        Assert.That(e!.Message, Does.Contain(point));
    }

    private static void FailsInside(ResumeCase c, SesaBuildOptions options, string trace, bool sa)
    {
        SesaBuildOptions failing = sa ? options with { SaTrace = FailAtTrace(trace) } : options with { LcpTrace = FailAtTrace(trace) };
        Assert.Throws<InjectedFailureException>(() => SesaBuilder.Build(c.Input, c.Target, failing));
    }

    private static int HookIndex(string point) => Array.IndexOf(HookNames, point);

    /// <summary>The stages that are complete when the build is stopped at <paramref name="point"/> (a build with PSI, or without).</summary>
    internal static string[] DoneAt(string point, bool psi)
    {
        int h = HookIndex(point);
        if (h < HookIndex(SesaHooks.Packed)) return [];
        if (h < HookIndex(SesaHooks.SaBuilt)) return [SesaStages.Packed];
        if (h < HookIndex(SesaHooks.LcpChildBuilt)) return [SesaStages.Packed, SesaStages.Sa];
        if (!psi) return [SesaStages.Packed, SesaStages.Sa, SesaStages.Lcp];
        if (h < HookIndex(SesaHooks.BwtBuilt)) return [SesaStages.Packed, SesaStages.Sa, SesaStages.Lcp];
        if (h < HookIndex(SesaHooks.PsiBuilt)) return [SesaStages.Packed, SesaStages.Sa, SesaStages.Lcp, SesaStages.Bwt];
        return [SesaStages.Packed, SesaStages.Sa, SesaStages.Lcp, SesaStages.Psi];
    }

    private static void AssertOnlyTarget(ResumeCase c, string label)
    {
        Assert.That(c.Entries, Is.EqualTo(new[] { "x.sesa" }), $"{label}: nothing but the target is left (no temporary file, scratch, state or lock)");
    }

    private static void AssertSame(ResumeCase c, byte[] reference, string label)
    {
        byte[] got = Read(c.Target);
        Assert.That(got.Length, Is.EqualTo(reference.Length), $"{label}: length");
        Assert.That(SHA256.HashData(got), Is.EqualTo(SHA256.HashData(reference)), $"{label}: the bytes of an uninterrupted build");
    }

    private static void FlipByte(string path, long offset)
    {
        using var fs = new FileStream(path, FileMode.Open, FileAccess.ReadWrite, FileShare.ReadWrite);
        fs.Position = offset;
        int b = fs.ReadByte();
        fs.Position = offset;
        fs.WriteByte((byte)(b ^ 0x5A));
    }

    private static SesaPlan PlanOf(ResumeCase c, bool psi) => SesaLayout.Plan(PackedTextWriter.Scan(c.Input), psi);

    // ------------------------------------------------------------------ the default does not change

    [Test]
    public void TheDefaultBuild_WritesNoStateFile_TakesNoLock_AndFlushesTwice()
    {
        using var c = new ResumeCase(s_text);
        var counters = new IoCounters();
        var seen = new List<string>();
        SesaBuildResult r = SesaBuilder.Build(c.Input, c.Target, new SesaBuildOptions
        {
            Threads = 2,
            Counters = counters,
            MeasureRss = false,
            Hook = p =>
            {
                seen.Add(p);
                Assert.That(File.Exists(c.State), Is.False, $"{p}: a state file");
                Assert.That(File.Exists(c.Lock), Is.False, $"{p}: a lock file");
            },
        });
        Assert.That(seen, Is.EqualTo(HookNames));
        Assert.That(r.Resume, Is.Null);
        Assert.That(counters.Flushes, Is.EqualTo(2), "the two fsyncs of the publication, none more");
        AssertOnlyTarget(c, "default");
    }

    [Test]
    public void TheDefaultBuild_WhenItFails_RemovesEverything_AsBefore()
    {
        foreach (string point in HookNames)
        {
            using var c = new ResumeCase(s_text);
            Assert.Throws<InjectedFailureException>(() => SesaBuilder.Build(c.Input, c.Target, new SesaBuildOptions { Threads = 2, MeasureRss = false, Hook = FailAt(point) }), point);
            Assert.That(c.Entries, Is.Empty, $"{point}: nothing is left by a failed build without Resume");
        }
    }

    // ------------------------------------------------------------------ a fresh build with Resume

    [Test]
    public void AFreshBuildWithResume_HasTheBytesOfTheDefault_RecordsEveryStageBeforeItsHook_AndLeavesNothing([Values(true, false)] bool psi, [Values(false, true)] bool packedChild)
    {
        ChildCodec codec = packedChild ? ChildCodec.Packed : ChildCodec.Tagged;
        using var c = new ResumeCase(s_text);
        var states = new List<(string Point, string[] Stages)>();
        var seen = new List<string>();
        SesaBuildResult r = SesaBuilder.Build(c.Input, c.Target, Opts(3, 2, psi, codec) with
        {
            Hook = p =>
            {
                seen.Add(p);
                if (File.Exists(c.State))
                {
                    Assert.That(SesaResumeFormat.TryParse(File.ReadAllBytes(c.State), out SesaResumeState? s), Is.True, $"{p}: the state file is valid");
                    states.Add((p, [.. s!.Stages.Select(x => x.Name)]));
                }
            },
        });
        AssertSame(c, Reference(s_text, psi, codec), $"psi {psi} {codec}");
        AssertOnlyTarget(c, "fresh");
        Assert.That(seen, Is.EqualTo(psi ? HookNames : HookNames.Where(h => h is not (SesaHooks.BwtBuilt or SesaHooks.PsiBuilt))));
        foreach ((string point, string[] stages) in states) Assert.That(stages, Is.EqualTo(DoneAt(point, psi)), $"the state at {point}");
        Assert.That(states.Select(s => s.Point), Is.EqualTo(seen.SkipWhile(h => h != SesaHooks.Packed).ToArray()), "a state exists from the boundary after PH0 on, and not before");
        Assert.That(r.Resume, Is.Not.Null);
        Assert.That((r.Resume!.StateFound, r.Resume.Accepted, r.Resume.Skipped), Is.EqualTo((false, false, Array.Empty<string>())));
        Assert.That(r.Resume.Checkpoints.Select(x => x.Stage), Is.EqualTo(psi ? new[] { "packed", "ranks", "sa", "lcp", "bwt", "psi" } : new[] { "packed", "ranks", "sa", "lcp" }));
        Assert.That(r.Resume.StateWrites, Is.EqualTo(psi ? 6 : 4));
        Assert.That(r.Resume.StateBytes, Is.GreaterThan(0));
    }

    [Test]
    public void TheStateFile_IsTheOnlyThingThatResumeAddsToTheCounters_AndFsyncOffMeansNoFlush()
    {
        foreach (bool psi in new[] { true, false })
        {
            using var plain = new ResumeCase(s_text);
            using var rec = new ResumeCase(s_text);
            var c0 = new IoCounters();
            var c1 = new IoCounters();
            SesaBuildOptions baseOptions = SesaRun.Tiny(3, 2, psi) with { MeasureRss = false };
            SesaBuildResult a = SesaBuilder.Build(plain.Input, plain.Target, baseOptions with { Counters = c0 });
            SesaBuildResult b = SesaBuilder.Build(rec.Input, rec.Target, baseOptions with { Counters = c1, Resume = true });
            Assert.That(c1.ReadBytes, Is.EqualTo(c0.ReadBytes), "the CRC pass reads through a mapping, as the final one does: no explicit read is added");
            Assert.That(c1.MappedWriteBytes, Is.EqualTo(c0.MappedWriteBytes));
            Assert.That(c1.ScratchPeakBytes, Is.EqualTo(c0.ScratchPeakBytes));
            Assert.That(c1.WriteBytes - c0.WriteBytes, Is.EqualTo(b.Resume!.StateBytes), "the extra explicit writes are the bytes of the state files");
            Assert.That(c1.WriteCalls - c0.WriteCalls, Is.EqualTo(b.Resume.StateWrites));
            Assert.That(c1.Flushes, Is.EqualTo(0), "Fsync is off: the state is recorded without a flush");
            Assert.That(c0.Flushes, Is.EqualTo(0));
            Assert.That(a.Resume, Is.Null);
        }
    }

    [Test]
    public void WithFsync_EveryStageIsForcedToStorageBeforeItIsMarked_TwoFsyncsPerStage()
    {
        foreach ((bool psi, int stages) in new[] { (true, 6), (false, 4) })
        {
            using var c = new ResumeCase(s_text);
            var counters = new IoCounters();
            SesaBuildResult r = SesaBuilder.Build(c.Input, c.Target, Opts(3, 2, psi) with { Fsync = true, Counters = counters });
            Assert.That(counters.Flushes, Is.EqualTo((2 * stages) + 2), "one fsync of the outputs and one of the state file per stage, and the two of the publication");
            Assert.That(r.Resume!.Checkpoints.Count, Is.EqualTo(stages));
            Assert.That(r.Resume.Checkpoints.All(x => x.FsyncSeconds >= 0 && x.CrcSeconds >= 0 && x.StateSeconds >= 0 && x.CrcBytes > 0), Is.True);
        }
    }

    // ------------------------------------------------------------------ a failure at every boundary, then a resume

    [Test]
    public void AFailureAtEveryBoundary_ThenAResume_IsByteIdentical_AndSkipsExactlyTheCompletedStages(
        [ValueSource(nameof(HookNames))] string point,
        [Values(true, false)] bool psi,
        [Values(0, 1, 2)] int threadsCase)
    {
        if (!psi && point is SesaHooks.BwtBuilt or SesaHooks.PsiBuilt) Assert.Ignore("a build without PSI has no such boundary");
        (int first, int second) = threadsCase switch { 0 => (1, 4), 1 => (4, 1), _ => (4, 4) };
        using var c = new ResumeCase(s_text);
        Fails(c, Opts(3, first, psi), point);
        string[] done = DoneAt(point, psi);
        Assert.That(File.Exists(c.Target), Is.False, "no target after a failure");
        Assert.That(File.Exists(c.State), Is.EqualTo(done.Length > 0), $"{point}: the state of a failed build stays when a stage is complete, and only then");
        Assert.That(File.Exists(c.Tmp), Is.EqualTo(done.Length > 0), $"{point}: the temporary file stays with its state");
        Assert.That(File.Exists(c.Lock), Is.False, "the lock is given up");

        var counters = new IoCounters();
        var seen = new List<string>();
        SesaBuildResult r = SesaBuilder.Build(c.Input, c.Target, Opts(8, second, psi) with { Counters = counters, Hook = seen.Add });
        AssertSame(c, Reference(s_text, psi), $"{point} psi {psi} {first}->{second} threads");
        AssertOnlyTarget(c, point);
        Assert.That(r.Resume!.Skipped, Is.EqualTo(done), $"{point}: the stages that were not run again");
        Assert.That(r.Resume.Accepted, Is.EqualTo(done.Length > 0));
        Assert.That(seen, Is.EqualTo(psi ? HookNames : HookNames.Where(h => h is not (SesaHooks.BwtBuilt or SesaHooks.PsiBuilt))), "every boundary is reached, by running the stage or by finding it done");
        Assert.That(counters.ScratchLiveBytes, Is.Zero, "the scratch on the account is back to zero");
    }

    [Test]
    public void AFailureInsideAPhase_ThenAResume_IsByteIdentical([Values("ph1", "rank-written", "ph23", "partition", "leaf", "handover", "scan")] string trace, [Values(1, 4)] int threads)
    {
        using var c = new ResumeCase(s_text);
        SesaBuildOptions o = Opts(3, threads);
        FailsInside(c, o, trace, sa: true);
        string[] done = trace is "ph1" or "rank-written" ? [SesaStages.Packed] : [SesaStages.Packed, SesaStages.Ranks];
        Assert.That(File.Exists(c.State), Is.True);
        SesaBuildResult r = SesaBuilder.Build(c.Input, c.Target, Opts(3, threads == 1 ? 4 : 1));
        AssertSame(c, Reference(s_text), trace);
        AssertOnlyTarget(c, trace);
        Assert.That(r.Resume!.Skipped, Is.EqualTo(done), trace);
    }

    [Test]
    public void AFailureInsidePh4_ThenAResume_RunsPh4Again_FromTheSuffixArray([Values("b1", "b2", "b3")] string trace)
    {
        using var c = new ResumeCase(s_text);
        FailsInside(c, Opts(5, 2), trace, sa: false);
        SesaBuildResult r = SesaBuilder.Build(c.Input, c.Target, Opts(5, 2));
        AssertSame(c, Reference(s_text), trace);
        AssertOnlyTarget(c, trace);
        Assert.That(r.Resume!.Skipped, Is.EqualTo(new[] { SesaStages.Packed, SesaStages.Sa }));
        Assert.That(r.Sa.Ph23Seconds, Is.Zero, "PH2 and PH3 did not run");
        Assert.That(r.Lcp.Seconds, Is.GreaterThan(0), "PH4 and PH5 did");
    }

    [Test]
    public void AResumeThatIsItselfStopped_AtAnotherBoundary_ResumesAgain([Values(0, 1, 2, 3, 4, 5)] int variant)
    {
        (string a, string b) = variant switch
        {
            0 => (SesaHooks.Packed, SesaHooks.SaBuilt),
            1 => (SesaHooks.SaBuilt, SesaHooks.LcpChildBuilt),
            2 => (SesaHooks.LcpChildBuilt, SesaHooks.PsiBuilt),
            3 => (SesaHooks.PsiBuilt, SesaHooks.Checksummed),
            4 => (SesaHooks.SaBuilt, SesaHooks.Packed), // the second run stops at the boundary of a stage that it found done: nothing changes
            _ => (SesaHooks.MagicWritten, SesaHooks.HeaderWritten),
        };
        using var c = new ResumeCase(s_text);
        Fails(c, Opts(1, 4), a);
        Fails(c, Opts(2, 1), b);
        SesaBuildResult r = SesaBuilder.Build(c.Input, c.Target, Opts(9, 2));
        AssertSame(c, Reference(s_text), $"{a} then {b}");
        AssertOnlyTarget(c, $"{a} then {b}");
        string[] expected = HookIndex(b) >= HookIndex(a) ? DoneAt(b, true) : DoneAt(a, true);
        Assert.That(r.Resume!.Skipped, Is.EqualTo(expected));
    }

    [Test]
    public void ABuildThatIsStoppedAtEveryBoundaryInTurn_IsFinishedByTheLastResume()
    {
        using var c = new ResumeCase(s_text);
        for (int i = 0; i < HookNames.Length; i++) Fails(c, Opts(i, 1 + (i % 4)), HookNames[i]);
        SesaBuildResult r = SesaBuilder.Build(c.Input, c.Target, Opts(4, 3));
        AssertSame(c, Reference(s_text), "stopped at all 14 boundaries in turn");
        AssertOnlyTarget(c, "stopped at all 14 boundaries in turn");
        Assert.That(r.Resume!.Skipped, Is.EqualTo(DoneAt(SesaHooks.MagicWritten, true)));
    }

    // ------------------------------------------------------------------ what the key is and what it is not

    [Test]
    public void ParametersThatCannotChangeAByte_DoNotInvalidateTheState_ButThePsiOptionAndTheChildCodecDo()
    {
        using var c = new ResumeCase(s_text);
        Fails(c, Opts(3, 1), SesaHooks.LcpChildBuilt);
        SesaBuildOptions other = Opts(7, 4) with
        {
            ScaleN = 10_000,
            RuleN = 2_700,
            Seed = 99,
            TextBlockSymbols = 16,
            TextWindowBytes = 77,
            ArenaCapBytes = 300L << 20,
            Fsync = true,
        };
        SesaBuildResult r = SesaBuilder.Build(c.Input, c.Target, other);
        AssertSame(c, Reference(s_text), "other parameters");
        Assert.That(r.Resume!.Skipped, Is.EqualTo(new[] { "packed", "sa", "lcp" }), "the state of the stages that are done is used with every other value of the parameters that the output does not depend on");

        foreach ((string name, SesaBuildOptions changed, byte[] reference) in new (string, SesaBuildOptions, byte[])[]
        {
            ("without PSI", Opts(3, 1, psi: false), Reference(s_text, psi: false)),
            ("the packed child codec", Opts(3, 1, codec: ChildCodec.Packed), Reference(s_text, codec: ChildCodec.Packed)),
        })
        {
            using var d = new ResumeCase(s_text);
            Fails(d, Opts(3, 1), SesaHooks.LcpChildBuilt);
            SesaBuildResult q = SesaBuilder.Build(d.Input, d.Target, changed);
            AssertSame(d, reference, name);
            Assert.That(q.Resume!.Skipped, Is.Empty, $"{name}: another layout starts clean");
            Assert.That((q.Resume.StateFound, q.Resume.Accepted), Is.EqualTo((true, false)), name);
            Assert.That(q.Resume.Reason, Is.Not.Empty, name);
            AssertOnlyTarget(d, name);
        }
    }

    [Test]
    public void TheRankFile_IsReusedForTheSameCoverPeriod_AndNotForAnother()
    {
        int pa = SesaRun.Tiny(3).Period, pb = SesaRun.Tiny(4).Period;
        Assert.That(pa, Is.Not.EqualTo(pb), "the two variants use different periods");
        foreach ((int v2, bool same) in new[] { (3, true), (4, false) })
        {
            using var c = new ResumeCase(s_text);
            FailsInside(c, Opts(3, 2), "ph23", sa: true);
            SesaBuildResult r = SesaBuilder.Build(c.Input, c.Target, Opts(v2, 2));
            AssertSame(c, Reference(s_text), $"period {SesaRun.Tiny(v2).Period}");
            Assert.That(r.Resume!.Skipped, Is.EqualTo(same ? new[] { "packed", "ranks" } : new[] { "packed" }), $"variant {v2}");
            AssertOnlyTarget(c, "period");
        }
    }

    // ------------------------------------------------------------------ what is not trusted

    public enum Damage
    {
        TextByte,
        SaByte,
        ChildByte,
        LcpByte,
        PsiByte,
        BwtByte,
        RanksByte,
        ContainerShort,
        ContainerLong,
        ContainerGone,
        BwtGone,
        RanksGone,
        ScratchDeleted,
    }

    [Test]
    public void AStageWhoseOutputIsMissingShortOrChanged_IsNotTrusted_AndEverythingAfterItIsDoneAgain(
        [Values(Damage.TextByte, Damage.SaByte, Damage.ChildByte, Damage.LcpByte, Damage.PsiByte, Damage.BwtByte, Damage.RanksByte, Damage.ContainerShort, Damage.ContainerLong, Damage.ContainerGone, Damage.BwtGone, Damage.RanksGone, Damage.ScratchDeleted)] Damage damage)
    {
        string point = damage switch
        {
            Damage.BwtByte or Damage.BwtGone => SesaHooks.BwtBuilt,
            Damage.RanksByte or Damage.RanksGone => "ph23",
            _ => SesaHooks.MagicWritten,
        };
        using var c = new ResumeCase(s_text);
        if (point == "ph23") FailsInside(c, Opts(3, 2), "ph23", sa: true);
        else Fails(c, Opts(3, 2), point);
        SesaPlan plan = PlanOf(c, true);
        string[] expected;
        switch (damage)
        {
            case Damage.TextByte: FlipByte(c.Tmp, plan.TextOffset + (plan.TextLength / 2)); expected = []; break;
            case Damage.SaByte: FlipByte(c.Tmp, plan.SaOffset + (plan.SaLength / 2)); expected = ["packed"]; break;
            case Damage.ChildByte: FlipByte(c.Tmp, plan.ChildOffset + (plan.ChildLength / 2)); expected = ["packed", "sa"]; break;
            case Damage.LcpByte: FlipByte(c.Tmp, plan.LcpOffset + 100); expected = ["packed", "sa"]; break;
            case Damage.PsiByte: FlipByte(c.Tmp, plan.PsiOffset + (plan.PsiLength / 2)); expected = ["packed", "sa", "lcp"]; break;
            case Damage.BwtByte: FlipByte(Path.Combine(c.Scratch, "ph6", "bwt.bin"), 1000); expected = ["packed", "sa", "lcp"]; break;
            case Damage.RanksByte: FlipByte(Path.Combine(c.Scratch, "ranks.bin"), 40); expected = ["packed"]; break;
            case Damage.ContainerShort:
                using (var fs = new FileStream(c.Tmp, FileMode.Open, FileAccess.Write, FileShare.ReadWrite)) fs.SetLength(plan.ReservedLength - 1);
                expected = [];
                break;
            case Damage.ContainerLong:
                using (var fs = new FileStream(c.Tmp, FileMode.Open, FileAccess.Write, FileShare.ReadWrite)) fs.SetLength(new FileInfo(c.Tmp).Length + 4096);
                expected = ["packed", "sa"];
                break;
            case Damage.ContainerGone: File.Delete(c.Tmp); expected = []; break;
            case Damage.BwtGone: File.Delete(Path.Combine(c.Scratch, "ph6", "bwt.bin")); expected = ["packed", "sa", "lcp"]; break;
            case Damage.RanksGone: File.Delete(Path.Combine(c.Scratch, "ranks.bin")); expected = ["packed"]; break;
            default: Directory.Delete(c.Scratch, true); expected = []; break;
        }
        SesaBuildResult r = SesaBuilder.Build(c.Input, c.Target, Opts(3, 2));
        AssertSame(c, Reference(s_text), damage.ToString());
        AssertOnlyTarget(c, damage.ToString());
        Assert.That(r.Resume!.Skipped, Is.EqualTo(expected), damage.ToString());
        if (expected.Length < DoneAt(point, true).Length) Assert.That(r.Resume.Reason, Is.Not.Empty, "the reason is reported");
    }

    [Test]
    public void BytesThatTheUnfinishedStagesLeftBehind_DoNotChangeTheResult_EveryStageWritesItsWholeRange([ValueSource(nameof(HookNames))] string point)
    {
        if (HookIndex(point) < HookIndex(SesaHooks.Packed)) Assert.Ignore("no state: a clean start");
        using var c = new ResumeCase(s_text);
        Fails(c, Opts(3, 2), point);
        SesaPlan plan = PlanOf(c, true);
        string[] done = DoneAt(point, true);
        using (var fs = new FileStream(c.Tmp, FileMode.Open, FileAccess.ReadWrite, FileShare.ReadWrite))
        {
            void Poison(long offset, long length)
            {
                fs.Position = offset;
                fs.Write(Enumerable.Repeat((byte)0xCD, (int)length).ToArray());
            }
            if (!done.Contains("sa")) Poison(plan.SaOffset, plan.SaLength);
            if (!done.Contains("lcp"))
            {
                Poison(plan.ChildOffset, plan.ChildLength);
                Poison(plan.LcpOffset, fs.Length - plan.LcpOffset);
            }
            if (!done.Contains("psi")) Poison(plan.PsiOffset, plan.PsiLength);
        }
        SesaBuildResult r = SesaBuilder.Build(c.Input, c.Target, Opts(3, 2));
        AssertSame(c, Reference(s_text), $"poison after {point}");
        Assert.That(r.Resume!.Skipped, Is.EqualTo(done));
        AssertOnlyTarget(c, point);
    }

    // ------------------------------------------------------------------ the state file

    private static void KillAtPsiBuilt(ResumeCase c) => Fails(c, Opts(3, 2), SesaHooks.PsiBuilt);

    private static void AssertStartsClean(ResumeCase c, bool stateFound, string label)
    {
        SesaBuildResult r = SesaBuilder.Build(c.Input, c.Target, Opts(3, 2));
        AssertSame(c, Reference(s_text), label);
        AssertOnlyTarget(c, label);
        Assert.That(r.Resume!.Skipped, Is.Empty, label);
        Assert.That(r.Resume.Accepted, Is.False, label);
        Assert.That(r.Resume.StateFound, Is.EqualTo(stateFound), label);
        if (stateFound) Assert.That(r.Resume.Reason, Is.Not.Empty, label);
    }

    [Test]
    public void AStateFile_ThatIsMissingEmptyCutDamagedOrNotAFile_StartsClean_AndNeverThrows([Values(0, 1, 2, 3, 4, 5, 6, 7, 8)] int variant)
    {
        using var c = new ResumeCase(s_text);
        KillAtPsiBuilt(c);
        byte[] good = File.ReadAllBytes(c.State);
        bool found = true;
        switch (variant)
        {
            case 0: File.Delete(c.State); found = false; break;
            case 1: File.WriteAllBytes(c.State, []); break;
            case 2: File.WriteAllBytes(c.State, good[..(good.Length / 2)]); break;
            case 3: File.WriteAllBytes(c.State, good[..^2]); break;
            case 4: { byte[] bad = (byte[])good.Clone(); bad[good.Length / 3] ^= 1; File.WriteAllBytes(c.State, bad); break; }
            case 5: { byte[] bad = (byte[])good.Clone(); bad[^4] = (byte)(bad[^4] == '0' ? '1' : '0'); File.WriteAllBytes(c.State, bad); break; }
            case 6: File.Delete(c.State); Directory.CreateDirectory(c.State); break;
            case 7: File.WriteAllBytes(c.State, new byte[5000]); break;
            default: File.WriteAllBytes(c.State, "{\"format\":\"sesa-resume\"}"u8.ToArray()); break;
        }
        AssertStartsClean(c, found, $"variant {variant}");
    }

    [Test]
    public void AStateFile_OfAnotherVersionWriterOrInput_StartsClean([Values(0, 1, 2, 3, 4, 5, 6, 7, 8, 9)] int variant)
    {
        using var c = new ResumeCase(s_text);
        KillAtPsiBuilt(c);
        Assert.That(SesaResumeFormat.TryParse(File.ReadAllBytes(c.State), out SesaResumeState? s), Is.True);
        SesaResumeState edited = variant switch
        {
            0 => s! with { Version = 2 },
            1 => s! with { Writer = "SuffixTree.Esa.External.SesaBuilder 2" },
            2 => s! with { N = s.N + 1 },
            3 => s! with { Sha256 = new string('f', 64) },
            4 => s! with { Psi = false, Stages = [.. s.Stages.Where(x => x.Name is "packed" or "sa" or "lcp")] },
            5 => s! with { ChildCodec = "packed" },
            6 => s! with { ReservedLength = s.ReservedLength + 4096 },
            7 => s! with { Counts = [.. s.Counts.Select((x, i) => i == 0 ? x + 1 : x)] },
            8 => s! with { Letters = [.. s.Letters.Select((x, i) => i == 0 ? x + 1 : x)] },
            _ => s! with { Bits = s.Bits + 1 },
        };
        File.WriteAllBytes(c.State, SesaResumeFormat.Serialize(edited));
        bool parses = SesaResumeFormat.TryParse(File.ReadAllBytes(c.State), out _);
        Assert.That(parses, Is.EqualTo(variant != 0), "the state of another version does not even parse");
        AssertStartsClean(c, true, $"edited {variant}");
    }

    [Test]
    public void AStageRecord_WithAWrongRange_ALongerLengthOrAWrongCrc_IsNotTrusted([Values(0, 1, 2, 3, 4)] int variant)
    {
        using var c = new ResumeCase(s_text);
        KillAtPsiBuilt(c);
        Assert.That(SesaResumeFormat.TryParse(File.ReadAllBytes(c.State), out SesaResumeState? s), Is.True);
        SesaStageRecord Edit(string stage, Func<SesaFileCheck, SesaFileCheck> f) => s!.Find(stage)! with { Checks = [.. s.Find(stage)!.Checks.Select(f)] };
        (SesaStageRecord record, string[] expected) = variant switch
        {
            0 => (Edit("sa", x => x with { Offset = x.Offset + 4096 }), new[] { "packed" }),
            1 => (Edit("sa", x => x with { Length = x.Length + 4 }), new[] { "packed" }),
            2 => (Edit("lcp", x => x with { Crc = x.Crc ^ 1 }), new[] { "packed", "sa" }),
            3 => (Edit("packed", x => x with { Length = x.Length + 1 }), Array.Empty<string>()),
            _ => (Edit("psi", x => x with { Offset = 0 }), new[] { "packed", "sa", "lcp" }),
        };
        File.WriteAllBytes(c.State, SesaResumeFormat.Serialize(s! with { Stages = [.. s.Stages.Select(x => x.Name == record.Name ? record : x)] }));
        SesaBuildResult r = SesaBuilder.Build(c.Input, c.Target, Opts(3, 2));
        AssertSame(c, Reference(s_text), $"record {variant}");
        Assert.That(r.Resume!.Skipped, Is.EqualTo(expected));
        AssertOnlyTarget(c, "record");
    }

    [Test]
    public void ALcpRecordThatDisagreesWithTheLayoutOrWithItself_IsNotTrusted([Values("kind", "overflow", "final_length", "node_count", "internal_count")] string value)
    {
        using var c = new ResumeCase(s_text);
        KillAtPsiBuilt(c);
        Assert.That(SesaResumeFormat.TryParse(File.ReadAllBytes(c.State), out SesaResumeState? s), Is.True);
        SesaStageRecord lcp = s!.Find("lcp")!;
        Assert.That(lcp.Values.Keys, Does.Contain(value));
        var values = new Dictionary<string, long>(lcp.Values) { [value] = lcp.Values[value] + 1 };
        File.WriteAllBytes(c.State, SesaResumeFormat.Serialize(s with { Stages = [.. s.Stages.Select(x => x.Name == "lcp" ? lcp with { Values = values } : x)] }));
        SesaBuildResult r = SesaBuilder.Build(c.Input, c.Target, Opts(3, 2));
        AssertSame(c, Reference(s_text), value);
        AssertOnlyTarget(c, value);
        // the width, the escapes and the length must agree with the layout and with each other, and NODE_COUNT = 1 + I + n + 1; the other statistics (the depth and the rank of the longest repeat) have no
        // cheap check and rest on the CRC of the state file: a state edited with a valid CRC is a forgery, not damage
        Assert.That(r.Resume!.Skipped, Is.EqualTo(new[] { "packed", "sa" }), value);
    }

    [Test]
    public void ATempFileOfTheState_AndForeignFilesInTheScratch_AreRemoved_NotTrusted()
    {
        using var c = new ResumeCase(s_text);
        Fails(c, Opts(3, 2), SesaHooks.LcpChildBuilt);
        File.WriteAllText(Path.Combine(c.Scratch, SesaResumeFormat.TempName), "{ half a state");
        File.WriteAllText(Path.Combine(c.Scratch, "junk.bin"), "x");
        Directory.CreateDirectory(Path.Combine(c.Scratch, "ph4"));
        File.WriteAllText(Path.Combine(c.Scratch, "ph4", "ph4phi-0.tmp"), "left by PH4");
        SesaBuildResult r = SesaBuilder.Build(c.Input, c.Target, Opts(3, 2));
        AssertSame(c, Reference(s_text), "foreign");
        AssertOnlyTarget(c, "foreign");
        Assert.That(r.Resume!.Skipped, Is.EqualTo(new[] { "packed", "sa", "lcp" }));
    }

    [Test]
    public void AGrownContainer_IsCutBackAndGrownAgainOnce_AndAGrowthThatWasCutOffInItsZerosIsReservedWhole([Values("overflow", "u16", "u32")] string kind, [Values(0, 1)] int variant)
    {
        string text = kind switch { "overflow" => TextOverflow, "u16" => TextU16, _ => TextU32 };
        using var c = new ResumeCase(text);
        SesaBuildOptions o = new() { Threads = 2, Fsync = false, MeasureRss = false, Resume = true };
        Fails(c, o, SesaHooks.LcpPlaced);
        PackedTextInfo info = PackedTextWriter.Scan(c.Input);
        SesaPlan plan = SesaLayout.Plan(info, true);
        long grownLength = new FileInfo(c.Tmp).Length;
        Assert.That(grownLength, Is.GreaterThan(plan.ReservedLength), "the failed run had grown the container");
        if (variant == 1)
        {
            // a growth that was killed after the length was set and before its zeros were written: the added range is a hole
            using var fs = new FileStream(c.Tmp, FileMode.Open, FileAccess.Write, FileShare.ReadWrite);
            fs.SetLength(plan.ReservedLength);
            fs.SetLength(grownLength);
        }
        var counters = new IoCounters();
        SesaBuildResult r = SesaBuilder.Build(c.Input, c.Target, o with { Counters = counters });
        AssertSame(c, Reference(text), $"{kind} {variant}");
        Assert.That(r.Resume!.Skipped, Is.EqualTo(new[] { "packed", "sa" }));
        Assert.That(r.GrownBytes, Is.EqualTo(grownLength - plan.ReservedLength), "the growth is done once, by this run");
        Assert.That(r.GrowthReservation, Is.EqualTo(ReservationOutcome.Filled), "and reserved whole, by zeros");
        Assert.That(counters.ReservedWriteBytes, Is.EqualTo(grownLength - plan.ReservedLength), "every byte of the added range was written as a zero: the reservation is not lost");
        AssertOnlyTarget(c, kind);
    }

    [Test]
    public void AContainerThatIsGrownAndWhoseLcpStageIsDone_IsNotGrownAgain([Values("overflow", "u16", "u32")] string kind)
    {
        string text = kind switch { "overflow" => TextOverflow, "u16" => TextU16, _ => TextU32 };
        using var c = new ResumeCase(text);
        SesaBuildOptions o = new() { Threads = 2, Fsync = false, MeasureRss = false, Resume = true };
        Fails(c, o, SesaHooks.LcpChildBuilt);
        long length = new FileInfo(c.Tmp).Length;
        var counters = new IoCounters();
        SesaBuildResult r = SesaBuilder.Build(c.Input, c.Target, o with { Counters = counters });
        AssertSame(c, Reference(text), kind);
        Assert.That(r.Resume!.Skipped, Is.EqualTo(new[] { "packed", "sa", "lcp" }));
        Assert.That((r.GrownBytes, r.GrowthReservation, counters.ReservedWriteBytes), Is.EqualTo((0L, (ReservationOutcome?)null, 0L)), "no second growth");
        Assert.That(r.FileSize, Is.EqualTo(length));
        Assert.That(r.LcpKind != LcpKind.U8 || r.OverflowCount > 0, Is.True, "the numbers of the skipped stage reach the result");
    }

    // ------------------------------------------------------------------ the numbers of a skipped stage reach the result and the header

    [Test]
    public void TheResultOfAResumedBuild_CarriesTheNumbersOfTheSkippedStages_AsAnUninterruptedBuildDoes()
    {
        foreach (string text in new[] { s_text, TextOverflow })
        {
            using var whole = new ResumeCase(text);
            SesaBuildResult a = SesaBuilder.Build(whole.Input, whole.Target, new SesaBuildOptions { Threads = 2, MeasureRss = false, Fsync = false });
            using var c = new ResumeCase(text);
            Fails(c, new SesaBuildOptions { Threads = 2, MeasureRss = false, Fsync = false, Resume = true }, SesaHooks.PsiBuilt);
            SesaBuildResult b = SesaBuilder.Build(c.Input, c.Target, new SesaBuildOptions { Threads = 3, MeasureRss = false, Fsync = false, Resume = true });
            Assert.That(b.Resume!.Skipped, Is.EqualTo(new[] { "packed", "sa", "lcp", "psi" }));
            Assert.That((b.N, b.Sigma, b.Bits, b.FileSize, b.LcpKind, b.OverflowCount, b.NodeCount, b.InternalCount, b.LrsDepth, b.LrsRank, b.MaxLcp, b.HasPsi, b.PrimaryRow),
                Is.EqualTo((a.N, a.Sigma, a.Bits, a.FileSize, a.LcpKind, a.OverflowCount, a.NodeCount, a.InternalCount, a.LrsDepth, a.LrsRank, a.MaxLcp, a.HasPsi, a.PrimaryRow)));
            Assert.That(b.Header.NodeCount, Is.EqualTo(a.Header.NodeCount));
            Assert.That(b.Sections.Select(s => (s.Kind, s.Offset, s.Length, s.Crc)), Is.EqualTo(a.Sections.Select(s => (s.Kind, s.Offset, s.Length, s.Crc))));
            Assert.That(b.CrcRows.Count, Is.EqualTo(a.CrcRows.Count));
            Assert.That(b.Sa.Ph1Seconds + b.Sa.Ph23Seconds, Is.Zero, "the suffix array was not built by this run");
            Assert.That(b.Seconds.Ph45 + b.Seconds.Bwt + b.Seconds.Psi, Is.Zero, "nor the rest");
            Assert.That(b.Seconds.Total, Is.GreaterThan(0));
            AssertSame(c, Read(whole.Target), "numbers");
        }
    }

    // ------------------------------------------------------------------ failures keep what is worth keeping, and only then

    [Test]
    public void AFailureBeforeTheFirstCheckpoint_LeavesNothing_AsABuildWithoutResumeDoes([Values("start", "scanned", "created")] string point)
    {
        using var c = new ResumeCase(s_text);
        Fails(c, Opts(3, 2), point);
        Assert.That(c.Entries, Is.Empty);
    }

    [Test]
    public void AFailedResumeBuild_LeavesAnExistingTargetAsItWas_AndAPlainBuildAfterItCleansUp()
    {
        using var c = new ResumeCase(s_text);
        File.WriteAllBytes(c.Target, [1, 2, 3, 4]);
        Fails(c, Opts(3, 2), SesaHooks.SaBuilt);
        Assert.That(File.ReadAllBytes(c.Target), Is.EqualTo(new byte[] { 1, 2, 3, 4 }), "the target of an earlier build is not touched by a build that fails");
        Assert.That(c.Entries, Is.EqualTo(new[] { "x.sesa", "x.sesa.tmp", "x.sesa.tmp.d" }));
        SesaBuilder.Build(c.Input, c.Target, new SesaBuildOptions { Threads = 2, MeasureRss = false, Fsync = false });
        AssertSame(c, Reference(s_text), "a plain build over the state of a resumable one");
        AssertOnlyTarget(c, "plain after resume");
    }

    [Test]
    public void ABuildThatFinished_AndALeftoverScratchAfterTheRename_StartCleanWithResume()
    {
        using var c = new ResumeCase(s_text);
        Fails(c, Opts(3, 2), SesaHooks.MagicWritten);
        File.Move(c.Tmp, c.Target); // the kill came after the rename and before the scratch was removed
        Assert.That(Directory.Exists(c.Scratch), Is.True);
        SesaBuildResult r = SesaBuilder.Build(c.Input, c.Target, Opts(3, 2));
        AssertSame(c, Reference(s_text), "after the rename");
        AssertOnlyTarget(c, "after the rename");
        Assert.That(r.Resume!.Skipped, Is.Empty);
    }

    [Test]
    public void AResumeBuildOfAnInputThatChanged_StartsCleanAndTheChangeIsInTheContainer()
    {
        string changed = s_text[..1000] + "T" + s_text[1001..];
        Assert.That(changed, Is.Not.EqualTo(s_text));
        using var c = new ResumeCase(s_text);
        Fails(c, Opts(3, 2), SesaHooks.LcpChildBuilt);
        File.WriteAllBytes(c.Input, System.Text.Encoding.Latin1.GetBytes(changed));
        SesaBuildResult r = SesaBuilder.Build(c.Input, c.Target, Opts(3, 2));
        AssertSame(c, Reference(changed), "input changed, same length");
        Assert.That(r.Resume!.Skipped, Is.Empty);
        Assert.That(r.Resume.StateFound, Is.True);
        string longer = changed + "ACGT";
        Fails(c, Opts(3, 2), SesaHooks.LcpChildBuilt);
        File.WriteAllBytes(c.Input, System.Text.Encoding.Latin1.GetBytes(longer));
        File.Delete(c.Target);
        r = SesaBuilder.Build(c.Input, c.Target, Opts(3, 2));
        AssertSame(c, Reference(longer), "input longer");
        Assert.That(r.Resume!.Skipped, Is.Empty);
        AssertOnlyTarget(c, "input changed");
    }

    // ------------------------------------------------------------------ two builds on one target

    [Test]
    public void TwoResumableBuildsOfOneTarget_TheSecondFailsFast_AndTouchesNothing_TheFirstIsUnharmed()
    {
        using var c = new ResumeCase(s_text);
        using var reached = new ManualResetEventSlim();
        using var release = new ManualResetEventSlim();
        SesaBuildResult? first = null;
        Exception? firstError = null;
        var thread = new Thread(() =>
        {
            try
            {
                first = SesaBuilder.Build(c.Input, c.Target, Opts(3, 2) with
                {
                    Hook = p =>
                    {
                        if (p != SesaHooks.SaBuilt) return;
                        reached.Set();
                        release.Wait();
                    },
                });
            }
#pragma warning disable CA1031 // the test thread reports whatever happened
            catch (Exception e)
            {
                firstError = e;
            }
#pragma warning restore CA1031
        });
        thread.Start();
        Assert.That(reached.Wait(TimeSpan.FromSeconds(60)), Is.True);
        string[] before = c.Entries;
        long tmpLength = new FileInfo(c.Tmp).Length;
        IOException? e = Assert.Throws<IOException>(() => SesaBuilder.Build(c.Input, c.Target, Opts(5, 1)));
        Assert.That(e!.Message, Does.Contain("lock"));
        Assert.That(c.Entries, Is.EqualTo(before), "the refused build changed nothing: it did not remove the temporary file or the scratch of the running one");
        Assert.That(new FileInfo(c.Tmp).Length, Is.EqualTo(tmpLength));
        Assert.That(File.Exists(c.Lock), Is.True);
        release.Set();
        thread.Join();
        Assert.That(firstError, Is.Null);
        AssertSame(c, Reference(s_text), "the first build");
        AssertOnlyTarget(c, "after both");
        SesaBuildResult again = SesaBuilder.Build(c.Input, c.Target, Opts(5, 1));
        AssertSame(c, Reference(s_text), "the second build, afterwards");
        Assert.That(first, Is.Not.Null);
        Assert.That(again.Resume!.Skipped, Is.Empty);
    }

    [Test]
    public void AStaleLockFile_FromAKilledBuild_DoesNotBlock_AndIsRemoved()
    {
        using var c = new ResumeCase(s_text);
        File.WriteAllText(c.Lock, "stale");
        SesaBuilder.Build(c.Input, c.Target, Opts(3, 2));
        AssertSame(c, Reference(s_text), "stale lock");
        AssertOnlyTarget(c, "stale lock");
    }

    [Test]
    public void ABuildWithoutResume_TakesNoLock_SoItDoesNotCollideWithTheHolderOfOne()
    {
        using var c = new ResumeCase(s_text);
        var seen = new List<bool>();
        SesaBuilder.Build(c.Input, c.Target, new SesaBuildOptions { Threads = 2, MeasureRss = false, Fsync = false, Hook = p => seen.Add(File.Exists(c.Lock)) });
        Assert.That(seen.Any(x => x), Is.False);
        using SesaBuildLock held = SesaBuildLock.Acquire(c.Lock);
        SesaBuilder.Build(c.Input, c.Target, new SesaBuildOptions { Threads = 2, MeasureRss = false, Fsync = false });
        Assert.That(File.Exists(c.Lock), Is.True, "the holder's file is still there");
    }

    // ------------------------------------------------------------------ the free-space pre-flight of a resume

    [Test]
    public void TheResumeNeed_IsWhatTheRemainingStagesWillWrite_NotTheWholeBuild()
    {
        const long n = 1000;
        long full = SesaBuilder.ScratchBytesPerChar * n;
        Assert.That(SesaBuilder.ResumeNeedBytes(n, true, [], 0), Is.EqualTo(full));
        Assert.That(SesaBuilder.ResumeNeedBytes(n, true, ["packed"], 0), Is.EqualTo(full));
        Assert.That(SesaBuilder.ResumeNeedBytes(n, true, ["packed", "ranks"], 170), Is.EqualTo(full - 170), "the rank file is already on the disk");
        Assert.That(SesaBuilder.ResumeNeedBytes(n, true, ["packed", "sa"], 0), Is.EqualTo(full));
        Assert.That(SesaBuilder.ResumeNeedBytes(n, true, ["packed", "sa", "lcp"], 0), Is.EqualTo(n + 1), "the BWT file is the only thing left to create");
        Assert.That(SesaBuilder.ResumeNeedBytes(n, true, ["packed", "sa", "lcp", "bwt"], n + 1), Is.Zero);
        Assert.That(SesaBuilder.ResumeNeedBytes(n, true, ["packed", "sa", "lcp", "psi"], 0), Is.Zero);
        Assert.That(SesaBuilder.ResumeNeedBytes(n, false, ["packed", "sa", "lcp"], 0), Is.Zero, "without PSI there is nothing after the LCP");
        Assert.That(SesaBuilder.ResumeNeedBytes(n, true, ["packed", "ranks"], 13 * n), Is.Zero, "never negative");
    }

    [Test]
    public void AResumeNeedsLessFreeSpaceThanAFreshBuild_AndRefusesBelowItsOwnNeed_LosingNothing()
    {
        using var c = new ResumeCase(s_text);
        Fails(c, Opts(3, 2), SesaHooks.LcpChildBuilt);
        long n = s_text.Length;
        long need = SesaBuilder.ResumeNeedBytes(n, true, ["packed", "sa", "lcp"], 0);
        SesaPlan plan = PlanOf(c, true);
        Assert.That(need, Is.LessThan(plan.ReservedLength + (SesaBuilder.ScratchBytesPerChar * n)), "a resume at this stage needs less than a fresh build");
        string[] before = c.Entries;
        IOException? e = Assert.Throws<IOException>(() => SesaBuilder.Build(c.Input, c.Target, Opts(3, 2) with { Files = new MappedFileOptions { Source = new FixedFreeSource(need - 1) } }));
        Assert.That(e!.Message, Does.Contain("not enough free space"));
        Assert.That(c.Entries, Is.EqualTo(before), "the pre-flight refused before it deleted anything: the work of the earlier run is still there");
        SesaBuildResult r = SesaBuilder.Build(c.Input, c.Target, Opts(3, 2) with { Files = new MappedFileOptions { Source = new FixedFreeSource(need) } });
        AssertSame(c, Reference(s_text), "exactly the need");
        Assert.That(r.Resume!.Skipped, Is.EqualTo(new[] { "packed", "sa", "lcp" }));
        AssertOnlyTarget(c, "need");
    }

    [Test]
    public void AFreshBuildWithResume_KeepsThePreFlightOfTheBuildWithout()
    {
        using var c = new ResumeCase(s_text);
        PackedTextInfo info = PackedTextWriter.Scan(c.Input);
        long need = SesaLayout.Plan(info, true).ReservedLength + (SesaBuilder.ScratchBytesPerChar * info.N);
        Assert.Throws<IOException>(() => SesaBuilder.Build(c.Input, c.Target, Opts(3, 2) with { Files = new MappedFileOptions { Source = new FixedFreeSource(need - 1) } }));
        Assert.That(c.Entries, Is.Empty);
        SesaBuilder.Build(c.Input, c.Target, Opts(3, 2) with { Files = new MappedFileOptions { Source = new FixedFreeSource(need) } });
        AssertSame(c, Reference(s_text), "fresh at the need");
    }

    [Test]
    public void TheScratchOfAResumedBuild_IsOnTheLedger_AndBackAtZero([Values("ph23", SesaHooks.LcpChildBuilt)] string point)
    {
        using var c = new ResumeCase(s_text);
        if (point == "ph23") FailsInside(c, Opts(3, 2), "ph23", sa: true);
        else Fails(c, Opts(3, 2), point);
        var counters = new IoCounters();
        SesaBuildResult r = SesaBuilder.Build(c.Input, c.Target, Opts(3, 2) with
        {
            Counters = counters,
            Hook = p =>
            {
                if (p == SesaHooks.BwtBuilt) Assert.That(counters.ScratchLiveBytes, Is.GreaterThanOrEqualTo(s_text.Length + 1), "the BWT file of the resumed run is on the account");
            },
        });
        Assert.That(r.Resume!.Skipped, Is.EqualTo(point == "ph23" ? new[] { "packed", "ranks" } : new[] { "packed", "sa", "lcp" }));
        Assert.That(counters.ScratchLiveBytes, Is.Zero);
        Assert.That(counters.ScratchPeakBytes, Is.GreaterThan(0));
        if (point == "ph23") Assert.That(r.Sa.RankFileBytes, Is.GreaterThan(0), "the rank file of the earlier run was counted again");
    }

    [Test]
    public void ABwtFileOfAnEarlierRun_IsOnTheLedgerOfTheResumedRun_AndItsCountsAreTheTextsCounts()
    {
        using var c = new ResumeCase(s_text);
        Fails(c, Opts(3, 2), SesaHooks.BwtBuilt);
        var counters = new IoCounters();
        long liveAtPsi = -1;
        SesaBuildResult r = SesaBuilder.Build(c.Input, c.Target, Opts(3, 2) with { Counters = counters, Hook = p => { if (p == SesaHooks.PsiBuilt) liveAtPsi = counters.ScratchLiveBytes; } });
        Assert.That(r.Resume!.Skipped, Is.EqualTo(new[] { "packed", "sa", "lcp", "bwt" }));
        Assert.That(r.Bwt, Is.Not.Null);
        Assert.That(r.Bwt!.PrimaryRow, Is.EqualTo(r.PrimaryRow));
        Assert.That(r.Bwt.Counts, Is.EqualTo(PackedTextWriter.Scan(c.Input).Counts));
        Assert.That(counters.ScratchLiveBytes, Is.Zero);
        Assert.That(liveAtPsi, Is.EqualTo(s_text.Length + 1), "at the boundary after PSI the BWT file is still on the account (it is removed with the directory)");
        AssertSame(c, Reference(s_text), "bwt");
    }

    // ------------------------------------------------------------------ determinism

    [Test]
    public void AResumeAtFourThreads_GivesTheBytesOfAnUninterruptedBuildAtFourThreadsAndAtOne([Values(SesaHooks.Packed, SesaHooks.SaBuilt, SesaHooks.LcpChildBuilt, SesaHooks.PsiBuilt)] string point)
    {
        byte[][] plain = new byte[2][];
        foreach ((int i, int threads) in new[] { (0, 1), (1, 4) })
        {
            using var w = new ResumeCase(s_text);
            SesaBuilder.Build(w.Input, w.Target, SesaRun.Tiny(3, threads) with { MeasureRss = false });
            plain[i] = Read(w.Target);
        }
        Assert.That(plain[0], Is.EqualTo(plain[1]), "the build is deterministic without Resume");
        using var c = new ResumeCase(s_text);
        Fails(c, Opts(3, 4), point);
        SesaBuildResult r = SesaBuilder.Build(c.Input, c.Target, Opts(3, 4));
        Assert.That(r.Resume!.Skipped, Is.EqualTo(DoneAt(point, true)));
        Assert.That(Read(c.Target), Is.EqualTo(plain[1]), "resumed at four threads: the bytes of the build at four threads");
        Assert.That(Read(c.Target), Is.EqualTo(plain[0]), "and of the build at one thread");
    }
}
