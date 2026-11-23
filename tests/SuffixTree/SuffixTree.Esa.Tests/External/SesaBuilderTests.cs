using System.Security.Cryptography;
using System.Text;
using SuffixTree.Esa.External;
using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.Esa.Tests.External;

/// <summary>The container writer end to end (gate G1 of work package 1b.7e): every section and the header against the in-RAM index and against ADR-003, written by the independent reader <see cref="SesaOracle"/>.</summary>
[TestFixture]
[HangLimit(30)]
public sealed class SesaBuilderTests
{
    private static readonly string[] Samples =
    [
        "", "a", "ab", "aa", "ba", "abc", "mississippi", "abracadabra", "banana", "aaaaaaaaaaaaaaaaaaaa", "ACGTACGTTGCAACGT", "abababababab", "ZZZZ", "xéÿéx", "the quick brown fox jumps over the lazy dog",
    ];

    [TestCaseSource(nameof(Samples))]
    public void EverySection_EqualsTheInMemoryIndex_WithAndWithoutPsi_DefaultAndTiny(string text)
    {
        using var dir = new TempDir();
        int v = 0;
        foreach (bool psi in new[] { true, false })
        {
            SesaRun.Outcome o = SesaRun.Build(dir, text, new SesaBuildOptions { IncludePsi = psi, Fsync = false }, $"-d{psi}");
            SesaRun.Compare(o, $"'{text}' default psi={psi}", psi);
            v++;
            o = SesaRun.Build(dir, text, SesaRun.Tiny(v + 5, psi: psi), $"-t{psi}");
            SesaRun.Compare(o, $"'{text}' tiny psi={psi}", psi);
        }
    }

    [Test]
    [HangLimit(400)]
    public void EveryStringOfSmallAlphabets_GivesTheInMemoryIndex()
    {
        using var dir = new TempDir();
        int v = 0;
        foreach ((int sigma, int maxLength) in new[] { (1, 5), (2, 6), (3, 4), (4, 3) })
        {
            foreach (byte[] s in TextOracle.AllTexts(sigma, maxLength))
            {
                v++;
                string text = new([.. s.Select(b => (char)('a' + b))]);
                SesaRun.Outcome o = SesaRun.Build(dir, text, SesaRun.Tiny(v, 1 + (v % 3), psi: v % 4 != 0), "-x");
                SesaRun.Compare(o, $"'{text}' variant {v}", v % 4 != 0);
                Directory.Delete(o.OutDir, true);
            }
        }
        Assert.That(v, Is.GreaterThan(300));
    }

    [Test]
    [HangLimit(120)]
    public void TheFamilies_GiveTheInMemoryIndex_AtTheDefaultAndATinyGeometry()
    {
        using var dir = new TempDir();
        int v = 0;
        foreach (int n in new[] { 600, 3000 })
        {
            foreach ((string name, string text) in LcpChildRun.Families(n))
            {
                v++;
                SesaRun.Outcome o = SesaRun.Build(dir, text, new SesaBuildOptions { Fsync = false }, "-a");
                SesaRun.Compare(o, $"{name} {n} defaults", true);
                Directory.Delete(o.OutDir, true);
                o = SesaRun.Build(dir, text, SesaRun.Tiny(v, 2), "-b");
                SesaRun.Compare(o, $"{name} {n} tiny {v}", true);
                Directory.Delete(o.OutDir, true);
            }
        }
    }

    [Test]
    public void TheLcpWidths_U8WithoutOverflow_U8WithOverflow_AndU16_AreWrittenAsTheirSections()
    {
        using var dir = new TempDir();
        string random = StringFamilies.RandomDna(3000, 1);
        string longRun = StringFamilies.RandomDna(2800, 2) + new string('A', 400) + StringFamilies.RandomDna(2800, 3);
        string u16 = new string('A', 500) + StringFamilies.RandomDna(100, 4);
        (string Label, string Text, LcpKind Kind, bool Overflow)[] cases =
        [
            ("U8, no escapes", random, LcpKind.U8, false),
            ("U8 with an overflow list", longRun, LcpKind.U8, true),
            ("U16 without escapes", u16, LcpKind.U16, false),
        ];
        foreach ((string label, string text, LcpKind kind, bool overflow) in cases)
        {
            SesaRun.Outcome o = SesaRun.Build(dir, text, new SesaBuildOptions { Fsync = false }, "-w");
            SesaRun.Compare(o, label, true);
            Assert.That(o.Result.LcpKind, Is.EqualTo(kind), label);
            Assert.That(o.Result.OverflowCount > 0, Is.EqualTo(overflow), label);
            Assert.That(o.Oracle.Find(5) is not null, Is.EqualTo(overflow), $"{label}: LCP_OVERFLOW iff PARAM0 > 0");
            Directory.Delete(o.OutDir, true);
        }
    }

    [Test]
    public void TheBuildInfo_IsTheFixedTextOfTheOutput_WithNoTimestampAndNoParameterOfTheRun()
    {
        using var dir = new TempDir();
        SesaRun.Outcome o = SesaRun.Build(dir, "mississippi", new SesaBuildOptions { Fsync = false });
        string expected = "format=SESA 1.0\nwriter=SuffixTree.Esa.External.SesaBuilder 1\nalgorithm=B-DCS (ADR-004)\nn=11\nsigma=4\ntext=codec 4, 2 bits per symbol\nlcp=U8, 0 overflow entries\nchild=layout A, tagged\npsi=yes\n";
        Assert.That(Encoding.UTF8.GetString(o.Oracle.Payload(11)), Is.EqualTo(expected));
        o = SesaRun.Build(dir, "mississippi", SesaRun.Tiny(3, 3, false), "-2");
        Assert.That(Encoding.UTF8.GetString(o.Oracle.Payload(11)), Is.EqualTo(expected.Replace("psi=yes", "psi=no", StringComparison.Ordinal)), "the threads, the sizes and the seed are not recorded");
        o = SesaRun.Build(dir, "mississippi", new SesaBuildOptions { Fsync = false, ChildCodec = ChildCodec.Packed }, "-3");
        Assert.That(Encoding.UTF8.GetString(o.Oracle.Payload(11)), Is.EqualTo(expected.Replace("tagged", "packed", StringComparison.Ordinal)), "the layout of the child table is a fact about the output");
    }

    [Test]
    [HangLimit(490)]
    public void TheBytes_DoNotDependOnThreadsSizesSeedsOrArenaCaps()
    {
        using var dir = new TempDir();
        foreach ((string name, string text) in LcpChildRun.Families(1500))
        {
            string reference = "";
            for (int v = 0; v < 8; v++)
            {
                SesaBuildOptions options = SesaRun.Tiny(v, new[] { 1, 2, 4 }[v % 3]);
                if (v == 7) options = new SesaBuildOptions { Threads = 4, Fsync = false, ArenaCapBytes = 96L << 20, ChildArenaCapBytes = 32L << 20 };
                SesaRun.Outcome o = SesaRun.Build(dir, text, options, $"-{v}");
                string sha = SesaRun.Sha(o.Target);
                if (reference.Length == 0) reference = sha;
                Assert.That(sha, Is.EqualTo(reference), $"{name}: variant {v} writes other bytes");
            }
        }
    }

    [Test]
    public void TwoBuildsOfTheSameInput_AtDifferentTimes_AreByteIdentical_WithAndWithoutPsi()
    {
        using var dir = new TempDir();
        foreach (bool psi in new[] { true, false })
        {
            SesaBuildOptions options = SesaRun.Tiny(4, 2, psi);
            string a = SesaRun.Sha(SesaRun.Build(dir, "abracadabra abracadabra", options, "-a").Target);
            SpinWait.SpinUntil(() => false, 1100); // a second apart: nothing of the clock may be in the file
            string b = SesaRun.Sha(SesaRun.Build(dir, "abracadabra abracadabra", options, "-b").Target);
            Assert.That(b, Is.EqualTo(a), $"psi={psi}");
        }
    }

    [Test]
    public void TheHooks_AreCalledInOrder_AndASkippedPsiSkipsTwoOfThem()
    {
        using var dir = new TempDir();
        foreach (bool psi in new[] { true, false })
        {
            var seen = new List<string>();
            SesaBuildOptions options = SesaRun.Tiny(2, 1, psi) with { Hook = seen.Add };
            SesaRun.Build(dir, "mississippi", options, $"-{psi}");
            string[] expected = psi ? SesaHooks.All : [.. SesaHooks.All.Where(h => h is not (SesaHooks.BwtBuilt or SesaHooks.PsiBuilt))];
            Assert.That(seen, Is.EqualTo(expected), $"psi={psi}");
        }
    }

    [Test]
    public void TheResult_CarriesTheCountersThePhasesAndThePeaks()
    {
        using var dir = new TempDir();
        string text = StringFamilies.RandomDna(5000, 7);
        var counters = new IoCounters();
        SesaRun.Outcome o = SesaRun.Build(dir, text, SesaRun.Tiny(1, 2) with { Counters = counters });
        SesaBuildResult r = o.Result;
        Assert.That(r.Io, Is.EqualTo(counters.Snapshot()), "the counters passed in are the ones of the build");
        Assert.That(r.Io.ReadBytes, Is.GreaterThan(5000));
        Assert.That(r.Io.MappedWriteBytes, Is.GreaterThan(4L * 5000 * 2), "SA, CHILD and PSI go through windows");
        Assert.That(r.Io.ScratchLiveBytes, Is.Zero);
        Assert.That(r.Io.ScratchPeakBytes, Is.GreaterThan(5000));
        Assert.That(r.ScratchPeakBytes, Is.EqualTo(r.Io.ScratchPeakBytes));
        Assert.That(r.Seconds.Total, Is.GreaterThan(0));
        Assert.That(new[] { r.Seconds.Scan, r.Seconds.Create, r.Seconds.Pack, r.Seconds.Ph1, r.Seconds.Ph23, r.Seconds.Ph45, r.Seconds.Bwt, r.Seconds.Psi, r.Seconds.Crc, r.Seconds.Seal }.Sum(), Is.LessThanOrEqualTo(r.Seconds.Total + 1e-6));
        Assert.That(r.Seconds.Ph23, Is.GreaterThan(0));
        Assert.That(r.ArenaPeakBytes, Is.GreaterThan(0));
        Assert.That(r.Sa.Ph23ArenaPeakBytes, Is.GreaterThan(0));
        Assert.That(r.Lcp.ArenaPeakBytes, Is.GreaterThan(0));
        Assert.That((r.NodeCount, r.InternalCount, r.LrsDepth, r.LrsRank, r.MaxLcp), Is.EqualTo((o.Case.Statistics.NodeCount, (long)o.Case.Statistics.InternalCount, (long)o.Case.Statistics.LrsDepth, (long)o.Case.Statistics.LrsRank, (long)o.Case.MaxLcp)));
        Assert.That(r.PrimaryRow, Is.EqualTo(SesaRun.PrimaryRow(o.Case.Sa)));
        Assert.That(r.Psi!.PsiBytesWritten, Is.EqualTo(4L * 5001));
        Assert.That(r.Bwt!.PrimaryRow, Is.EqualTo(r.PrimaryRow));
        Assert.That(r.CrcRows.Select(c => c.Kind).OrderBy(k => k), Is.EqualTo(r.Sections.Select(s => s.Kind).OrderBy(k => k)), "one checksum row per section");
        Assert.That(r.CrcRows.Sum(c => c.Bytes), Is.EqualTo(r.Sections.Sum(s => s.Length)));
        Assert.That(r.Reservation, Is.EqualTo(ReservationOutcome.Preallocated).Or.EqualTo(ReservationOutcome.FreeSpaceChecked).Or.EqualTo(ReservationOutcome.Unchecked));
        Assert.That(r.GrowthReservation, Is.Null, "an LCP of one byte per entry without escapes does not grow the file");
        Assert.That(r.GrownBytes, Is.Zero);
        Assert.That(r.TargetPath, Is.EqualTo(o.Target));
    }

    [Test]
    public void AnOverflowList_GrowsTheFileBehindTheReservedLength_AndTheResultSaysSo()
    {
        using var dir = new TempDir();
        string text = StringFamilies.RandomDna(2800, 2) + new string('A', 400) + StringFamilies.RandomDna(2800, 3);
        SesaRun.Outcome o = SesaRun.Build(dir, text, SesaRun.Tiny(2));
        SesaRun.Compare(o, "overflow", true);
        SesaPlan plan = SesaLayout.Plan(new PackedTextInfo(text.Length, 2, [0, 1, 2, 3], [1, 1, 1, 1], new string('0', 64)), true);
        Assert.That(o.Result.OverflowCount, Is.GreaterThan(0));
        Assert.That(o.Result.GrownBytes, Is.EqualTo(o.Result.FileSize - plan.ReservedLength));
        Assert.That(o.Result.GrownBytes, Is.GreaterThan(8L * o.Result.OverflowCount - 64));
        Assert.That(o.Result.GrowthReservation, Is.EqualTo(ReservationOutcome.Filled), "the range that was added has its blocks (work package 1b.9d)");
    }

    [Test]
    public void OnAPlatformWithoutProcAndWithoutFallocate_TheBuildRuns_TheFiguresAreNull_AndTheVerdictIsUnavailable()
    {
        using var dir = new TempDir();
        foreach (long free in new long[] { -1, 1L << 40 })
        {
            foreach (IProcSource? proc in new IProcSource?[] { new FakeProcSource(), new FakeProcSource { Throws = new IOException("hostile") }, new FakeProcSource { Throws = new UnauthorizedAccessException("hostile") } })
            {
                var source = new FakeReservationSource { CanPreallocate = false, Free = free };
                SesaBuildOptions options = SesaRun.Tiny(3) with { Files = new MappedFileOptions { Source = source }, ProcSource = proc };
                SesaRun.Outcome o = SesaRun.Build(dir, "mississippi river", options, $"-{free}-{proc!.GetHashCode()}");
                SesaRun.Compare(o, $"free {free}", true);
                Assert.That(o.Result.RssVerdict, Is.EqualTo(RssVerdict.Unavailable));
                Assert.That(o.Result.RssAnonPeakMb, Is.Null);
                Assert.That((o.Result.Rss.Ph0, o.Result.Rss.Ph1, o.Result.Rss.Ph23, o.Result.Rss.Ph45, o.Result.Rss.Ph6, o.Result.Rss.Seal), Is.EqualTo(((double?)null, (double?)null, (double?)null, (double?)null, (double?)null, (double?)null)));
                Assert.That(o.Result.Reservation, Is.EqualTo(free < 0 ? ReservationOutcome.Unchecked : ReservationOutcome.FreeSpaceChecked));
                Assert.That(source.PreallocateCalls, Is.Zero);
                Assert.That(source.FreeCalls, Is.GreaterThan(0));
                Directory.Delete(o.OutDir, true);
            }
        }
    }

    [Test]
    public void OnLinux_EveryPhaseHasAFigure_AndTheVerdictFollowsTheFigure()
    {
        if (!OperatingSystem.IsLinux() || !File.Exists("/proc/self/status")) Assert.Ignore("no /proc");
        using var dir = new TempDir();
        SesaRun.Outcome o = SesaRun.Build(dir, StringFamilies.RandomDna(20000, 3), new SesaBuildOptions { Fsync = false });
        Assert.That(o.Result.RssAnonPeakMb, Is.GreaterThan(0));
        Assert.That(o.Result.RssVerdict, Is.EqualTo(o.Result.RssAnonPeakMb <= 320 ? RssVerdict.Within : RssVerdict.Exceeded), "the verdict is the figure against the default cap of 320 MB (the test host holds what the tests before this one left, so the figure itself is not asserted to be below the cap here: G3 measures it in a process of its own)");
        Assert.That(SesaRun.Build(dir, StringFamilies.RandomDna(2000, 3), new SesaBuildOptions { Fsync = false, RssAnonCapMb = 1e6 }, "-big").Result.RssVerdict, Is.EqualTo(RssVerdict.Within), "a cap above any process is always met");
        SesaPhaseRss r = o.Result.Rss;
        Assert.That(new[] { r.Ph0, r.Ph1, r.Ph23, r.Ph45, r.Ph6, r.Seal }.All(x => x is > 0), Is.True);
        Assert.That(o.Result.RssAnonPeakMb, Is.GreaterThanOrEqualTo(new[] { r.Ph0, r.Ph1, r.Ph23, r.Ph45, r.Ph6, r.Seal }.Max()!.Value - 5.0), "the whole-build sampler and the phase samplers read every 10 ms at different instants: a process that allocates (the test host) lets them differ by the change of one interval");
        o = SesaRun.Build(dir, StringFamilies.RandomDna(2000, 3), new SesaBuildOptions { Fsync = false, RssAnonCapMb = 0.001 }, "-cap");
        Assert.That(o.Result.RssVerdict, Is.EqualTo(RssVerdict.Exceeded), "the gate is a verdict, never an exception");
        o = SesaRun.Build(dir, StringFamilies.RandomDna(2000, 3), new SesaBuildOptions { Fsync = false, MeasureRss = false }, "-off");
        Assert.That((o.Result.RssVerdict, o.Result.RssAnonPeakMb), Is.EqualTo((RssVerdict.Unavailable, (double?)null)));
    }

    [Test]
    public void TheDefaults_AreTheDecisionsOfTheOwner()
    {
        var o = new SesaBuildOptions();
        Assert.That((o.Threads, o.ArenaCapBytes, o.ChildArenaCapBytes, o.Ph6CapBytes), Is.EqualTo((4, 160L << 20, 64L << 20, 160L << 20)));
        Assert.That((o.RssAnonCapMb, o.Fsync, o.MeasureRss, o.Period, o.ScaleN), Is.EqualTo((320.0, true, true, 1024, 0L)));
        Assert.That(o.TextWindowBytes, Is.GreaterThanOrEqualTo(1 << 20));
        Assert.That(o.TextBlockSymbols, Is.EqualTo(1 << 20));
        Assert.That(o.Files, Is.Null, "null = the Auto reservation");
    }

    [Test]
    public void ASymbolAboveU00FF_IsNotInTheInput_AndABinaryInputOf256ByteValuesIsBuilt()
    {
        using var dir = new TempDir();
        string input = dir.File("bytes.bin");
        byte[] bytes = [.. Enumerable.Range(0, 256).Select(i => (byte)i), .. Enumerable.Range(0, 256).Select(i => (byte)(255 - i)), 0, 0, 255];
        File.WriteAllBytes(input, bytes);
        string target = dir.File("bytes.sesa");
        SesaBuildResult r = SesaBuilder.Build(input, target, SesaRun.Tiny(2));
        SesaOracle f = SesaOracle.Read(target);
        f.AssertConforms("256 symbols", true);
        Assert.That((r.Sigma, r.Bits), Is.EqualTo((256, 8)));
        Assert.That(f.TextSha256, Is.EqualTo(SHA256.HashData(bytes)));
        Assert.That(MemoryMarshal_Cast(f.Payload(2)), Is.EqualTo(Enumerable.Range(0, 256).Select(i => (ushort)i).ToArray()));
    }

    private static ushort[] MemoryMarshal_Cast(byte[] bytes) => System.Runtime.InteropServices.MemoryMarshal.Cast<byte, ushort>(bytes).ToArray();


    [Test]
    public void EveryOption_ReachesItsStage()
    {
        using var dir = new TempDir();
        string text = StringFamilies.RandomDna(4000, 5);
        // the sizes of the stages
        SesaRun.Outcome o = SesaRun.Build(dir, text, new SesaBuildOptions
        {
            Threads = 3,
            Fsync = false,
            Period = 12,
            ExactPeriod = true,
            ScaleN = 3_000_000_000,
            RuleN = 1_500_000_000,
            SaTuning = new SaBuildTuning { BucketCapacity = 777, DcsBucketCapacity = 555 },
            LcpTuning = new LcpChildTuning { LcpChunk = 321, WindowBytes = 4096 },
            Ph6Tuning = new SesaPh6Tuning { SymbolsPerPass = 3, WindowBytes = 40 },
        }, "-a");
        SesaBuildResult r = o.Result;
        Assert.That((r.Sa.Period, r.Sa.CoverKind), Is.EqualTo((12, CoverKind.Pruned)), "ExactPeriod uses the period itself");
        Assert.That((r.Sa.Parameters.Scaled, r.Sa.Parameters.Factor), Is.EqualTo((true, 0.5)), "ScaleN and RuleN reach PH1-PH3");
        Assert.That((r.Lcp.Parameters.Scaled, r.Lcp.Parameters.Factor), Is.EqualTo((true, 0.5)), "ScaleN and RuleN reach PH4-PH5");
        Assert.That((r.Sa.Parameters.BucketCapacity, r.Sa.Parameters.DcsBucketCapacity), Is.EqualTo((777L, 555L)), "the tuning of PH1-PH3");
        Assert.That((r.Lcp.Parameters.LcpChunk, r.Lcp.Parameters.WindowBytes), Is.EqualTo((321L, 4096)), "the tuning of PH4-PH5");
        Assert.That(r.Lcp.Parameters.LcpWorkers, Is.EqualTo(2), "three threads, two workers of PH4 under the cap");
        Assert.That((r.Psi!.Passes, r.Psi.SymbolsPerPass, r.Psi.MappedBoundBytes), Is.EqualTo((2, 3, 3 * 10 * 4L)), "the tuning of PH6");
        SesaRun.Compare(o, "every option", true);
        o = SesaRun.Build(dir, text, new SesaBuildOptions { Fsync = false, Period = 12, Threads = 1, SaTuning = new SaBuildTuning { BucketCapacity = 777 } }, "-b");
        Assert.That((o.Result.Sa.Period, o.Result.Sa.CoverKind), Is.EqualTo((13, CoverKind.Singer)), "without ExactPeriod the nearest Singer period is used");
        Assert.That(o.Result.Lcp.Parameters.LcpWorkers, Is.EqualTo(1));
        Assert.That(o.Result.Sa.Parameters.Scaled, Is.False);
        // the options that only show in what the stages do
        var steps = new List<string>();
        SesaRun.Build(dir, text, new SesaBuildOptions { Fsync = false, VerifyRanks = true, SaTrace = steps.Add, LcpTrace = steps.Add }, "-c");
        Assert.That(steps, Does.Contain("verify").And.Contain("ph1").And.Contain("ph23").And.Contain("b1").And.Contain("b3"), "VerifyRanks and the trace hooks reach their stages");
        steps.Clear();
        SesaRun.Build(dir, text, new SesaBuildOptions { Fsync = false, SaTrace = steps.Add }, "-d");
        Assert.That(steps, Does.Not.Contain("verify"), "the optional check is off by default");
        // the windows and the blocks of PH0 are the ones asked for
        var small = new IoCounters();
        var large = new IoCounters();
        SesaRun.Build(dir, text, new SesaBuildOptions { Fsync = false, Counters = small, TextBlockSymbols = 8, TextWindowBytes = 16 }, "-e");
        SesaRun.Build(dir, text, new SesaBuildOptions { Fsync = false, Counters = large }, "-f");
        Assert.That(small.ReadCalls, Is.GreaterThan(large.ReadCalls + 500), "blocks of 8 symbols read the input in many more calls");
        Assert.That(small.MappedWindows, Is.GreaterThan(large.MappedWindows + 50), "windows of 16 bytes map many more windows");
        // the verdicts of the stages are the gate of the owner and the cap that was passed
        o = SesaRun.Build(dir, text, new SesaBuildOptions { Fsync = false, RssAnonCapMb = 1e6 }, "-g");
        if (OperatingSystem.IsLinux() && File.Exists("/proc/self/status"))
        {
            Assert.That((o.Result.Sa.RssVerdict, o.Result.Lcp.RssVerdict), Is.EqualTo((RssVerdict.Within, RssVerdict.Within)));
            o = SesaRun.Build(dir, text, new SesaBuildOptions { Fsync = false, RssAnonCapMb = 0.001 }, "-h");
            Assert.That((o.Result.Sa.RssVerdict, o.Result.Lcp.RssVerdict), Is.EqualTo((RssVerdict.Exceeded, RssVerdict.Exceeded)));
        }
        o = SesaRun.Build(dir, text, new SesaBuildOptions { Fsync = false, ProcSource = new FakeProcSource() }, "-i");
        Assert.That((o.Result.Sa.RssVerdict, o.Result.Lcp.RssVerdict), Is.EqualTo((RssVerdict.Unavailable, RssVerdict.Unavailable)));
        o = SesaRun.Build(dir, text, new SesaBuildOptions { Fsync = false, MeasureRss = false }, "-j");
        Assert.That((o.Result.Sa.RssVerdict, o.Result.Lcp.RssVerdict), Is.EqualTo((RssVerdict.Unavailable, RssVerdict.Unavailable)));
        // the largest peak of the native buffers is the largest of the stages (PH4 and PH5 live together)
        Assert.That(r.ArenaPeakBytes, Is.GreaterThanOrEqualTo(Math.Max(r.Sa.Ph1ArenaPeakBytes, r.Sa.Ph23ArenaPeakBytes)));
        Assert.That(r.ArenaPeakBytes, Is.GreaterThanOrEqualTo(r.Lcp.ArenaPeakBytes + r.Lcp.ChildArenaPeakBytes));
        Assert.That(r.ArenaPeakBytes, Is.GreaterThanOrEqualTo(r.Psi.MappedBoundBytes));
        Assert.That(r.Lcp.Reservation, Is.EqualTo(r.Reservation), "the LCP is placed in the container, whose outcome it reports");
    }

    [Test]
    public void OptionsThatCannotBeUsed_AreRefusedBeforeAnythingIsCreated()
    {
        using var dir = new TempDir();
        string input = SesaRun.WriteInput(dir, "mississippi");
        string target = dir.File("x.sesa");
        var calls = new List<string>();
        foreach ((string label, SesaBuildOptions options) in new[]
        {
            ("no thread", new SesaBuildOptions { Threads = 0 }),
            ("negative arena cap", new SesaBuildOptions { ArenaCapBytes = -1 }),
            ("negative child cap", new SesaBuildOptions { ChildArenaCapBytes = -1 }),
            ("negative PH6 cap", new SesaBuildOptions { Ph6CapBytes = -1 }),
            ("negative RSS cap", new SesaBuildOptions { RssAnonCapMb = -1 }),
        })
        {
            SesaBuildOptions hooked = options with { Hook = calls.Add };
            Assert.Throws<ArgumentOutOfRangeException>(() => SesaBuilder.Build(input, target, hooked), label);
            Assert.That(calls, Is.Empty, $"{label}: the build had started (a hook point was reached) before the option was refused");
            Assert.That(Directory.GetFileSystemEntries(dir.Path).Select(Path.GetFileName), Is.EquivalentTo(new[] { "input.txt" }), $"{label}: something was created");
        }
    }
}
