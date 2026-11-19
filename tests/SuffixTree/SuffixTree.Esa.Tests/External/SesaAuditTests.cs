using SuffixTree.Esa.External;
using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.Esa.Tests.External;

/// <summary>
/// Tests that close the gaps which the mutation audit of work package 1b.7e showed (rounds r1 and r2 of <c>analysis/mutation</c>): the exact entry of every section in the table, the widest LCP, the accounting of the buffers of
/// PH6 against the shared budget, and the limits of the arguments of the stages.
/// </summary>
[TestFixture]
public class SesaAuditTests
{
    private TempDir _dir = null!;

    [SetUp]
    public void SetUp() => _dir = new TempDir();

    [TearDown]
    public void TearDown() => _dir.Dispose();

    /// <summary>Every field of every section entry that this writer chooses, from the bytes of the file: flags, element width, codec and the two parameters.</summary>
    private static void AssertEveryEntry(SesaRun.Outcome o, bool psi, string label)
    {
        SesaOracle f = o.Oracle;
        const bool req = true, opt = false;
        foreach (OracleSection s in f.Sections)
        {
            (bool required, int width, int codec, long p0, long p1) expected = s.Kind switch
            {
                1 => (req, 0, 4, o.Case.Bits, 0L),
                2 => (req, 2, 1, 0L, 0L),
                13 => (opt, 8, 1, 0L, 0L),
                3 => (req, 4, 1, 0L, 0L),
                6 => (req, 4, 1, 0L, 0L),
                12 => (opt, 4, 1, SesaRun.PrimaryRow(o.Case.Sa), 0L),
                11 => (opt, 0, 1, 0L, 0L),
                4 => (req, o.Case.Lcp.Kind switch { LcpKind.U8 => 1, LcpKind.U16 => 2, _ => 4 }, o.Case.Lcp.Kind switch { LcpKind.U8 => 1, LcpKind.U16 => 2, _ => 3 }, o.Case.Lcp.OverflowCount, 0L),
                5 => (req, 8, 1, 0L, 0L),
                _ => throw new InvalidOperationException($"{label}: unknown kind {s.Kind}"),
            };
            Assert.That((s.Required, s.ElemWidth, s.Codec, s.Param0, s.Param1), Is.EqualTo(expected), $"{label}: entry of the section {s.Kind}");
            Assert.That(s.HasCrc, Is.True, $"{label}: section {s.Kind} carries HAS_CRC");
        }
        Assert.That(f.Sections.Select(s => s.Kind), Is.EqualTo(f.Sections.OrderBy(s => s.Offset).Select(s => s.Kind)), $"{label}: the table is in file order");
        Assert.That(f.Find(12) is not null, Is.EqualTo(psi), $"{label}: PSI");
    }

    [Test]
    public void EveryEntryOfTheTable_IsAsChosen_ForEveryWidthOfTheLcp()
    {
        // U8 without an overflow list, U8 with one, U16 without, U16 with one (a run above 65,534 symbols), and the 32-bit width (more than one escape in eight entries at both narrower widths: a long run of one symbol)
        string[] texts =
        [
            StringFamilies.RandomDna(3000, 1),
            StringFamilies.RandomDna(2800, 2) + new string('A', 400) + StringFamilies.RandomDna(2800, 3),
            new string('A', 500) + StringFamilies.RandomDna(100, 4),
            new string('A', 70_000),
            new string('A', 80_000),
        ];
        LcpKind[] kinds = [LcpKind.U8, LcpKind.U8, LcpKind.U16, LcpKind.U16, LcpKind.U32];
        bool[] overflow = [false, true, false, true, false];
        for (int k = 0; k < texts.Length; k++)
        {
            foreach (bool psi in k < 3 ? new[] { true, false } : new[] { true })
            {
                SesaRun.Outcome o = SesaRun.Build(_dir, texts[k], new SesaBuildOptions { Threads = 2, Fsync = false, IncludePsi = psi, MeasureRss = false }, $"-{k}{(psi ? "p" : "n")}");
                Assert.That(o.Case.Lcp.Kind, Is.EqualTo(kinds[k]), $"text {k}: the width of the LCP");
                Assert.That(o.Result.LcpKind, Is.EqualTo(kinds[k]));
                SesaRun.Compare(o, $"text {k}, psi {psi}", psi);
                AssertEveryEntry(o, psi, $"text {k}, psi {psi}");
                Assert.That(o.Oracle.Find(5) is not null, Is.EqualTo(overflow[k]), $"text {k}: the overflow list exists where the LCP has escapes");
                Directory.Delete(o.OutDir, true);
            }
        }
    }

    // ---- what the result reports

    [Test]
    public void TheResultFigures_AreInTheirUnits_AndEverySecondIsMeasured()
    {
        string text = StringFamilies.RandomDna(3000, 21);
        SesaRun.Outcome o = SesaRun.Build(_dir, text, new SesaBuildOptions { Threads = 2, Fsync = false }, "-u");
        SesaBuildResult r = o.Result;
        SesaPhaseSeconds t = r.Seconds;
        foreach ((string name, double v) in new[] { ("Scan", t.Scan), ("Create", t.Create), ("Pack", t.Pack), ("Ph1", t.Ph1), ("Ph23", t.Ph23), ("Ph45", t.Ph45), ("Bwt", t.Bwt), ("Psi", t.Psi), ("Crc", t.Crc), ("Seal", t.Seal), ("Total", t.Total) })
        {
            Assert.That(v, Is.GreaterThan(0).And.LessThan(120), $"{name}: seconds, not ticks, and measured");
        }
        Assert.That(t.Total, Is.GreaterThanOrEqualTo(t.Scan + t.Create + t.Pack + t.Ph1 + t.Ph23 + t.Ph45 + t.Bwt + t.Psi + t.Crc + t.Seal - 1e-3), "the phases lie inside the total");
        if (OperatingSystem.IsLinux() && File.Exists("/proc/self/status"))
        {
            foreach ((string name, double? v) in new (string, double?)[] { ("Ph0", r.Rss.Ph0), ("Ph1", r.Rss.Ph1), ("Ph23", r.Rss.Ph23), ("Ph45", r.Rss.Ph45), ("Ph6", r.Rss.Ph6), ("Seal", r.Rss.Seal), ("whole", r.RssAnonPeakMb) })
            {
                Assert.That(v, Is.Not.Null.And.InRange(5.0, 2000.0), $"{name}: RssAnon in MB (10^6 bytes) of a process that holds the runtime, not a fraction of a megabyte and not kilobytes");
            }
        }
        Assert.That((r.HasPsi, r.PrimaryRow), Is.EqualTo((true, SesaRun.PrimaryRow(o.Case.Sa))));
        Assert.That(r.CrcRows.Select(c => c.Kind), Is.EqualTo(r.Sections.Select(s => s.Kind)), "one checksum row per section, in the order of the table");
        Assert.That(r.CrcRows.Select(c => c.Bytes), Is.EqualTo(r.Sections.Select(s => s.Length)));
        Assert.That(r.CrcRows.All(c => c.Seconds >= 0), Is.True);
        Assert.That((r.Reservation, r.GrownBytes, r.GrowthReservation), Is.EqualTo((ReservationOutcome.Preallocated, 0L, (ReservationOutcome?)null)), "a U8 LCP without escapes: the file never grows");
        Assert.That(r.ScratchPeakBytes, Is.GreaterThanOrEqualTo(12L * text.Length), "the scratch peak is the one of PH4 (12 B per character and more) as the counters saw it");

        // without PSI: no BWT, no Psi, no primary row, no time for them
        o = SesaRun.Build(_dir, text, new SesaBuildOptions { Threads = 2, Fsync = false, IncludePsi = false }, "-v");
        r = o.Result;
        Assert.That((r.HasPsi, r.PrimaryRow, r.Psi, r.Bwt), Is.EqualTo((false, -1L, (PsiStats?)null, (BwtResult?)null)));
        Assert.That((r.Seconds.Bwt, r.Seconds.Psi), Is.EqualTo((0.0, 0.0)));
        Assert.That(r.Rss.Ph6, Is.Null, "the phase that did not run has no figure");
    }

    [Test]
    public void TheArenaPeak_IsTheLargestPeakOfAnyStage_WhicheverStageThatIs()
    {
        string text = StringFamilies.RandomDna(900, 22);
        var dominating = new HashSet<string>();
        var configs = new List<(string Label, SesaBuildOptions Options)>
        {
            ("defaults", new SesaBuildOptions { Threads = 2, Fsync = false }),
            ("no PSI", new SesaBuildOptions { Threads = 2, Fsync = false, IncludePsi = false }),
        };
        for (int v = 0; v < 4; v++) configs.Add(($"tiny {v}", SesaRun.Tiny(v, 2, v % 3 != 0)));
        configs.Add(("windows of PSI of 32 MiB", new SesaBuildOptions { Threads = 2, Fsync = false, Ph6Tuning = new SesaPh6Tuning { WindowBytes = 32 << 20 } }));
        foreach ((string label, SesaBuildOptions options) in configs)
        {
            SesaBuildResult r = SesaRun.Build(_dir, text, options, "-w").Result;
            long ph1 = r.Sa.Ph1ArenaPeakBytes, ph23 = r.Sa.Ph23ArenaPeakBytes, ph45 = r.Lcp.ArenaPeakBytes + r.Lcp.ChildArenaPeakBytes, ph6 = r.Psi?.MappedBoundBytes ?? 0;
            Assert.That(r.ArenaPeakBytes, Is.EqualTo(new[] { ph1, ph23, ph45, ph6 }.Max()), $"{label}: ph1 {ph1}, ph23 {ph23}, ph45 {ph45}, ph6 {ph6}");
            if (ph1 == r.ArenaPeakBytes) dominating.Add("ph1");
            if (ph23 == r.ArenaPeakBytes) dominating.Add("ph23");
            if (ph45 == r.ArenaPeakBytes) dominating.Add("ph45");
            if (ph6 == r.ArenaPeakBytes) dominating.Add("ph6");
            TestContext.Out.WriteLine($"{label}: ph1 {ph1}, ph23 {ph23}, ph45 {ph45}, ph6 {ph6}");
        }
        // PH1 is never the largest at a size that a test can build (a few KB against the 25 MB of the buckets of PH2 / PH3 at the defaults, and PH4 / PH5 hold at least 1.5 KB under the tiniest tuning); its term is checked in the equality only
        Assert.That(dominating, Is.EquivalentTo(new[] { "ph23", "ph45", "ph6" }), "the three other terms of the maximum are each the largest in one of the configurations");
    }

    [Test]
    public void TheVerdicts_AreWithinAtTheCapExactly_AndExceededJustAbove_AtEveryLevel()
    {
        string text = StringFamilies.RandomDna(400, 41);
        // a kernel that always says RssAnon = 1,000 kB: 1,000 x 1,024 / 10^6 = 1.024 MB
        FakeProcSource Proc() => new FakeProcSource().With("/proc/self/status", "Name:\tx\nRssAnon:\t    1000 kB\nRssFile:\t      10 kB\nVmHWM:\t    2000 kB\nThreads:\t4\n");
        SesaBuildResult within = SesaRun.Build(_dir, text, new SesaBuildOptions { Threads = 1, Fsync = false, ProcSource = Proc(), RssAnonCapMb = 1.024 }, "-a").Result;
        Assert.That((within.RssVerdict, within.Sa.RssVerdict, within.Lcp.RssVerdict), Is.EqualTo((RssVerdict.Within, RssVerdict.Within, RssVerdict.Within)), "at the cap exactly");
        foreach ((string name, double? v) in new (string, double?)[] { ("Ph0", within.Rss.Ph0), ("Ph1", within.Rss.Ph1), ("Ph23", within.Rss.Ph23), ("Ph45", within.Rss.Ph45), ("Ph6", within.Rss.Ph6), ("Seal", within.Rss.Seal), ("whole", within.RssAnonPeakMb) })
        {
            Assert.That(v, Is.EqualTo(1.024), $"{name}: 1,000 kB in MB (10^6 bytes)");
        }
        SesaBuildResult above = SesaRun.Build(_dir, text, new SesaBuildOptions { Threads = 1, Fsync = false, ProcSource = Proc(), RssAnonCapMb = 1.0239 }, "-b").Result;
        Assert.That((above.RssVerdict, above.Sa.RssVerdict, above.Lcp.RssVerdict), Is.EqualTo((RssVerdict.Exceeded, RssVerdict.Exceeded, RssVerdict.Exceeded)), "just above the cap");
    }

    [Test]
    public void AnArenaCapOfOneByte_StopsTheBuildInPh1_BeforePh4Starts_AndAnUnreachablePh4CapInPh4()
    {
        string text = StringFamilies.RandomDna(400, 42);
        string input = SesaRun.WriteInput(_dir, text);
        var steps = new List<string>();
        var options = new SesaBuildOptions { Threads = 1, Fsync = false, ArenaCapBytes = 1, SaTrace = steps.Add, LcpTrace = steps.Add };
        Assert.Throws<ArenaBudgetExceededException>(() => SesaBuilder.Build(input, _dir.File("a.sesa"), options));
        Assert.That(steps, Does.Not.Contain("b1"), "PH4 never started: the cap reached PH1 (and PH2 / PH3), not only PH4");
        steps.Clear();
        options = new SesaBuildOptions { Threads = 1, Fsync = false, ChildArenaCapBytes = 1, SaTrace = steps.Add, LcpTrace = steps.Add };
        Assert.Throws<ArenaBudgetExceededException>(() => SesaBuilder.Build(input, _dir.File("b.sesa"), options));
        Assert.That(steps, Does.Contain("ph23"), "PH1-PH3 ran under their own cap");
        steps.Clear();
        options = new SesaBuildOptions { Threads = 1, Fsync = false, Ph6CapBytes = 1, SaTrace = steps.Add, LcpTrace = steps.Add };
        var hooks = new List<string>();
        options = options with { Hook = hooks.Add };
        Assert.Throws<ArenaBudgetExceededException>(() => SesaBuilder.Build(input, _dir.File("c.sesa"), options));
        Assert.That(hooks, Does.Contain(SesaHooks.LcpChildBuilt).And.Not.Contain(SesaHooks.BwtBuilt), "PH6 failed under its own cap after PH4 / PH5 had run");
    }

    // ---- the ranges of the stages: refused at the entry, to the byte, before any work is done

    [Test]
    public void ARangeThatDoesNotFitTheFile_ByOneToFourBytes_IsRefusedAtTheEntryOfEveryStage_BeforeAnyWorkIsDone()
    {
        string text = StringFamilies.RandomDna(300, 51);
        LcpChildCase c = LcpChildCase.Of(text);
        SaCase sc = SaCase.Of(text);
        long arrayBytes = 4L * (c.N + 1);
        string textPath = TextOracle.WriteFile(_dir, "t.pk", c.Symbols, c.Bits);
        using MappedTextReader reader = MappedTextReader.Open(textPath, c.N, c.Bits);
        string satPath = SaRun.WriteText(_dir, sc, "-r");
        using MappedTextReader satReader = MappedTextReader.Open(satPath, sc.N, sc.Bits);
        var sparse = new MappedFileOptions { Reservation = ReservationMode.Sparse };
        for (int late = 1; late <= 4; late++)
        {
            string label = $"{late} byte(s) too late";
            // the suffix array stage writes a range of the container
            var steps = new List<string>();
            using (MappedOutputFile file = MappedOutputFile.Create(_dir.File($"sa{late}.bin"), 64 + arrayBytes, sparse))
            {
                Assert.Throws<ArgumentOutOfRangeException>(() => SaBuilder.Build(satReader, file, 64 + late, SaTiny.Build(_dir.Sub($"s{late}"), 2) with { Trace = steps.Add }), label);
                Assert.That(steps, Is.Empty, $"{label}: PH1 never started");
                Assert.That(Directory.Exists(_dir.Sub($"s{late}")), Is.False.Or.Property("Length").Zero, $"{label}: no scratch");
            }
            // the child range of the section stage
            steps.Clear();
            using (MappedOutputFile file = MappedOutputFile.Create(_dir.File($"ch{late}.bin"), 128 + arrayBytes + arrayBytes, sparse))
            {
                LcpChildOptions o = LcpChildRun.Options(_dir.Sub($"c{late}"), LcpChildRun.TinyTuning(2)) with { Trace = steps.Add };
                Assert.Throws<ArgumentOutOfRangeException>(() => SectionLcpChildBuilder.Build(_dir.File($"ch{late}.bin"), 128, reader, file, 128 + arrayBytes + late, _ => throw new InvalidOperationException("never asked"), o), label);
                Assert.That(steps, Is.Empty, $"{label}: PH4 never started");
                Assert.That(Directory.Exists(_dir.Sub($"c{late}")), Is.False, $"{label}: no scratch directory was made");
            }
            // PSI
            using (MappedOutputFile file = MappedOutputFile.Create(_dir.File($"ps{late}.bin"), 32 + arrayBytes, sparse))
            {
                var budget = new ArenaBudget(long.MaxValue);
                Assert.Throws<ArgumentOutOfRangeException>(() => PsiStream.Build(_dir.File("no-bwt.bin"), c.N, 0, [c.N], file, 32 + late, new PsiOptions { Budget = budget }), label);
                Assert.That(budget.PeakBytes, Is.Zero, $"{label}: nothing was reserved");
            }
            // the BWT pass reads a suffix array range that ends past the end of the file
            string container = _dir.File($"bw{late}.bin");
            File.WriteAllBytes(container, new byte[128 + arrayBytes - late]);
            string bwt = _dir.File($"bw{late}.bwt");
            Assert.Throws<InvalidDataException>(() => BwtPass.Build(container, 128, reader, c.Sigma, bwt), label);
            Assert.That(File.Exists(bwt), Is.False, $"{label}: no BWT file");
        }
    }

    [Test]
    public void ANullArgument_IsAnArgumentNullException_AndANegativeOffset_IsRefusedBeforePh1()
    {
        string text = StringFamilies.RandomDna(100, 52);
        SaCase sc = SaCase.Of(text);
        string satPath = SaRun.WriteText(_dir, sc, "-n");
        using MappedTextReader reader = MappedTextReader.Open(satPath, sc.N, sc.Bits);
        using MappedOutputFile file = MappedOutputFile.Create(_dir.File("n.bin"), 64 + (4L * (sc.N + 1)), new MappedFileOptions { Reservation = ReservationMode.Sparse });
        SaBuildOptions options = SaTiny.Build(_dir.Sub("sn"), 2);
        Assert.Throws<ArgumentNullException>(() => SaBuilder.Build(null!, file, 64, options));
        Assert.Throws<ArgumentNullException>(() => SaBuilder.Build(reader, null!, 64, options));
        Assert.Throws<ArgumentNullException>(() => SaBuilder.Build(reader, file, 64, null!));
        var steps = new List<string>();
        Assert.Throws<ArgumentOutOfRangeException>(() => SaBuilder.Build(reader, file, -4, options with { Trace = steps.Add }));
        Assert.That(steps, Is.Empty, "a negative offset is refused before PH1");
    }

    [Test]
    public void ThePsiStats_CountTheWindowsOfEveryPass_AndTheCountersSeeThem()
    {
        byte[] s = [.. StringFamilies.RandomDna(500, 53).Select(ch => (byte)"ACGT".IndexOf(ch))];
        var counters = new IoCounters();
        // 12-byte windows = 3 entries, 2 symbols per pass: 2 passes over 4 symbols
        Ph6Run.Outcome o = Ph6Run.Run(_dir, s, 4, null, new PsiOptions { WindowBytes = 12, SymbolsPerPass = 2, Counters = counters });
        Ph6Run.Compare(s, 4, o, "two passes");
        long expected = o.Bwt.Counts.Sum(count => (count + 2) / 3);
        Assert.That(o.Psi.Passes, Is.EqualTo(2));
        Assert.That(o.Psi.Windows, Is.EqualTo(expected), "every symbol opens ceil(count / 3) windows, summed over the passes");
        Assert.That(counters.MappedWindows, Is.EqualTo(expected + 1), "the counters saw the windows of the cursors and the one of PSI[0]");
        Assert.That(counters.MappedWriteBytes, Is.EqualTo((4L * (s.Length + 1)) - 4 + 4 + 0), "every byte of Psi is a mapped write, once");
    }

    [Test]
    public void ASymbolThatOccursLessOftenThanCounted_IsFound_AlsoWhenItsLastWindowIsUsedUp()
    {
        byte[] s = [.. StringFamilies.RandomDna(200, 54).Select(ch => (byte)"ACGT".IndexOf(ch))];
        string dir = _dir.Sub("p");
        Directory.CreateDirectory(dir);
        int[] sa = Ph6Run.Sa(s);
        string container = Ph6Run.WriteContainer(_dir, sa, 64, 64 + (4L * (s.Length + 1)) + 4 + 64 + (4L * (s.Length + 1)));
        string textPath = TextOracle.WriteFile(_dir, "tt.pk", s, 2);
        using MappedTextReader reader = MappedTextReader.Open(textPath, s.Length, 2);
        string bwtPath = _dir.File("w.bwt");
        BwtResult b = BwtPass.Build(container, 64, reader, 4, bwtPath);
        using MappedOutputFile file = MappedOutputFile.Create(_dir.File("psi.bin"), 4L * (s.Length + 1) + 64, new MappedFileOptions { Reservation = ReservationMode.Sparse });
        // a byte of the BWT that is no symbol of the alphabet is skipped by the pass, so the symbol it replaced is one short of its count; with windows of one entry the cursor of that symbol is
        // always at the end of its window when the BWT runs out, and only the comparison of the entries written with the count finds it
        byte[] bwtBytes = File.ReadAllBytes(bwtPath);
        int at = Array.IndexOf(bwtBytes, (byte)1, (int)b.PrimaryRow + 1);
        bwtBytes[at] = 9;
        File.WriteAllBytes(bwtPath, bwtBytes);
        Assert.Throws<InvalidDataException>(() => PsiStream.Build(bwtPath, s.Length, b.PrimaryRow, b.Counts, file, 0, new PsiOptions { WindowBytes = 4 }));
        Assert.Throws<InvalidDataException>(() => PsiStream.Build(bwtPath, s.Length, b.PrimaryRow, b.Counts, file, 0, new PsiOptions { WindowBytes = 4, SymbolsPerPass = 1 }));
        Assert.Throws<InvalidDataException>(() => PsiStream.Build(bwtPath, s.Length, b.PrimaryRow, b.Counts, file, 0, new PsiOptions { WindowBytes = 1 << 20 }), "and with a window that is not used up");
    }

    // ---- the counted calls of a fixed geometry (the sizes of every stage reach their stage: a size that is dropped or changed on the way changes a count)

    private static string Snapshot(IoSnapshot s) => $"{s.ReadBytes}/{s.WriteBytes}/{s.ReadCalls}/{s.WriteCalls}/{s.MappedWriteBytes}/{s.MappedWindows}/{s.MappedMaxWindowBytes}/{s.ScratchPeakBytes}";

    /// <summary>The figures of a whole build that do not depend on how the workers of PH4 interleave (the number of explicit write calls does).</summary>
    private static string Stable(IoSnapshot s) => $"{s.ReadBytes}/{s.WriteBytes}/{s.ReadCalls}/{s.MappedWriteBytes}/{s.MappedWindows}/{s.MappedMaxWindowBytes}/{s.ScratchPeakBytes}";

    [Test]
    public void TheCountedCalls_OfTwoTinyGeometries_ArePinned()
    {
        string text = StringFamilies.RandomDna(700, 31);
        var seen = new List<string>();
        foreach (int v in new[] { 3, 6, 11 })
        {
            SesaBuildResult r = SesaRun.Build(_dir, text, SesaRun.Tiny(v, 1, v != 6), $"-g{v}").Result;
            seen.Add($"v{v} total {Stable(r.Io)} ph1 {Snapshot(r.Sa.Ph1Io)} ph23 {Snapshot(r.Sa.Ph23Io)} arena {r.Sa.Ph1ArenaPeakBytes}/{r.Sa.Ph23ArenaPeakBytes}/{r.Lcp.ArenaPeakBytes}/{r.Lcp.ChildArenaPeakBytes}/{r.Psi?.MappedBoundBytes ?? 0}");
        }
        // read bytes / written bytes / read calls / (calls) / mapped bytes / windows / largest window / scratch peak, then the arena peaks of PH1, PH2 / PH3, PH4, PH5 and the window bound of PH6; measured once on the head of 1b.7e, identical in five runs; the phases of the suffix array have all eight figures,
        // the build as a whole the six that do not depend on the interleaving of the workers of PH4 (the 12 B per character of the scratch peak is the one of PH4)
        Assert.That(seen, Is.EqualTo(new[]
        {
            "v3 total 38385/32796/1282/10668/653/160/8400 ph1 0/0/0/0/860/8/112/864 ph23 16676/14692/558/263/2804/26/112/5476 arena 220/144/1452/160/64",
            "v6 total 36264/35012/733/7106/563/64/8400 ph1 0/0/0/0/804/35/55/816 ph23 18012/16860/360/301/2804/122/55/5328 arena 136/216/1868/96/0",
            "v11 total 31434/26800/1000/11008/749/224/8400 ph1 464/464/20/3/1200/150/80/1212 ph23 8512/8184/111/143/2804/351/80/5628 arena 480/1104/4548/96/40",
        }), "the counted I/O of three fixed geometries on a fixed text: a size that is not passed on to its stage changes a count");
    }

    // ---- the buffers of PH6 and the shared budget

    [Test]
    public void ThePh6Buffers_AreChargedToTheSharedBudget_ByTheSizeOfTheInput_AndGivenBack()
    {
        byte[] s = [.. StringFamilies.RandomDna(200, 11).Select(c => (byte)"ACGT".IndexOf(c))];
        int n = s.Length;
        var bwtBudget = new ArenaBudget(long.MaxValue);
        var psiBudget = new ArenaBudget(long.MaxValue);
        Ph6Run.Outcome o = Ph6Run.Run(_dir, s, 4, new BwtOptions { Budget = bwtBudget }, new PsiOptions { Budget = psiBudget });
        Ph6Run.Compare(s, 4, o, "default sizes");
        Assert.That(bwtBudget.PeakBytes, Is.EqualTo((4L * (n + 1)) + (n + 1)), "the BWT pass holds the suffix-array block (4 B per row) and the output block (1 B per row), not a megabyte each");
        Assert.That(bwtBudget.LiveBytes, Is.Zero, "and gives them back");
        long windows = 4L * (1 << 20) / 4 * 4; // four symbols, one window of 1 MiB each
        Assert.That(o.Psi.MappedBoundBytes, Is.EqualTo(windows));
        Assert.That(psiBudget.PeakBytes, Is.EqualTo(windows + n + 1), "the windows of the pass and one read block of the BWT (one byte per row)");
        Assert.That(psiBudget.LiveBytes, Is.Zero, "and gives them back");
    }

    [Test]
    public void ThePsiStream_RefusesATextAtTheLimitOfThePositionsAsARangeError_AndOneAboveItAsUnsupported()
    {
        using MappedOutputFile file = MappedOutputFile.Create(_dir.File("c.bin"), 256, new MappedFileOptions { Reservation = ReservationMode.Sparse });
        long limit = PackedTextFormat.MaxSymbols;
        string bwt = _dir.File("never-opened.bin");
        Assert.Throws<ArgumentOutOfRangeException>(() => PsiStream.Build(bwt, limit, 0, [limit], file, 0), "at the limit the text is allowed, the range of 4 (n + 1) bytes does not fit the file");
        Assert.Throws<NotSupportedException>(() => PsiStream.Build(bwt, limit + 1, 0, [limit + 1], file, 0), "above it the positions do not fit 32 bits");
    }
}
