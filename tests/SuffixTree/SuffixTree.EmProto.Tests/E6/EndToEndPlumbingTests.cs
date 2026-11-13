using System.Globalization;
using SuffixTree.EmProto.Build;
using SuffixTree.EmProto.Dcs;
using SuffixTree.EmProto.Experiments;
using SuffixTree.EmProto.Harness;
using SuffixTree.EmProto.Ph4;
using SuffixTree.EmProto.Ph5;
using SuffixTree.EmProto.Pipeline;
using SuffixTree.EmProto.Tests.Ph5;
using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.EmProto.Tests.WholePipeline;

/// <summary>The glue of the driver: which parameter reaches which stage, the scale-down run, the command line, the additive changes to the harness.</summary>
[TestFixture]
public class EndToEndPlumbingTests : TempDirFixture
{
    [Test]
    public void EveryStage_IsGivenTheResolvedParameters_TheThreadCount_TheModeAndTheCodec()
    {
        string text = StringFamilies.RandomDna(40_000, 21);
        string d = Run("o");
        EndToEndOptions o = PipelineChecks.Options(d, text, 9, 3) with
        {
            ScaleN = 3_000_000_000,
            Tuning = new EndToEndTuning(),
            Period = 64,
            ArenaCapBytes = 100L << 20,
            ChildArenaCapBytes = 30L << 20,
        };
        EndToEndResult res = EndToEnd.Run(o);
        EndToEndParameters p = ScaleRules.Resolve(40_000, o);
        Assert.That(res.Parameters, Is.EqualTo(p));
        Assert.That(p.Scaled, Is.True);

        Assert.That(res.DcsOptions.BucketCapacity, Is.EqualTo(p.DcsBucketCapacity));
        Assert.That(res.DcsOptions.BaseCase, Is.EqualTo(p.DcsBaseCase));
        Assert.That(res.DcsOptions.ChunkSlots, Is.EqualTo(p.DcsChunkSlots));
        Assert.That(res.DcsOptions.SortLeafCapacity, Is.EqualTo(p.DcsSortCapacity));
        Assert.That(res.DcsOptions.RunSpillBytes, Is.EqualTo(p.DcsRunSpillBytes));
        Assert.That(res.DcsOptions.FanOut, Is.EqualTo(p.DcsFanOut));
        Assert.That(res.DcsOptions.ReservoirSize, Is.EqualTo(p.DcsReservoir));
        Assert.That(res.DcsOptions.IoBufferBytes, Is.EqualTo(p.DcsIoBytes));
        Assert.That(res.DcsOptions.MergeBufferBudget, Is.EqualTo(p.DcsMergeBudget));
        Assert.That(res.DcsOptions.WindowBytes, Is.EqualTo(p.WindowBytes));
        Assert.That(res.DcsOptions.Threads, Is.EqualTo(3));
        Assert.That(res.DcsOptions.ScratchDir, Is.EqualTo(Path.Combine(Path.GetFullPath(o.ScratchDir), "ph1")));

        Assert.That(res.BuildOptions.LeafCapacity, Is.EqualTo(p.BucketCapacity));
        Assert.That(res.BuildOptions.FanOut, Is.EqualTo(p.FanOut));
        Assert.That(res.BuildOptions.Oversample, Is.EqualTo(p.Oversample));
        Assert.That(res.BuildOptions.SpillBytes, Is.EqualTo(p.SpillBytes));
        Assert.That(res.BuildOptions.Threads, Is.EqualTo(3));
        Assert.That(res.BuildOptions.ScratchDir, Is.EqualTo(Path.Combine(Path.GetFullPath(o.ScratchDir), "ph3")));

        Assert.That(res.LcpOptions.ChunkPositions, Is.EqualTo(p.LcpChunk));
        Assert.That(res.LcpOptions.Workers, Is.EqualTo(p.LcpWorkers));
        Assert.That(res.LcpOptions.Threads, Is.EqualTo(3));
        Assert.That(res.LcpOptions.Mode, Is.EqualTo(LcpMode.Distributed), "PH4 in base mode");
        Assert.That(res.LcpOptions.ReadBufferBytes, Is.EqualTo(p.LcpReadBytes));
        Assert.That(res.LcpOptions.TupleBufferBytes, Is.EqualTo(p.LcpTupleBytes));
        Assert.That(res.LcpOptions.WindowBytes, Is.EqualTo(p.WindowBytes));
        Assert.That(res.LcpOptions.ScratchDir, Is.EqualTo(Path.Combine(Path.GetFullPath(o.ScratchDir), "ph4")));
        Assert.That(res.Lcp.Workers, Is.EqualTo(Math.Min(p.LcpWorkers, (int)res.Lcp.TextChunks)));
        Assert.That(res.Lcp.TextChunks, Is.EqualTo((40_000 + p.LcpChunk - 1) / p.LcpChunk));

        Assert.That(res.ChildOptions.Codec, Is.EqualTo(ChildCodec.Packed));
        Assert.That(res.ChildOptions.BlockEntries, Is.EqualTo(p.BlockEntries));
        Assert.That(res.ChildOptions.WindowBytes, Is.EqualTo(p.ChildWindowBytes));
        Assert.That(res.ChildOptions.WindowSlots, Is.EqualTo(p.ChildSlots));
        Assert.That(res.ChildOptions.SpillFlushBytes, Is.EqualTo(p.SpillFlushBytes));
        Assert.That(res.ChildOptions.FlushWindows, Is.True);
        Assert.That(res.ChildOptions.ScratchDir, Is.EqualTo(Path.Combine(Path.GetFullPath(o.ScratchDir), "ph5")));

        Assert.That(res.DcsOptions.Budget!.CapBytes, Is.EqualTo(100L << 20));
        Assert.That(res.BuildOptions.Budget!.CapBytes, Is.EqualTo(100L << 20));
        Assert.That(res.LcpOptions.Budget!.CapBytes, Is.EqualTo(100L << 20));
        Assert.That(res.ChildOptions.Budget!.CapBytes, Is.EqualTo(30L << 20));
        Assert.That(res.DcsOptions.Budget, Is.Not.SameAs(res.BuildOptions.Budget), "a budget per stage: the stages run one after the other");
        Assert.That(res.LcpOptions.Budget, Is.Not.SameAs(res.ChildOptions.Budget), "PH5 has its own budget next to PH4");
    }

    [Test]
    public void TheScaledRun_At100000Symbols_EqualsTheInRamStructures_WithTheScaleDownParameters()
    {
        string text = StringFamilies.RandomDna(100_000, 31);
        string d = Run("s");
        EndToEndOptions o = PipelineChecks.Options(d, text, 0, 4) with { ScaleN = 3_000_000_000, Tuning = new EndToEndTuning(), Period = 1024, ExactPeriod = false };
        EndToEndResult res = PipelineChecks.RunAndCompare(ChildReference.Of(text), o, "scale-down 100k");
        Assert.That(res.CoverPeriod, Is.EqualTo(993), "the nominal 1024 becomes the Singer period 993");
        Assert.That(res.Parameters.BucketCapacity, Is.EqualTo(33), "K = n / 3000");
        Assert.That(res.Build.Leaves, Is.GreaterThan(1000), "n / K buckets, as at 3x10^9");
        Assert.That(res.Lcp.TextChunks, Is.EqualTo((100_000 + res.Parameters.LcpChunk - 1) / res.Parameters.LcpChunk));
    }

    [Test]
    public void TheRealBlockSizes_OnASmallText_AreEqualTooAndLeaveNoScratch()
    {
        string text = StringFamilies.RandomDna(200_000, 41);
        string d = Run("u");
        EndToEndOptions o = PipelineChecks.Options(d, text, 0, 4) with { ScaleN = 0, Tuning = new EndToEndTuning(), Period = 1024, ExactPeriod = false, ArenaGuards = null, SyncOutputs = true };
        EndToEndResult res = PipelineChecks.RunAndCompare(ChildReference.Of(text), o, "unscaled real sizes");
        Assert.That(res.Build.Leaves, Is.EqualTo(1), "one bucket: K = 2^20");
        Assert.That(res.Parameters.LcpChunk, Is.EqualTo(8_000_000));
        Assert.That(res.Phases.Select(x => x.Name), Does.Contain(EndToEnd.PhaseSync));
    }

    [Test]
    public void TheCommandLine_BecomesTheOptions()
    {
        EndToEndOptions o = E6.ParseOptions(Options.Parse(["--input", "i.txt", "--out", "o", "--scratch", "s", "--threads", "3", "--scale-n", "3000000000", "--period", "64", "--exact", "--codec", "tagged", "--budget-mb", "100", "--ph5-budget-mb", "20", "--no-sync", "--touch-buffers"]));
        Assert.That((o.InputPath, o.OutPrefix, o.ScratchDir, o.Threads, o.ScaleN, o.Period, o.ExactPeriod, o.Codec), Is.EqualTo(("i.txt", "o", "s", 3, 3_000_000_000L, 64, true, ChildCodec.Tagged)));
        Assert.That((o.ArenaCapBytes, o.ChildArenaCapBytes, o.SyncOutputs, o.TouchBuffers), Is.EqualTo((100L << 20, 20L << 20, false, true)));
        EndToEndOptions d = E6.ParseOptions(Options.Parse(["--input", "i.txt", "--out", "o", "--scratch", "s"]));
        Assert.That((d.Threads, d.ScaleN, d.Period, d.ExactPeriod, d.Codec, d.SyncOutputs, d.TouchBuffers, d.ArenaCapBytes, d.ChildArenaCapBytes), Is.EqualTo((4, 0L, 1024, false, ChildCodec.Packed, true, false, 160L << 20, 64L << 20)));
        Assert.That(d.Tuning, Is.EqualTo(new EndToEndTuning()), "no explicit value unless given");
        EndToEndTuning t = E6.ParseOptions(Options.Parse(["--input", "i", "--out", "o", "--scratch", "s", "--k", "11", "--fan", "3", "--k-s", "12", "--nbase", "13", "--dcs-chunk", "14", "--sort-k", "15", "--chunk", "16", "--b2-workers", "2", "--block-entries", "17", "--slots", "3", "--window-kb", "4", "--spill-kb", "5", "--io-kb", "6", "--tuple-kb", "7", "--child-window-kb", "8", "--prepare-block", "24", "--spill-flush-mb", "2"])).Tuning;
        Assert.That((t.BucketCapacity, t.FanOut, t.DcsBucketCapacity, t.DcsBaseCase, t.DcsChunkSlots, t.DcsSortCapacity, t.LcpChunk, t.LcpWorkers, t.BlockEntries, t.ChildSlots), Is.EqualTo((11L, 3, 12L, 13L, 14L, 15L, 16L, 2, 17, 3)));
        Assert.That((t.WindowBytes, t.SpillBytes, t.LcpReadBytes, t.LcpTupleBytes, t.ChildWindowBytes, t.PrepareBlockSymbols, t.SpillFlushBytes), Is.EqualTo((4096, 5120, 6144, 7168, 8192, 24, 2L << 20)));
        Assert.Throws<ArgumentException>(() => E6.ParseOptions(Options.Parse(["--input", "i", "--out", "o", "--scratch", "s", "--codec", "zip"])));
        Assert.Throws<ArgumentException>(() => E6.ParseOptions(Options.Parse(["--out", "o", "--scratch", "s"])));
    }

    [Test]
    public void TheCommands_MakeTheText_RunThePipeline_AndCompareWithTheReferenceFilesOfTheOtherStages()
    {
        if (!File.Exists("/proc/self/stat")) Assert.Ignore("the E6 report reads /proc (Linux)");
        string text = Path.Combine(Dir, "in.txt");
        var report = new Report(null);
        Assert.That(E6.MakeText(Options.Parse(["--family", "genome:150000", "--out", text]), report), Is.EqualTo(0));
        Assert.That(new FileInfo(text).Length, Is.EqualTo(150_000));
        string refPk = Path.Combine(Dir, "ref.pk");
        string refSa = Path.Combine(Dir, "ref.sa");
        string refLcp = Path.Combine(Dir, "ref.lcp");
        Assert.That(E3.Pack(Options.Parse(["--family", "genome:150000", "--out", refPk]), report), Is.EqualTo(0));
        Assert.That(E3.OracleHash(Options.Parse(["--family", "genome:150000", "--out", refSa]), report), Is.EqualTo(0));
        Assert.That(new FileInfo(refSa).Length, Is.EqualTo(4L * 150_001));
        Assert.That(E4.Oracle(Options.Parse(["--text", refPk, "--sa", refSa, "--out", Path.Combine(Dir, "ref")]), report), Is.EqualTo(0));
        Assert.That(E5.Oracle(Options.Parse(["--text", refPk, "--lcp", Path.Combine(Dir, "ref"), "--out", Path.Combine(Dir, "ref")]), report), Is.EqualTo(0));
        string[] args = ["--input", text, "--out", Path.Combine(Dir, "o"), "--scratch", Path.Combine(Dir, "scratch"), "--threads", "2", "--scale-n", "3000000000", "--rss-cap-mb", "100000", "--verify-pk", refPk, "--verify-sa", refSa, "--verify-lcp", Path.Combine(Dir, "ref"), "--verify-child", Path.Combine(Dir, "ref") + ".cldP", "--verify-stats", Path.Combine(Dir, "ref") + ".stats"];
        Assert.That(E6.Run(Options.Parse(args), report), Is.EqualTo(0), "every output equals its reference");
        Assert.That(Directory.Exists(Path.Combine(Dir, "scratch")), Is.False);

        // a corrupted reference is a failure (exit code 1), for each of the comparisons
        foreach (string victim in new[] { refPk, refSa, Path.Combine(Dir, "ref.lcp"), Path.Combine(Dir, "ref.ovf"), Path.Combine(Dir, "ref.cldP") })
        {
            byte[] bytes = File.ReadAllBytes(victim);
            if (bytes.Length == 0) continue;
            byte[] bad = [.. bytes];
            bad[bytes.Length / 2] ^= 0x55;
            File.WriteAllBytes(victim, bad);
            Assert.That(E6.Run(Options.Parse(args), report), Is.EqualTo(1), $"a corrupted {Path.GetFileName(victim)} must fail the run");
            File.WriteAllBytes(victim, bytes);
        }
        string stats = Path.Combine(Dir, "ref.stats");
        string[] lines = File.ReadAllLines(stats);
        File.WriteAllLines(stats, lines.Select(l => l.StartsWith("lrsdepth=", StringComparison.Ordinal) ? "lrsdepth=1" : l));
        Assert.That(E6.Run(Options.Parse(args), report), Is.EqualTo(1), "different statistics fail the run");
        File.WriteAllLines(stats, lines);
        Assert.That(E6.Run(Options.Parse(args), report), Is.EqualTo(0));
    }

    [Test]
    public void TheRssCap_IsPartOfTheExitCode()
    {
        if (!File.Exists("/proc/self/stat")) Assert.Ignore("the E6 report reads /proc (Linux)");
        string text = Path.Combine(Dir, "in.txt");
        File.WriteAllBytes(text, PipelineChecks.Latin1(StringFamilies.RandomDna(5000, 3)));
        string[] common = ["--input", text, "--out", Path.Combine(Dir, "o"), "--scratch", Path.Combine(Dir, "scratch"), "--scale-n", "3000000000"];
        Assert.That(E6.Run(Options.Parse([.. common, "--rss-cap-mb", "100000"]), new Report(null)), Is.EqualTo(0));
        Assert.That(E6.Run(Options.Parse([.. common, "--rss-cap-mb", "0.001"]), new Report(null)), Is.EqualTo(1), "a cap below the process's RssAnon fails the run");
    }

    [Test]
    public void TheReport_NamesEveryPhaseAndTheTotals()
    {
        if (!File.Exists("/proc/self/stat")) Assert.Ignore("the E6 report reads /proc (Linux)");
        string text = Path.Combine(Dir, "in.txt");
        File.WriteAllBytes(text, PipelineChecks.Latin1(StringFamilies.RandomDna(20_000, 3)));
        string rows = Path.Combine(Dir, "rows.tsv");
        int rc = E6.Run(Options.Parse(["--input", text, "--out", Path.Combine(Dir, "o"), "--scratch", Path.Combine(Dir, "scratch"), "--scale-n", "3000000000", "--rss-cap-mb", "100000", "--case", "rep"]), new Report(rows));
        Assert.That(rc, Is.EqualTo(0));
        string[][] lines = [.. File.ReadAllLines(rows).Select(l => l.Split('\t'))];
        Dictionary<string, string> v = lines.Where(l => l[2] == "rep").ToDictionary(l => l[3], l => l[4]);
        foreach (string ph in new[] { "ph0_prepare", "ph1_dcs", "ph23_sa", "ph45_lcp_child", "sync" })
        {
            foreach (string metric in new[] { "_s", "_io_explicit_read", "_io_explicit_write", "_io_mapped_write", "_proc_rchar", "_proc_wchar", "_rssanon_peak", "_scratch_peak_cum", "_dirty_peak", "_proc_dirty_peak" })
                Assert.That(v, Does.ContainKey(ph + metric));
        }
        foreach (string metric in new[] { "total_s", "sa_stage_s", "wall_s", "ns_per_symbol", "throughput", "rssanon_peak", "io_total_counted", "proc_rw_total", "scratch_peak", "output_bytes", "dirty_peak", "proc_dirty_peak", "node_count", "lrs", "ph4_width", "clean", "sha256", "parameters" })
            Assert.That(v, Does.ContainKey(metric));
        Assert.That(double.Parse(v["total_s"], CultureInfo.InvariantCulture), Is.EqualTo(EndToEnd.StagePhases.Sum(ph => double.Parse(v[ph + "_s"], CultureInfo.InvariantCulture))).Within(1e-4));
        Assert.That(double.Parse(v["sa_stage_s"], CultureInfo.InvariantCulture), Is.EqualTo(new[] { "ph0_prepare", "ph1_dcs", "ph23_sa" }.Sum(ph => double.Parse(v[ph + "_s"], CultureInfo.InvariantCulture))).Within(1e-4));
        long counted = long.Parse(v["io_total_counted"], CultureInfo.InvariantCulture);
        Assert.That(counted, Is.GreaterThan(30L * 20_000), "SA + LCP + CHILD + the tuples of PH4 and the buckets of PH3 are counted");
        Assert.That(counted, Is.LessThan(100L * 20_000), "gate G3: at most 100 B/char");
        Assert.That(v["clean"], Is.EqualTo("yes"));
    }

    [Test]
    public void TheDirtySamplerAddition_TakesAndResetsThePhaseMaxima()
    {
        if (!File.Exists("/proc/meminfo")) Assert.Ignore("Linux only");
        using var sampler = new DirtySampler(5);
        sampler.Start();
        Assert.That(sampler.BaselineKb, Is.GreaterThanOrEqualTo(0));
        SpinWait.SpinUntil(static () => false, 80);
        sampler.Stop();
        (long dirty, long process) = sampler.TakePhaseMax();
        Assert.That(dirty, Is.GreaterThanOrEqualTo(0));
        Assert.That(process, Is.GreaterThanOrEqualTo(0));
        (long dirty2, long process2) = sampler.TakePhaseMax();
        Assert.That((dirty2, process2), Is.EqualTo((0L, 0L)), "the maxima are reset by the call");
    }

    [Test]
    public void TheDirtySamplerAddition_SeesTheDirtyPagesOfAMappedFile_ThenResets()
    {
        if (!File.Exists("/proc/self/smaps")) Assert.Ignore("Linux only");
        string path = Path.Combine(Dir, "mapped.bin");
        using var sampler = new DirtySampler(2, filePrefix: path);
        sampler.Start();
        sampler.TakePhaseMax();
        using (var w = new MappedWindowWriter(path, 8 << 20, 8 << 20))
        {
            w.WriteAt(0, new byte[8 << 20].Select(_ => (byte)1).ToArray());
            SpinWait.SpinUntil(static () => false, 60);
            (long dirty, long process) = sampler.TakePhaseMax();
            Assert.That(process, Is.GreaterThanOrEqualTo(4096), "the window is mapped and dirty: at least 4 MiB of the 8 MiB are seen");
            Assert.That(dirty, Is.GreaterThanOrEqualTo(0));
        }
        sampler.Stop();
    }

    [Test]
    public void TheMappedSaSink_CountsItsWrites_InTheMappedIoCounter()
    {
        MappedIo.Reset();
        using (var sink = new MappedSaSink(Path.Combine(Dir, "sa.bin"), 1000))
        {
            sink.Write(0, new uint[10]);
            sink.Write(500, new uint[100]);
            sink.Write(999, new uint[1]);
            sink.Write(3, ReadOnlySpan<uint>.Empty);
        }
        Assert.That(MappedIo.WriteBytes, Is.EqualTo(4 * 111));
    }
}
