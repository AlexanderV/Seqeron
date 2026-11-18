using System.Diagnostics;
using System.Globalization;
using SuffixTree.EmProto.Build;
using SuffixTree.EmProto.Dcs;
using SuffixTree.EmProto.Harness;
using SuffixTree.EmProto.Text;
using Prod = SuffixTree.Esa.External;
using ProtoPipe = SuffixTree.EmProto.Pipeline;

namespace SuffixTree.EmProto.Tests.External;

/// <summary>
/// Gates G2 (the whole chain at 10^8) and G3 / G5 (cost parity of the production PH1 and PH2 / PH3 with the prototype's) of work package 1b.7d: one run per process, driven by
/// <c>docs/refactoring/suffix-tree-upgrade/analysis/em/run_1b7d_g2.sh</c> and <c>run_1b7d_g3.sh</c> (cold or warm cache, processes one after the other, environment variables name the files). Not part of any run
/// (<c>[Explicit]</c>). Variants: <c>proto</c> (PH1 with <c>ExternalDcs.BuildRankFile</c> and PH2 / PH3 with <c>EmBuilder.Run</c> and the mapped sink, the parameters of <see cref="ProtoPipe.ScaleRules"/>, phase by phase as
/// <c>EndToEnd.Run</c> does, the rank file deleted when PH3 is done), <c>prod</c> (<c>SaBuilder.Build</c>), <c>prodnr</c> (without the RSS sampler of the builder), <c>prodsparse</c> (the files created sparse instead of
/// reserved: G5), <c>prodverify</c> (with the optional check of the rank file between PH1 and PH2 / PH3: G5) and <c>e2e</c> (PH0 <c>PackedTextWriter</c>, then <c>SaBuilder</c>, then <c>LcpChildBuilder</c>, one shared <c>IoCounters</c>: G2). The timed region of the stage variants is the build call including what it
/// removes at its end (the rank file, the scratch directories): the prototype's side does the same.
/// </summary>
[TestFixture]
[Explicit("gates G2 / G3 / G5 of 1b.7d: run by analysis/em/run_1b7d_g2.sh and run_1b7d_g3.sh, which set the G3_* variables")]
[Category("G3Cost")]
public unsafe class SaBuildCostTests
{
    private static string Env(string name) => Environment.GetEnvironmentVariable(name) ?? throw new InvalidOperationException($"{name} is not set");

    private static string EnvOr(string name, string fallback) => Environment.GetEnvironmentVariable(name) ?? fallback;

    private static string F(double v) => v.ToString("F4", CultureInfo.InvariantCulture);

    [Test]
    public void RunOnce()
    {
        string variant = Env("G3_VARIANT");
        string rows = Env("G3_ROWS");
        string label = Env("G3_LABEL");
        string scratch = Env("G3_SCRATCH");
        long scale = long.Parse(EnvOr("G3_SCALE", "3000000000"), CultureInfo.InvariantCulture);
        int threads = int.Parse(EnvOr("G3_THREADS", "4"), CultureInfo.InvariantCulture);
        Directory.CreateDirectory(scratch);
        var rss = new Prod.RssSampler(new Prod.RssReader(), 10);
        long before = new Prod.RssReader().Read().RssAnonKb;
        rss.Start();
        string row;
        switch (variant)
        {
            case "proto":
                row = RunProto(scratch, scale, threads);
                Prod.RssStats ps = rss.Stop();
                File.AppendAllText(rows, string.Join('\t', label, variant, row, before, ps.MaxRssAnonKb, ps.Samples) + "\n");
                Assert.That(File.Exists(Env("G3_SA")), Is.True);
                return;
            case "prod":
            case "prodnr":
            case "prodsparse":
            case "prodverify":
                row = RunProd(variant, scratch, scale, threads, out Prod.SaBuildResult? result);
                Prod.RssStats s = rss.Stop();
                File.AppendAllText(rows, string.Join('\t', label, variant, row, before, s.MaxRssAnonKb, s.Samples) + "\n");
                TestContext.Out.WriteLine(row);
                Assert.That(result, Is.Not.Null);
                return;
            case "e2e":
                RunEndToEnd(label, scratch, scale, threads, rows, rss, before);
                return;
            default:
                throw new InvalidOperationException("G3_VARIANT must be proto, prod, prodnr, prodsparse, prodverify or e2e");
        }
    }

    private static string RunProto(string scratch, long scale, int threads)
    {
        string pk = Env("G3_PK");
        long n = long.Parse(Env("G3_N"), CultureInfo.InvariantCulture);
        int bits = int.Parse(Env("G3_BITS"), CultureInfo.InvariantCulture);
        string saPath = Env("G3_SA");
        ProtoPipe.EndToEndParameters p = ProtoPipe.ScaleRules.Resolve(n, new ProtoPipe.EndToEndOptions { InputPath = "x", OutPrefix = "x", ScratchDir = "x", Threads = threads, ScaleN = scale });
        DifferenceCover cover = DifferenceCover.Create(DifferenceCover.NearestSingerPeriod(1024));
        using var mapped = new MappedFile(pk, writable: false);
        using var packed = new PackedText(mapped.Pointer, n, bits, mapped, ownsOwner: false);
        string dirPh1 = Path.Combine(scratch, "ph1"), dirPh3 = Path.Combine(scratch, "ph3"), ranksPath = Path.Combine(scratch, "ranks.bin");
        var dcsBudget = new ArenaBudget(160L << 20);
        var saBudget = new ArenaBudget(160L << 20);
        IoCounters.Reset();
        MappedIo.Reset();
        ScratchTracker.Reset();
        long t0 = Stopwatch.GetTimestamp();
        ExtDcsStats dcs = ExternalDcs.BuildRankFile(packed, cover, ranksPath, new ExtDcsOptions
        {
            ScratchDir = dirPh1,
            BucketCapacity = p.DcsBucketCapacity,
            BaseCase = p.DcsBaseCase,
            ChunkSlots = p.DcsChunkSlots,
            SortLeafCapacity = p.DcsSortCapacity,
            Threads = threads,
            Seed = 1,
            FanOut = p.DcsFanOut,
            ReservoirSize = p.DcsReservoir,
            RunSpillBytes = p.DcsRunSpillBytes,
            IoBufferBytes = p.DcsIoBytes,
            MergeBufferBudget = p.DcsMergeBudget,
            WindowBytes = p.WindowBytes,
            Budget = dcsBudget,
        });
        long rankBytes = new FileInfo(ranksPath).Length;
        long scratchPh1 = ScratchTracker.Peak; // the working files of PH1 (the rank file is added to the account only from here on, as the pipeline of E6 does)
        ScratchTracker.Add(rankBytes);
        long t1 = Stopwatch.GetTimestamp();
        long r1 = IoCounters.ReadBytes, w1 = IoCounters.WriteBytes, m1 = MappedIo.WriteBytes;
        Directory.CreateDirectory(dirPh3);
        using (var ranks = new DcsRankTable(cover, new MappedFile(ranksPath, writable: false), n))
        using (var sink = new MappedSaSink(saPath, n + 1))
        {
            _ = EmBuilder.Run(packed, new BuildOptions { LeafCapacity = p.BucketCapacity, FanOut = p.FanOut, Threads = threads, Seed = 1, Oversample = p.Oversample, ScratchDir = dirPh3, SpillBytes = p.SpillBytes, Budget = saBudget }, sink, ranks);
        }
        File.Delete(ranksPath);
        ScratchTracker.Remove(rankBytes);
        long t2 = Stopwatch.GetTimestamp();
        foreach (string d in new[] { dirPh1, dirPh3 }) if (Directory.Exists(d)) Directory.Delete(d, true);
        long t3 = Stopwatch.GetTimestamp();
        double ph1 = (t1 - t0) / (double)Stopwatch.Frequency, ph23 = (t2 - t1) / (double)Stopwatch.Frequency, total = (t3 - t0) / (double)Stopwatch.Frequency;
        long r2 = IoCounters.ReadBytes, w2 = IoCounters.WriteBytes, m2 = MappedIo.WriteBytes;
        string row = string.Join('\t', F(ph1), F(ph23), F(total), r1, w1, m1, r2 - r1, w2 - w1, m2 - m1, scratchPh1, ScratchTracker.Peak, dcsBudget.PeakBytes, saBudget.PeakBytes, dcs.MaxLevel, "-");
        TestContext.Out.WriteLine(row);
        return row;
    }

    private static string RunProd(string variant, string scratch, long scale, int threads, out Prod.SaBuildResult result)
    {
        string pk = Env("G3_PK");
        long n = long.Parse(Env("G3_N"), CultureInfo.InvariantCulture);
        int bits = int.Parse(Env("G3_BITS"), CultureInfo.InvariantCulture);
        string saPath = Env("G3_SA");
        using Prod.MappedTextReader text = Prod.MappedTextReader.Open(pk, n, bits);
        var options = new Prod.SaBuildOptions
        {
            ScratchDir = scratch,
            Threads = threads,
            ScaleN = scale,
            MeasureRss = variant != "prodnr",
            Files = variant == "prodsparse" ? new Prod.MappedFileOptions { Reservation = Prod.ReservationMode.Sparse } : null,
            VerifyRanks = variant == "prodverify" || EnvOr("G3_VERIFY", "0") == "1",
        };
        long t0 = Stopwatch.GetTimestamp();
        result = Prod.SaBuilder.Build(text, saPath, options);
        double total = (Stopwatch.GetTimestamp() - t0) / (double)Stopwatch.Frequency;
        string row = string.Join('\t', F(result.Ph1Seconds), F(result.Ph23Seconds), F(total), result.Ph1Io.ReadBytes, result.Ph1Io.WriteBytes, result.Ph1Io.MappedWriteBytes, result.Ph23Io.ReadBytes, result.Ph23Io.WriteBytes, result.Ph23Io.MappedWriteBytes, result.Ph1ScratchPeakBytes, result.Ph23ScratchPeakBytes, result.Ph1ArenaPeakBytes, result.Ph23ArenaPeakBytes, result.RecursionLevels, result.RssVerdict);
        return row;
    }

    private static void RunEndToEnd(string label, string scratch, long scale, int threads, string rows, Prod.RssSampler rss, long before)
    {
        string input = Env("G2_TXT");
        string outPrefix = Env("G2_OUT");
        string sigmaText = Env("G2_SIGMA");
        var counters = new Prod.IoCounters();
        long t0 = Stopwatch.GetTimestamp();
        Prod.PackedTextInfo info = Prod.PackedTextWriter.WriteFile(input, outPrefix + ".pk", new Prod.PackedTextWriterOptions { Counters = counters });
        double ph0 = (Stopwatch.GetTimestamp() - t0) / (double)Stopwatch.Frequency;
        File.WriteAllText(outPrefix + ".pk.meta", string.Create(CultureInfo.InvariantCulture, $"n={info.N}\nbits={info.Bits}\nsigma={info.Sigma}\nletters_hex={Convert.ToHexString(info.Letters)}\ncounts={string.Join(',', info.Counts)}\nsha256={info.Sha256}\nname={Path.GetFileName(input)}\n"));
        Assert.That(info.Sigma.ToString(CultureInfo.InvariantCulture), Is.EqualTo(sigmaText));
        using Prod.MappedTextReader text = Prod.MappedTextReader.Open(outPrefix + ".pk", info.N, info.Bits);
        Prod.SaBuildResult sa = Prod.SaBuilder.Build(text, outPrefix + ".sa", new Prod.SaBuildOptions { ScratchDir = scratch, Threads = threads, ScaleN = scale, Counters = counters });
        var files = new Prod.LcpChildOutputs(outPrefix + ".lcp", outPrefix + ".ovf", outPrefix + ".cld", outPrefix + ".stats");
        Prod.LcpChildResult lcp = Prod.LcpChildBuilder.Build(outPrefix + ".sa", text, files, new Prod.LcpChildOptions { ScratchDir = scratch, Threads = threads, ScaleN = scale, Sigma = info.Sigma, Counters = counters });
        Prod.RssStats s = rss.Stop();
        string row = string.Join('\t', label, "e2e", F(ph0), F(sa.Ph1Seconds), F(sa.Ph23Seconds), F(lcp.Seconds), before, s.MaxRssAnonKb, s.Samples, sa.RssAnonPeakMb?.ToString("F1", CultureInfo.InvariantCulture) ?? "na", sa.RssVerdict, lcp.RssVerdict, counters.ReadBytes, counters.WriteBytes, counters.MappedWriteBytes, counters.ScratchPeakBytes);
        File.AppendAllText(rows, row + "\n");
        TestContext.Out.WriteLine(row);
        Assert.That(File.Exists(outPrefix + ".cld"), Is.True);
    }
}
