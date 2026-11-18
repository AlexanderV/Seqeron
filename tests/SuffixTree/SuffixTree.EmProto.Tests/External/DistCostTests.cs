using System.Diagnostics;
using System.Globalization;
using SuffixTree.Esa;
using Prod = SuffixTree.Esa.External;
using Proto = SuffixTree.EmProto.Ph6;
using ProtoBuild = SuffixTree.EmProto.Build;
using ProtoDcs = SuffixTree.EmProto.Dcs;
using ProtoHarness = SuffixTree.EmProto.Harness;

namespace SuffixTree.EmProto.Tests.External;

/// <summary>
/// Gate G3 of work package 1b.7a (cost parity of the production <c>Dist</c> with the prototype <c>IsaByDist</c>): one timed run of ISA by DIST on a suffix-array file, driven by
/// <c>docs/refactoring/suffix-tree-upgrade/analysis/em/run_1b7a_g3.sh</c> (cold cache, processes one after the other, environment variables name the files). Not part of any run (<c>[Explicit]</c>).
/// Variants: <c>proto</c> (<c>IsaByDist</c>), <c>prod</c> (<c>Dist.BuildIsa</c> with <c>IoCounters</c>), <c>prodnc</c> (the same without counters, for the counter overhead).
/// </summary>
[TestFixture]
[Explicit("gate G3 of 1b.7a: run by analysis/em/run_1b7a_g3.sh, which sets the G3_* variables")]
[Category("G3Cost")]
public class DistCostTests
{
    private static string Env(string name) => Environment.GetEnvironmentVariable(name) ?? throw new InvalidOperationException($"{name} is not set");

    [Test]
    public void RunOnce()
    {
        string variant = Env("G3_VARIANT");
        string saPath = Env("G3_SA");
        long n = long.Parse(Env("G3_N"), CultureInfo.InvariantCulture);
        long r = long.Parse(Env("G3_R"), CultureInfo.InvariantCulture);
        int io = int.Parse(Env("G3_IO"), CultureInfo.InvariantCulture);
        string outPath = Env("G3_OUT");
        string scratch = Env("G3_SCRATCH");
        string rows = Env("G3_ROWS");
        string label = Env("G3_LABEL");

        var rss = new Prod.RssSampler(new Prod.RssReader(), 10);
        long before = new Prod.RssReader().Read().RssAnonKb;
        rss.Start();
        long arenaPeak;
        long read;
        long write;
        long mapped;
        long scratchPeak;
        long t0 = Stopwatch.GetTimestamp();
        switch (variant)
        {
            case "proto":
            {
                ProtoHarness.IoCounters.Reset();
                ProtoBuild.ScratchTracker.Reset();
                ProtoDcs.MappedIo.Reset();
                var budget = new ProtoHarness.ArenaBudget(long.MaxValue);
                Proto.IsaByDist.Build(saPath, n, outPath, new Proto.DistOptions { ScratchDir = scratch, ChunkPositions = r, IoBufferBytes = io, Budget = budget });
                arenaPeak = budget.PeakBytes;
                read = ProtoHarness.IoCounters.ReadBytes;
                write = ProtoHarness.IoCounters.WriteBytes;
                mapped = ProtoDcs.MappedIo.WriteBytes;
                scratchPeak = ProtoBuild.ScratchTracker.Peak;
                break;
            }
            case "prod":
            case "prodnc":
            {
                var counters = variant == "prod" ? new Prod.IoCounters() : null;
                var budget = new Prod.ArenaBudget(long.MaxValue);
                Prod.Dist.BuildIsa(saPath, n, outPath, new Prod.DistOptions { ScratchDir = scratch, ChunkEntries = r, IoBufferBytes = io, Budget = budget, Counters = counters });
                arenaPeak = budget.PeakBytes;
                read = counters?.ReadBytes ?? -1;
                write = counters?.WriteBytes ?? -1;
                mapped = counters?.MappedWriteBytes ?? -1;
                scratchPeak = counters?.ScratchPeakBytes ?? -1;
                break;
            }
            default:
                throw new InvalidOperationException("G3_VARIANT must be proto, prod or prodnc");
        }
        double seconds = (Stopwatch.GetTimestamp() - t0) / (double)Stopwatch.Frequency;
        Prod.RssStats stats = rss.Stop();
        string row = string.Join('\t', label, variant, r, io, seconds.ToString("F4", CultureInfo.InvariantCulture), before, stats.MaxRssAnonKb, stats.Samples, arenaPeak, read, write, mapped, scratchPeak);
        File.AppendAllText(rows, row + "\n");
        TestContext.Out.WriteLine(row);
        Assert.That(File.Exists(outPath), Is.True);
    }

    /// <summary>Cost of one counter update (the counters overhead of gate G3 is this times the number of updates of a run: about 50,000 for ISA of 10^8 entries).</summary>
    [Test]
    public void CounterOpCost()
    {
        var c = new Prod.IoCounters();
        const int n = 20_000_000;
        for (int i = 0; i < 1_000_000; i++) c.AddRead(1); // warm up the JIT
        long t0 = Stopwatch.GetTimestamp();
        for (int i = 0; i < n; i++)
        {
            c.AddRead(65536);
            c.AddWrite(65536);
            c.AddScratch(65536);
            c.RemoveScratch(65536);
        }
        double ns = (Stopwatch.GetTimestamp() - t0) * 1e9 / Stopwatch.Frequency / (4.0 * n);
        string row = string.Join('\t', "counter-op", "prod", "-", "-", ns.ToString("F2", CultureInfo.InvariantCulture), "ns per counter update, one thread", c.ReadBytes, c.ScratchPeakBytes);
        TestContext.Out.WriteLine(row);
        Assert.That(ns, Is.GreaterThan(0));
    }
}
