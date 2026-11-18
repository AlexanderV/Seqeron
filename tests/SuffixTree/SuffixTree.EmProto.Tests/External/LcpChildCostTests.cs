using System.Diagnostics;
using System.Globalization;
using SuffixTree.EmProto.Experiments;
using SuffixTree.EmProto.Harness;
using SuffixTree.EmProto.Ph4;
using SuffixTree.EmProto.Ph5;
using SuffixTree.EmProto.Text;
using Prod = SuffixTree.Esa.External;
using ProtoBuild = SuffixTree.EmProto.Build;
using ProtoDcs = SuffixTree.EmProto.Dcs;
using ProtoPipe = SuffixTree.EmProto.Pipeline;

namespace SuffixTree.EmProto.Tests.External;

/// <summary>
/// Gate G3 of work package 1b.7c (cost parity of the production PH4 + PH5 with the prototype's, fused through the B3 hook as in E6): one timed run on a packed text and its suffix array file, driven by
/// <c>docs/refactoring/suffix-tree-upgrade/analysis/em/run_1b7c_g3.sh</c> (cold or warm cache, processes one after the other, environment variables name the files). Not part of any run (<c>[Explicit]</c>).
/// Variants: <c>proto</c> (<c>ExternalLcp.Build</c> with <c>E5.LazyScan</c>, the parameters of <see cref="ProtoPipe.ScaleRules"/>), <c>prod</c> (<c>LcpChildBuilder.Build</c>), <c>prodnr</c> (the same without the RSS sampler of the builder).
/// The timed region is the build call alone (the text is mapped before it) including what it removes at its end (the spill file, the scratch directories): the prototype's side does the same.
/// </summary>
[TestFixture]
[Explicit("gate G3 of 1b.7c: run by analysis/em/run_1b7c_g3.sh, which sets the G3_* variables")]
[Category("G3Cost")]
public unsafe class LcpChildCostTests
{
    private static string Env(string name) => Environment.GetEnvironmentVariable(name) ?? throw new InvalidOperationException($"{name} is not set");

    private static string EnvOr(string name, string fallback) => Environment.GetEnvironmentVariable(name) ?? fallback;

    [Test]
    public void RunOnce()
    {
        string variant = Env("G3_VARIANT");
        string pk = Env("G3_PK");
        long n = long.Parse(Env("G3_N"), CultureInfo.InvariantCulture);
        int bits = int.Parse(Env("G3_BITS"), CultureInfo.InvariantCulture);
        string saPath = Env("G3_SA");
        string outPrefix = Env("G3_OUT");
        string scratch = Env("G3_SCRATCH");
        string rows = Env("G3_ROWS");
        string label = Env("G3_LABEL");
        long scale = long.Parse(EnvOr("G3_SCALE", "3000000000"), CultureInfo.InvariantCulture);
        int threads = int.Parse(EnvOr("G3_THREADS", "4"), CultureInfo.InvariantCulture);
        Directory.CreateDirectory(scratch);

        var rss = new Prod.RssSampler(new Prod.RssReader(), 10);
        long before = new Prod.RssReader().Read().RssAnonKb;
        rss.Start();
        string extra;
        double seconds;
        long read, write, mapped, childBytes, arenaPeak, childArenaPeak;
        switch (variant)
        {
            case "proto":
            {
                ProtoPipe.EndToEndParameters p = ProtoPipe.ScaleRules.Resolve(n, new ProtoPipe.EndToEndOptions { InputPath = "x", OutPrefix = "x", ScratchDir = "x", Threads = threads, ScaleN = scale });
                var lcpBudget = new ArenaBudget(160L << 20);
                var childBudget = new ArenaBudget(64L << 20);
                using var mappedFile = new MappedFile(pk, writable: false);
                using var packed = new PackedText(mappedFile.Pointer, n, bits, mappedFile, ownsOwner: false);
                var lcp = new ExtLcpOptions { ScratchDir = Path.Combine(scratch, "ph4"), ChunkPositions = p.LcpChunk, Threads = threads, Workers = p.LcpWorkers, ReadBufferBytes = p.LcpReadBytes, TupleBufferBytes = p.LcpTupleBytes, WindowBytes = p.WindowBytes, Mode = LcpMode.Distributed, Budget = lcpBudget };
                var child = new ChildBuildOptions { ScratchDir = Path.Combine(scratch, "ph5"), BlockEntries = p.BlockEntries, WindowBytes = p.ChildWindowBytes, WindowSlots = p.ChildSlots, FlushWindows = true, SpillFlushBytes = p.SpillFlushBytes, ReadBufferBytes = p.LcpReadBytes, Codec = ChildCodec.Packed, Budget = childBudget };
                IoCounters.Reset();
                ProtoDcs.MappedIo.Reset();
                ProtoBuild.ScratchTracker.Reset();
                long t0 = Stopwatch.GetTimestamp();
                ExtLcpStats st;
                ChildScanResult r;
                using (var scan = new E5.LazyScan(n, new ChildScanOutputs { ChildPath = outPrefix + ".cld" }, child, null))
                {
                    st = ExternalLcp.Build(saPath, packed, outPrefix + ".lcp", outPrefix + ".ovf", lcp, scan);
                    r = scan.Complete();
                } // the disposal removes the spill file: part of the call, as it is for the production builder
                foreach (string d in new[] { lcp.ScratchDir, child.ScratchDir }) if (Directory.Exists(d)) Directory.Delete(d, true);
                seconds = (Stopwatch.GetTimestamp() - t0) / (double)Stopwatch.Frequency;
                read = IoCounters.ReadBytes;
                write = IoCounters.WriteBytes;
                mapped = ProtoDcs.MappedIo.WriteBytes;
                childBytes = 4 * r.ChildWindows.Writes;
                arenaPeak = lcpBudget.PeakBytes;
                childArenaPeak = childBudget.PeakBytes;
                extra = string.Create(CultureInfo.InvariantCulture, $"{r.Spill.BytesWritten}\t{r.Spill.BytesRead}\t{r.ChildWindows.FlushedBytes}\t{r.ChildWindows.Reopened}\t{r.StackPeak}\t{st.Workers}\t{p.LcpChunk}\t{ProtoBuild.ScratchTracker.Peak}\t{st.B1Seconds:F3}\t{st.B2Seconds:F3}\t{st.B3Seconds:F3}\t{st.CompareSteps}");
                break;
            }
            case "prod":
            case "prodnr":
            {
                using Prod.MappedTextReader text = Prod.MappedTextReader.Open(pk, n, bits);
                string? sigma = Environment.GetEnvironmentVariable("G3_SIGMA"); // set by run_1b7c_g2.sh: the statistics file is written too
                var options = new Prod.LcpChildOptions { ScratchDir = scratch, Threads = threads, ScaleN = scale, MeasureRss = variant == "prod", Sigma = sigma is null ? null : int.Parse(sigma, CultureInfo.InvariantCulture) };
                var files = new Prod.LcpChildOutputs(outPrefix + ".lcp", outPrefix + ".ovf", outPrefix + ".cld", sigma is null ? null : outPrefix + ".stats");
                long t0 = Stopwatch.GetTimestamp();
                Prod.LcpChildResult r = Prod.LcpChildBuilder.Build(saPath, text, files, options);
                seconds = (Stopwatch.GetTimestamp() - t0) / (double)Stopwatch.Frequency;
                read = r.Io.ReadBytes;
                write = r.Io.WriteBytes;
                childBytes = 4 * r.Scan.ChildWindows.Writes;
                mapped = r.Io.MappedWriteBytes - childBytes;
                arenaPeak = r.ArenaPeakBytes;
                childArenaPeak = r.ChildArenaPeakBytes;
                extra = string.Create(CultureInfo.InvariantCulture, $"{r.Scan.Spill.BytesWritten}\t{r.Scan.Spill.BytesRead}\t{r.Scan.ChildWindows.FlushedBytes}\t{r.Scan.ChildWindows.Reopened}\t{r.Scan.StackPeak}\t{r.Lcp.Workers}\t{r.Parameters.LcpChunk}\t{r.Io.ScratchPeakBytes}\t{r.Lcp.B1Seconds:F3}\t{r.Lcp.B2Seconds:F3}\t{r.Lcp.B3Seconds:F3}\t{r.Lcp.CompareSteps}");
                break;
            }
            default:
                throw new InvalidOperationException("G3_VARIANT must be proto, prod or prodnr");
        }
        Prod.RssStats stats = rss.Stop();
        string row = string.Join('\t', label, variant, seconds.ToString("F4", CultureInfo.InvariantCulture), before, stats.MaxRssAnonKb, stats.Samples, arenaPeak, childArenaPeak, read, write, mapped, childBytes, extra);
        File.AppendAllText(rows, row + "\n");
        TestContext.Out.WriteLine(row);
        Assert.That(File.Exists(outPrefix + ".cld"), Is.True);
    }
}
