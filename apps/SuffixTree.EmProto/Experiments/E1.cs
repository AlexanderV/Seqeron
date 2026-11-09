using System.Diagnostics;
using SuffixTree.EmProto.Build;
using SuffixTree.EmProto.Gen;
using SuffixTree.EmProto.Harness;
using SuffixTree.EmProto.Text;

namespace SuffixTree.EmProto.Experiments;

/// <summary>
/// E1 of ADR-004 section 5 (and, with <c>--dcs X</c>, E2): in-RAM simulation of the bucket partition and the bucket sort, suffix array only, one family per process. Bucket capacity K = n / 1000
/// (the 3x10^9 / 3 M ratio) unless <c>--k</c> is given, fan-out 64. Without a rank table the comparator is the direct suffix comparison with a vectorised LCE (E1); with one it is the
/// difference-cover comparator (E2). Correctness gate: the suffix array equals the committed libsais hash where there is one, otherwise the suffix array of <c>SaisKernel</c> on the same text.
/// </summary>
internal static class E1
{
    internal static int Run(Options o, Report r, string experiment)
    {
        string spec = o.Require("family");
        string data = o.Str("data", string.Empty);
        double limit = o.Double("limit", 600);
        int threads = o.Int("threads", 1);
        SymbolText text = FamilyLoader.Load(spec, data);
        long n = text.Length;
        long k = o.Long("k", Math.Max(2, n / o.Long("k-div", 1000)));
        int fan = o.Int("fan", 64);
        string caseName = o.Str("case", spec);

        long t0 = Stopwatch.GetTimestamp();
        int bits = Math.Max(1, PackedText.BitsFor(text.Sigma));
        using PackedText packed = PackedText.Pack(text.Symbols, bits);
        double prep = Elapsed(t0);

        using var dog = new Watchdog(limit, r, experiment, caseName);
        double dcsSeconds = 0;
        int dcsPeriod = o.Int("dcs", 0);
        DcsRankTable? ranks = null;
        var counters = new WorkCounters();
        Dcs.DifferenceCover? cover = null;
        if (dcsPeriod > 0)
        {
            cover = o.Flag("exact") ? Dcs.DifferenceCover.Create(dcsPeriod) : Dcs.DifferenceCover.Create(Dcs.DifferenceCover.NearestSingerPeriod(dcsPeriod));
            long d0 = Stopwatch.GetTimestamp();
            ranks = Dcs.DcsBuilder.Build(packed, cover, counters, o.Int("dcs-base", 16_000_000), out Dcs.DcsStats dst);
            dcsSeconds = Elapsed(d0);
            r.Row(experiment, caseName, "dcs_s", dcsSeconds, "s", $"period {cover.Period} ({cover.Kind}), |D|={cover.Residues.Length}, sample={dst.SampleSize}, unique names={dst.NamesUnique}, levels={dst.Levels}, sorted tie runs={dst.TieRuns}");
        }

        var options = new BuildOptions { LeafCapacity = k, FanOut = fan, Threads = threads, Seed = 1, Oversample = o.Int("oversample", 16), TieSort = o.Str("tie-sort", "merge") == "intro" ? TieSort.Intro : TieSort.Merge };
        using var sink = new MemorySaSink(n + 1);
        long b0 = Stopwatch.GetTimestamp();
        BuildStats st = EmBuilder.Run(packed, options, sink, ranks);
        double build = Elapsed(b0);
        counters.Add(st.Work);
        double total = prep + dcsSeconds + build;

        r.Row(experiment, caseName, "n", n, "symbols", $"sigma={text.Sigma}, b={bits} bits, K={k}, fan-out={fan}, threads={threads}, tie sort={options.TieSort}");
        r.Row(experiment, caseName, "prep_s", prep, "s", "pack the text (excluded from sa_s)");
        r.Row(experiment, caseName, "sa_s", build, "s", "PH2-PH3: partition + bucket sort, suffix array only");
        r.Row(experiment, caseName, "sa_plus_dcs_s", build + dcsSeconds, "s", "bucket build plus the rank table construction");
        r.Row(experiment, caseName, "total_s", total, "s", "pack + rank table + build");
        r.Row(experiment, caseName, "ns_per_symbol", build * 1e9 / n, "ns", "sa_s / n");
        r.Row(experiment, caseName, "buckets", st.Buckets, "n", $"leaves={st.Leaves}, partitions={st.Partitions}, max depth={st.MaxDepth}, largest leaf={st.MaxLeaf}, degenerate retries={st.Retries}");
        r.Row(experiment, caseName, "tie_share", 100.0 * st.Work.TiePositions / Math.Max(1, st.Work.LeafPositions), "%", $"positions of leaf buckets whose 64-bit key equals another's in the same bucket: {st.Work.TiePositions} of {st.Work.LeafPositions}");
        r.Row(experiment, caseName, "compares", st.Work.Compares, "n", $"comparator calls; from key ties: {st.Work.TiedCompares}; scan positions (all levels): {st.Work.ScanPositions}");
        r.Row(experiment, caseName, "lce_symbols", st.Work.LceSymbols, "symbols", $"LCE work: symbols matched beyond the key; LCE calls {st.Work.LceCalls}, loop steps {st.Work.LceSteps}, rank lookups {st.Work.RankLookups}");
        r.Row(experiment, caseName, "lce_symbols_per_n", (double)st.Work.LceSymbols / n, "symbols/symbol", "LCE work per text symbol");
        r.Row(experiment, caseName, "phase_cpu_s", st.PartitionSeconds, "s", $"summed over workers: partition {st.PartitionSeconds:F3}, leaf sort {st.LeafSeconds:F3}");
        r.Row(experiment, caseName, "vmhwm", ReadStatusKb("VmHWM:") * 1024.0 / 1e6, "MB", "peak resident set of the process (text, symbols, packed text, SA, buckets)");

        bool ok = true;
        if (!o.Flag("no-verify"))
        {
            string? vname = FamilyLoader.VectorName(spec);
            string source;
            if (vname is not null && o.Get("vectors") is { } vpath && VectorTable.Load(vpath).TryGetValue(vname, out string? expected))
            {
                string got = sink.HashLikeLibsais();
                ok = got == expected;
                source = $"libsais hash of {vname}";
            }
            else
            {
                int[] reference = Oracle.SaisKernelSa(text, out double refSeconds);
                ok = Oracle.Equal(sink, reference);
                source = "SaisKernel on the same text";
                r.Row(experiment, caseName, "ref_saiskernel_s", refSeconds, "s", "SaisKernel SA time on the same text (cold JIT, 1 thread), for scale");
            }
            r.Row(experiment, caseName, "sa_correct", ok ? "yes" : "NO", "", $"compared with {source}");
        }
        return ok ? 0 : 1;
    }

    private static double Elapsed(long t0) => (Stopwatch.GetTimestamp() - t0) / (double)Stopwatch.Frequency;

    internal static long ReadStatusKb(string key)
    {
        foreach (string line in File.ReadLines("/proc/self/status"))
        {
            if (!line.StartsWith(key, StringComparison.Ordinal)) continue;
            string[] parts = line.Split([' ', '\t'], StringSplitOptions.RemoveEmptyEntries);
            return long.Parse(parts[1], System.Globalization.CultureInfo.InvariantCulture);
        }
        return 0;
    }
}
