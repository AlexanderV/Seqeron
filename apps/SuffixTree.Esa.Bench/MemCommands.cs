using System.Diagnostics;
using System.Globalization;
using SuffixTree.Esa.Bench.Links;
using UkkonenTree = global::SuffixTree.SuffixTree;

namespace SuffixTree.Esa.Bench;

/// <summary>
/// Work package 1b.4a, prototype P2: end-to-end MEM of the 100 kbp query with 1 % substitutions (<c>queries/&lt;name&gt;.memq100k.txt</c> of the benchmark generator, the query of
/// <c>BENCHMARK.md</c> section 4) against the text, minLength 20, with the throwaway matching-statistics driver on the ESA. The first call is the cold-JIT time, the warm time is the
/// median of the following calls. With <c>--ukkonen true</c> the in-memory Ukkonen tree of the same text answers the same query in the same process (reference time) and its MEM set
/// is compared with the driver's.
/// </summary>
internal static class MemCommands
{
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;

    private static double Seconds(long start) => (Stopwatch.GetTimestamp() - start) / (double)Stopwatch.Frequency;

    internal static void Run(string[] args)
    {
        string? input = null, name = null, queries = null, synthetic = null;
        LinkKernel kernel = LinkKernel.RmqPsvNsv;
        EndpointSource source = EndpointSource.Psi;
        int loci = 2, calls = 10, minLength = 20, fanout = 16;
        bool ukkonen = true, timeLinks = true;
        for (int i = 0; i < args.Length; i += 2)
        {
            string value = i + 1 < args.Length ? args[i + 1] : throw new ArgumentException($"option {args[i]} needs a value");
            switch (args[i])
            {
                case "--input": input = value; break;
                case "--synthetic": synthetic = value; break;
                case "--name": name = value; break;
                case "--queries": queries = value; break;
                case "--kernel": kernel = LinksCommands.ParseKernel(value); break;
                case "--source": source = Enum.Parse<EndpointSource>(value, ignoreCase: true); break;
                case "--loci": loci = int.Parse(value, Inv); break;
                case "--calls": calls = int.Parse(value, Inv); break;
                case "--min-length": minLength = int.Parse(value, Inv); break;
                case "--fanout": fanout = int.Parse(value, Inv); break;
                case "--ukkonen": ukkonen = bool.Parse(value); break;
                case "--time-links": timeLinks = bool.Parse(value); break;
                default: throw new ArgumentException($"unknown option '{args[i]}'");
            }
        }
        if (input is null || queries is null) throw new ArgumentException("--input and --queries are required");
        if (calls < 6) throw new ArgumentException("--calls must be at least 6 (one cold call and five warm ones)");
        _ = synthetic;
        name ??= Path.GetFileNameWithoutExtension(input);

        string text = Pipeline.LoadText(input);
        string query = File.ReadAllText(Path.Combine(queries, "queries", $"{name}.memq100k.txt")).Trim();
        string label = $"{LinksCommands.KernelLabel(kernel)}_{source.ToString().ToLowerInvariant()}_loci{loci}";
        Emit.Row("mem", "text_chars", text.Length, "chars", name);
        Emit.Row("mem", "query_chars", query.Length, "chars", $"{name}.memq100k, minLength {minLength}, variant {label}");

        long t0 = Stopwatch.GetTimestamp();
        EsaIndex index = EsaIndex.Build(text, ChildLayoutKind.Tagged);
        Emit.Row("mem", "build_esa_s", Seconds(t0), "s", "EsaIndex.Build (tagged child table), cold JIT");
        int[]? isa = null, psi = null;
        long t1 = Stopwatch.GetTimestamp();
        isa = LinkTables.BuildIsa(index.SuffixArray);
        double isaSeconds = Seconds(t1);
        t1 = Stopwatch.GetTimestamp();
        psi = LinkTables.BuildPsi(index.SuffixArray, isa);
        double psiSeconds = Seconds(t1);
        t1 = Stopwatch.GetTimestamp();
        var helper = new LcpRangeHelper(index.Lcp, fanout);
        Emit.Row("mem", "build_nav_s", isaSeconds + psiSeconds + Seconds(t1), "s", $"ISA {isaSeconds.ToString("R", Inv)} s + Psi {psiSeconds.ToString("R", Inv)} s + pyramid {Seconds(t1).ToString("R", Inv)} s");
        var ctx = new LinkContext(index, source == EndpointSource.Isa ? isa : null, source == EndpointSource.Psi ? psi : null, helper);
        var driver = new MemDriver(index, ctx, kernel, source, loci);
        GC.Collect();

        IReadOnlyList<(int PositionInText, int PositionInQuery, int Length)>? result = null;
        var seconds = new List<double>();
        for (int call = 0; call < calls; call++)
        {
            long c0 = Stopwatch.GetTimestamp();
            result = driver.Run(query, minLength);
            double s = Seconds(c0);
            seconds.Add(s);
            Emit.Row("mem_esa", call == 0 ? "first_call_s" : $"call{call + 1}_s", s, "s", call == 0 ? "cold JIT included" : "warm");
        }
        double[] warm = seconds.Skip(1).Order().ToArray();
        double[] steady = seconds.Skip(calls - 5).Order().ToArray();
        Emit.Row("mem_esa", "warm_median_s", warm[warm.Length / 2], "s", $"median of {warm.Length} warm calls (2..{calls}) [{warm[0].ToString("R", Inv)} .. {warm[^1].ToString("R", Inv)}]");
        Emit.Row("mem_esa", "steady_median_s", steady[2], "s", $"median of the last 5 calls (after the tiered JIT settled) [{steady[0].ToString("R", Inv)} .. {steady[^1].ToString("R", Inv)}]");
        MemStats st = driver.Stats;
        Emit.Row("mem_esa", "matches", result!.Count, "count");
        Emit.Row("mem_esa", "checksum", Checksum(result), "count", "sum of r*31 + q*17 + len, as in BENCHMARK.md");
        Emit.Row("mem_esa", "starts", st.Starts, "count");
        Emit.Row("mem_esa", "links_followed", st.Links, "count", $"kernel calls ({loci} loci)");
        Emit.Row("mem_esa", "extend_steps", st.ExtendSteps, "count", "query symbols consumed by the extension");
        Emit.Row("mem_esa", "child_lookups", st.ChildLookups, "count", "extensions from a node");
        Emit.Row("mem_esa", "rescan_lookups", st.RescanLookups, "count", "child lookups of the rescans after the links");
        Emit.Row("mem_esa", "emitted", st.Emitted, "count");

        if (timeLinks)
        {
            driver.TimeLinks = true;
            long c0 = Stopwatch.GetTimestamp();
            driver.Run(query, minLength);
            double total = Seconds(c0);
            driver.TimeLinks = false;
            Emit.Row("mem_esa", "link_kernel_share", driver.LinkSeconds / total, "fraction", $"one extra warm call with a timer around every kernel call: {driver.LinkSeconds.ToString("R", Inv)} s of {total.ToString("R", Inv)} s (the 2 timestamps per link are inside both numbers)");
            Emit.Row("mem_esa", "link_kernel_ns_per_link", driver.LinkSeconds * 1e9 / Math.Max(1, driver.Stats.Links), "ns", "includes the timer floor");
        }
        Emit.Row("mem_esa", "peak_rss_hwm", Process.StatusBytes("VmHWM"), "bytes", "VmHWM before the Ukkonen reference is built");

        if (ukkonen)
        {
            long u0 = Stopwatch.GetTimestamp();
            UkkonenTree tree = UkkonenTree.Build(text);
            Emit.Row("mem_ukkonen", "build_s", Seconds(u0), "s", "in-memory Ukkonen tree, same process");
            IReadOnlyList<(int PositionInText, int PositionInQuery, int Length)>? reference = null;
            var useconds = new List<double>();
            for (int call = 0; call < calls; call++)
            {
                long c0 = Stopwatch.GetTimestamp();
                reference = tree.FindMaximalExactMatches(query, minLength);
                useconds.Add(Seconds(c0));
            }
            double[] uwarm = useconds.Skip(1).Order().ToArray();
            double[] usteady = useconds.Skip(calls - 5).Order().ToArray();
            Emit.Row("mem_ukkonen", "steady_median_s", usteady[2], "s", $"median of the last 5 calls [{usteady[0].ToString("R", Inv)} .. {usteady[^1].ToString("R", Inv)}]");
            Emit.Row("mem_ukkonen", "first_call_s", useconds[0], "s", "cold JIT of the Ukkonen MEM path in this process (the ESA driver ran first)");
            Emit.Row("mem_ukkonen", "warm_median_s", uwarm[uwarm.Length / 2], "s", $"median of {uwarm.Length} warm calls [{uwarm[0].ToString("R", Inv)} .. {uwarm[^1].ToString("R", Inv)}]");
            Emit.Row("mem_ukkonen", "matches", reference!.Count, "count");
            Emit.Row("mem_ukkonen", "checksum", Checksum(reference), "count");
            bool equal = SameSet(result, reference);
            Emit.Row("mem_check", "mem_set_equal", equal ? 1 : 0, "bool", "set of (position in text, position in query, length): ESA driver == Ukkonen FindMaximalExactMatches");
            if (!equal) Environment.ExitCode = 3;
            Emit.Row("mem", "peak_rss_hwm_final", Process.StatusBytes("VmHWM"), "bytes", "VmHWM after the Ukkonen tree");
            GC.KeepAlive(tree);
        }
    }

    private static long Checksum(IReadOnlyList<(int PositionInText, int PositionInQuery, int Length)> mems)
    {
        long ck = 0;
        foreach (var m in mems) ck += (m.PositionInText * 31L) + (m.PositionInQuery * 17L) + m.Length;
        return ck;
    }

    private static bool SameSet(IReadOnlyList<(int PositionInText, int PositionInQuery, int Length)> a, IReadOnlyList<(int PositionInText, int PositionInQuery, int Length)> b)
    {
        if (a.Count != b.Count) return false;
        var x = a.Order().ToArray();
        var y = b.Order().ToArray();
        return x.AsSpan().SequenceEqual(y);
    }
}
