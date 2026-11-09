using System.Diagnostics;
using System.Globalization;
using SuffixTree.Esa.Bench.Links;

namespace SuffixTree.Esa.Bench;

/// <summary>
/// Work package 1b.4a, prototype P1: cost of the suffix-link kernels (a) outward LCP scan, (b) child-table descent, (c) RMQ + PSV/NSV, (c2) PSV/NSV on the ESA, for the
/// endpoint sources ISA and Psi. One process = one text: build the index, the navigation structures (timed, with bytes per character), then per (kernel, source) a deterministic
/// sample of internal nodes (steps with a budget, ns per link one by one) and sweeps of the worst nodes.
/// </summary>
internal static class LinksCommands
{
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;

    private sealed class LinkOptions
    {
        internal string Input { get; init; } = "";

        internal string? Synthetic { get; init; }

        internal string Name { get; init; } = "";

        internal int Sample { get; init; } = 200_000;

        internal ulong Seed { get; init; } = 1;

        internal int Worst { get; init; } = 1000;

        internal int Fanout { get; init; } = 16;

        internal LinkKernel[] Kernels { get; init; } = Enum.GetValues<LinkKernel>();

        internal EndpointSource[] Sources { get; init; } = Enum.GetValues<EndpointSource>();

        internal int Rotate { get; init; }

        internal long LinkBudget { get; init; } = 3_000_000;

        internal long TotalBudget { get; init; } = 400_000_000;

        internal double TimeCap { get; init; } = 4.0;

        internal int SearchCheck { get; init; } = 300;

        internal bool Sweeps { get; init; } = true;

        internal double Warmup { get; init; } = 1.0;

        internal static LinkOptions Parse(string[] args)
        {
            string? input = null, synthetic = null, name = null;
            int sample = 200_000, worst = 1000, fanout = 16, rotate = 0, searchCheck = 300;
            ulong seed = 1;
            long linkBudget = 3_000_000, totalBudget = 400_000_000;
            double timeCap = 4.0, warmup = 1.0;
            bool sweeps = true;
            LinkKernel[] kernels = Enum.GetValues<LinkKernel>();
            EndpointSource[] sources = Enum.GetValues<EndpointSource>();
            for (int i = 0; i < args.Length; i += 2)
            {
                string value = i + 1 < args.Length ? args[i + 1] : throw new ArgumentException($"option {args[i]} needs a value");
                switch (args[i])
                {
                    case "--input": input = value; break;
                    case "--synthetic": synthetic = value; break;
                    case "--name": name = value; break;
                    case "--sample": sample = int.Parse(value, Inv); break;
                    case "--seed": seed = ulong.Parse(value, Inv); break;
                    case "--worst": worst = int.Parse(value, Inv); break;
                    case "--fanout": fanout = int.Parse(value, Inv); break;
                    case "--rotate": rotate = int.Parse(value, Inv); break;
                    case "--link-budget": linkBudget = long.Parse(value, Inv); break;
                    case "--total-budget": totalBudget = long.Parse(value, Inv); break;
                    case "--time-cap": timeCap = double.Parse(value, Inv); break;
                    case "--search-check": searchCheck = int.Parse(value, Inv); break;
                    case "--sweeps": sweeps = bool.Parse(value); break;
                    case "--warmup": warmup = double.Parse(value, Inv); break;
                    case "--kernels": kernels = value.Equals("all", StringComparison.OrdinalIgnoreCase) ? kernels : value.Split(',').Select(ParseKernel).ToArray(); break;
                    case "--sources": sources = value.Equals("all", StringComparison.OrdinalIgnoreCase) ? sources : value.Split(',').Select(v => Enum.Parse<EndpointSource>(v, ignoreCase: true)).ToArray(); break;
                    default: throw new ArgumentException($"unknown option '{args[i]}'");
                }
            }
            if (input is null && synthetic is null) throw new ArgumentException("--input or --synthetic is required");
            return new LinkOptions
            {
                Input = input ?? "",
                Synthetic = synthetic,
                Name = name ?? (synthetic is not null ? synthetic.Split(':')[0] : Path.GetFileNameWithoutExtension(input!)),
                Sample = sample,
                Seed = seed,
                Worst = worst,
                Fanout = fanout,
                Kernels = kernels,
                Sources = sources,
                Rotate = rotate,
                LinkBudget = linkBudget,
                TotalBudget = totalBudget,
                TimeCap = timeCap,
                SearchCheck = searchCheck,
                Sweeps = sweeps,
                Warmup = warmup,
            };
        }
    }

    /// <summary>The short names of the kernels on the command line: scan, descent, rmq, psv (or the enum names).</summary>
    internal static LinkKernel ParseKernel(string name) => name.ToLowerInvariant() switch
    {
        "a" or "scan" => LinkKernel.Scan,
        "b" or "descent" => LinkKernel.Descent,
        "c" or "rmq" or "rmqpsvnsv" => LinkKernel.RmqPsvNsv,
        "c2" or "psv" or "psvnsv" => LinkKernel.PsvNsv,
        _ => throw new ArgumentException($"unknown kernel '{name}' (scan, descent, rmq, psv)"),
    };

    internal static string KernelLabel(LinkKernel k) => k switch
    {
        LinkKernel.Scan => "a_scan",
        LinkKernel.Descent => "b_descent",
        LinkKernel.RmqPsvNsv => "c_rmq",
        _ => "c2_psv",
    };

    private static double Seconds(long start) => (Stopwatch.GetTimestamp() - start) / (double)Stopwatch.Frequency;

    internal static void Run(string[] args)
    {
        LinkOptions o = LinkOptions.Parse(args);
        string text = o.Synthetic is not null ? QueryCommands.Synthetic(o.Synthetic) : Pipeline.LoadText(o.Input);
        int n = text.Length;
        Emit.Row("links", "input_chars", n, "chars", $"{o.Name}; sample {o.Sample} seed {o.Seed}; fanout {o.Fanout}; rotate {o.Rotate}");

        long t0 = Stopwatch.GetTimestamp();
        EsaIndex index = EsaIndex.Build(text, ChildLayoutKind.Tagged);
        Emit.Row("links", "build_total_s", Seconds(t0), "s", "EsaIndex.Build (tagged child table), cold JIT");
        int[] sa = index.SuffixArray;

        // navigation structures: build time and bytes per character are part of the cost of a link
        long t1 = Stopwatch.GetTimestamp();
        int[] isa = LinkTables.BuildIsa(sa);
        double isaSeconds = Seconds(t1);
        long t2 = Stopwatch.GetTimestamp();
        int[] psi = LinkTables.BuildPsi(sa, isa);
        double psiSeconds = Seconds(t2);
        long t3 = Stopwatch.GetTimestamp();
        var helper = new LcpRangeHelper(index.Lcp, o.Fanout);
        double rmqSeconds = Seconds(t3);
        Emit.Row("nav", "build_isa_s", isaSeconds, "s", "ISA from SA: one scattered pass (cold JIT)");
        Emit.Row("nav", "build_psi_s", psiSeconds, "s", "Psi from SA and ISA: one gather pass (add build_isa_s for the cost from SA alone)");
        Emit.Row("nav", "build_rmq_s", rmqSeconds, "s", $"block-minima pyramid, fanout {o.Fanout}, {helper.Levels} levels: one sequential pass over the compact LCP");
        Emit.Row("nav", "isa_bytes_per_char", 4.0 * isa.Length / n, "bytes/char");
        Emit.Row("nav", "psi_bytes_per_char", 4.0 * psi.Length / n, "bytes/char");
        Emit.Row("nav", "rmq_bytes_per_char", helper.PayloadBytes / (double)n, "bytes/char", $"pyramid only; the LCP array ({index.Lcp.PayloadBytes / (double)n:F3} B/char, kind {index.Lcp.Kind}) is shared with the index");
        Emit.Row("nav", "lcp_kind", index.Lcp.Kind.ToString(), "kind", $"escapes {index.Lcp.OverflowCount}");

        var ctx = new LinkContext(index, isa, psi, helper);
        t1 = Stopwatch.GetTimestamp();
        LinkNodeSet nodes = LinkNodeSet.Enumerate(index.Lcp);
        Emit.Row("nodes", "internal_nodes", nodes.Count, "count", $"enumerate_s {Seconds(t1).ToString("R", Inv)}; header statistic {index.Statistics.InternalCount}");
        if (nodes.Count != index.Statistics.InternalCount) Environment.ExitCode = 3;
        if (nodes.Count == 0) return;

        double timerFloor = CalibrateTimer();
        Emit.Row("timer", "per_link_timer_floor_ns", timerFloor, "ns", "median of back-to-back Stopwatch.GetTimestamp pairs; the per-link medians below are raw (include it) and the *_net rows subtract it");

        // samples: A (per-link timing) and B (batch time), disjoint, plus a small warm-up set
        int[] drawn = nodes.Sample(2 * o.Sample + 4000, o.Seed);
        int half = Math.Min(o.Sample, Math.Max(1, (drawn.Length - 4000) / 2));
        int[] sampleA = drawn[..half], sampleB = drawn[half..(2 * half)], warm = drawn[^Math.Min(4000, drawn.Length)..];
        Emit.Row("nodes", "sample_nodes", sampleA.Length, "count", $"sample A (per-link) and B (batch) of {sampleA.Length} each, warm-up {warm.Length}");

        var combos = new List<(LinkKernel, EndpointSource)>();
        foreach (LinkKernel k in o.Kernels)
            foreach (EndpointSource s in o.Sources)
                combos.Add((k, s));
        if (combos.Count > 0)
        {
            int rot = o.Rotate % combos.Count;
            combos = [.. combos.Skip(rot), .. combos.Take(rot)];
        }

        bool ok = true;

        // 0. cold: the very first link of each kernel in this process (JIT of the kernel, the delegate table and the shared readers), before anything else ran
        foreach (var (kernel, source) in combos)
        {
            LinkNode first = nodes.Node(warm[0], sa);
            var none = new NoCount();
            LinkFunc<NoCount> f = LinkKernels.Get<NoCount>(kernel, source);
            long c0 = Stopwatch.GetTimestamp();
            f(ctx, first, ref none, out _);
            Emit.Row($"sample_{KernelLabel(kernel)}_{source.ToString().ToLowerInvariant()}", "first_link_ns", (Stopwatch.GetTimestamp() - c0) * 1e9 / Stopwatch.Frequency, "ns", "very first call of this kernel in the process (JIT + first touch); the first kernel also pays the JIT of the shared readers");
        }

        // warm-up: all kernels round robin on nodes that every kernel finishes cheaply, long enough for the tiered JIT to promote the shared readers (CompactLcp) to their final code
        var warmNodes = new Dictionary<(LinkKernel, EndpointSource), LinkNode[]>();
        foreach (var combo in combos)
        {
            var list = new List<LinkNode>();
            LinkFunc<BudgetSteps> probe = LinkKernels.Get<BudgetSteps>(combo.Item1, combo.Item2);
            foreach (int id in warm)
            {
                if (list.Count >= 400) break;
                LinkNode v = nodes.Node(id, sa);
                var counter = new BudgetSteps(20_000);
                if (probe(ctx, v, ref counter, out _)) list.Add(v);
            }
            warmNodes[combo] = [.. list];
        }
        long w0 = Stopwatch.GetTimestamp();
        long warmLinks = 0;
        while (Seconds(w0) < o.Warmup)
        {
            foreach (var combo in combos)
            {
                LinkFunc<NoCount> f = LinkKernels.Get<NoCount>(combo.Item1, combo.Item2);
                var none = new NoCount();
                foreach (LinkNode v in warmNodes[combo])
                {
                    f(ctx, v, ref none, out _);
                    warmLinks++;
                }
            }
            if (combos.Count == 0 || warmNodes.Values.All(a => a.Length == 0)) break;
        }
        Emit.Row("links", "warmup_links", warmLinks, "count", $"round robin over all kernels for at least {o.Warmup} s before the first timed pass");

        // 1. the sample, every combination
        var reference = new Dictionary<string, (int[] Lb, int[] Rb)>();
        foreach (var (kernel, source) in combos)
            ok &= Measure($"sample_{KernelLabel(kernel)}_{source.ToString().ToLowerInvariant()}", kernel, source, ctx, nodes, sa, sampleA, sampleB, warmNodes[(kernel, source)], o, timerFloor, reference);

        // 1b. drift check: the first combination again, after all the others
        if (combos.Count > 0)
        {
            var (k0, s0) = combos[0];
            var scratch = new Dictionary<string, (int[] Lb, int[] Rb)>();
            ok &= Measure($"repeat_{KernelLabel(k0)}_{s0.ToString().ToLowerInvariant()}", k0, s0, ctx, nodes, sa, sampleA, sampleB, warmNodes[(k0, s0)], o, timerFloor, scratch);
        }

        // 2. definition check of a part of the sample by the pattern search of the index
        if (o.SearchCheck > 0 && reference.Count > 0)
            ok &= SearchCheck(text, index, nodes, sampleA, reference.Values.First(), o.SearchCheck);

        // 3. worst nodes
        if (o.Sweeps && o.Worst > 0)
            ok &= Sweeps(ctx, nodes, sa, isa, psi, combos, o, timerFloor);

        Emit.Row("links", "peak_rss_hwm", Process.StatusBytes("VmHWM"), "bytes", "VmHWM of the process");
        if (!ok) Environment.ExitCode = 3;
    }

    private static double CalibrateTimer()
    {
        const int N = 400_000;
        var d = new long[N];
        for (int round = 0; round < 2; round++)
            for (int i = 0; i < N; i++)
            {
                long a = Stopwatch.GetTimestamp();
                long b = Stopwatch.GetTimestamp();
                d[i] = b - a;
            }
        Array.Sort(d);
        return d[N / 2] * 1e9 / Stopwatch.Frequency;
    }

    /// <summary>The nearest-rank percentile of an ascending array.</summary>
    private static double Percentile(long[] sorted, double p)
    {
        if (sorted.Length == 0) return double.NaN;
        int rank = (int)Math.Ceiling(p * sorted.Length) - 1;
        return sorted[Math.Clamp(rank, 0, sorted.Length - 1)];
    }

    /// <summary>Counting pass (steps, results), then cold first call, warm-up, per-link timing on A and batch timing on B for one combination.</summary>
    private static bool Measure(
        string op, LinkKernel kernel, EndpointSource source, LinkContext ctx, LinkNodeSet nodes, int[] sa, int[] sampleA, int[] sampleB, LinkNode[] warm, LinkOptions o, double timerFloor,
        Dictionary<string, (int[] Lb, int[] Rb)> reference)
    {
        LinkFunc<NoCount> timed = LinkKernels.Get<NoCount>(kernel, source);
        LinkFunc<BudgetSteps> counting = LinkKernels.Get<BudgetSteps>(kernel, source);

        // counting pass: steps per link, results for the cross-check; the budget turns an unbounded kernel into a measured "budget exceeded"
        int count = sampleA.Length;
        var steps = new long[count];
        var lbs = new int[count];
        var rbs = new int[count];
        int done = 0, exceeded = 0;
        long totalSteps = 0;
        long cp = Stopwatch.GetTimestamp();
        for (int i = 0; i < count; i++)
        {
            LinkNode v = nodes.Node(sampleA[i], sa);
            var counter = new BudgetSteps(o.LinkBudget);
            bool okLink = counting(ctx, v, ref counter, out LinkTarget t);
            steps[i] = counter.Steps;
            totalSteps += counter.Steps;
            lbs[i] = okLink ? t.Lb : -1;
            rbs[i] = okLink ? t.Rb : -1;
            if (!okLink) exceeded++;
            done = i + 1;
            if (totalSteps > o.TotalBudget) break;
        }
        Emit.Row(op, "counting_pass_s", Seconds(cp), "s", "steps and results of the sample (BudgetSteps instantiation)");
        EmitSteps(op, steps.AsSpan(0, done), done, count, exceeded, o);

        // cross-check of this combination with the first one that ran
        bool ok = CrossCheck(op, lbs, rbs, done, reference);

        // warm-up on a separate set, then the per-link timing
        for (int pass = 0; pass < 2; pass++)
            foreach (LinkNode v in warm)
            {
                var w = new NoCount();
                timed(ctx, v, ref w, out _);
            }

        int limit = done;           // do not time past what the counting pass could finish within its budget
        var ns = new long[limit];
        int timedCount = 0;
        long sink = 0;
        long start = Stopwatch.GetTimestamp();
        for (int i = 0; i < limit; i++)
        {
            if (steps[i] > o.LinkBudget) continue;      // budget exceeded: not timed
            LinkNode v = nodes.Node(sampleA[i], sa);
            var c = new NoCount();
            long a = Stopwatch.GetTimestamp();
            timed(ctx, v, ref c, out LinkTarget t);
            long b = Stopwatch.GetTimestamp();
            ns[timedCount++] = b - a;
            sink += t.Lb;
            if ((i & 1023) == 0 && Seconds(start) > o.TimeCap) break;
        }
        long[] raw = ns.AsSpan(0, timedCount).ToArray();
        Array.Sort(raw);
        double toNs = 1e9 / Stopwatch.Frequency;
        Emit.Row(op, "timed_links", timedCount, "count", $"of {count} sampled; time cap {o.TimeCap} s");
        if (timedCount > 0)
        {
            Emit.Row(op, "median_ns", Percentile(raw, 0.5) * toNs, "ns", "raw per-link timing (includes the timer floor)");
            Emit.Row(op, "p99_ns", Percentile(raw, 0.99) * toNs, "ns", "raw");
            Emit.Row(op, "max_ns", raw[^1] * toNs, "ns", "raw");
            Emit.Row(op, "median_net_ns", Math.Max(0, (Percentile(raw, 0.5) * toNs) - timerFloor), "ns", "median minus the timer floor");
            Emit.Row(op, "p99_net_ns", Math.Max(0, (Percentile(raw, 0.99) * toNs) - timerFloor), "ns", "p99 minus the timer floor");
        }

        // batch: two timestamps around the whole loop on the other sample (no per-link timer, no cache reuse of sample A)
        int bCount = 0;
        long bStart = Stopwatch.GetTimestamp();
        for (int i = 0; i < sampleB.Length; i++)
        {
            LinkNode v = nodes.Node(sampleB[i], sa);
            var c = new NoCount();
            timed(ctx, v, ref c, out LinkTarget t);
            sink += t.Rb;
            bCount++;
            if ((i & 1023) == 0 && Seconds(bStart) > o.TimeCap) break;
        }
        Emit.Row(op, "mean_batch_ns", Seconds(bStart) * 1e9 / bCount, "ns", $"mean over {bCount} links of sample B, no per-link timer");
        GC.KeepAlive(sink);
        return ok;
    }

    private static void EmitSteps(string op, ReadOnlySpan<long> steps, int done, int count, int exceeded, LinkOptions o)
    {
        long[] sorted = steps.ToArray();
        Array.Sort(sorted);
        Emit.Row(op, "steps_links", done, "count", $"counted of {count} sampled (total step budget {o.TotalBudget})");
        Emit.Row(op, "steps_budget_exceeded", exceeded, "count", $"links that needed more than {o.LinkBudget} steps; the sample stops at the total budget");
        if (sorted.Length == 0) return;
        Emit.Row(op, "steps_median", Percentile(sorted, 0.5), "steps", "array entries read");
        Emit.Row(op, "steps_p99", Percentile(sorted, 0.99), "steps");
        Emit.Row(op, "steps_max", sorted[^1], "steps");
        Emit.Row(op, "steps_mean", sorted.Average(), "steps");
    }

    private static bool CrossCheck(string op, int[] lbs, int[] rbs, int done, Dictionary<string, (int[] Lb, int[] Rb)> reference)
    {
        bool ok = true;
        if (reference.Count == 0)
        {
            reference[op] = (lbs, rbs);
            Emit.Row(op, "cross_check", "reference", "label", "first combination of this process: the others are compared with it");
            return true;
        }
        var (rl, rr) = reference.Values.First();
        int compared = 0, different = 0;
        for (int i = 0; i < done; i++)
        {
            if (lbs[i] < 0 || rl[i] < 0) continue;
            compared++;
            if (lbs[i] != rl[i] || rbs[i] != rr[i]) different++;
        }
        ok = different == 0;
        Emit.Row(op, "cross_check_different", different, "count", $"of {compared} comparable links (both finished) against {reference.Keys.First()}");
        return ok;
    }

    /// <summary>The links of a part of the sample against the pattern search of the index: the target is the rank interval of the label without its first symbol.</summary>
    private static bool SearchCheck(string text, EsaIndex index, LinkNodeSet nodes, int[] sampleA, (int[] Lb, int[] Rb) reference, int count)
    {
        int bad = 0, checkedNodes = 0;
        for (int i = 0; i < sampleA.Length && checkedNodes < count; i++)
        {
            if (reference.Lb[i] < 0) continue;
            int id = sampleA[i];
            int depth = nodes.Depth[id];
            int pos = index.SuffixArray[nodes.Lb[id]];
            int lo, hi;
            if (depth == 1)
            {
                lo = 0;
                hi = index.Length;
            }
            else if (!index.TryMatch(text.AsSpan(pos + 1, depth - 1), out lo, out hi))
            {
                bad++;
                checkedNodes++;
                continue;
            }
            checkedNodes++;
            if (lo != reference.Lb[i] || hi != reference.Rb[i]) bad++;
        }
        Emit.Row("search_check", "nodes_checked", checkedNodes, "count", "link target == TryMatch(label without its first symbol)");
        Emit.Row("search_check", "mismatches", bad, "count");
        return bad == 0;
    }

    /// <summary>
    /// One untimed pass of the cheapest kernel over all nodes gives the exact scan cost and the target width of every node; the worst nodes by scan cost, target width, string
    /// depth and node depth are then run through every combination (steps with the budget, per-link timing).
    /// </summary>
    private static bool Sweeps(LinkContext ctx, LinkNodeSet nodes, int[] sa, int[] isa, int[] psi, List<(LinkKernel, EndpointSource)> combos, LinkOptions o, double timerFloor)
    {
        _ = isa;
        long p0 = Stopwatch.GetTimestamp();
        var scanCost = new int[nodes.Count];
        var width = new int[nodes.Count];
        LinkFunc<NoCount> psv = LinkKernels.Get<NoCount>(LinkKernel.PsvNsv, EndpointSource.Psi);
        var none = new NoCount();
        for (int i = 0; i < nodes.Count; i++)
        {
            LinkNode v = nodes.Node(i, sa);
            psv(ctx, v, ref none, out LinkTarget t);
            if (v.Depth <= 1) continue;     // a depth-1 node links to the root: no search at all
            width[i] = t.Rb - t.Lb + 1;
            scanCost[i] = (psi[v.Lb] - t.Lb) + (t.Rb - psi[v.Rb]);
        }
        Emit.Row("sweeps", "full_pass_s", Seconds(p0), "s", "kernel c2 over every internal node (untimed per link): target width and exact scan cost of each node");
        Emit.Row("sweeps", "max_scan_cost", scanCost.Max(), "steps", "largest number of LCP entries the outward scan has to read beyond [Psi(l), Psi(r)] on any node of this text");
        Emit.Row("sweeps", "max_target_width", width.Max(), "ranks");
        Emit.Row("sweeps", "max_string_depth", nodes.Depth.Max(), "chars");
        Emit.Row("sweeps", "max_node_depth", nodes.NodeDepth.Max(), "edges");

        var lists = new (string Name, int[] Ids)[]
        {
            ("widest_scan", LinkNodeSet.Largest(scanCost, o.Worst)),
            ("widest_target", LinkNodeSet.Largest(width, o.Worst)),
            ("deepest_string", LinkNodeSet.Largest(nodes.Depth, o.Worst)),
            ("deepest_node", LinkNodeSet.Largest(nodes.NodeDepth, o.Worst)),
        };
        var empty = new Dictionary<string, (int[] Lb, int[] Rb)>();
        bool ok = true;
        foreach (var (name, ids) in lists)
        {
            Emit.Row($"worst_{name}", "nodes", ids.Length, "count", "the nodes with the largest key (see sweeps rows for the maxima)");
            var reference = new Dictionary<string, (int[] Lb, int[] Rb)>();
            foreach (var (kernel, source) in combos)
                ok &= Measure($"worst_{name}_{KernelLabel(kernel)}_{source.ToString().ToLowerInvariant()}", kernel, source, ctx, nodes, sa, ids, ids[..Math.Min(ids.Length, 200)], [], o, timerFloor, reference);
        }
        _ = empty;
        return ok;
    }
}
