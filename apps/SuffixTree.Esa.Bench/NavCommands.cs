using System.Diagnostics;
using System.Globalization;
using SuffixTree.Esa.Bench.Links;
using UkkonenTree = global::SuffixTree.SuffixTree;

namespace SuffixTree.Esa.Bench;

using Triple = (int PositionInText, int PositionInQuery, int Length);

/// <summary>
/// Work package 1b.4: the production navigator behind the generic algorithms. One process = one text and one operation (<c>--op</c>): <c>mem</c>, <c>mum</c>, <c>anchors</c> (the benchmark's
/// <c>&lt;name&gt;.memq100k</c> query, 100 kbp with 1 % substitutions, <c>--min-length</c> 20), <c>lcs</c> (the LCS family against a second text), <c>links</c> (ns per production link against
/// the prototype kernel (c2) of 1b.4a on the same internal nodes) or <c>split</c> (MEM through a counting and timing navigator: where the time goes). The first call of an operation is the cold-JIT
/// time, <c>warm_median_s</c> is the median over calls 2..N, <c>steady_median_s</c> over the last five; with <c>--ukkonen true</c> (the default) the in-memory Ukkonen tree of the same text
/// answers the same query in the same process (reference time) and its answer is compared with the ESA's (exit code 3 on a difference).
/// </summary>
internal static class NavCommands
{
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;

    /// <summary>Always 0 at run time, but a plain static field: the JIT cannot fold it, so a value combined with it is a real data dependency.</summary>
    private static int s_zero = int.Parse("0", CultureInfo.InvariantCulture);

    private static double Seconds(long start) => (Stopwatch.GetTimestamp() - start) / (double)Stopwatch.Frequency;

    internal static void Run(string[] args)
    {
        string? input = null, name = null, queries = null, op = null, other = null;
        string otherKind = "mutated";
        int calls = 10, minLength = 20, sample = 200_000;
        bool ukkonen = true;
        ulong seed = 1;
        for (int i = 0; i < args.Length; i += 2)
        {
            string value = i + 1 < args.Length ? args[i + 1] : throw new ArgumentException($"option {args[i]} needs a value");
            switch (args[i])
            {
                case "--input": input = value; break;
                case "--name": name = value; break;
                case "--queries": queries = value; break;
                case "--op": op = value; break;
                case "--other": other = value; break;
                case "--other-kind": otherKind = value; break;
                case "--calls": calls = int.Parse(value, Inv); break;
                case "--min-length": minLength = int.Parse(value, Inv); break;
                case "--sample": sample = int.Parse(value, Inv); break;
                case "--seed": seed = ulong.Parse(value, Inv); break;
                case "--ukkonen": ukkonen = bool.Parse(value); break;
                default: throw new ArgumentException($"unknown option '{args[i]}'");
            }
        }
        if (input is null || op is null) throw new ArgumentException("--input and --op are required");
        if (calls < 6) throw new ArgumentException("--calls must be at least 6 (one cold call and five warm ones)");
        name ??= Path.GetFileNameWithoutExtension(input);
        string text = Pipeline.LoadText(input);
        Emit.Row("nav", "text_chars", text.Length, "chars", $"{name}, op {op}");

        long t0 = Stopwatch.GetTimestamp();
        EsaIndex index = EsaIndex.Build(text);
        Emit.Row("nav", "build_esa_s", Seconds(t0), "s", $"EsaIndex.Build (child layout {index.Child.Kind}), cold JIT");
        Emit.Row("nav", "pack_built_before", index.NavigationBuilt ? 1 : 0, "bool", "the navigation pack must not exist before the first matching call");

        if (op == "links")
        {
            Links(index, name, sample, seed);
            return;
        }

        string query = "";
        if (op is "mem" or "mum" or "anchors" or "split")
        {
            if (queries is null) throw new ArgumentException("--queries <bench dir> is required for this operation");
            query = File.ReadAllText(Path.Combine(queries, "queries", $"{name}.memq100k.txt")).Trim();
            Emit.Row("nav", "query_chars", query.Length, "chars", $"{name}.memq100k, minLength {minLength}");
        }

        // the pack is built by the first operation; time it on its own so that the first call is not charged with it
        long p0 = Stopwatch.GetTimestamp();
        EsaNavigationPack pack = index.Navigation;
        double packSeconds = Seconds(p0);
        Emit.Row("nav", "pack_build_s", packSeconds, "s", "ISA + LCP pyramid (fanout 32) + pinned views, cold JIT; paid once by the first suffix-link user");
        Emit.Row("nav", "pack_bytes_per_char", pack.PayloadBytes / (double)text.Length, "B/char", $"ISA {4.0 * (index.Length + 1) / text.Length:F3} + pyramid {pack.Pyramid.LevelBytes / (double)text.Length:F3} ({pack.Pyramid.Levels} levels)");
        using (LcpPyramid again = LcpPyramid.Build(index.Lcp))
        {
            long q0 = Stopwatch.GetTimestamp();
            using LcpPyramid timed = LcpPyramid.Build(index.Lcp);
            Emit.Row("nav", "pyramid_build_s", Seconds(q0), "s", $"second build, warm ({again.Levels} levels)");
        }
        GC.Collect();

        // the ESA runs first, with nothing else on the heap; the Ukkonen tree (about 1.9 GB of small objects) is built afterwards: with it alive the gen-2 collections that the large result lists
        // provoke (LOH allocations) walk its whole object graph and would be charged to the ESA
        var steps = new List<Step>();
        switch (op)
        {
            case "mem":
                steps.Add(new Step("mem", () => index.FindMaximalExactMatches(query, minLength), t => t.FindMaximalExactMatches(query, minLength), (a, b) => SameTriples((IReadOnlyList<Triple>)a, (IReadOnlyList<Triple>)b)));
                break;
            case "mum":
                steps.Add(new Step("mum", () => index.FindMaximalUniqueMatches(query, minLength), t => t.FindMaximalUniqueMatches(query, minLength), (a, b) => SameTriples((IReadOnlyList<Triple>)a, (IReadOnlyList<Triple>)b)));
                break;
            case "anchors":
                steps.Add(new Step("anchors", () => index.FindExactMatchAnchors(query, minLength), t => t.FindExactMatchAnchors(query, minLength), (a, b) => SameAnchors((IReadOnlyList<Triple>)a, (IReadOnlyList<Triple>)b)));
                break;
            case "lcs":
                string second = other is not null ? Pipeline.LoadText(other) : MakeOther(text, otherKind, seed);
                Emit.Row("nav", "other_chars", second.Length, "chars", $"the second text ({(other is not null ? other : otherKind)})");
                steps.Add(new Step("lcs_info", () => index.LongestCommonSubstringInfo(second), t => t.LongestCommonSubstringInfo(second), (a, b) => SameInfo(((string, int, int))a, ((string, int, int))b)));
                steps.Add(new Step("lcs_all", () => index.FindAllLongestCommonSubstrings(second), t => t.FindAllLongestCommonSubstrings(second), (a, b) => SameAll(((string, IReadOnlyList<int>, IReadOnlyList<int>))a, ((string, IReadOnlyList<int>, IReadOnlyList<int>))b)));
                steps.Add(new Step("lcs_distinct", () => index.FindAllDistinctLongestCommonSubstrings(second), t => t.FindAllDistinctLongestCommonSubstrings(second), (a, b) => SameDistinct((IReadOnlyList<(string Substring, IReadOnlyList<int> PositionsInText, IReadOnlyList<int> PositionsInOther)>)a, (IReadOnlyList<(string Substring, IReadOnlyList<int> PositionsInText, IReadOnlyList<int> PositionsInOther)>)b)));
                break;
            case "split":
                Split(index, query, minLength, calls);
                break;
            default:
                throw new ArgumentException($"unknown --op '{op}' (mem, mum, anchors, lcs, links, split)");
        }

        var esaResults = new List<object>();
        foreach (Step step in steps) esaResults.Add(RunCalls("nav_esa", step.Label, calls, step.Esa));
        UkkonenTree? tree = null;
        if (ukkonen && steps.Count > 0)
        {
            tree = Build(text);
            for (int i = 0; i < steps.Count; i++)
            {
                Step step = steps[i];
                object expected = RunCalls("nav_ukkonen", step.Label, calls, () => step.Ukkonen(tree));
                bool equal = step.Same(esaResults[i], expected);
                Emit.Row("nav_check", $"{step.Label}_equal", equal ? 1 : 0, "bool", "the ESA answer equals the Ukkonen tree's");
                if (!equal) Environment.ExitCode = 3;
            }
        }
        Emit.Row("nav", "peak_rss_hwm", Process.StatusBytes("VmHWM"), "bytes", "VmHWM of the process at the end");
        GC.KeepAlive(tree);
        GC.KeepAlive(pack);
    }

    private static UkkonenTree Build(string text)
    {
        long u0 = Stopwatch.GetTimestamp();
        UkkonenTree tree = UkkonenTree.Build(text);
        Emit.Row("nav_ukkonen", "build_s", Seconds(u0), "s", "in-memory Ukkonen tree, same process");
        return tree;
    }

    private sealed record Step(string Label, Func<object> Esa, Func<UkkonenTree, object> Ukkonen, Func<object, object, bool> Same);

    /// <summary>Calls <paramref name="call"/> <paramref name="calls"/> times, reports the cold first call, the warm median and the steady median, and returns the last answer.</summary>
    private static object RunCalls(string group, string label, int calls, Func<object> call)
    {
        object? result = null;
        var seconds = new List<double>();
        for (int i = 0; i < calls; i++)
        {
            long c0 = Stopwatch.GetTimestamp();
            result = call();
            seconds.Add(Seconds(c0));
        }
        Report(group, label, seconds, result);
        return result!;
    }

    private static void Report(string group, string label, List<double> seconds, object? result)
    {
        double[] warm = seconds.Skip(1).Order().ToArray();
        double[] steady = seconds.Skip(seconds.Count - 5).Order().ToArray();
        Emit.Row(group, $"{label}_first_call_s", seconds[0], "s", "cold JIT included");
        Emit.Row(group, $"{label}_warm_median_s", warm[warm.Length / 2], "s", $"median of {warm.Length} warm calls (2..{seconds.Count}) [{warm[0].ToString("R", Inv)} .. {warm[^1].ToString("R", Inv)}]");
        Emit.Row(group, $"{label}_steady_median_s", steady[2], "s", $"median of the last 5 calls [{steady[0].ToString("R", Inv)} .. {steady[^1].ToString("R", Inv)}]");
        switch (result)
        {
            case IReadOnlyList<Triple> triples:
                Emit.Row(group, $"{label}_results", triples.Count, "count");
                long ck = 0;
                foreach (Triple m in triples) ck += (m.PositionInText * 31L) + (m.PositionInQuery * 17L) + m.Length;
                Emit.Row(group, $"{label}_checksum", ck, "count", "sum of r*31 + q*17 + len, as in BENCHMARK.md");
                break;
            case ValueTuple<string, int, int> info:
                Emit.Row(group, $"{label}_length", info.Item1.Length, "chars");
                break;
            case ValueTuple<string, IReadOnlyList<int>, IReadOnlyList<int>> all:
                Emit.Row(group, $"{label}_length", all.Item1.Length, "chars", $"{all.Item2.Count} positions in the text, {all.Item3.Count} in the other text");
                break;
            case IReadOnlyList<(string Substring, IReadOnlyList<int> PositionsInText, IReadOnlyList<int> PositionsInOther)> distinct:
                Emit.Row(group, $"{label}_ties", distinct.Count, "count", distinct.Count > 0 ? $"length {distinct[0].Substring.Length}" : "no common character");
                break;
            default:
                break;
        }
    }

    private static bool SameTriples(IReadOnlyList<Triple> a, IReadOnlyList<Triple> b) => a.Count == b.Count && a.SequenceEqual(b);

    private static bool SameAnchors(IReadOnlyList<Triple> a, IReadOnlyList<Triple> b)
        => a.Count == b.Count && a.Select(x => (x.PositionInQuery, x.Length)).SequenceEqual(b.Select(x => (x.PositionInQuery, x.Length)));

    private static bool SameInfo((string Substring, int PositionInText, int PositionInOther) a, (string Substring, int PositionInText, int PositionInOther) b)
        => a.Substring == b.Substring && a.PositionInOther == b.PositionInOther;

    private static bool SameAll((string Substring, IReadOnlyList<int> PositionsInText, IReadOnlyList<int> PositionsInOther) a, (string Substring, IReadOnlyList<int> PositionsInText, IReadOnlyList<int> PositionsInOther) b)
        => a.Substring == b.Substring && a.PositionsInText.SequenceEqual(b.PositionsInText) && a.PositionsInOther.SequenceEqual(b.PositionsInOther);

    private static bool SameDistinct(IReadOnlyList<(string Substring, IReadOnlyList<int> PositionsInText, IReadOnlyList<int> PositionsInOther)> a, IReadOnlyList<(string Substring, IReadOnlyList<int> PositionsInText, IReadOnlyList<int> PositionsInOther)> b)
        => a.Count == b.Count && a.Zip(b).All(p => SameAll(p.First, p.Second));

    /// <summary>The second text of the LCS measurement, deterministic: <c>random</c> (independent, same letters) or <c>mutated</c> (a copy of the text with 1 % substitutions and 0.1 % single-letter insertions and deletions).</summary>
    private static string MakeOther(string text, string kind, ulong seed)
    {
        ulong state = seed * 0x9E3779B97F4A7C15UL;
        ulong Next()
        {
            unchecked
            {
                state += 0x9E3779B97F4A7C15UL;
                ulong z = state;
                z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9UL;
                z = (z ^ (z >> 27)) * 0x94D049BB133111EBUL;
                return z ^ (z >> 31);
            }
        }
        char[] letters = text.Distinct().Order().ToArray();
        var sb = new System.Text.StringBuilder(text.Length + 16);
        if (kind == "random")
        {
            for (int i = 0; i < text.Length; i++) sb.Append(letters[(int)(Next() % (ulong)letters.Length)]);
            return sb.ToString();
        }
        if (kind != "mutated") throw new ArgumentException($"unknown --other-kind '{kind}' (random, mutated)");
        foreach (char c in text)
        {
            ulong r = Next() % 10_000;
            if (r < 100) sb.Append(letters[(int)(Next() % (ulong)letters.Length)]);
            else if (r >= 105 && r < 110) sb.Append(c).Append(letters[(int)(Next() % (ulong)letters.Length)]);
            else if (r >= 110) sb.Append(c);
        }
        return sb.ToString();
    }

    // ---------------------------------------------------------------- links: production kernel against the prototype

    private static void Links(EsaIndex index, string name, int sampleSize, ulong seed)
    {
        EsaNavigationPack pack = index.Navigation;
        int[] isa = LinkTables.BuildIsa(index.SuffixArray);
        var helper = new LcpRangeHelper(index.Lcp, 32);
        var ctx = new LinkContext(index, isa, null, helper);
        var set = LinkNodeSet.Enumerate(index.Lcp);
        int[] sa = index.SuffixArray;
        int[] sample = set.Sample(sampleSize, seed);
        Emit.Row("nav_links", "internal_nodes", set.Count, "count", $"{name}; the prototype is timed WITH the extra SA[Lb] read that completes the production handle (Pos), which the 1b.4a figures did not include");
        Emit.Row("nav_links", "sample", sample.Length, "count", $"seed {seed}");

        // timer floor: two back-to-back timestamps, as in P1
        long floorTicks = 0;
        for (int i = 0; i < 1_000_000; i++)
        {
            long a = Stopwatch.GetTimestamp();
            long b = Stopwatch.GetTimestamp();
            floorTicks += b - a;
        }
        double floorNs = floorTicks * 1e9 / Stopwatch.Frequency / 1_000_000;
        Emit.Row("nav_links", "timer_floor_ns", floorNs, "ns", "mean of an empty timed interval");

        var nodes = new EsaNode[sample.Length];
        var protoNodes = new LinkNode[sample.Length];
        for (int k = 0; k < sample.Length; k++)
        {
            int i = sample[k];
            nodes[k] = new EsaNode(set.Lb[i], set.Rb[i], set.Depth[i], 0, index.SuffixArray[set.Lb[i]], false);
            protoNodes[k] = set.Node(i, index.SuffixArray);
        }

        // agreement first (also warms both)
        int differences = 0;
        for (int k = 0; k < sample.Length; k++)
        {
            EsaNode actual = pack.SuffixLink(nodes[k]);
            LinkTarget expected = LinkKernels.Link(LinkKernel.PsvNsv, EndpointSource.Isa, ctx, protoNodes[k]);
            if (set.Depth[sample[k]] == 1) continue;
            if (actual.Lb != expected.Lb || actual.Rb != expected.Rb || actual.ParentDepth != expected.ParentDepth || actual.LastChild != expected.IsLast) differences++;
        }
        Emit.Row("nav_links", "differences", differences, "count", "production kernel against prototype (c2) with ISA, same nodes");
        if (differences != 0) Environment.ExitCode = 3;

        // paired design: the sample is split into two interleaved halves; in every round the production kernel gets one half and the prototype the other, the halves swap and the order alternates
        // between rounds, so that neither kernel profits from the cache lines that the other one has just touched
        var ratios = new List<double>();
        var allProduction = new List<double>();
        var allPrototype = new List<double>();
        for (int round = 0; round < 6; round++)
        {
            int parityProduction = round % 2, parityPrototype = 1 - parityProduction;
            double[] production, prototype;
            if (round % 4 < 2)
            {
                production = TimeLinks(sample.Length, parityProduction, k => pack.SuffixLink(nodes[k]).Lb, floorNs);
                prototype = TimeLinks(sample.Length, parityPrototype, k => sa[LinkKernels.Link(LinkKernel.PsvNsv, EndpointSource.Isa, ctx, protoNodes[k]).Lb], floorNs);
            }
            else
            {
                prototype = TimeLinks(sample.Length, parityPrototype, k => sa[LinkKernels.Link(LinkKernel.PsvNsv, EndpointSource.Isa, ctx, protoNodes[k]).Lb], floorNs);
                production = TimeLinks(sample.Length, parityProduction, k => pack.SuffixLink(nodes[k]).Lb, floorNs);
            }
            double mp = production[production.Length / 2], mq = prototype[prototype.Length / 2];
            ratios.Add(mp / mq);
            allProduction.AddRange(production);
            allPrototype.AddRange(prototype);
            Emit.Row("nav_links", $"production_median_ns_r{round}", mp, "ns", $"p99 {production[(int)(production.Length * 0.99)].ToString("F0", Inv)} ns; half {parityProduction}, {(round % 4 < 2 ? "first" : "second")} in the round");
            Emit.Row("nav_links", $"prototype_c2_isa_f32_pos_median_ns_r{round}", mq, "ns", $"p99 {prototype[(int)(prototype.Length * 0.99)].ToString("F0", Inv)} ns; half {parityPrototype}");
        }
        allProduction.Sort();
        allPrototype.Sort();
        ratios.Sort();
        Emit.Row("nav_links", "production_median_ns", allProduction[allProduction.Count / 2], "ns", $"median of all {allProduction.Count} timed links; p99 {allProduction[(int)(allProduction.Count * 0.99)].ToString("F0", Inv)} ns");
        Emit.Row("nav_links", "prototype_c2_isa_f32_pos_median_ns", allPrototype[allPrototype.Count / 2], "ns", $"median of all {allPrototype.Count} timed links; p99 {allPrototype[(int)(allPrototype.Count * 0.99)].ToString("F0", Inv)} ns");
        Emit.Row("nav_links", "ratio_production_over_prototype", (ratios[2] + ratios[3]) / 2, "ratio", $"median of the six per-round ratios [{ratios[0].ToString("F3", Inv)} .. {ratios[^1].ToString("F3", Inv)}]");
        Emit.Row("nav_links", "peak_rss_hwm", Process.StatusBytes("VmHWM"), "bytes", "VmHWM of the process");
        GC.KeepAlive(helper);
    }

    private static double[] TimeLinks(int count, int parity, Func<int, int> link, double floorNs)
    {
        var ns = new double[count / 2];
        long sink = 0;
        for (int j = 0; j < ns.Length; j++)
        {
            int k = (2 * j) + parity;
            long a = Stopwatch.GetTimestamp();
            sink += link(k);
            long b = Stopwatch.GetTimestamp();
            ns[j] = Math.Max(0, ((b - a) * 1e9 / Stopwatch.Frequency) - floorNs);
        }
        GC.KeepAlive(sink);
        Array.Sort(ns);
        return ns;
    }

    // ---------------------------------------------------------------- split: where the MEM time goes

    private static void Split(EsaIndex index, string query, int minLength, int calls)
    {
        // untimed reference: the plain navigator
        var seconds = new List<double>();
        IReadOnlyList<Triple>? plain = null;
        long allocated0 = GC.GetAllocatedBytesForCurrentThread();
        int gen0 = GC.CollectionCount(0);
        for (int call = 0; call < calls; call++)
        {
            long c0 = Stopwatch.GetTimestamp();
            plain = index.FindMaximalExactMatches(query, minLength);
            seconds.Add(Seconds(c0));
        }
        long allocatedPerCall = (GC.GetAllocatedBytesForCurrentThread() - allocated0) / calls;
        Report("split_plain", "mem", seconds, plain);
        double[] steady = seconds.Skip(calls - 5).Order().ToArray();
        double plainSeconds = steady[2];
        Emit.Row("split_plain", "allocated_bytes_per_call", allocatedPerCall, "bytes", $"{(GC.CollectionCount(0) - gen0)} gen-0 collections in {calls} calls; the result list and the lists of the generic algorithm");

        // counting pass (no clock): exact call counts
        var counts = new NavProfile(timed: false);
        var countingNav = new ProfilingNavigator(index.CreateNavigator(), counts);
        IReadOnlyList<Triple> viaCounting = SuffixTreeAlgorithms.FindMaximalExactMatches<EsaNode, ProfilingNavigator>(ref countingNav, query, minLength);
        bool same = viaCounting.SequenceEqual(plain!);
        Emit.Row("split", "counting_navigator_same_result", same ? 1 : 0, "bool", "the counting wrapper does not change the answer");
        if (!same) Environment.ExitCode = 3;

        // timed passes: heavy members are timed one by one (the timer floor is subtracted); the median of three
        double floorNs = Calibrate();
        Emit.Row("split", "timer_floor_ns", floorNs, "ns", "mean of an empty timed interval");
        NavProfile best = new(timed: true);
        double bestTotal = double.MaxValue;
        for (int round = 0; round < 3; round++)
        {
            var timed = new NavProfile(timed: true);
            var nav = new ProfilingNavigator(index.CreateNavigator(), timed);
            long c0 = Stopwatch.GetTimestamp();
            SuffixTreeAlgorithms.FindMaximalExactMatches<EsaNode, ProfilingNavigator>(ref nav, query, minLength);
            double total = Seconds(c0);
            if (total < bestTotal)
            {
                bestTotal = total;
                best = timed;
            }
        }
        double heavy = 0;
        string[] names = NavProfile.Names;
        for (int m = 0; m < names.Length; m++)
        {
            long calls0 = counts.Calls[m];
            Emit.Row("split_calls", names[m], calls0, "count", names[m] is "TryGetChildAfterLink" or "TryGetChildOther" ? "TryGetChild: after an IsRoot call (a rescan, or the one extension lookup that follows it) / otherwise" : "");
            if (!NavProfile.Timed[m]) continue;
            double netSeconds = Math.Max(0, (best.Ticks[m] / (double)Stopwatch.Frequency) - (best.Calls[m] * floorNs * 1e-9));
            heavy += netSeconds;
            Emit.Row("split_time", names[m], netSeconds, "s", $"{(netSeconds / plainSeconds * 100).ToString("F1", Inv)} % of the plain steady call ({plainSeconds.ToString("R", Inv)} s); {(netSeconds * 1e9 / Math.Max(1, best.Calls[m])).ToString("F0", Inv)} ns per call");
        }
        Emit.Row("split_time", "timed_members_total", heavy, "s", $"{(heavy / plainSeconds * 100).ToString("F1", Inv)} % of the plain steady call");
        Emit.Row("split_time", "remainder", plainSeconds - heavy, "s", $"{((plainSeconds - heavy) / plainSeconds * 100).ToString("F1", Inv)} %: GetEdgeSymbol / LengthOf / IsNull / IsRoot ({counts.Calls[NavProfile.GetEdgeSymbol] + counts.Calls[NavProfile.LengthOf] + counts.Calls[NavProfile.IsNull] + counts.Calls[NavProfile.IsRoot]} calls), the generic algorithm itself and its allocations");
        Emit.Row("split_time", "timed_run_total", bestTotal, "s", "the best of three timed runs (inflated by the timers); the shares of THIS total are the better estimate of the inline numbers");

        // replay: the recorded arguments of each member are played back on their own, one after the other with a data dependency (like the real chain), one timer around the whole loop
        var rec = new NavProfile(timed: false) { Record = true };
        var recNav = new ProfilingNavigator(index.CreateNavigator(), rec);
        SuffixTreeAlgorithms.FindMaximalExactMatches<EsaNode, ProfilingNavigator>(ref recNav, query, minLength);
        EsaNavigator nav2 = index.CreateNavigator();
        var linkNodes = rec.LinkNodes.ToArray();
        var childCalls = rec.ChildCalls.ToArray();
        var edgeCalls = rec.EdgeCalls.ToArray();
        var lengthNodes = rec.LengthNodes.ToArray();
        var collectNodes = rec.CollectNodes.ToArray();
        var childrenNodes = rec.ChildrenNodes.ToArray();
        double linkNs = ReplayNs(linkNodes.Length, (k, d) => nav2.GetSuffixLink(linkNodes[k + d]).Lb);
        double childNs = ReplayNs(childCalls.Length, (k, d) => { nav2.TryGetChild(childCalls[k + d].Node, childCalls[k + d].Key, out EsaNode c); return c.Lb + 1; });
        double edgeNs = ReplayNs(edgeCalls.Length, (k, d) => nav2.GetEdgeSymbol(edgeCalls[k + d].Node, edgeCalls[k + d].Offset) + 1);
        double lengthNs = ReplayNs(lengthNodes.Length, (k, d) => nav2.LengthOf(lengthNodes[k + d]));
        var scratch = new List<int>();
        double collectNs = ReplayNs(collectNodes.Length, (k, d) => { scratch.Clear(); nav2.CollectLeaves(collectNodes[k + d], 0, scratch); return scratch.Count; });
        var childList = new List<EsaNode>();
        double childrenNs = ReplayNs(childrenNodes.Length, (k, d) => { nav2.GetChildren(childrenNodes[k + d], childList); return childList.Count; });
        double emptyNs = ReplayNs(1_000_000, (k, d) => d);
        Emit.Row("replay", "loop_overhead_ns", emptyNs, "ns", "an empty replay loop with the dependency chain; subtract from the figures below when the member is cheap");
        double childrenSeconds = childrenNodes.Length * childrenNs * 1e-9, linkSeconds = linkNodes.Length * linkNs * 1e-9, childSeconds = childCalls.Length * childNs * 1e-9, edgeSeconds = edgeCalls.Length * edgeNs * 1e-9, lengthSeconds = lengthNodes.Length * lengthNs * 1e-9, collectSeconds = collectNodes.Length * collectNs * 1e-9;
        double lightSeconds = (counts.Calls[NavProfile.IsNull] + counts.Calls[NavProfile.IsRoot]) * 1e-9 * 1.0;
        double sum = linkSeconds + childSeconds + childrenSeconds + edgeSeconds + lengthSeconds + collectSeconds + lightSeconds;
        string Share(double x) => $"{(x / plainSeconds * 100).ToString("F1", Inv)} % of the plain steady call";
        Emit.Row("replay", "GetSuffixLink_ns", linkNs, "ns", $"{linkNodes.Length} calls -> {linkSeconds.ToString("R", Inv)} s, {Share(linkSeconds)}");
        Emit.Row("replay", "TryGetChild_ns", childNs, "ns", $"{childCalls.Length} calls -> {childSeconds.ToString("R", Inv)} s, {Share(childSeconds)}");
        Emit.Row("replay", "GetChildren_ns", childrenNs, "ns", $"{childrenNodes.Length} calls -> {childrenSeconds.ToString("R", Inv)} s, {Share(childrenSeconds)}");
        Emit.Row("replay", "GetEdgeSymbol_ns", edgeNs, "ns", $"{edgeCalls.Length} calls -> {edgeSeconds.ToString("R", Inv)} s, {Share(edgeSeconds)}");
        Emit.Row("replay", "LengthOf_ns", lengthNs, "ns", $"{lengthNodes.Length} calls -> {lengthSeconds.ToString("R", Inv)} s, {Share(lengthSeconds)}");
        Emit.Row("replay", "CollectLeaves_ns", collectNs, "ns", $"{collectNodes.Length} calls -> {collectSeconds.ToString("R", Inv)} s, {Share(collectSeconds)}");
        Emit.Row("replay", "IsNull_IsRoot_estimate", lightSeconds, "s", "1 ns per call (inlined compares)");
        Emit.Row("replay", "members_total", sum, "s", Share(sum));
        Emit.Row("replay", "remainder", plainSeconds - sum, "s", $"{Share(plainSeconds - sum)}: the generic algorithm itself (loop control, list handling, left-maximality tests, the result list), allocations, and what the replay cannot see (cache warmth between members)");
    }

    /// <summary>Plays <paramref name="step"/> for k = 0..count-1 three times and returns the best ns per call; the delegate returns a non-negative value that is added to the next index, a data dependency the compiler cannot remove.</summary>
    private static double ReplayNs(int count, Func<int, int, int> step)
    {
        if (count == 0) return 0;
        double best = double.MaxValue;
        for (int rep = 0; rep < 3; rep++)
        {
            long t0 = Stopwatch.GetTimestamp();
            int dependency = 0;
            for (int k = 0; k < count; k++) dependency = step(k, dependency & s_zero) & s_zero;
            double ns = (Stopwatch.GetTimestamp() - t0) * 1e9 / Stopwatch.Frequency / count;
            if (dependency != 0) Console.Error.WriteLine("unreachable");
            best = Math.Min(best, ns);
        }
        return best;
    }


    private static double Calibrate()
    {
        long ticks = 0;
        for (int i = 0; i < 2_000_000; i++)
        {
            long a = Stopwatch.GetTimestamp();
            long b = Stopwatch.GetTimestamp();
            ticks += b - a;
        }
        return ticks * 1e9 / Stopwatch.Frequency / 2_000_000;
    }
}

/// <summary>Call counters (and, for the heavy members, tick sums) of one <see cref="ProfilingNavigator"/> run.</summary>
internal sealed class NavProfile(bool timed)
{
    internal const int GetEdgeSymbol = 0, LengthOf = 1, IsNull = 2, IsRoot = 3, GetSuffixLink = 4, TryGetChildAfterLink = 5, TryGetChildOther = 6, GetChildren = 7, CollectLeaves = 8, FindAnyLeafPosition = 9, ChildrenReturned = 10, LeavesCollected = 11;

    internal static readonly string[] Names =
    [
        "GetEdgeSymbol", "LengthOf", "IsNull", "IsRoot", "GetSuffixLink", "TryGetChildAfterLink", "TryGetChildOther", "GetChildren", "CollectLeaves", "FindAnyLeafPosition", "children_returned", "leaves_collected",
    ];

    internal static readonly bool[] Timed = [false, false, false, false, true, true, true, true, true, true, false, false];

    internal long[] Calls { get; } = new long[Names.Length];

    internal long[] Ticks { get; } = new long[Names.Length];

    internal bool TimeIt { get; } = timed;

    /// <summary>Set by <c>IsRoot</c> (called only at the start of the link-and-rescan step of the generic algorithms) and cleared by the next member that the extension step uses.</summary>
    internal bool AfterLink { get; set; }

    /// <summary>When set, the arguments of the heavy calls are recorded for the replay.</summary>
    internal bool Record { get; set; }

    internal List<EsaNode> LinkNodes { get; } = [];

    internal List<(EsaNode Node, int Key)> ChildCalls { get; } = [];

    internal List<(EsaNode Node, int Offset)> EdgeCalls { get; } = [];

    internal List<EsaNode> LengthNodes { get; } = [];

    internal List<EsaNode> CollectNodes { get; } = [];

    internal List<EsaNode> ChildrenNodes { get; } = [];
}

/// <summary>A navigator that forwards to <see cref="EsaNavigator"/> and counts (and optionally times) the calls of the generic algorithm.</summary>
internal readonly struct ProfilingNavigator(EsaNavigator inner, NavProfile profile) : ISuffixTreeNavigator<EsaNode>
{
    public ITextSource Text => inner.Text;

    public EsaNode Root => inner.Root;

    public EsaNode NullNode => inner.NullNode;

    public bool IsNull(EsaNode node)
    {
        profile.Calls[NavProfile.IsNull]++;
        return inner.IsNull(node);
    }

    public bool IsRoot(EsaNode node)
    {
        profile.Calls[NavProfile.IsRoot]++;
        profile.AfterLink = true;
        return inner.IsRoot(node);
    }

    public int GetEdgeSymbol(EsaNode node, int offset)
    {
        profile.Calls[NavProfile.GetEdgeSymbol]++;
        if (profile.Record) profile.EdgeCalls.Add((node, offset));
        profile.AfterLink = false;
        return inner.GetEdgeSymbol(node, offset);
    }

    public int LengthOf(EsaNode node)
    {
        profile.Calls[NavProfile.LengthOf]++;
        if (profile.Record) profile.LengthNodes.Add(node);
        return inner.LengthOf(node);
    }

    public EsaNode GetSuffixLink(EsaNode node)
    {
        profile.Calls[NavProfile.GetSuffixLink]++;
        if (profile.Record) profile.LinkNodes.Add(node);
        if (!profile.TimeIt) return inner.GetSuffixLink(node);
        long t = Stopwatch.GetTimestamp();
        EsaNode link = inner.GetSuffixLink(node);
        profile.Ticks[NavProfile.GetSuffixLink] += Stopwatch.GetTimestamp() - t;
        return link;
    }

    public bool TryGetChild(EsaNode node, int key, out EsaNode child)
    {
        int slot = profile.AfterLink ? NavProfile.TryGetChildAfterLink : NavProfile.TryGetChildOther;
        profile.Calls[slot]++;
        if (profile.Record) profile.ChildCalls.Add((node, key));
        if (!profile.TimeIt) return inner.TryGetChild(node, key, out child);
        long t = Stopwatch.GetTimestamp();
        bool found = inner.TryGetChild(node, key, out child);
        profile.Ticks[slot] += Stopwatch.GetTimestamp() - t;
        return found;
    }

    public void CollectLeaves(EsaNode node, int depth, List<int> results)
    {
        profile.Calls[NavProfile.CollectLeaves]++;
        profile.AfterLink = false;
        if (profile.Record) profile.CollectNodes.Add(node);
        int before = results.Count;
        if (!profile.TimeIt)
        {
            inner.CollectLeaves(node, depth, results);
        }
        else
        {
            long t = Stopwatch.GetTimestamp();
            inner.CollectLeaves(node, depth, results);
            profile.Ticks[NavProfile.CollectLeaves] += Stopwatch.GetTimestamp() - t;
        }
        profile.Calls[NavProfile.LeavesCollected] += results.Count - before;
    }

    public int FindAnyLeafPosition(EsaNode node, int depthFromRoot)
    {
        profile.Calls[NavProfile.FindAnyLeafPosition]++;
        profile.AfterLink = false;
        return inner.FindAnyLeafPosition(node, depthFromRoot);
    }

    public void GetChildren(EsaNode node, List<EsaNode> children)
    {
        profile.Calls[NavProfile.GetChildren]++;
        profile.AfterLink = false;
        if (profile.Record) profile.ChildrenNodes.Add(node);
        if (!profile.TimeIt)
        {
            inner.GetChildren(node, children);
        }
        else
        {
            long t = Stopwatch.GetTimestamp();
            inner.GetChildren(node, children);
            profile.Ticks[NavProfile.GetChildren] += Stopwatch.GetTimestamp() - t;
        }
        profile.Calls[NavProfile.ChildrenReturned] += children.Count;
    }
}
