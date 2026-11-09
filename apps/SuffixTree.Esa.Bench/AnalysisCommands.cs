using System.Diagnostics;
using System.Globalization;
using SuffixTree.Esa;

namespace SuffixTree.Esa.Bench;

/// <summary>
/// Work package 1b.3: time of the native analysis operations of the in-memory ESA (longest repeated substring and its ties, distinct-substring counts, suffix
/// enumeration) and of the operations that run on <c>Traverse</c> through the generic algorithms (maximal repeated pairs, k-common substrings), on one text per process.
/// With <c>--check true</c> the native results are also compared with the <c>Traverse</c>-default path of the same index (an equality failure sets exit code 3).
/// </summary>
internal static class AnalysisCommands
{
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;

    /// <summary>The ESA behind <see cref="ISuffixTree"/> with the default members left to the interface (so they run on <c>Traverse</c>); the members that need other machinery throw.</summary>
    private sealed class TraverseView(EsaIndex index) : ISuffixTree
    {
        public ITextSource Text { get; } = new StringTextSource(index.Text);

        public int NodeCount => checked((int)index.Statistics.NodeCount);

        public int LeafCount => index.Statistics.LeafCount;

        public int MaxDepth => index.Statistics.MaxDepth;

        public bool IsEmpty => index.Length == 0;

        public void Traverse(ISuffixTreeVisitor visitor) => index.Traverse(visitor);

        public string PrintTree() => throw Unsupported();

        public bool Contains(string value) => index.Contains(value);

        public bool Contains(ReadOnlySpan<char> value) => index.Contains(value);

        public IReadOnlyList<int> FindAllOccurrences(string pattern) => index.FindAllOccurrences(pattern);

        public IReadOnlyList<int> FindAllOccurrences(ReadOnlySpan<char> pattern) => index.FindAllOccurrences(pattern);

        public int CountOccurrences(string pattern) => index.CountOccurrences(pattern);

        public int CountOccurrences(ReadOnlySpan<char> pattern) => index.CountOccurrences(pattern);

        public string LongestRepeatedSubstring() => index.LongestRepeatedSubstring();

        public IReadOnlyList<string> GetAllSuffixes() => index.GetAllSuffixes();

        public IEnumerable<string> EnumerateSuffixes() => index.EnumerateSuffixes();

        public string LongestCommonSubstring(string other) => throw Unsupported();

        public string LongestCommonSubstring(ReadOnlySpan<char> other) => throw Unsupported();

        public (string Substring, int PositionInText, int PositionInOther) LongestCommonSubstringInfo(string other) => throw Unsupported();

        public (string Substring, IReadOnlyList<int> PositionsInText, IReadOnlyList<int> PositionsInOther) FindAllLongestCommonSubstrings(string other) => throw Unsupported();

        public IReadOnlyList<(int PositionInText, int PositionInQuery, int Length)> FindExactMatchAnchors(string query, int minLength) => throw Unsupported();

        private static NotSupportedException Unsupported() => new("the benchmark view serves Traverse-based operations only");
    }

    internal static void Run(Options o)
    {
        string text = o.Synthetic is not null ? QueryCommands.Synthetic(o.Synthetic) : Pipeline.LoadText(o.Input);
        int n = text.Length;
        Emit.Row("analysis", "input_chars", n, "chars");

        long t0 = Stopwatch.GetTimestamp();
        EsaIndex index = EsaIndex.Build(text);
        Emit.Row("analysis", "build_total_s", Seconds(t0), "s", "EsaIndex.Build: prepare + SA-IS + LCP + child table + statistics (cold JIT)");
        EsaStatistics st = index.Statistics;
        Emit.Row("analysis", "lrs_depth", st.LrsDepth, "chars", "header statistic");

        bool ok = true;

        // 1. longest repeated substring (O(L) copy of the header's deepest node)
        string? lrs = null;
        Time("lrs", "LongestRepeatedSubstring: Substring of the header's deepest interval", () => lrs = index.LongestRepeatedSubstring());
        Emit.Row("lrs", "length", lrs!.Length, "chars");
        Time("lrs_memory", "LongestRepeatedSubstringMemory: a slice of the text, no copy", () => _ = index.LongestRepeatedSubstringMemory().Length);

        // 2. all ties: one pass over the LCP array
        IReadOnlyList<(string Substring, IReadOnlyList<int> Positions)>? all = null;
        Time("find_all_lrs", "FindAllLongestRepeatedSubstrings: one pass over the LCP array + suffix-array slices", () => all = index.FindAllLongestRepeatedSubstrings());
        Emit.Row("find_all_lrs", "ties", all!.Count, "count");
        Emit.Row("find_all_lrs", "positions", all.Sum(e => (long)e.Positions.Count), "count", "sum over the ties (disjoint intervals: at most n)");
        Emit.Row("find_all_lrs", "output_chars", all.Sum(e => (long)e.Substring.Length), "chars", "sum of the substring lengths (the output of the operation)");

        // 3. distinct substrings
        long total = 0;
        Time("distinct_total", "CountDistinctSubstrings: n(n+1)/2 - sum LCP", () => total = index.CountDistinctSubstrings());
        Emit.Row("distinct_total", "value", total, "count");
        var byLength = new Dictionary<int, long[]>();
        foreach (int m in new[] { 100, 10_000, n })
        {
            int mm = Math.Min(m, n);
            Time($"distinct_by_length_{(m == n ? "n" : m.ToString(Inv))}", $"CountDistinctSubstringsByLength({m}): LCP histogram + suffix sums, result {8L * (mm + 1)} bytes", () => byLength[m] = index.CountDistinctSubstringsByLength(m));
        }

        // 4. suffix enumeration: the SA walk is O(1) per suffix, the cost is the copy of the string
        const int Taken = 20;
        long a0 = GC.GetAllocatedBytesForCurrentThread();
        long s0 = Stopwatch.GetTimestamp();
        long copied = 0;
        foreach (string suffix in index.EnumerateSuffixes().Take(Taken)) copied += suffix.Length;
        double enumerate = Seconds(s0);
        Emit.Row("enumerate_suffixes", "first_20_s", enumerate, "s", $"EnumerateSuffixes().Take({Taken}): lazy, one string per step, cold JIT");
        Emit.Row("enumerate_suffixes", "first_20_chars", copied, "chars", "characters materialised (shortest suffixes first in rank order, so the sizes are the ranks' lengths)");
        Emit.Row("enumerate_suffixes", "first_20_allocated_bytes", GC.GetAllocatedBytesForCurrentThread() - a0, "bytes");
        Emit.Row("enumerate_suffixes", "ns_per_char", enumerate * 1e9 / Math.Max(1, copied), "ns", "includes the walk and the allocation");

        // 5. GetAllSuffixes materialises Theta(n^2) characters by definition: a prefix of the text
        foreach (int prefix in new[] { 2_000, 10_000 })
        {
            if (prefix > n) continue;
            EsaIndex small = EsaIndex.Build(text[..prefix]);
            long a1 = GC.GetAllocatedBytesForCurrentThread();
            long g0 = Stopwatch.GetTimestamp();
            IReadOnlyList<string> suffixes = small.GetAllSuffixes();
            double elapsed = Seconds(g0);
            Emit.Row($"get_all_suffixes_{prefix}", "time_s", elapsed, "s", $"GetAllSuffixes of the first {prefix} characters ({(long)prefix * (prefix + 1) / 2} characters of output)");
            Emit.Row($"get_all_suffixes_{prefix}", "allocated_bytes", GC.GetAllocatedBytesForCurrentThread() - a1, "bytes");
            GC.KeepAlive(suffixes);
        }

        // 6. the operations that run on Traverse (generic algorithms): maximal repeated pairs and k-common substrings
        if (o.Generic)
        {
            ISuffixTree view = new TraverseView(index);
            IReadOnlyList<(int FirstPosition, int SecondPosition, int Length)>? pairs = null;
            long p0 = Stopwatch.GetTimestamp();
            pairs = view.FindMaximalRepeatedPairs(o.MinLength);
            Emit.Row("maximal_repeated_pairs", "time_s", Seconds(p0), "s", $"FindMaximalRepeatedPairs({o.MinLength}) through the interface default on Traverse (Traverse itself includes the 4 B/char leftmost-occurrence pre-pass), cold");
            Emit.Row("maximal_repeated_pairs", "pairs", pairs.Count, "count");
            Emit.Row("maximal_repeated_pairs", "peak_rss_after", Process.StatusBytes("VmHWM"), "bytes", "VmHWM so far");

            int half = n / 2;
            string[] halves = [text[..half], text[half..]];
            long k0 = Stopwatch.GetTimestamp();
            IReadOnlyList<string> common = SuffixTreeAlgorithms.FindLongestCommonSubstrings(halves, 2, static s => new TraverseView(EsaIndex.Build(s)));
            Emit.Row("k_common", "time_s", Seconds(k0), "s", "FindLongestCommonSubstrings of the two halves of the text, support 2: generalized text + EsaIndex.Build + Traverse + colour counting, cold");
            Emit.Row("k_common", "results", common.Count, "count");
            Emit.Row("k_common", "length", common.Count == 0 ? 0 : common[0].Length, "chars");
        }

        // 7. the Traverse-default path as the oracle of the native results
        if (o.Check)
        {
            ISuffixTree view = new TraverseView(index);
            long c0 = Stopwatch.GetTimestamp();
            long defTotal = SuffixTreeAlgorithms.CountDistinctSubstrings(view);
            Emit.Row("check", "default_total_s", Seconds(c0), "s", "CountDistinctSubstrings through Traverse (the interface default), for comparison with the native time");
            ok &= Compare("distinct_total", defTotal == total);
            foreach (var (m, native) in byLength)
            {
                long c1 = Stopwatch.GetTimestamp();
                long[] def = SuffixTreeAlgorithms.CountDistinctSubstringsByLength(view, m);
                Emit.Row("check", $"default_by_length_{(m == n ? "n" : m.ToString(Inv))}_s", Seconds(c1), "s", "CountDistinctSubstringsByLength through Traverse");
                ok &= Compare($"distinct_by_length_{m}", def.AsSpan().SequenceEqual(native));
            }
            long c2 = Stopwatch.GetTimestamp();
            var defAll = SuffixTreeAlgorithms.FindAllLongestRepeatedSubstrings(view);
            Emit.Row("check", "default_find_all_lrs_s", Seconds(c2), "s", "FindAllLongestRepeatedSubstrings through Traverse + a search per tie");
            bool same = defAll.Count == all.Count;
            for (int i = 0; same && i < all.Count; i++)
                same = defAll[i].Substring == all[i].Substring && defAll[i].Positions.SequenceEqual(all[i].Positions);
            ok &= Compare("find_all_lrs", same);
            ok &= Compare("lrs_is_a_tie", all.Count == 0 ? lrs.Length == 0 : all.Any(e => e.Substring == lrs));
            ok &= Compare("lrs_length", lrs.Length == st.LrsDepth);
            // the first suffixes against the definition (comparison of the suffix array neighbours)
            bool ordered = true;
            string? previous = null;
            foreach (string suffix in index.EnumerateSuffixes().Take(Taken))
            {
                ordered &= previous is null || string.CompareOrdinal(previous, suffix) < 0;
                ordered &= text.EndsWith(suffix, StringComparison.Ordinal);
                previous = suffix;
            }
            ok &= Compare("suffix_order", ordered);
        }

        Emit.Row("analysis", "peak_rss_hwm", Process.StatusBytes("VmHWM"), "bytes", "VmHWM of the process after everything above");
        if (!ok) Environment.ExitCode = 3;
    }

    private static bool Compare(string what, bool equal)
    {
        Emit.Row("check", what, equal ? 1 : 0, "bool", "native result equals the Traverse-default / definition-level result");
        return equal;
    }

    private static double Seconds(long start) => (Stopwatch.GetTimestamp() - start) / (double)Stopwatch.Frequency;

    /// <summary>The first call (cold JIT included) and the median of three further calls; the bytes allocated by the first call.</summary>
    private static void Time(string op, string note, Action run)
    {
        long a0 = GC.GetAllocatedBytesForCurrentThread();
        long t0 = Stopwatch.GetTimestamp();
        run();
        double cold = Seconds(t0);
        long allocated = GC.GetAllocatedBytesForCurrentThread() - a0;
        var warm = new double[3];
        for (int i = 0; i < warm.Length; i++)
        {
            long t1 = Stopwatch.GetTimestamp();
            run();
            warm[i] = Seconds(t1);
        }
        Array.Sort(warm);
        Emit.Row(op, "cold_s", cold, "s", note);
        Emit.Row(op, "warm_median_s", warm[1], "s", $"median of 3 further calls in the same process [{warm[0].ToString("R", Inv)} .. {warm[2].ToString("R", Inv)}]");
        Emit.Row(op, "allocated_bytes", allocated, "bytes", "first call");
    }
}
