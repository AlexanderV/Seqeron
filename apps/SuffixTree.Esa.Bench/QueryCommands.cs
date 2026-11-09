using System.Diagnostics;
using System.Globalization;
using SuffixTree.Esa;

namespace SuffixTree.Esa.Bench;

/// <summary>
/// Work package 1b.1: the child-table layout comparison of prototype P3 (size, scan time, write pattern). One process measures one text; the child-table scan is timed
/// on its own on the same SA and LCP.
/// </summary>
internal static class QueryCommands
{
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;

    /// <summary>The texts that are generated in the process (no file): the adversarial rows of P3. <c>name:n</c>.</summary>
    internal static string Synthetic(string spec)
    {
        string[] f = spec.Split(':');
        int n = f.Length > 1 ? int.Parse(f[1], Inv) : 1_000_000;
        switch (f[0])
        {
            case "a^n":
                return new string('a', n);
            case "abab":
                return string.Create(n, 0, static (span, _) => { for (int i = 0; i < span.Length; i++) span[i] = (i & 1) == 0 ? 'a' : 'b'; });
            case "fibonacci":
                {
                    string prev = "a", cur = "ab";
                    while (cur.Length < n) (prev, cur) = (cur, cur + prev);
                    return cur[..n];
                }
            case "thue-morse":
                return string.Create(n, 0, static (span, _) => { for (int i = 0; i < span.Length; i++) span[i] = (char)('a' + (int.PopCount(i) & 1)); });
            case "n-runs":
                {
                    // a genome-like text with long N runs: 4 runs of n/8 'N' between random ACGT blocks
                    var rng = new Random(12345);
                    var chars = new char[n];
                    for (int i = 0; i < n; i++) chars[i] = "ACGT"[rng.Next(4)];
                    int run = n / 8;
                    for (int r = 0; r < 4; r++) Array.Fill(chars, 'N', ((2 * r) + 1) * run - (run / 2), run);
                    return new string(chars);
                }
            default:
                throw new ArgumentException($"unknown synthetic text '{f[0]}' (a^n, abab, fibonacci, thue-morse, n-runs)");
        }
    }

    internal static void Query(Options o)
    {
        string text = o.Synthetic is not null ? Synthetic(o.Synthetic) : Pipeline.LoadText(o.Input);
        int n = text.Length;
        Emit.Row("query", "input_chars", n, "chars");

        long t0 = Stopwatch.GetTimestamp();
        PreparedText prepared = TextPreparation.Prepare(text);
        int[] sa = SaisKernel.Build(SaisKernelKind.Pointer, prepared);
        CompactLcp lcp = LcpKernel.Build(LcpKernelKind.Direct, prepared, sa);
        double saLcpSeconds = (Stopwatch.GetTimestamp() - t0) / (double)Stopwatch.Frequency;
        Emit.Row("query", "build_sa_lcp_s", saLcpSeconds, "s", "prepare + SA-IS pointer kernel + direct LCP kernel (cold JIT)");

        foreach (ChildLayoutKind layout in o.Layouts)
        {
            string op = $"layout_{layout}";
            long c0 = Stopwatch.GetTimestamp();
            EsaIndex index = EsaIndex.FromArrays(text, prepared, sa, lcp, layout);
            double scan = (Stopwatch.GetTimestamp() - c0) / (double)Stopwatch.Frequency;
            Emit.Row(op, "child_scan_s", scan, "s", "ChildTable.Build: one scan over the LCP stream (cold JIT of the scan on the first layout)");
            Emit.Row(op, "child_bytes", index.Child.PayloadBytes, "bytes");
            Emit.Row(op, "child_bytes_per_char", index.Child.PayloadBytes / (double)n, "bytes/char");
            Emit.Row(op, "internal_nodes", index.Statistics.InternalCount, "count");
            Emit.Row(op, "stack_peak", index.Statistics.StackPeak, "count", "peak height of the scan stack (16 bytes per entry in the in-memory scan)");

            if (o.Trace)
            {
                foreach (int cachePages in new[] { 16, 64, 1024 })
                {
                    ChildWriteTrace t = ChildWriteTrace.Run(layout, n, lcp, 1024, cachePages, trackSlots: false);
                    string note = $"LRU of {cachePages} pages of 4,096 bytes";
                    if (cachePages == 16)
                    {
                        Emit.Row(op, "trace_writes", t.Writes, "count");
                        Emit.Row(op, "trace_at_cursor", t.AtCursor, "count", "writes at rank cursor-1 or cursor");
                        Emit.Row(op, "trace_behind", t.Behind, "count", "writes behind the cursor (left boundary of a closing interval, previous l-index of an open one)");
                        Emit.Row(op, "trace_max_back_slots", t.MaxBackDistance, "slots", "largest cursor - target of a behind write");
                        Emit.Row(op, "trace_backward_runs", t.BackwardRuns, "count", "maximal runs of consecutive behind writes with descending targets");
                    }
                    Emit.Row(op, $"trace_page_misses_{cachePages}", t.PageMisses, "count", note);
                }
            }
        }
    }
}
