using System.Globalization;

namespace SuffixTree.PerfHarness;

/// <summary>How a metric is judged against its baseline.</summary>
internal enum Direction
{
    /// <summary>Smaller is better (time, memory, bytes per character).</summary>
    LowerIsBetter,

    /// <summary>Larger is better (throughput).</summary>
    HigherIsBetter,

    /// <summary>Must be exactly equal (oracle: counts and checksums of the query results).</summary>
    Exact,

    /// <summary>Informational ceiling: the value must not exceed the baseline times <see cref="Expectation.Tolerance"/>.</summary>
    CeilingFactor,
}

/// <summary>One expected value of the 10 Mbp row.</summary>
/// <param name="Target">Measured subject: <c>inmem</c> or <c>persistent</c>.</param>
/// <param name="Op">Operation name of the line protocol (<c>build</c>, <c>count_present20</c>, ...).</param>
/// <param name="Metric">Metric name of the line protocol.</param>
/// <param name="Value">Baseline value from <c>analysis/bench/bench.tsv</c> (the numbers of BENCHMARK.md).</param>
/// <param name="Direction">How the measurement is compared.</param>
/// <param name="Unit">Unit shown in the report.</param>
/// <param name="Tolerance">Relative band (0.15 = +-15 %), or the allowed factor for <see cref="Direction.CeilingFactor"/>.</param>
internal sealed record Expectation(string Target, string Op, string Metric, double Value, Direction Direction, string Unit, double Tolerance)
{
    public string Key => $"{Target}/{Op}/{Metric}";
}

/// <summary>
/// Baselines of the <c>dna_10m</c> row of BENCHMARK.md section 1 (median of 3 runs, 4 vCPU Intel Xeon under KVM, .NET 10.0.403,
/// library commit 242318429a5e), copied from <c>docs/refactoring/suffix-tree-upgrade/analysis/bench/bench.tsv</c>.
/// </summary>
internal static class Baseline
{
    public const double DefaultTolerance = 0.15;

    /// <summary>Returns the expectations for <paramref name="input"/>, or an empty list when the input has no recorded baseline.</summary>
    public static IReadOnlyList<Expectation> For(string input, double tolerance)
        => input == "dna_10m" ? Dna10M(tolerance) : Array.Empty<Expectation>();

    private static List<Expectation> Dna10M(double tol)
    {
        const string im = "inmem", ps = "persistent";
        var list = new List<Expectation>
        {
            // ── in-memory tree (BENCHMARK.md section 1, column "Seqeron in-memory") ──
            new(im, "build", "time", 14.640838, Direction.LowerIsBetter, "s", tol),
            new(im, "build", "bytes_per_char", 168.790601, Direction.LowerIsBetter, "B/char", tol),
            new(im, "build", "peak_rss_hwm_after_build", 1760518144, Direction.LowerIsBetter, "bytes", tol),
            new(im, "count_present20", "throughput", 589810.434926, Direction.HigherIsBetter, "queries/s", tol),
            new(im, "count_present20", "latency_p50", 1969, Direction.LowerIsBetter, "ns", tol),
            new(im, "count_present20_4t", "throughput", 2561619.8, Direction.HigherIsBetter, "queries/s", tol),
            new(im, "locate8", "throughput", 49653.310654, Direction.HigherIsBetter, "patterns/s", tol),

            // ── persistent tree (column "Seqeron persistent (mmap)") ──
            new(ps, "build", "time", 14.905865, Direction.LowerIsBetter, "s", tol),
            new(ps, "build", "bytes_per_char", 56.181149, Direction.LowerIsBetter, "B/char", tol),
            new(ps, "build_process", "peak_rss_hwm", 1178689536, Direction.LowerIsBetter, "bytes", tol),
            new(ps, "count_present20", "throughput", 388370.705275, Direction.HigherIsBetter, "queries/s", tol),
            new(ps, "count_present20", "latency_p50", 2580, Direction.LowerIsBetter, "ns", tol),
            new(ps, "count_present20_4t", "throughput", 1716118.0, Direction.HigherIsBetter, "queries/s", tol),
            new(ps, "locate8", "throughput", 39942.861736, Direction.HigherIsBetter, "patterns/s", tol),
            // Load is a few milliseconds: relative jitter is large, so only a generous ceiling (x5) is enforced.
            new(ps, "warm_load", "time", 0.00393, Direction.CeilingFactor, "s", 5.0),
        };

        // Oracles: results of the query sets are identical for every correct tree (all tools in the benchmark agree).
        foreach (string t in new[] { im, ps })
        {
            list.Add(new(t, "build", "node_count", 16229747, Direction.Exact, "count", 0));
            list.Add(new(t, "build", "leaf_count", 10000000, Direction.Exact, "count", 0));
            list.Add(new(t, "count_present20", "checksum", 100000, Direction.Exact, "count", 0));
            list.Add(new(t, "count_absent20", "checksum", 1, Direction.Exact, "count", 0));
            list.Add(new(t, "locate8", "total_hits", 152362, Direction.Exact, "count", 0));
            list.Add(new(t, "locate8", "checksum", 761250397446, Direction.Exact, "count", 0));
            list.Add(new(t, "lrs", "length", 25, Direction.Exact, "chars", 0));
            list.Add(new(t, "mem_q100k_min20", "matches", 817, Direction.Exact, "count", 0));
            list.Add(new(t, "mem_q100k_min20", "checksum", 77915707478, Direction.Exact, "count", 0));
            list.Add(new(t, "mum_q100k_min20", "matches", 816, Direction.Exact, "count", 0));
            list.Add(new(t, "anchors_q100k_min20", "matches", 816, Direction.Exact, "count", 0));
            list.Add(new(t, "count_present20_4t", "checksum_ok", 1, Direction.Exact, "bool", 0));
        }
        return list;
    }

    /// <summary>Verdict of one comparison.</summary>
    public enum Verdict
    {
        Ok,
        Regression,
        Improved,
        Mismatch,
        Missing,
    }

    /// <summary>Compares a measured value with an expectation.</summary>
    public static Verdict Judge(Expectation e, double? measured)
    {
        if (measured is not { } m)
            return Verdict.Missing;

        switch (e.Direction)
        {
            case Direction.Exact:
                return m == e.Value ? Verdict.Ok : Verdict.Mismatch;
            case Direction.CeilingFactor:
                return m <= e.Value * e.Tolerance ? Verdict.Ok : Verdict.Regression;
            case Direction.LowerIsBetter:
                if (m > e.Value * (1 + e.Tolerance)) return Verdict.Regression;
                return m < e.Value * (1 - e.Tolerance) ? Verdict.Improved : Verdict.Ok;
            default: // HigherIsBetter
                if (m < e.Value * (1 - e.Tolerance)) return Verdict.Regression;
                return m > e.Value * (1 + e.Tolerance) ? Verdict.Improved : Verdict.Ok;
        }
    }

    public static string Format(double v)
    {
        double a = Math.Abs(v);
        if (a >= 1000)
            return v.ToString("N0", CultureInfo.InvariantCulture);
        if (a >= 10)
            return v.ToString("N1", CultureInfo.InvariantCulture);
        return v.ToString("0.####", CultureInfo.InvariantCulture);
    }
}
