using System.Diagnostics;
using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.Esa.Tests.Fuzz;

/// <summary>A size band of the differential fuzz: lengths in <c>[MinN, MaxN]</c>, uniformly or log-uniformly, and the number of cases of the acceptance run.</summary>
internal sealed record FuzzTier(string Name, int MinN, int MaxN, bool LogUniform, long Cases, ulong Salt)
{
    /// <summary>n = 0..64 (acceptance of work package 1a.7: 10^6 cases).</summary>
    internal static readonly FuzzTier Tiny = new("n<=64", 0, 64, false, 1_000_000, 1);

    /// <summary>n = 65..2000, naive oracle (acceptance: 10^5 cases).</summary>
    internal static readonly FuzzTier Small = new("n<=2000", 65, 2000, true, 100_000, 2);

    /// <summary>n = 2001..10^5, verifier and Kasai oracles (acceptance: 10^3 cases).</summary>
    internal static readonly FuzzTier Large = new("n<=1e5", 2001, 100_000, true, 1_000, 3);
}

/// <summary>Parameters of one fuzz run.</summary>
internal sealed class FuzzOptions
{
    /// <summary>The run seed; case i of a tier is a pure function of (MasterSeed, tier salt, i).</summary>
    internal ulong MasterSeed { get; init; } = 1;

    /// <summary>The size band.</summary>
    internal required FuzzTier Tier { get; init; }

    /// <summary>Number of cases to run (at most); defaults to the tier's acceptance count.</summary>
    internal long? Cases { get; init; }

    /// <summary>Stop starting new cases after this time (the run is then reported as budget-limited).</summary>
    internal TimeSpan? TimeBudget { get; init; }

    /// <summary>The system under test.</summary>
    internal EsaBuilder Candidate { get; init; } = FuzzOracle.ReferencePipeline;

    /// <summary>Optional reference engine for entry-by-entry equality (the optimized kernels are compared with <see cref="SaisReference"/>).</summary>
    internal EsaBuilder? Reference { get; init; }

    /// <summary>The run stops after this many failures.</summary>
    internal int MaxFailures { get; init; } = 5;

    /// <summary>Optional progress/diagnostic sink.</summary>
    internal Action<string>? Log { get; init; }
}

/// <summary>A failing case with its minimized text.</summary>
internal sealed record FuzzFailure(FuzzCase Case, string Op, string Detail, string MinimizedText);

/// <summary>Outcome of a fuzz run.</summary>
internal sealed class FuzzSummary
{
    internal required ulong MasterSeed { get; init; }

    internal required string TierName { get; init; }

    internal long Cases { get; set; }

    internal long TotalUnits { get; set; }

    internal TimeSpan Elapsed { get; set; }

    internal bool BudgetExhausted { get; set; }

    internal Dictionary<string, long> PerFamily { get; } = [];

    internal List<FuzzFailure> Failures { get; } = [];

    internal string FirstCaseId { get; set; } = "";

    internal string LastCaseId { get; set; } = "";

    public override string ToString() =>
        $"tier {TierName}, master seed {MasterSeed}, {Cases} cases ({TotalUnits} code units) in {Elapsed.TotalSeconds:F1} s{(BudgetExhausted ? " [time budget reached]" : "")}, " +
        $"{Failures.Count} failure(s); first case {FirstCaseId}, last case {LastCaseId}; families: " +
        string.Join(", ", PerFamily.OrderBy(kv => kv.Key, StringComparer.Ordinal).Select(kv => $"{kv.Key}={kv.Value}"));
}

/// <summary>The differential fuzz driver (work package 1a.7): seeded families x alphabets x sizes, a time budget, a minimizer for every failure.</summary>
internal static class FuzzRunner
{
    /// <summary>The case with index <paramref name="index"/> of a run: a pure function of the master seed, the tier and the index.</summary>
    internal static FuzzCase CaseAt(ulong masterSeed, FuzzTier tier, long index)
    {
        ulong h = masterSeed ^ (tier.Salt * 0xD1B54A32D192ED03UL) ^ ((ulong)index * 0xD6E8FEB86659FD93UL);
        var rng = new SplitMix64(new SplitMix64(h).Next()); // one hashing round: consecutive indices must not give shifted copies of one stream
        return FuzzTexts.Pick(rng, SampleLength(rng, tier));
    }

    /// <summary>Runs the fuzz and returns the summary; never throws for a mismatch (the failures are in the summary).</summary>
    internal static FuzzSummary Run(FuzzOptions options)
    {
        FuzzTier tier = options.Tier;
        long cases = options.Cases ?? tier.Cases;
        var summary = new FuzzSummary { MasterSeed = options.MasterSeed, TierName = tier.Name };
        long start = Stopwatch.GetTimestamp();
        long nextReport = 0;
        for (long i = 0; i < cases; i++)
        {
            if (options.TimeBudget is { } budget && Stopwatch.GetElapsedTime(start) >= budget)
            {
                summary.BudgetExhausted = true;
                break;
            }

            FuzzCase c = CaseAt(options.MasterSeed, tier, i);
            string text = c.Text();
            if (i == 0) summary.FirstCaseId = c.Id;
            summary.LastCaseId = c.Id;
            summary.Cases++;
            summary.TotalUnits += text.Length;
            summary.PerFamily[c.Family] = summary.PerFamily.GetValueOrDefault(c.Family) + 1;

            FuzzMismatch? mismatch = FuzzOracle.Check(text, options.Candidate, options.Reference);
            if (mismatch is { } m)
            {
                string minimized = FuzzMinimizer.Minimize(text,
                    t => FuzzOracle.Check(t, options.Candidate, options.Reference) is { } mm && mm.Class == m.Class);
                var failure = new FuzzFailure(c, m.Op, m.Detail, minimized);
                summary.Failures.Add(failure);
                options.Log?.Invoke($"FAIL {c.Id} op {m.Op}: {m.Detail}; minimized to {minimized.Length} unit(s): {FuzzRegressions.EncodeUnits(minimized)}");
                if (summary.Failures.Count >= options.MaxFailures) break;
            }

            if (options.Log is not null && Stopwatch.GetElapsedTime(start).Ticks >= nextReport)
            {
                options.Log($"[{tier.Name}] {summary.Cases}/{cases} cases, {summary.Failures.Count} failure(s), {Stopwatch.GetElapsedTime(start).TotalSeconds:F0} s");
                nextReport = Stopwatch.GetElapsedTime(start).Ticks + TimeSpan.FromSeconds(30).Ticks;
            }
        }
        summary.Elapsed = Stopwatch.GetElapsedTime(start);
        return summary;
    }

    private static int SampleLength(SplitMix64 rng, FuzzTier tier)
    {
        if (!tier.LogUniform) return tier.MinN + rng.NextInt(tier.MaxN - tier.MinN + 1);
        double u = (rng.Next() >> 11) * (1.0 / (1UL << 53));
        double v = Math.Exp(Math.Log(tier.MinN) + (u * (Math.Log(tier.MaxN + 1.0) - Math.Log(tier.MinN))));
        return Math.Clamp((int)v, tier.MinN, tier.MaxN);
    }
}
