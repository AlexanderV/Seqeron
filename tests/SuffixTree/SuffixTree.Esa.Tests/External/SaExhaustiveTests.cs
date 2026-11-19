using SuffixTree.Esa.External;
using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.Esa.Tests.External;

/// <summary>
/// Heavy checks of work package 1b.7d (gate G1). Category <c>Exhaustive</c>, <c>[Explicit]</c>:
/// <c>TMPDIR=/dev/shm dotnet test tests/SuffixTree/SuffixTree.Esa.Tests/SuffixTree.Esa.Tests.csproj -c Release --filter FullyQualifiedName~SaExhaustiveTests</c> (the runs create and map tens of small files each: a memory
/// file system makes them cheap). Every string over sigma 1 up to length 14, sigma 2 up to 12, sigma 3 up to 7, sigma 4 up to 5 goes through PH1 and through the whole stage with rotating tiny parameters (the period, K_s, N_base, K,
/// the fan-out, the oversample, the windows, the threads); the families go through the whole stage at 10^5 and 10^6 symbols with several bucket capacities and thread counts. Every output is compared with the oracles.
/// </summary>
[TestFixture]
[Category("Exhaustive")]
[Explicit("Heavy enumeration; run with --filter FullyQualifiedName~SaExhaustiveTests (see PROGRESS-S1.md, Log 1b.7d)")]
public class SaExhaustiveTests
{
    private static readonly (int Sigma, int Max)[] Tiers = [(1, 14), (2, 12), (3, 7), (4, 5)];

    [Test]
    [CancelAfter(3_600_000)]
    public void EveryString_ThroughPh1_IsTheRankFileOfTheDefinition()
    {
        int count = 0;
        foreach ((int sigma, int max) in Tiers)
        {
            foreach (string s in ExhaustiveStrings.All(sigma, max))
            {
                if (s.Length == 0) continue;
                SaCase c = SaCase.Of(s);
                using var dir = new TempDir();
                int v = count++;
                DifferenceCover cover = DifferenceCover.Create(SaTiny.Periods[v % SaTiny.Periods.Length]);
                ExternalDcsTests.RunAndCompare(dir, c, cover, SaTiny.Dcs(dir.Sub("scratch"), v, 1 + (v % 4)), $"'{s}' variant {v}", out _); // strings of up to 14 symbols do not recurse: the levels are exercised by the families, ExternalDcsTests and SaPinnedTests
            }
        }
        Assert.That(count, Is.EqualTo(14 + 8190 + 3279 + 1364), "unary 1-14, binary 1-12, ternary 1-7, 4-ary 1-5");
    }

    [Test]
    [CancelAfter(3_600_000)]
    public void EveryString_ThroughTheWholeStage_IsTheSuffixArray()
    {
        int count = 0;
        foreach ((int sigma, int max) in Tiers)
        {
            foreach (string s in ExhaustiveStrings.All(sigma, max))
            {
                if (s.Length == 0) continue;
                SaCase c = SaCase.Of(s);
                using var dir = new TempDir();
                int v = count++;
                SaRun.Outcome o = SaRun.Run(dir, c, SaTiny.Build(dir.Sub("scratch"), v, 1 + (v % 4)));
                SaRun.Compare(c, o, $"'{s}' variant {v}");
            }
        }
        Assert.That(count, Is.EqualTo(14 + 8190 + 3279 + 1364));
    }

    [Test]
    [CancelAfter(3_600_000)]
    public void TheFamilies_AtTenToTheFifth_AndTenToTheSixth_WithSeveralBucketCapacitiesAndThreadCounts()
    {
        int v = 0;
        foreach (int n in new[] { 100_000, 1_000_000 })
        {
            foreach ((string name, string text) in SaRun.Families(n))
            {
                SaCase c = SaCase.Of(text);
                string? reference = null;
                foreach ((long k, int threads, long scale) in new[] { (n / 7L, 4, 0L), (1L << 20, 1, 0L), (n / 3L + 1, 2, 0L), (0L, 4, 3_000_000_000L) })
                {
                    using var dir = new TempDir();
                    var o = new SaBuildOptions { ScratchDir = dir.Sub("scratch"), Threads = threads, ScaleN = scale, Tuning = k > 0 ? new SaBuildTuning { BucketCapacity = k } : new SaBuildTuning() };
                    SaRun.Outcome r = SaRun.Run(dir, c, o);
                    string sha = SaRun.Sha(r.SaPath);
                    if (reference is null) SaRun.Compare(c, r, $"{name} {n} K={k}");
                    reference ??= sha;
                    Assert.That(sha, Is.EqualTo(reference), $"{name} {n}: K {k}, {threads} threads, scale {scale}");
                    v++;
                }
            }
        }
        Assert.That(v, Is.EqualTo(2 * 9 * 4));
    }
}
