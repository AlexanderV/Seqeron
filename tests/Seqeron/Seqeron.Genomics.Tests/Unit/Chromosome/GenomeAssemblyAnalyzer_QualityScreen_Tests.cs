using Seqeron.Genomics.Analysis;
using static Seqeron.Genomics.Chromosome.GenomeAssemblyAnalyzer;

namespace Seqeron.Genomics.Tests.Unit.Chromosome;

/// <summary>
/// ASSEMBLY-QC-001: <see cref="CalculateLocalQuality"/> (windowed GC / N count / linguistic complexity)
/// and <see cref="FindSuspiciousRegions"/> (N-rich windows, GC-outlier windows, SDUST intervals).
/// Expected values come from tools/review/oracles/B10/qc_oracle.py: Biopython <c>gc_fraction</c> windows,
/// the Orlov-Potapov/Rosalind LING linguistic complexity (ΣU_i/ΣUmax_i, m = 10) and a port of lh3/sdust <c>sdust_core</c>.
/// </summary>
[TestFixture]
public class GenomeAssemblyAnalyzer_QualityScreen_Tests
{
    private const double Tol = 1e-9;

    /// <summary>Deterministic generator mirrored by qc_oracle.py <c>rand_dna</c>.</summary>
    private static string Rand(long seed, int n, double gc = 0.5)
    {
        long x = seed;
        var chars = new char[n];
        for (int i = 0; i < n; i++)
        {
            x = (x * 1103515245L + 12345L) & 0x7FFFFFFFL;
            double u = ((x >> 8) & 0xFFFF) / 65536.0;
            x = (x * 1103515245L + 12345L) & 0x7FFFFFFFL;
            bool bit = ((x >> 16) & 1) == 1;
            chars[i] = u < gc ? (bit ? 'G' : 'C') : (bit ? 'A' : 'T');
        }
        return new string(chars);
    }

    #region CalculateLocalQuality

    [Test]
    public void CalculateLocalQuality_MatchesBiopythonGcAndLinguisticComplexityOracle()
    {
        (int Pos, int Len, double Gc, int N, double Lc)[] expected =
        {
            (0, 1000, 0.529, 0, 0.913347087763847),
            (500, 1000, 0.516, 0, 0.9154102523408982),
            (1000, 500, 0.49, 0, 0.9509239624356256),
        };
        var got = CalculateLocalQuality(new[] { ("s", Rand(12, 1500)) }, 1000).ToList();
        Assert.That(got, Has.Count.EqualTo(expected.Length));
        for (int i = 0; i < expected.Length; i++)
        {
            Assert.That(got[i].Position, Is.EqualTo(expected[i].Pos));
            Assert.That(got[i].WindowSize, Is.EqualTo(expected[i].Len));
            Assert.That(got[i].GcContent, Is.EqualTo(expected[i].Gc).Within(Tol));
            Assert.That(got[i].NCount, Is.EqualTo(expected[i].N));
            Assert.That(got[i].Complexity, Is.EqualTo(expected[i].Lc).Within(Tol));
        }
    }

    [Test]
    public void CalculateLocalQuality_NRunAndExplicitStep_MatchOracle()
    {
        (int Pos, int Len, double Gc, int N, double Lc)[] expected =
        {
            (0, 300, 0.4866666666666667, 0, 0.9500237981913374),
            (100, 300, 0.47, 0, 0.9504997620180866),
            (200, 300, 0.48, 50, 0.7844554902846814),
            (300, 300, 0.504, 50, 0.7830998644374153),
            (400, 300, 0.544, 50, 0.7677361048350655),
            (500, 250, 0.508, 0, 0.9140625),
            (600, 150, 0.5133333333333333, 0, 0.9386446886446886),
            (700, 50, 0.44, 0, 0.9574468085106383),
        };
        string seq = Rand(13, 400) + new string('N', 50) + Rand(14, 300);
        var got = CalculateLocalQuality(new[] { ("s", seq) }, 300, stepSize: 100).ToList();
        Assert.That(got, Has.Count.EqualTo(expected.Length));
        for (int i = 0; i < expected.Length; i++)
        {
            Assert.That(got[i].Position, Is.EqualTo(expected[i].Pos), $"window {i}");
            Assert.That(got[i].WindowSize, Is.EqualTo(expected[i].Len), $"window {i}");
            Assert.That(got[i].GcContent, Is.EqualTo(expected[i].Gc).Within(Tol), $"window {i}");
            Assert.That(got[i].NCount, Is.EqualTo(expected[i].N), $"window {i}");
            Assert.That(got[i].Complexity, Is.EqualTo(expected[i].Lc).Within(Tol), $"window {i}");
        }
    }

    [Test]
    public void CalculateLocalQuality_GcDenominatorExcludesNAndIupac_AsQuast()
    {
        // F10 / QUAST: GC = (G+C)/(A+C+G+T). "ACGTRYSWNN": G,C = 2 of A,C,G,T = 4 → 0.5
        // (Biopython 'remove' would give 3/6 because it counts S and W; the all-IUPAC reading 2/10).
        var q = CalculateLocalQuality(new[] { ("s", "ACGTRYSWNN") }).Single();
        Assert.That(q.GcContent, Is.EqualTo(0.5).Within(Tol));
        Assert.That(q.NCount, Is.EqualTo(2));
    }

    [Test]
    public void CalculateLocalQuality_LowercaseNAndBases_AreCounted()
    {
        var q = CalculateLocalQuality(new[] { ("s", "ggccnnaa") }).Single();
        Assert.That(q.GcContent, Is.EqualTo(4.0 / 6.0).Within(Tol));
        Assert.That(q.NCount, Is.EqualTo(2));
    }

    [Test]
    public void CalculateLocalQuality_WindowOfOne_TerminatesWithOneWindowPerBase()
    {
        // windowSize 1 used to give step 1/2 = 0 and never advance.
        var got = CalculateLocalQuality(new[] { ("s", "ACGN") }, windowSize: 1).ToList();
        Assert.That(got.Select(w => w.Position), Is.EqualTo(new[] { 0, 1, 2, 3 }));
        Assert.That(got.Select(w => w.GcContent), Is.EqualTo(new[] { 0.0, 1.0, 1.0, 0.0 }));
        Assert.That(got.Select(w => w.NCount), Is.EqualTo(new[] { 0, 0, 0, 1 }));
    }

    [Test]
    public void CalculateLocalQuality_EmptySequence_YieldsNothing()
    {
        Assert.That(CalculateLocalQuality(new[] { ("s", "") }).ToList(), Is.Empty);
    }

    [Test]
    public void CalculateLocalQuality_InvalidArguments_Throw()
    {
        var seqs = new[] { ("s", "ACGT") };
        Assert.Throws<ArgumentNullException>(() => CalculateLocalQuality(null!));
        Assert.Throws<ArgumentOutOfRangeException>(() => CalculateLocalQuality(seqs, 0));
        Assert.Throws<ArgumentOutOfRangeException>(() => CalculateLocalQuality(seqs, -5));
        Assert.Throws<ArgumentOutOfRangeException>(() => CalculateLocalQuality(seqs, 10, 0));
        Assert.Throws<ArgumentException>(() => CalculateLocalQuality(new[] { ("s", (string)null!) }).ToList());
    }

    [Test]
    public void CalculateLocalQuality_GcIsTheCanonicalWindowedGcProfile()
    {
        // No private GC loop: equals GcSkewCalculator.CalculateWindowedGcContent(includePartialWindow: true).
        string seq = Rand(21, 777, 0.6);
        var q = CalculateLocalQuality(new[] { ("s", seq) }, 100, 30).ToList();
        var core = GcSkewCalculator.CalculateWindowedGcContent(seq, 100, 30, fraction: true, includePartialWindow: true).ToList();
        Assert.That(q.Select(w => w.GcContent), Is.EqualTo(core.Select(c => c.GcContent)));
        Assert.That(q.Select(w => w.Position), Is.EqualTo(core.Select(c => c.WindowStart)));
    }

    #endregion

    #region GcSkewCalculator.CalculateWindowedGcContent(includePartialWindow)

    [Test]
    public void WindowedGcContent_IncludePartialWindow_TruncatesTrailingWindows()
    {
        // Biopython-style slices seq[i:i+w] for every i = 0, step, … < len.
        const string seq = "GGGGAAAACC"; // len 10
        var full = GcSkewCalculator.CalculateWindowedGcContent(seq, 4, 3, fraction: true, includePartialWindow: false).ToList();
        var part = GcSkewCalculator.CalculateWindowedGcContent(seq, 4, 3, fraction: true, includePartialWindow: true).ToList();

        Assert.That(full.Select(p => p.WindowStart), Is.EqualTo(new[] { 0, 3, 6 }));
        Assert.That(part.Select(p => p.WindowStart), Is.EqualTo(new[] { 0, 3, 6, 9 }));
        Assert.That(part.Select(p => p.WindowEnd), Is.EqualTo(new[] { 3, 6, 9, 9 }));
        Assert.That(part.Select(p => p.Position), Is.EqualTo(new[] { 2, 5, 8, 9 }));
        Assert.That(part.Select(p => p.GcContent), Is.EqualTo(new[] { 1.0, 0.25, 0.5, 1.0 }).Within(Tol));
        Assert.That(part.Take(3).Select(p => p.GcContent), Is.EqualTo(full.Select(p => p.GcContent)));
    }

    [Test]
    public void WindowedGcContent_IncludePartialWindow_ShortSequenceGivesOneWindow()
    {
        var p = GcSkewCalculator.CalculateWindowedGcContent("GCAT", 1000, 500, fraction: false, includePartialWindow: true).Single();
        Assert.That((p.WindowStart, p.WindowEnd, p.Position), Is.EqualTo((0, 3, 2)));
        Assert.That(p.GcContent, Is.EqualTo(50.0).Within(Tol));
        Assert.That(GcSkewCalculator.CalculateWindowedGcContent("GCAT", 1000, 500, false, false), Is.Empty);
    }

    #endregion

    #region FindSuspiciousRegions

    private static readonly (string Id, int Start, int End, string Reason, double Gc, double N, double Dust)[] PolyAInRandom =
    {
        ("s", 0, 749, "GC deviation", 0.16883333333333334, 0.0, 0.0),
        ("s", 1000, 4000, "GC deviation; Low complexity", 0.33316666666666667, 0.0, 999.5),
        ("s", 4750, 5999, "GC deviation", 0.21483333333333338, 0.0, 0.0),
    };

    private static void AssertRegions(
        IEnumerable<SuspiciousRegion> got, (string Id, int Start, int End, string Reason, double Gc, double N, double Dust)[] expected)
    {
        var list = got.ToList();
        Assert.That(list, Has.Count.EqualTo(expected.Length));
        for (int i = 0; i < expected.Length; i++)
        {
            Assert.That(list[i].SequenceId, Is.EqualTo(expected[i].Id), $"region {i}");
            Assert.That(list[i].Start, Is.EqualTo(expected[i].Start), $"region {i}");
            Assert.That(list[i].End, Is.EqualTo(expected[i].End), $"region {i}");
            Assert.That(list[i].Reason, Is.EqualTo(expected[i].Reason), $"region {i}");
            Assert.That(list[i].MaxGcDeviation, Is.EqualTo(expected[i].Gc).Within(Tol), $"region {i}");
            Assert.That(list[i].MaxNFraction, Is.EqualTo(expected[i].N).Within(Tol), $"region {i}");
            Assert.That(list[i].MaxDustScore, Is.EqualTo(expected[i].Dust).Within(Tol), $"region {i}");
        }
    }

    [Test]
    public void FindSuspiciousRegions_PolyAInRandom_MatchesOracle()
        => AssertRegions(FindSuspiciousRegions(new[] { ("s", Rand(1, 2000) + new string('A', 2000) + Rand(2, 2000)) }), PolyAInRandom);

    [Test]
    public void FindSuspiciousRegions_SmallPolyA_MatchesOracle()
        => AssertRegions(
            FindSuspiciousRegions(new[] { ("s", Rand(15, 3000) + new string('A', 300) + Rand(16, 3000)) }),
            new[]
            {
                ("s", 1711, 1717, "Low complexity", 0.0, 0.0, 2.5),
                ("s", 2750, 3499, "GC deviation; Low complexity", 0.2568571428571429, 0.0, 149.0),
                ("s", 5372, 5379, "Low complexity", 0.0, 0.0, 3.0),
                ("s", 6237, 6244, "Low complexity", 0.0, 0.0, 3.0),
            });

    [Test]
    public void FindSuspiciousRegions_GcIsland_MatchesOracle()
        => AssertRegions(
            FindSuspiciousRegions(new[] { ("s", Rand(3, 3000, 0.4) + Rand(4, 1000, 0.95) + Rand(5, 3000, 0.4)) }),
            new[]
            {
                ("s", 2473, 2479, "Low complexity", 0.0, 0.0, 2.5),
                ("s", 2750, 4249, "GC deviation; Low complexity", 0.476, 0.0, 48.545271629778675),
                ("s", 4851, 4858, "Low complexity", 0.0, 0.0, 3.0),
            });

    [Test]
    public void FindSuspiciousRegions_NRun_MatchesOracle()
        => AssertRegions(
            FindSuspiciousRegions(new[] { ("s", Rand(6, 1500) + new string('N', 400) + Rand(7, 1500)) }),
            new[] { ("s", 1250, 2249, "High N content; Low complexity", 0.0, 0.8, 2.1) });

    [Test]
    public void FindSuspiciousRegions_MultipleSequences_PoolsGcAndReportsPerSequence()
        => AssertRegions(
            FindSuspiciousRegions(new[]
            {
                ("a", Rand(8, 1200, 0.45)),
                ("b", Rand(9, 1100, 0.45) + string.Concat(Enumerable.Repeat("AC", 150)) + Rand(10, 900, 0.45)),
            }),
            new[] { ("b", 1100, 1399, "Low complexity", 0.0, 0.0, 74.24915824915826) });

    [Test]
    public void FindSuspiciousRegions_PolyA600_SingleDustIntervalWithExactScore()
    {
        // 598 identical triplets: DUST = c(c−1)/2/(ℓ−1) = 598·597/2/597 = 299 (Morgulis 2006).
        var r = FindSuspiciousRegions(new[] { ("s", new string('A', 600)) }).Single();
        Assert.That((r.Start, r.End, r.Reason), Is.EqualTo((0, 599, "Low complexity")));
        Assert.That(r.MaxDustScore, Is.EqualTo(299.0).Within(Tol));
        Assert.That(r.MaxGcDeviation, Is.EqualTo(0.0));
        Assert.That(r.MaxNFraction, Is.EqualTo(0.0));
    }

    [Test]
    public void FindSuspiciousRegions_RegionEndIsLastFlaggedBase_NotNextCleanWindow()
    {
        // Old code ended a region at the start of the first clean window − 1, i.e. inside the last
        // flagged (overlapping) window. N run [250,349] ⇒ windows 0,250 (and 100.. etc.) flagged → union end.
        string seq = Rand(31, 250) + new string('N', 100) + Rand(32, 650);
        var regions = FindSuspiciousRegions(new[] { ("s", seq) }, windowSize: 200).ToList();
        var n = regions.Single(r => r.Reason.Contains("High N content"));
        // windows (step 100): [100,299]=0 N… [200,399] 50 N… [300,499] 50 N... ; 100/200=0.5 [250..]:
        Assert.That(n.Start, Is.LessThanOrEqualTo(250));
        Assert.That(n.End, Is.GreaterThanOrEqualTo(349));
        Assert.That(n.MaxNFraction, Is.GreaterThan(0.1));
    }

    [Test]
    public void FindSuspiciousRegions_GcDeviationUsesPooledGcOverAllSequences()
    {
        // 1000 bp at 100 % GC pooled with 3000 bp at 0 % GC ⇒ pooled 0.25; every window deviates ≥ 0.25 > 0.15.
        var got = FindSuspiciousRegions(new[] { ("gc", new string('G', 500) + new string('C', 500)), ("at", Rand(41, 3000, 0.0)) },
            dustThreshold: 1e6).ToList();
        Assert.That(got.Select(r => (r.SequenceId, r.Start, r.End)), Is.EqualTo(new[] { ("gc", 0, 999), ("at", 0, 2999) }));
        Assert.That(got.All(r => r.Reason == "GC deviation"), Is.True);
        Assert.That(got[0].MaxGcDeviation, Is.EqualTo(0.75).Within(Tol));
        Assert.That(got[1].MaxGcDeviation, Is.EqualTo(0.25).Within(Tol));
    }

    [Test]
    public void FindSuspiciousRegions_WindowWithoutAcgt_IsNotAGcOutlier()
    {
        // An all-N window has no GC (0 by the canonical helper); it must be flagged only for N content.
        string seq = Rand(51, 1000) + new string('N', 600) + Rand(52, 1000);
        var regions = FindSuspiciousRegions(new[] { ("s", seq) }, dustThreshold: 1e6).ToList();
        Assert.That(regions.All(r => r.Reason == "High N content"), Is.True);
    }

    [Test]
    public void FindSuspiciousRegions_TruncatedTrailingWindowCoveredByPreviousWindow_IsNotTested()
    {
        // len 1100, window 500, step 250: starts 0,250,500,750,1000; the window at 1000 (100 bp) is inside
        // the one at 750 [750,1249]→[750,1099] and is skipped, so a GC-extreme 100-bp tail that does not
        // shift the 350-bp-diluted full window is not flagged on its own.
        string seq = Rand(61, 1000, 0.5) + new string('G', 100);
        var regions = FindSuspiciousRegions(new[] { ("s", seq) }, gcDeviation: 0.2, dustThreshold: 1e6).ToList();
        Assert.That(regions.All(r => r.End != 1099 || r.Start <= 750), Is.True);
        Assert.That(regions.Where(r => r.Start >= 1000), Is.Empty);
    }

    [Test]
    public void FindSuspiciousRegions_SequenceShorterThanWindow_IsScoredAsOneWindow()
    {
        var regions = FindSuspiciousRegions(new[] { ("s", new string('N', 50)) }, windowSize: 500, dustThreshold: 1e6).ToList();
        var r = regions.Single();
        Assert.That((r.Start, r.End, r.Reason, r.MaxNFraction), Is.EqualTo((0, 49, "High N content", 1.0)));
    }

    [Test]
    public void FindSuspiciousRegions_NoValidBase_NoGcTest_AndEmptyInputYieldsNothing()
    {
        Assert.That(FindSuspiciousRegions(new[] { ("s", "") }), Is.Empty);
        Assert.That(FindSuspiciousRegions(Array.Empty<(string, string)>()), Is.Empty);
        Assert.That(FindSuspiciousRegions(new[] { ("s", "RYKM") }), Is.Empty);
    }

    [Test]
    public void FindSuspiciousRegions_DustThreshold_IsPassedToSdust()
    {
        // (ACGT)n scores DUST = 1.0: not > 2.0 (dustmasker level 20) but > 0.5.
        string seq = string.Concat(Enumerable.Repeat("ACGTTGCA", 100));
        Assert.That(SequenceComplexity.CalculateDustScore(seq), Is.GreaterThan(0.5));
        var strict = FindSuspiciousRegions(new[] { ("s", seq) }, dustThreshold: 1000).ToList();
        var loose = FindSuspiciousRegions(new[] { ("s", seq) }, dustThreshold: 0.5).ToList();
        Assert.That(strict.Any(r => r.Reason.Contains("Low complexity")), Is.False);
        Assert.That(loose.Any(r => r.Reason.Contains("Low complexity")), Is.True);
    }

    [Test]
    public void FindSuspiciousRegions_ReasonsAreJoinedInFixedOrder()
    {
        // N-rich AND GC-deviant AND low-complexity windows overlap in one region.
        string seq = Rand(71, 2000, 0.5) + new string('A', 300) + new string('N', 300) + Rand(72, 2000, 0.5);
        var r = FindSuspiciousRegions(new[] { ("s", seq) }, windowSize: 200).ToList();
        Assert.That(r.Any(x => x.Reason == "High N content; GC deviation; Low complexity" || x.Reason == "GC deviation; Low complexity"
            || x.Reason == "High N content; Low complexity"), Is.True);
        Assert.That(r.All(x => !x.Reason.Contains("Low complexity; ")), Is.True, "order is N, GC, Low complexity");
    }

    [Test]
    public void FindSuspiciousRegions_InvalidArguments_Throw()
    {
        var s = new[] { ("s", "ACGT") };
        Assert.Throws<ArgumentNullException>(() => FindSuspiciousRegions(null!));
        Assert.Throws<ArgumentOutOfRangeException>(() => FindSuspiciousRegions(s, gcDeviation: -0.1));
        Assert.Throws<ArgumentOutOfRangeException>(() => FindSuspiciousRegions(s, gcDeviation: double.NaN));
        Assert.Throws<ArgumentOutOfRangeException>(() => FindSuspiciousRegions(s, dustThreshold: -1));
        Assert.Throws<ArgumentOutOfRangeException>(() => FindSuspiciousRegions(s, maxNFraction: 1.5));
        Assert.Throws<ArgumentOutOfRangeException>(() => FindSuspiciousRegions(s, windowSize: 0));
        Assert.Throws<ArgumentException>(() => FindSuspiciousRegions(new[] { ("s", (string)null!) }));
    }

    #endregion

    #region CalculateLengthDistribution (F13 re-check)

    [Test]
    public void CalculateLengthDistribution_HalfOpenBinsAndOverflow()
    {
        // numpy.histogram-style half-open bins [prev, b); lengths equal to a bound go to the NEXT bin.
        var d = CalculateLengthDistribution(new[] { 0, 99, 100, 499, 500, 999, 1000 }, 100, 500, 1000);
        Assert.That(d["<100"], Is.EqualTo(2));
        Assert.That(d["<500"], Is.EqualTo(2));
        Assert.That(d["<1000"], Is.EqualTo(2));
        Assert.That(d[">=1000"], Is.EqualTo(1));
        Assert.That(d.Values.Sum(), Is.EqualTo(7));
    }

    #endregion
}
