// GENOMIC-SIMILARITY-001 (B09 A3-11) — sourmash ANI confidence intervals (containment_to_distance / jaccard_to_distance, Hera et al. 2023) and the `mash dist` comparison.
// Oracles: scripts/analysis/generate_ani_ci_oracle.py (sourmash 4.9.4: containment_ani(estimate_ci=True) both ways, max_containment_ani, jaccard_ani on 280 random pairs = 1120 rows)
// and scripts/analysis/generate_mash_dist_oracle.py (mash 2.3 `sketch` + `dist`: distance, p-value, shared / denominator on 400 pairs).
// Evidence: docs/Evidence/GENOMIC-SIMILARITY-001-Evidence.md.

using System.Globalization;

namespace Seqeron.Genomics.Tests.Unit.Analysis;

[TestFixture]
public class GenomicAnalyzer_AniConfidence_Oracle_Tests
{
    private static List<string[]> ReadRows(string file)
    {
        string name = "Seqeron.Genomics.Tests.TestData.Analysis." + file;
        using var stream = TestDataResource.Open(name)
            ?? throw new InvalidOperationException($"Embedded oracle '{name}' not found.");
        using var reader = new StreamReader(stream);
        var rows = new List<string[]>();
        string? line;
        while ((line = reader.ReadLine()) is not null)
            if (line.Length > 0) rows.Add(line.Split('\t'));
        return rows;
    }

    private static double D(string s) => double.Parse(s, CultureInfo.InvariantCulture);

    private static double? Opt(string s) => s == "." ? null : D(s);

    private static void AssertClose(double? actual, double? expected, string label, double relative = 1e-9)
    {
        if (expected is null) { Assert.That(actual, Is.Null, label); return; }
        Assert.That(actual, Is.Not.Null, label);
        Assert.That(actual!.Value, Is.EqualTo(expected.Value).Within(relative * Math.Max(1.0, Math.Abs(expected.Value))), label);
    }

    [Test]
    public void AniEstimates_MatchSourmash_AllRows()
    {
        var rows = ReadRows("ani_ci_oracle.tsv");
        Assert.That(rows, Has.Count.GreaterThanOrEqualTo(1000));
        int withInterval = 0, nullInterval = 0, inaccurate = 0, exceeds = 0, jacExceeds = 0;
        foreach (var f in rows)
        {
            string seqA = f[0], seqB = f[1], kind = f[5];
            int k = int.Parse(f[2], CultureInfo.InvariantCulture);
            long scaled = long.Parse(f[3], CultureInfo.InvariantCulture);
            double confidence = D(f[4]);
            var est = GenomicAnalyzer.CompareSimilarity(seqA, seqB, k, canonical: true, scaled: scaled, aniOptions: new AniOptions(true, confidence));
            var all = est.AniEstimates!;
            var r = kind switch { "c_ab" => all.ContainmentAInB, "c_ba" => all.ContainmentBInA, "max" => all.MaxContainment, _ => all.Jaccard };
            string label = $"{kind} k={k} S={scaled} conf={confidence} {seqA.Length}/{seqB.Length}";
            AssertClose(r.Distance, D(f[6]), label + " dist");
            if (kind != "jac")
            {
                AssertClose(r.DistanceLow, Opt(f[7]), label + " dist_low");
                AssertClose(r.DistanceHigh, Opt(f[8]), label + " dist_high");
                if (r.DistanceLow is not null) withInterval++; else nullInterval++;
            }

            AssertClose(r.ProbabilityNothingInCommon, D(f[9]), label + " p_nothing", 1e-8);
            Assert.That(r.ProbabilityExceedsThreshold, Is.EqualTo(f[10] == "1"), label + " p_exceeds");
            Assert.That(r.SizeIsInaccurate, Is.EqualTo(f[11] == "1"), label + " size_is_inaccurate");
            if (kind == "jac")
            {
                AssertClose(r.JaccardError, Opt(f[12]), label + " jaccard_error", 1e-8);
                Assert.That(r.JaccardErrorExceedsThreshold, Is.EqualTo(f[13] == "1"), label + " je_exceeds");
                if (r.JaccardErrorExceedsThreshold) jacExceeds++;
            }

            if (r.SizeIsInaccurate) inaccurate++;
            if (r.ProbabilityExceedsThreshold) exceeds++;
        }

        Assert.That(withInterval, Is.GreaterThan(500));
        Assert.That(inaccurate, Is.GreaterThan(100));
        Assert.That(exceeds + jacExceeds, Is.GreaterThan(10));
        Assert.That(nullInterval, Is.GreaterThanOrEqualTo(0));
    }

    [Test]
    public void MashDist_DistancePValueAndSharedHashes_MatchMash_AllRows()
    {
        var rows = ReadRows("mash_dist_oracle.tsv");
        Assert.That(rows, Has.Count.GreaterThanOrEqualTo(200));
        int tinyP = 0, lowK = 0, noncanonical = 0;
        foreach (var f in rows)
        {
            int k = int.Parse(f[2], CultureInfo.InvariantCulture), s = int.Parse(f[3], CultureInfo.InvariantCulture);
            bool canonical = f[4] == "1";
            var est = GenomicAnalyzer.CompareSimilarity(f[0], f[1], k, canonical, sketchSize: s);
            string label = $"k={k} s={s} canonical={canonical} {f[0].Length}/{f[1].Length}";
            Assert.That(est.SharedHashes, Is.EqualTo(int.Parse(f[7], CultureInfo.InvariantCulture)), label + " shared");
            Assert.That(est.DenominatorHashes, Is.EqualTo(int.Parse(f[8], CultureInfo.InvariantCulture)), label + " denominator");
            double dist = D(f[5]), p = D(f[6]);
            Assert.That(est.MashDistance, Is.EqualTo(dist).Within(1e-5 * Math.Max(1e-3, dist)), label + " distance");
            Assert.That(est.MashPValue!.Value, Is.EqualTo(p).Within(1e-5 * p + 1e-300), label + " p-value");
            if (p < 1e-20) tinyP++;
            if (k <= 16) lowK++;
            if (!canonical) noncanonical++;
        }

        Assert.That(tinyP, Is.GreaterThan(50));
        Assert.That(lowK, Is.GreaterThan(100));
        Assert.That(noncanonical, Is.GreaterThan(60));
    }

    [Test]
    public void NormalQuantile_EqualsScipyNormPpf()
    {
        // scipy.stats.norm.ppf (4.9.4 environment)
        (double P, double Z)[] cases =
        {
            (0.5, 0.0), (0.975, 1.959963984540054), (0.995, 2.5758293035489004), (0.95, 1.6448536269514722), (1e-10, -6.361340902404056), (0.001, -3.090232306167813),
            (0.3, -0.5244005127080409), (0.9999, 3.719016485455709), (0.0001, -3.7190164854556804), (0.9, 1.2815515655446004), (0.99999999, 5.612001243305505),
            (1e-300, -37.0470962993612), (0.2, -0.8416212335729142), (0.7, 0.5244005127080407), (0.425, -0.18911842627279252), (0.075, -1.4395314709384563),
        };
        foreach (var (p, z) in cases)
            Assert.That(StatisticsHelper.NormalQuantile(p), Is.EqualTo(z).Within(1e-12 * Math.Max(1.0, Math.Abs(z))), $"p = {p}");
        Assert.That(StatisticsHelper.NormalQuantile(0.0), Is.EqualTo(double.NegativeInfinity));
        Assert.That(StatisticsHelper.NormalQuantile(1.0), Is.EqualTo(double.PositiveInfinity));
        Assert.Throws<ArgumentOutOfRangeException>(() => StatisticsHelper.NormalQuantile(1.5));
        Assert.Throws<ArgumentOutOfRangeException>(() => StatisticsHelper.NormalQuantile(double.NaN));
    }

    [Test]
    public void BrentRoot_FindsRoots_AndRefusesAnUnbracketedFunction()
    {
        // scipy.optimize.brentq(lambda x: x**3 - 2, 0, 3) = 1.2599210498948739; brentq(cos(x) - x, 0, 1) = 0.7390851332151559 (the same bracket gives the same digits)
        Assert.That(StatisticsHelper.BrentRoot(x => x * x * x - 2, 0, 3), Is.EqualTo(1.2599210498948739).Within(1e-14));
        Assert.That(StatisticsHelper.BrentRoot(x => Math.Cos(x) - x, 0, 1), Is.EqualTo(0.7390851332151559).Within(1e-14));
        Assert.That(StatisticsHelper.BrentRoot(x => x - 1, 1, 5), Is.EqualTo(1.0));
        Assert.Throws<ArgumentException>(() => StatisticsHelper.BrentRoot(x => x * x + 1, -1, 1));
    }

    [Test]
    public void BinomialCdf_EqualsScipy()
    {
        // scipy.stats.binom.cdf
        Assert.That(StatisticsHelper.BinomialCdf(3, 10, 0.3), Is.EqualTo(0.6496107184).Within(1e-9));
        Assert.That(StatisticsHelper.BinomialCdf(3.7, 10, 0.3), Is.EqualTo(0.6496107184).Within(1e-9));
        Assert.That(StatisticsHelper.BinomialCdf(-1, 10, 0.3), Is.EqualTo(0.0));
        Assert.That(StatisticsHelper.BinomialCdf(10, 10, 0.3), Is.EqualTo(1.0));
    }

    [Test]
    public void AniEstimate_WorkedExample_RunOnSourmash()
    {
        // sourmash.distance_utils.containment_to_distance(0.7, 21, 10, n_unique_kmers=50000, estimate_ci=True): dist 0.016841097303797903, low 0.015448819645416393, high 0.018329203382171456
        var r = GenomicAnalyzer.ContainmentToDistance(0.7, 21, 10, 50000, estimateConfidenceInterval: true);
        Assert.That(r.Distance, Is.EqualTo(0.016841097303797903).Within(1e-15));
        Assert.That(r.DistanceLow, Is.EqualTo(0.015448819645416393).Within(1e-11));
        Assert.That(r.DistanceHigh, Is.EqualTo(0.018329203382171456).Within(1e-11));
        Assert.That(r.AniLow, Is.EqualTo(1.0 - r.DistanceHigh));
        Assert.That(r.Ani, Is.EqualTo(1.0 - r.Distance));
        // C = 0 / C = 1
        var zero = GenomicAnalyzer.ContainmentToDistance(0, 21, 10, 5000, estimateConfidenceInterval: true);
        Assert.That((zero.Distance, zero.DistanceLow, zero.DistanceHigh, zero.ProbabilityNothingInCommon), Is.EqualTo((1.0, (double?)1.0, (double?)1.0, 1.0)));
        var one = GenomicAnalyzer.ContainmentToDistance(1, 21, 10, 5000, estimateConfidenceInterval: true);
        Assert.That((one.Distance, one.DistanceLow, one.DistanceHigh), Is.EqualTo((0.0, (double?)0.0, (double?)0.0)));
        Assert.Throws<ArgumentOutOfRangeException>(() => GenomicAnalyzer.ContainmentToDistance(1.5, 21, 10, 5000));
        Assert.Throws<ArgumentOutOfRangeException>(() => GenomicAnalyzer.ContainmentToDistance(0.5, 21, 10, 5000, confidence: 1.0));
    }
}
