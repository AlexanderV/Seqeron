namespace Seqeron.Genomics.Tests.Fuzzing;

/// <summary>
/// FIN-B24 heavy tier — fuzz tests for the input contracts added by the B24 finisher: R <c>NA</c> handling in
/// <c>ascat.aspcf</c> / <c>ascat.asmultipcf</c> (F59: NaN is accepted and follows R's NA path), the asmultipcf
/// one-sample / probe-set errors (F53: R fails with "incorrect number of subscripts on matrix"), the GISTIC2 ziggurat cap
/// and noise filter (F50–F54), the canonical copy-number guards (F28/F66) and the Student t / incomplete beta domain
/// (F39). Random inputs are locally seeded; every run either succeeds with a well-formed result or throws the documented
/// exception type.
/// </summary>
[TestFixture]
[Category("Fuzzing")]
[Category("Oncology")]
public class OncologyFinisherFuzzTests
{
    private static List<OncologyAnalyzer.AlleleSpecificLocus> NoisyLoci(Random rng, int chromosomes, int perChromosome, double naRate)
    {
        var loci = new List<OncologyAnalyzer.AlleleSpecificLocus>();
        for (int c = 1; c <= chromosomes; c++)
        {
            for (int i = 0; i < perChromosome; i++)
            {
                double logR = rng.NextDouble() < naRate ? double.NaN : (rng.NextDouble() - 0.5) * (rng.Next(10) == 0 ? 20 : 1);
                double baf = rng.NextDouble() < naRate ? double.NaN : rng.Next(8) switch { 0 => 0.0, 1 => 1.0, _ => rng.NextDouble() };
                loci.Add(new OncologyAnalyzer.AlleleSpecificLocus(c.ToString(), 1000L * (i + 1), logR, baf));
            }
        }

        return loci;
    }

    [Test]
    [Description("F59: NaN logR/BAF (R NA) anywhere in the germline-aware ascat.aspcf input never throws; every probe gets a finite segmented logR unless its whole level is missing, segmented BAF is NaN or in [0, 1], and segments tile the loci.")]
    public void SegmentAlleleSpecificAspcf_RandomMissingValues_WellFormed()
    {
        for (int seed = 0; seed < 60; seed++)
        {
            var rng = new Random(seed);
            var loci = NoisyLoci(rng, rng.Next(1, 4), rng.Next(1, 60), naRate: 0.15);
            var het = loci.Select(_ => rng.NextDouble() < 0.7).ToArray();

            var result = OncologyAnalyzer.SegmentAlleleSpecificAspcf(loci, het);

            Assert.Multiple(() =>
            {
                Assert.That(result.SegmentedLogR, Has.Count.EqualTo(loci.Count), $"seed {seed}");
                Assert.That(result.SegmentedBaf.All(b => double.IsNaN(b) || (b >= 0 && b <= 1)), Is.True, $"BAF range, seed {seed}");
                Assert.That(result.SegmentedLogR.All(v => !double.IsInfinity(v)), Is.True, $"no infinite logR, seed {seed}");
                Assert.That(result.Segments.Sum(s => s.LocusCount), Is.EqualTo(loci.Count), $"tiling, seed {seed}");
            });
        }
    }

    [Test]
    [Description("F53/F59: random multi-sample input with missing values either segments every sample (one result per sample, per-probe arrays aligned) or throws ArgumentException for a single-probe chromosome part.")]
    public void SegmentAlleleSpecificAsMultiPcf_RandomMissingValues_WellFormed()
    {
        for (int seed = 0; seed < 40; seed++)
        {
            var rng = new Random(1000 + seed);
            int chromosomes = rng.Next(1, 3), probes = rng.Next(2, 40), samplesCount = rng.Next(2, 4);
            var template = NoisyLoci(rng, chromosomes, probes, naRate: 0);
            var samples = Enumerable.Range(0, samplesCount)
                .Select(_ => (IReadOnlyList<OncologyAnalyzer.AlleleSpecificLocus>)template.Select(l => l with
                {
                    LogR = rng.NextDouble() < 0.1 ? double.NaN : l.LogR + (rng.NextDouble() - 0.5) * 0.1,
                    BAF = rng.NextDouble() < 0.1 ? double.NaN : l.BAF,
                }).ToList())
                .ToList();

            var result = OncologyAnalyzer.SegmentAlleleSpecificAsMultiPcf(samples);

            Assert.That(result, Has.Count.EqualTo(samplesCount), $"seed {seed}");
            foreach (var r in result)
            {
                Assert.Multiple(() =>
                {
                    Assert.That(r.SegmentedLogR, Has.Count.EqualTo(template.Count), $"seed {seed}");
                    Assert.That(r.SegmentedBaf.All(b => double.IsNaN(b) || (b >= 0.5 && b <= 1)), Is.True, $"mirrored BAF, seed {seed}");
                });
            }
        }
    }

    [Test]
    [Description("F53 (R: 'incorrect number of subscripts on matrix'): one sample, differing probe sets, a wrong genotype count and a single-probe chromosome part are ArgumentExceptions.")]
    public void SegmentAlleleSpecificAsMultiPcf_ContractViolations_Throw()
    {
        var rng = new Random(7);
        var loci = NoisyLoci(rng, 2, 10, naRate: 0);
        var shifted = loci.Select(l => l with { Position = l.Position + 1 }).ToList();
        var singleProbePart = loci.Take(10).Append(loci[10]).ToList();

        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentException>(() => OncologyAnalyzer.SegmentAlleleSpecificAsMultiPcf(new[] { loci }));
            Assert.Throws<ArgumentException>(() => OncologyAnalyzer.SegmentAlleleSpecificAsMultiPcf(new[] { loci, shifted }));
            Assert.Throws<ArgumentException>(() => OncologyAnalyzer.SegmentAlleleSpecificAsMultiPcf(new[] { loci, loci }, new bool[3]));
            Assert.Throws<ArgumentException>(() => OncologyAnalyzer.SegmentAlleleSpecificAsMultiPcf(new[] { singleProbePart, singleProbePart }));
            Assert.Throws<ArgumentNullException>(() => OncologyAnalyzer.SegmentAlleleSpecificAsMultiPcf(null!));
        });
    }

    [Test]
    [Description("F50–F54: extreme log2 values (±∞, ±1e6) are capped at ±1.5 by default, so the ziggurat deconstruction stays finite; without a cap +∞ is an ArgumentException; a noise filter that removes every sample is an ArgumentException.")]
    public void DeconstructZiggurat_ExtremeValues_CappedOrRejected()
    {
        var layout = new[] { new OncologyAnalyzer.ZigguratChromosome("1", 3, 4), new OncologyAnalyzer.ZigguratChromosome("2", 0, 5) };
        double[] extremes = { double.PositiveInfinity, double.NegativeInfinity, 1e6, -1e6, 0.0, 1e-300 };
        for (int seed = 0; seed < 40; seed++)
        {
            var rng = new Random(2000 + seed);
            var sample = new List<OncologyAnalyzer.ZigguratSegment>
            {
                new("1", 1, 3, extremes[rng.Next(extremes.Length)]),
                new("1", 4, 7, extremes[rng.Next(extremes.Length)]),
                new("2", 1, 5, extremes[rng.Next(extremes.Length)]),
            };

            var result = OncologyAnalyzer.DeconstructZiggurat(layout, new[] { sample });
            Assert.That(result.Events.All(e => double.IsFinite(e.Amplitude) && e.Amplitude > 0 && double.IsFinite(e.StartLevel) && double.IsFinite(e.EndLevel)),
                Is.True, $"seed {seed}");
        }

        var infinite = new[] { new OncologyAnalyzer.ZigguratSegment("1", 1, 7, double.PositiveInfinity), new OncologyAnalyzer.ZigguratSegment("2", 1, 5, 0.0) };
        Assert.Throws<ArgumentException>(() => OncologyAnalyzer.DeconstructZiggurat(layout, new[] { infinite },
            OncologyAnalyzer.ZigguratOptions.Default with { Cap = null }));
        var twoSegments = new[] { new OncologyAnalyzer.ZigguratSegment("1", 1, 7, 0.5), new OncologyAnalyzer.ZigguratSegment("2", 1, 5, -0.5) };
        Assert.Throws<ArgumentException>(() => OncologyAnalyzer.DeconstructZiggurat(layout, new[] { twoSegments },
            OncologyAnalyzer.ZigguratOptions.Default with { MaxSegmentsPerSample = 1 }));
    }

    [Test]
    [Description("F28/F29/F66: the shared copy-number guards reject purity ∉ (0, 1] and non-finite / non-positive ploidy (NaN included); NaN log2 propagates, −∞ gives 0 copies, and the purity path never returns a negative copy number.")]
    public void CopyNumberMath_BoundaryInputs()
    {
        double[] badPurity = { 0.0, -0.1, 1.0000001, double.NaN, double.PositiveInfinity };
        double[] badPloidy = { 0.0, -2.0, double.NaN, double.PositiveInfinity };
        Assert.Multiple(() =>
        {
            foreach (double p in badPurity)
            {
                Assert.Throws<ArgumentOutOfRangeException>(() => CopyNumberMath.ValidatePurity(p), $"purity {p}");
                Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.Log2RatioToCopyNumber(0.3, 2.0, p), $"purity {p}");
            }

            foreach (double q in badPloidy)
            {
                Assert.Throws<ArgumentOutOfRangeException>(() => CopyNumberMath.ValidatePloidy(q), $"ploidy {q}");
            }

            Assert.That(CopyNumberMath.Log2RatioToAbsolute(double.NaN, 2.0), Is.NaN);
            Assert.That(CopyNumberMath.Log2RatioToAbsolute(double.NegativeInfinity, 2.0), Is.EqualTo(0.0));
            Assert.That(CopyNumberMath.AbsoluteToLog2Ratio(0.0, 2.0), Is.EqualTo(Math.Log2(1e-3)));
        });

        var rng = new Random(11);
        for (int i = 0; i < 2000; i++)
        {
            double v = (rng.NextDouble() - 0.5) * 40;
            double purity = Math.Max(1e-6, rng.NextDouble());
            double n = CopyNumberMath.Log2RatioToAbsolute(v, rng.Next(1, 5), rng.Next(0, 5), purity);
            Assert.That(n, Is.GreaterThanOrEqualTo(0.0), $"v={v}, p={purity}");
        }
    }

    [Test]
    [Description("F39: the Student t CDF and the incomplete beta reject out-of-domain arguments and stay in [0, 1] (never NaN) for extreme in-domain arguments.")]
    public void StudentTAndIncompleteBeta_ExtremeArguments()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentOutOfRangeException>(() => StatisticsHelper.StudentTCdf(double.NaN, 3));
            Assert.Throws<ArgumentOutOfRangeException>(() => StatisticsHelper.StudentTCdf(1.0, 0));
            Assert.Throws<ArgumentOutOfRangeException>(() => StatisticsHelper.StudentTCdf(1.0, double.PositiveInfinity));
            Assert.Throws<ArgumentOutOfRangeException>(() => StatisticsHelper.RegularizedIncompleteBeta(-0.1, 1, 1));
            Assert.Throws<ArgumentOutOfRangeException>(() => StatisticsHelper.RegularizedIncompleteBeta(0.5, 0, 1));
            Assert.Throws<ArgumentOutOfRangeException>(() => StatisticsHelper.RegularizedIncompleteBeta(0.5, 1, double.NaN));
        });

        var rng = new Random(13);
        double[] ts = { 1e-300, 1e-8, 1, 30, 1e3, 1e8, 1e300, double.MaxValue };
        double[] nus = { 1e-3, 0.5, 1, 2.5, 30, 1e4, 1e6, 1e9 };
        foreach (double t in ts)
        {
            foreach (double nu in nus)
            {
                foreach (double sign in new[] { -1.0, 1.0 })
                {
                    double f = StatisticsHelper.StudentTCdf(sign * t, nu);
                    Assert.That(f, Is.InRange(0.0, 1.0), $"t={sign * t}, nu={nu}");
                }
            }
        }

        for (int i = 0; i < 2000; i++)
        {
            double x = rng.NextDouble();
            double a = Math.Exp((rng.NextDouble() - 0.3) * 20), b = Math.Exp((rng.NextDouble() - 0.3) * 20);
            Assert.That(StatisticsHelper.RegularizedIncompleteBeta(x, a, b), Is.InRange(0.0, 1.0), $"x={x}, a={a}, b={b}");
        }
    }
}
