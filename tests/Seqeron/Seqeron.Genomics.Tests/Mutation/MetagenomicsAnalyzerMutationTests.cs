using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Infrastructure;
using Seqeron.Genomics.Metagenomics;

namespace Seqeron.Genomics.Tests.Mutation;

/// <summary>
/// Targeted mutation-killing tests for MetagenomicsAnalyzer.cs (checklist 04 rows 53-57,
/// 194-197). The differential-abundance / Welch t-test path (DifferentialAbundance +
/// CalculateTTestPValue) was entirely untested, leaving a large block of Stryker survivors.
/// These pin the published Welch two-sample t-test (normal-approximation p-value) and the
/// log2 fold-change / significance rule against independent ground truth.
/// </summary>
[TestFixture]
public class MetagenomicsAnalyzerMutationTests
{
    // Independent Welch two-sample t-test with normal-approximation p-value (the theory):
    //   var_i = Σ(x−mean)²/(n−1);  se = √(var1/n1 + var2/n2);  t = |m1−m2|/se;
    //   p = 2·(1 − Φ(t)).  Φ = StatisticsHelper.NormalCDF (Infrastructure, not mutated).
    private static double ExpectedWelchP(double[] g1, double[] g2)
    {
        double m1 = g1.Average(), m2 = g2.Average();
        double v1 = g1.Sum(x => (x - m1) * (x - m1)) / (g1.Length - 1);
        double v2 = g2.Sum(x => (x - m2) * (x - m2)) / (g2.Length - 1);
        double se = Math.Sqrt(v1 / g1.Length + v2 / g2.Length);
        double t = Math.Abs(m1 - m2) / se;
        return 2 * (1 - StatisticsHelper.NormalCDF(t));
    }

    private static IReadOnlyDictionary<string, double>[] Samples(string taxon, params double[] values) =>
        values.Select(v => (IReadOnlyDictionary<string, double>)
            new Dictionary<string, double> { [taxon] = v }).ToArray();

    [Test]
    public void DifferentialAbundance_WelchTTest_PValueAndFoldChangeAreExact()
    {
        // taxon "T": condition1 = {1,2,3} (mean 2), condition2 = {4,5,6} (mean 5).
        var c1 = Samples("T", 1, 2, 3);
        var c2 = Samples("T", 4, 5, 6);

        var result = MetagenomicsAnalyzer.DifferentialAbundance(c1, c2).Single();

        result.Taxon.Should().Be("T");
        result.FoldChange.Should().BeApproximately(Math.Log2(5.0 / 2.0), 1e-9, "log2(mean2/mean1)");
        result.PValue.Should().BeApproximately(ExpectedWelchP(new[] { 1.0, 2, 3 }, new[] { 4.0, 5, 6 }), 1e-9);
        result.Significant.Should().BeTrue("p < 0.05 and |log2FC| > 1");
    }

    [Test]
    public void DifferentialAbundance_EqualMeans_NotSignificant_PValueOne()
    {
        // Equal means → t = 0 → p = 2·(1 − Φ(0)) = 1.0; fold change 1 → log2FC 0 → not significant.
        var c1 = Samples("T", 2, 4, 6); // mean 4
        var c2 = Samples("T", 3, 4, 5); // mean 4
        var result = MetagenomicsAnalyzer.DifferentialAbundance(c1, c2).Single();

        result.PValue.Should().BeApproximately(1.0, 1e-6); // t=0 ⇒ 2·(1−Φ(0)); Φ via Erf approx
        result.FoldChange.Should().BeApproximately(0.0, 1e-9);
        result.Significant.Should().BeFalse();
    }

    [Test]
    public void DifferentialAbundance_IdenticalConstantGroups_PValueOne()
    {
        // var1 == var2 == 0 and mean1 == mean2 → p = 1.0 (degenerate-variance guard).
        var c1 = Samples("T", 5, 5, 5);
        var c2 = Samples("T", 5, 5, 5);
        MetagenomicsAnalyzer.DifferentialAbundance(c1, c2).Single().PValue.Should().Be(1.0);
    }

    [Test]
    public void DifferentialAbundance_ConstantGroupsDifferentMeans_PValueZero()
    {
        // var1 == var2 == 0 but mean1 != mean2 → p = 0.0 (separable, zero-variance branch).
        var c1 = Samples("T", 2, 2, 2);
        var c2 = Samples("T", 9, 9, 9);
        MetagenomicsAnalyzer.DifferentialAbundance(c1, c2).Single().PValue.Should().Be(0.0);
    }

    [Test]
    public void DifferentialAbundance_FewerThanTwoSamples_PValueOne()
    {
        // group with <2 observations → t-test returns 1.0.
        var c1 = Samples("T", 5);        // single sample
        var c2 = Samples("T", 1, 9);
        MetagenomicsAnalyzer.DifferentialAbundance(c1, c2).Single().PValue.Should().Be(1.0);
    }

    [Test]
    public void DifferentialAbundance_ZeroBaselineMean_FoldChangeIsPositiveInfinityDirection()
    {
        // mean1 == 0, mean2 > 0 → fold change +∞ → log2 of (clamped) large value is positive/large.
        var c1 = Samples("T", 0, 0, 0);
        var c2 = Samples("T", 4, 5, 6);
        var result = MetagenomicsAnalyzer.DifferentialAbundance(c1, c2).Single();

        result.FoldChange.Should().Be(double.PositiveInfinity);
    }

    [Test]
    public void DifferentialAbundance_EmptyCondition_YieldsNothing()
    {
        var c2 = Samples("T", 1, 2, 3);
        MetagenomicsAnalyzer.DifferentialAbundance(System.Array.Empty<IReadOnlyDictionary<string, double>>(), c2)
            .Should().BeEmpty();
    }
}
