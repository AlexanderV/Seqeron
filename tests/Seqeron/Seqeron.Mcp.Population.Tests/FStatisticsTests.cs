using NUnit.Framework;
using Seqeron.Mcp.Population.Tools;

namespace Seqeron.Mcp.Population.Tests;

[TestFixture]
public class FStatisticsTests
{
    [Test]
    public void FStatistics_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => PopulationTools.FStatistics(
            "Pop1", "Pop2",
            new[] { new VariantDataItem(20, 50, 25, 50, 0.4, 0.5) }));

        // Empty variant data is valid input; returns zeros.
        Assert.DoesNotThrow(() => PopulationTools.FStatistics(
            "Pop1", "Pop2", Array.Empty<VariantDataItem>()));
    }

    [Test]
    public void FStatistics_Binding_InvokesSuccessfully()
    {
        // Two-locus hand calculation (POP-FST-001):
        //   H_I = 0.45, H_S = 0.475, H_T = 0.4875
        //   Fis = 1/19, Fit = 1/13, Fst = 1/39
        var stats = PopulationTools.FStatistics(
            "Pop1", "Pop2",
            new[]
            {
                new VariantDataItem(20, 50, 25, 50, 0.4, 0.5),
                new VariantDataItem(30, 50, 15, 50, 0.5, 0.3),
            });

        Assert.Multiple(() =>
        {
            Assert.That(stats.Population1, Is.EqualTo("Pop1"));
            Assert.That(stats.Population2, Is.EqualTo("Pop2"));
            Assert.That(stats.Fis, Is.EqualTo(1.0 / 19.0).Within(1e-10));
            Assert.That(stats.Fit, Is.EqualTo(1.0 / 13.0).Within(1e-10));
            Assert.That(stats.Fst, Is.EqualTo(1.0 / 39.0).Within(1e-10));

            // Wright partition identity holds exactly.
            Assert.That(1 - stats.Fit,
                Is.EqualTo((1 - stats.Fis) * (1 - stats.Fst)).Within(1e-10));
        });
    }

    [Test]
    public void FStatistics_Binding_ExcessHeterozygosity_NegativeFis()
    {
        // Single locus, HetObs 60/100 and 80/100, p=0.3 and 0.7.
        //   Fis = -2/3, Fit = -2/5, Fst = 4/25.
        var stats = PopulationTools.FStatistics(
            "Pop1", "Pop2",
            new[] { new VariantDataItem(60, 100, 80, 100, 0.3, 0.7) });

        Assert.Multiple(() =>
        {
            Assert.That(stats.Fis, Is.EqualTo(-2.0 / 3.0).Within(1e-10));
            Assert.That(stats.Fit, Is.EqualTo(-2.0 / 5.0).Within(1e-10));
            Assert.That(stats.Fst, Is.EqualTo(4.0 / 25.0).Within(1e-10));
        });
    }
}
