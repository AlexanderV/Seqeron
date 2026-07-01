using NUnit.Framework;
using Seqeron.Mcp.Population.Tools;

namespace Seqeron.Mcp.Population.Tests;

[TestFixture]
public class DiversityStatisticsTests
{
    [Test]
    public void DiversityStatistics_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => PopulationTools.DiversityStatistics(
            new[] { "ACGTACGT", "ACGTATGT", "ACGTACGA" }));

        // Single sequence is valid input; returns zeros (not an error).
        Assert.DoesNotThrow(() => PopulationTools.DiversityStatistics(new[] { "ACGT" }));
    }

    [Test]
    public void DiversityStatistics_Binding_InvokesSuccessfully()
    {
        // Wikipedia Tajima's D example: 5 sequences, length 20.
        // Source: https://en.wikipedia.org/wiki/Tajima%27s_D#Example
        //   π = 2.0/20 = 0.1 ; S = 4 ; θ_W = 4/(a₁·20), a₁ = 25/12 ≈ 0.096 ; D ≈ 0.273
        var wiki = PopulationTools.DiversityStatistics(new[]
        {
            "00000000000000000000",
            "00100000000010000010",
            "00000000000010000010",
            "00000010000000000010",
            "00000010000010000010",
        });

        Assert.Multiple(() =>
        {
            Assert.That(wiki.SampleSize, Is.EqualTo(5));
            Assert.That(wiki.SegregatingSites, Is.EqualTo(4));
            Assert.That(wiki.NucleotideDiversity, Is.EqualTo(0.1).Within(0.001));
            Assert.That(wiki.WattersonTheta, Is.EqualTo(0.096).Within(0.001));
            Assert.That(wiki.TajimasD, Is.EqualTo(0.273).Within(0.005));
        });
    }

    [Test]
    public void DiversityStatistics_Binding_ThreeSequences_ExactValues()
    {
        // Three sequences length 8, S = 2 (positions 5, 7).
        //   k̂ = 4/3, π = 4/(3×8) = 1/6
        //   a₁(3) = 3/2, θ_W = 2/(1.5×8) = 1/6
        //   k̂ = S/a₁ → D = 0
        //   H_exp = 1/9, H_obs = (3/2)·(1/9) = 1/6
        var stats = PopulationTools.DiversityStatistics(
            new[] { "ACGTACGT", "ACGTATGT", "ACGTACGA" });

        Assert.Multiple(() =>
        {
            Assert.That(stats.SampleSize, Is.EqualTo(3));
            Assert.That(stats.SegregatingSites, Is.EqualTo(2));
            Assert.That(stats.NucleotideDiversity, Is.EqualTo(1.0 / 6.0).Within(1e-4));
            Assert.That(stats.WattersonTheta, Is.EqualTo(1.0 / 6.0).Within(1e-4));
            Assert.That(stats.TajimasD, Is.EqualTo(0.0).Within(1e-3));
            Assert.That(stats.HeterozygosityObserved, Is.EqualTo(1.0 / 6.0).Within(1e-4));
            Assert.That(stats.HeterozygosityExpected, Is.EqualTo(1.0 / 9.0).Within(1e-4));
        });
    }
}
