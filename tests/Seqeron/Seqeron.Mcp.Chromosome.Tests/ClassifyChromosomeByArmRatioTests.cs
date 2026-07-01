using NUnit.Framework;
using Seqeron.Mcp.Chromosome.Tools;

namespace Seqeron.Mcp.Chromosome.Tests;

/// <summary>
/// Tests for <c>classify_chromosome_by_arm_ratio</c>. ChromosomeAnalyzer.ClassifyChromosomeByArmRatio
/// normalizes the ratio to r = long/short >= 1, then per Levan et al. (1964):
/// r <= 1.7 Metacentric, r <= 3.0 Submetacentric, r < 7.0 Subtelocentric, else Acrocentric;
/// ratio <= 0 -> Telocentric.
/// </summary>
[TestFixture]
public class ClassifyChromosomeByArmRatioTests
{
    [Test]
    public void ClassifyChromosomeByArmRatio_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => ChromosomeTools.ClassifyChromosomeByArmRatio(1.0));
        // Non-positive ratio is documented-valid (Telocentric), not an error.
        Assert.DoesNotThrow(() => ChromosomeTools.ClassifyChromosomeByArmRatio(0.0));
        Assert.Throws<ArgumentException>(() => ChromosomeTools.ClassifyChromosomeByArmRatio(double.NaN));
        Assert.Throws<ArgumentException>(() => ChromosomeTools.ClassifyChromosomeByArmRatio(double.PositiveInfinity));
    }

    [Test]
    public void ClassifyChromosomeByArmRatio_Binding_InvokesSuccessfully()
    {
        Assert.Multiple(() =>
        {
            Assert.That(ChromosomeTools.ClassifyChromosomeByArmRatio(1.0).Classification, Is.EqualTo("Metacentric"));
            Assert.That(ChromosomeTools.ClassifyChromosomeByArmRatio(1.7).Classification, Is.EqualTo("Metacentric"));
            Assert.That(ChromosomeTools.ClassifyChromosomeByArmRatio(2.5).Classification, Is.EqualTo("Submetacentric"));
            Assert.That(ChromosomeTools.ClassifyChromosomeByArmRatio(5.0).Classification, Is.EqualTo("Subtelocentric"));
            Assert.That(ChromosomeTools.ClassifyChromosomeByArmRatio(10.0).Classification, Is.EqualTo("Acrocentric"));
            Assert.That(ChromosomeTools.ClassifyChromosomeByArmRatio(0.0).Classification, Is.EqualTo("Telocentric"));
            // Reciprocal is normalized: 0.4 -> r = 2.5 -> Submetacentric.
            Assert.That(ChromosomeTools.ClassifyChromosomeByArmRatio(0.4).Classification, Is.EqualTo("Submetacentric"));
        });
    }
}
