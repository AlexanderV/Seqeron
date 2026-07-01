using NUnit.Framework;
using Seqeron.Mcp.Chromosome.Tools;

namespace Seqeron.Mcp.Chromosome.Tests;

/// <summary>
/// Tests for <c>arm_ratio</c>. ChromosomeAnalyzer.CalculateArmRatio returns p/q where
/// p = centromerePosition and q = chromosomeLength - centromerePosition.
/// </summary>
[TestFixture]
public class ArmRatioTests
{
    [Test]
    public void ArmRatio_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => ChromosomeTools.ArmRatio(40, 100));
        Assert.Throws<ArgumentException>(() => ChromosomeTools.ArmRatio(0, 100));
        Assert.Throws<ArgumentException>(() => ChromosomeTools.ArmRatio(-1, 100));
        Assert.Throws<ArgumentException>(() => ChromosomeTools.ArmRatio(40, 0));
        Assert.Throws<ArgumentException>(() => ChromosomeTools.ArmRatio(100, 100));
        Assert.Throws<ArgumentException>(() => ChromosomeTools.ArmRatio(120, 100));
    }

    [Test]
    public void ArmRatio_Binding_InvokesSuccessfully()
    {
        // p=40, q=60 -> 40/60 = 0.6667
        Assert.That(ChromosomeTools.ArmRatio(40, 100).ArmRatio, Is.EqualTo(40.0 / 60.0).Within(1e-9));
        // p=50, q=50 -> 1.0 (metacentric geometry)
        Assert.That(ChromosomeTools.ArmRatio(50, 100).ArmRatio, Is.EqualTo(1.0).Within(1e-9));
        // p=25, q=75 -> 1/3
        Assert.That(ChromosomeTools.ArmRatio(25, 100).ArmRatio, Is.EqualTo(25.0 / 75.0).Within(1e-9));
    }
}
