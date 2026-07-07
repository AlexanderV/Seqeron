using NUnit.Framework;
using Seqeron.Mcp.Chromosome.Tools;

namespace Seqeron.Mcp.Chromosome.Tests;

/// <summary>
/// Tests for <c>detect_ploidy</c>. Expected values follow ChromosomeAnalyzer.DetectPloidy:
/// ploidy = round(median/expectedDiploidDepth * 2) clamped to [1, 8]; confidence = 1 - 2*frac.
/// </summary>
[TestFixture]
public class DetectPloidyTests
{
    [Test]
    public void DetectPloidy_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => ChromosomeTools.DetectPloidy(new List<double> { 1.0 }));
        Assert.Throws<ArgumentException>(() => ChromosomeTools.DetectPloidy(null!));
        Assert.Throws<ArgumentException>(() => ChromosomeTools.DetectPloidy(new List<double> { 1.0 }, expectedDiploidDepth: 0));
        Assert.Throws<ArgumentException>(() => ChromosomeTools.DetectPloidy(new List<double> { 1.0 }, expectedDiploidDepth: -1));
    }

    [Test]
    public void DetectPloidy_Binding_InvokesSuccessfully()
    {
        var diploid = ChromosomeTools.DetectPloidy(new List<double> { 1.0, 1.0, 1.0 }, 1.0);
        var tetraploid = ChromosomeTools.DetectPloidy(new List<double> { 2.0, 2.0, 2.0 }, 1.0);

        Assert.Multiple(() =>
        {
            Assert.That(diploid.PloidyLevel, Is.EqualTo(2));
            Assert.That(diploid.Confidence, Is.EqualTo(1.0).Within(1e-9));
            Assert.That(tetraploid.PloidyLevel, Is.EqualTo(4));
            Assert.That(tetraploid.Confidence, Is.EqualTo(1.0).Within(1e-9));
        });
    }
}
