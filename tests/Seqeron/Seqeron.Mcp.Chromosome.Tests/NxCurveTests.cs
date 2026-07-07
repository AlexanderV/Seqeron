using NUnit.Framework;
using Seqeron.Mcp.Chromosome.Tools;

namespace Seqeron.Mcp.Chromosome.Tests;

/// <summary>
/// Tests for <c>nx_curve</c>. GenomeAssemblyAnalyzer.CalculateNxCurve computes Nx/Lx over the
/// requested thresholds (default 10..90 step 10), each via CalculateNx.
/// </summary>
[TestFixture]
public class NxCurveTests
{
    [Test]
    public void NxCurve_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => ChromosomeTools.NxCurve(new List<int> { 100, 50 }));
        Assert.Throws<ArgumentException>(() => ChromosomeTools.NxCurve(null!));
    }

    [Test]
    public void NxCurve_DefaultThresholds_ReturnsNineDeciles()
    {
        var result = ChromosomeTools.NxCurve(new List<int> { 100, 90, 80, 70, 60, 50, 40, 30, 20, 10 });

        Assert.That(result.Items, Has.Count.EqualTo(9));
        Assert.That(result.Items.Select(i => i.Threshold), Is.EqualTo(new[] { 10, 20, 30, 40, 50, 60, 70, 80, 90 }));

        // The N50 entry matches the single-threshold computation (Nx=70, Lx=4).
        var n50 = result.Items.Single(i => i.Threshold == 50);
        Assert.Multiple(() =>
        {
            Assert.That(n50.Nx, Is.EqualTo(70));
            Assert.That(n50.Lx, Is.EqualTo(4));
        });
    }

    [Test]
    public void NxCurve_ExplicitThresholds_AreSortedAscending()
    {
        var result = ChromosomeTools.NxCurve(new List<int> { 100, 50 }, new[] { 90, 10 });
        Assert.That(result.Items.Select(i => i.Threshold), Is.EqualTo(new[] { 10, 90 }));
    }
}
