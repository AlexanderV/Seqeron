using NUnit.Framework;
using Seqeron.Mcp.Population.Tools;

namespace Seqeron.Mcp.Population.Tests;

[TestFixture]
public class IntegratedHaplotypeScoreTests
{
    [Test]
    public void IntegratedHaplotypeScore_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => PopulationTools.IntegratedHaplotypeScore(
            new[] { 1.0, 0.5 }, new[] { 1.0, 1.0 }, new[] { 0, 1000 }));

        // Fewer than 2 points → 0 (not an error).
        Assert.DoesNotThrow(() => PopulationTools.IntegratedHaplotypeScore(
            new[] { 1.0 }, new[] { 1.0 }, new[] { 0 }));
    }

    [Test]
    public void IntegratedHaplotypeScore_Binding_InvokesSuccessfully()
    {
        // iHH0 = ½(1+0.5)·1000 = 750, iHH1 = ½(1+1)·1000 = 1000.
        // iHS = ln(1000/750) = ln(4/3) ≈ 0.28768207245178085.
        var result = PopulationTools.IntegratedHaplotypeScore(
            new[] { 1.0, 0.5 }, new[] { 1.0, 1.0 }, new[] { 0, 1000 });

        Assert.That(result.Ihs, Is.EqualTo(System.Math.Log(4.0 / 3.0)).Within(1e-12));
        Assert.That(result.Ihs, Is.GreaterThan(0));
    }

    [Test]
    public void IntegratedHaplotypeScore_Binding_BalancedAndInsufficient()
    {
        Assert.Multiple(() =>
        {
            // Balanced EHH decay → iHH0 = iHH1 → iHS = 0.
            var balanced = PopulationTools.IntegratedHaplotypeScore(
                new[] { 1.0, 0.8, 0.5, 0.2, 0.1 },
                new[] { 1.0, 0.8, 0.5, 0.2, 0.1 },
                new[] { 0, 1000, 2000, 3000, 4000 });
            Assert.That(balanced.Ihs, Is.EqualTo(0.0).Within(1e-12));

            // Insufficient data (single point) → 0.
            var insufficient = PopulationTools.IntegratedHaplotypeScore(
                new[] { 1.0 }, new[] { 1.0 }, new[] { 0 });
            Assert.That(insufficient.Ihs, Is.EqualTo(0.0));
        });
    }
}
