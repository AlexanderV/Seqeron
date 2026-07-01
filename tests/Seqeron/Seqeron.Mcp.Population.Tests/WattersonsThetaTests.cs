using NUnit.Framework;
using Seqeron.Mcp.Population.Tools;

namespace Seqeron.Mcp.Population.Tests;

[TestFixture]
public class WattersonsThetaTests
{
    [Test]
    public void WattersonsTheta_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => PopulationTools.WattersonsTheta(10, 10, 1000));
        // n < 2 or L ≤ 0 are valid inputs that return 0.
        Assert.DoesNotThrow(() => PopulationTools.WattersonsTheta(5, 1, 100));
        Assert.DoesNotThrow(() => PopulationTools.WattersonsTheta(5, 10, 0));
    }

    [Test]
    public void WattersonsTheta_Binding_InvokesSuccessfully()
    {
        // S=10, n=10, L=1000: a₁ = Σ_{i=1}^{9} 1/i ≈ 2.82897 → θ ≈ 0.00353.
        var theta = PopulationTools.WattersonsTheta(10, 10, 1000);
        Assert.That(theta.Theta, Is.EqualTo(0.00353).Within(0.0005));

        // Exact: n=2 → a₁ = 1 → θ = S/L = 5/100 = 0.05.
        var minN = PopulationTools.WattersonsTheta(5, 2, 100);
        Assert.That(minN.Theta, Is.EqualTo(0.05).Within(1e-10));
    }

    [Test]
    public void WattersonsTheta_Binding_ExactHarmonicAndDegenerate()
    {
        Assert.Multiple(() =>
        {
            // a₁(3) = 3/2 → θ = 10/(1.5×100) = 1/15.
            Assert.That(PopulationTools.WattersonsTheta(10, 3, 100).Theta,
                Is.EqualTo(10.0 / (1.5 * 100)).Within(1e-10));
            // a₁(5) = 25/12 → θ = 10/((25/12)×100).
            Assert.That(PopulationTools.WattersonsTheta(10, 5, 100).Theta,
                Is.EqualTo(10.0 / (25.0 / 12.0 * 100)).Within(1e-10));

            // Degenerate cases → 0.
            Assert.That(PopulationTools.WattersonsTheta(5, 1, 100).Theta, Is.EqualTo(0.0));
            Assert.That(PopulationTools.WattersonsTheta(5, 10, 0).Theta, Is.EqualTo(0.0));
            Assert.That(PopulationTools.WattersonsTheta(0, 10, 1000).Theta, Is.EqualTo(0.0));
        });
    }
}
