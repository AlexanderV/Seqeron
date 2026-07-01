using NUnit.Framework;
using Seqeron.Mcp.Population.Tools;

namespace Seqeron.Mcp.Population.Tests;

[TestFixture]
public class TajimasDTests
{
    [Test]
    public void TajimasD_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => PopulationTools.TajimasD(2.0, 4, 5));
        // n < 3 or S = 0 are valid inputs that return 0.
        Assert.DoesNotThrow(() => PopulationTools.TajimasD(5.0, 10, 2));
    }

    [Test]
    public void TajimasD_Binding_WikipediaExample()
    {
        // Wikipedia example: k̂ = 2.0, S = 4, n = 5 → D ≈ 0.273.
        var d = PopulationTools.TajimasD(2.0, 4, 5);
        Assert.That(d.TajimasD, Is.EqualTo(0.273).Within(0.005));
    }

    [Test]
    public void TajimasD_Binding_NeutralAndDegenerate()
    {
        Assert.Multiple(() =>
        {
            // Neutral: k̂ = S/a₁ → numerator 0 → D = 0.
            double a1 = 0;
            for (int i = 1; i < 50; i++) a1 += 1.0 / i;
            double kHat = 100.0 / a1;
            Assert.That(PopulationTools.TajimasD(kHat, 100, 50).TajimasD, Is.EqualTo(0.0).Within(1e-3));

            // n < 3 → 0.
            Assert.That(PopulationTools.TajimasD(5.0, 10, 2).TajimasD, Is.EqualTo(0.0));
            // S = 0 → 0.
            Assert.That(PopulationTools.TajimasD(0.0, 0, 50).TajimasD, Is.EqualTo(0.0));
        });
    }
}
