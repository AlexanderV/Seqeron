using NUnit.Framework;
using Seqeron.Mcp.Population.Tools;

namespace Seqeron.Mcp.Population.Tests;

[TestFixture]
public class HardyWeinbergTestTests
{
    [Test]
    public void HardyWeinbergTest_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => PopulationTools.HardyWeinbergTest("rs1", 25, 50, 25));
        // Zero samples is valid input; returns in-equilibrium with p-value 1.
        Assert.DoesNotThrow(() => PopulationTools.HardyWeinbergTest("rs0", 0, 0, 0));
    }

    [Test]
    public void HardyWeinbergTest_Binding_FordsMothData_InEquilibrium()
    {
        // Ford's scarlet tiger moth (Wikipedia): AA=1469, Aa=138, aa=5.
        //   E = (1467.40, 141.21, 3.40); χ² ≈ 0.8309; in equilibrium.
        var r = PopulationTools.HardyWeinbergTest("FORD_MOTH", 1469, 138, 5);

        Assert.Multiple(() =>
        {
            Assert.That(r.VariantId, Is.EqualTo("FORD_MOTH"));
            Assert.That(r.ExpectedAA, Is.EqualTo(1467.40).Within(0.1));
            Assert.That(r.ExpectedAa, Is.EqualTo(141.21).Within(0.1));
            Assert.That(r.Expectedaa, Is.EqualTo(3.40).Within(0.1));
            Assert.That(r.ChiSquare, Is.EqualTo(0.8309).Within(0.01));
            Assert.That(r.PValue, Is.GreaterThan(0.05));
            Assert.That(r.InEquilibrium, Is.True);
        });
    }

    [Test]
    public void HardyWeinbergTest_Binding_ExcessHeterozygotes_DeviatesFromHwe()
    {
        // AA=10, Aa=80, aa=10 → E=(25,50,25); χ² = 9+18+9 = 36; not in equilibrium.
        var r = PopulationTools.HardyWeinbergTest("EXCESS_HET", 10, 80, 10);

        Assert.Multiple(() =>
        {
            Assert.That(r.ExpectedAA, Is.EqualTo(25.0).Within(1e-2));
            Assert.That(r.ExpectedAa, Is.EqualTo(50.0).Within(1e-2));
            Assert.That(r.Expectedaa, Is.EqualTo(25.0).Within(1e-2));
            Assert.That(r.ChiSquare, Is.EqualTo(36.0).Within(1e-2));
            Assert.That(r.PValue, Is.LessThan(0.05));
            Assert.That(r.InEquilibrium, Is.False);
        });
    }
}
