using NUnit.Framework;
using Seqeron.Mcp.Population.Tools;

namespace Seqeron.Mcp.Population.Tests;

[TestFixture]
public class LinkageDisequilibriumTests
{
    private static GenotypePairItem G(int g1, int g2) => new(g1, g2);

    [Test]
    public void LinkageDisequilibrium_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => PopulationTools.LinkageDisequilibrium(
            "V1", "V2", new[] { G(0, 0), G(1, 1) }, 1000));

        // Empty genotypes is valid input; returns zero LD.
        Assert.DoesNotThrow(() => PopulationTools.LinkageDisequilibrium(
            "V1", "V2", Array.Empty<GenotypePairItem>(), 1000));
    }

    [Test]
    public void LinkageDisequilibrium_Binding_PerfectLD()
    {
        // Identical genotype vectors → r² = 1.0, D' = 1.0 (clamped, Lewontin 1964).
        var ld = PopulationTools.LinkageDisequilibrium(
            "V1", "V2",
            new[] { G(0, 0), G(0, 0), G(1, 1), G(1, 1), G(2, 2), G(2, 2) },
            1000);

        Assert.Multiple(() =>
        {
            Assert.That(ld.Variant1, Is.EqualTo("V1"));
            Assert.That(ld.Variant2, Is.EqualTo("V2"));
            Assert.That(ld.RSquared, Is.EqualTo(1.0).Within(1e-10));
            Assert.That(ld.DPrime, Is.EqualTo(1.0).Within(1e-10));
            Assert.That(ld.Distance, Is.EqualTo(1000));
        });
    }

    [Test]
    public void LinkageDisequilibrium_Binding_NoLD_BalancedDesign()
    {
        // Balanced 3×3 design → Cov = 0 → r² = 0, D' = 0.
        var ld = PopulationTools.LinkageDisequilibrium(
            "V1", "V2",
            new[]
            {
                G(0, 0), G(0, 1), G(0, 2),
                G(1, 0), G(1, 1), G(1, 2),
                G(2, 0), G(2, 1), G(2, 2),
            },
            1000);

        Assert.Multiple(() =>
        {
            Assert.That(ld.RSquared, Is.EqualTo(0.0).Within(1e-10));
            Assert.That(ld.DPrime, Is.EqualTo(0.0).Within(1e-10));
        });
    }
}
